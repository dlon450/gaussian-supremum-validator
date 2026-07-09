"""Unit tests for the non-solver infrastructure (run in a numpy/scipy env).

    /tmp/gsv_venv/bin/python tests/test_infra.py

No pytest dependency: plain asserts + a tiny runner. Solver-dependent behavior
(CCP_* solves, legacy-CSV equivalence) is covered by the regression fixtures,
which require Gurobi and run in the solver environment.
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gsv import rng, config, dgp as D, oracle as O, metrics as M
from gsv import validators as V
from gsv.config import DGP

TESTS = []
def test(fn): TESTS.append(fn); return fn


# ------------------------------- RNG ------------------------------- #
@test
def test_rng_deterministic_independent_order_free():
    a1 = rng.make_stream(2, "cfg", 7).normal(size=8)
    a2 = rng.make_stream(2, "cfg", 7).normal(size=8)
    assert np.array_equal(a1, a2), "same key must give identical draws"
    assert not np.allclose(a1, rng.make_stream(2, "cfg", 8).normal(size=8)), "rep id independence"
    assert not np.allclose(a1, rng.make_stream(2, "other", 7).normal(size=8)), "config id independence"
    _ = rng.make_stream(2, "cfg", 99).normal(size=1000)  # consume an unrelated stream
    assert np.array_equal(a1, rng.make_stream(2, "cfg", 7).normal(size=8)), "must be order-free"
    assert rng.stable_hash("x") == rng.stable_hash("x") and rng.stable_hash("x") != rng.stable_hash("y")


# ------------------------------ config ----------------------------- #
@test
def test_config_registry_values():
    leg, pap = config.get_config("legacy_so"), config.get_config("paper_so")
    assert (leg.alpha, leg.beta, leg.n_reps, leg.dgp.kind) == (0.05, 0.10, 100, "half_normal")
    assert (pap.alpha, pap.beta, pap.n_reps, pap.dgp.kind) == (0.10, 0.05, 1000, "gaussian")
    assert config.get_config("paper_dro_wasserstein").benchmark is None
    assert "SO_all" == config.get_config("paper_so").benchmark
    assert config.pilot("paper_so", n_reps=5).n_reps == 5

@test
def test_config_validation_rejects_bad():
    from dataclasses import replace
    for bad in [dict(alpha=1.5), dict(beta=0.0), dict(formulation="nope"),
                dict(validators=("ZZZ",)), dict(split_grid=(0.0, 0.5))]:
        try:
            replace(config.get_config("paper_so"), **bad).validate()
        except AssertionError:
            continue
        raise AssertionError(f"validation should reject {bad}")


# ------------------------------ dgp/oracle ------------------------- #
@test
def test_oracle_closedform_matches_large_sample():
    r = np.random.default_rng(0); d = 6
    g = DGP("gaussian", {"mu_scale": 0.8, "var_diag": 0.4, "corr": 0.15})
    mu, Sigma = D.moments(g, d); b = 2.0
    X = r.random((d, 4))
    cf = O.gaussian_feasibility(X, mu, Sigma, b)
    ls = O.large_sample_feasibility(X, D.sample(g, 2_000_000, d, r), b)
    assert np.max(np.abs(cf - ls)) < 3e-3, f"closed-form vs sample too far: {np.max(np.abs(cf-ls))}"

@test
def test_dgp_t_has_target_moments():
    r = np.random.default_rng(1); d = 5
    t = DGP("multivariate_t", {"mu_scale": 0.8, "var_diag": 0.4, "corr": 0.0, "df": 6})
    s = D.sample(t, 2_000_000, d, r)
    assert abs(s.mean() - 0.8) < 1e-2 and abs(s.var(axis=0).mean() - 0.4) < 2e-2

@test
def test_path_oracle():
    s = np.array([1., 2, 3, 4, 5]); obj = np.array([-5., -4, -3, -2, -1]); feas = np.array([.8, .88, .93, .97, .99])
    po = O.path_oracle(s, obj, feas, target=0.90)
    assert po["any_feasible"] and po["best_obj_idx"] == 2 and po["min_feasible_s"] == 3.0
    assert not O.path_oracle(s, obj, np.full(5, 0.5), 0.90)["any_feasible"]


# ------------------------------ metrics ---------------------------- #
@test
def test_metrics_coverage_and_paired():
    p, lo, hi = M.coverage_ci([1]*95 + [0]*5, level=0.95)
    assert abs(p - 0.95) < 1e-9 and 0.0 <= lo <= p <= hi <= 1.0 and lo > 0.88
    m, se = M.mean_se([1.0, 2.0, 3.0]); assert abs(m - 2.0) < 1e-9 and abs(se - (1.0/np.sqrt(3))) < 1e-9
    pd = M.paired_diff([3., 4, 5], [1., 1, 1]); assert abs(pd["mean_diff"] - 3.0) < 1e-9 and pd["lo"] < 3.0 < pd["hi"]
    s = M.summarize_method("UG", feasible_flags=[1]*96 + [0]*4, objectives=[-4.0]*100, target=0.95)
    assert s.meets_target and abs(s.coverage - 0.96) < 1e-9


# --------------------------- validators ---------------------------- #
@test
def test_select_index_branches():
    # some qualify -> cheapest qualified; folds K, threshold 1-beta on avg feasibility
    alpha, beta, K = 0.10, 0.05, 10
    # column feasibility rates across folds; col2 all feasible (qualifies), col0 cheaper but not qualifying
    V_mat = np.zeros((K, 3))
    V_mat[:, 0] = 0.5           # never >= 1-alpha=0.9 -> avg_ind 0
    V_mat[:, 1] = 1.0           # always feasible -> avg_ind 1 >= 1-beta
    V_mat[:, 2] = 1.0
    C = np.array([[-9.0, -3.0, -5.0]])   # col0 cheapest but unqualified; among qualified {1,2}, col2 cheaper
    assert V._select_index(V_mat, C, alpha, beta, K) == 2
    # none qualify -> argmax avg_ind (most feasible)
    V_none = np.zeros((K, 3)); V_none[:, 0] = 0.5; V_none[:, 1] = 0.8; V_none[:, 2] = 1.0
    V_none[:] = np.where(V_none >= 0.9, V_none, 0.5)  # keep < threshold counts low
    V_bad = np.full((K, 3), 0.5); V_bad[:3, 1] = 1.0; V_bad[:7, 2] = 1.0  # col2 most often feasible
    assert V._select_index(V_bad, C, alpha, beta, K) == 2

@test
def test_validator_objective_ordering():
    # Paper claim: c'X_UG <= c'X_NGS and <= c'X_UNGS; and c'X_NV <= c'X_UG (nested margins).
    # A guaranteed-feasible ultra-safe candidate (x=0) prevents the argmax fallback.
    d, p, n2, b = 4, 8, 400, 3.0
    c = -np.ones(d)
    worse = 0
    for seed in range(25):
        r = np.random.default_rng(seed)
        xx = r.random((d, p)); xx[:, 0] = 0.0            # column 0 always feasible (0 <= b)
        def mock(para, cc, bb, data, al, dd, n, _xx=xx):
            pa = np.atleast_1d(para); return _xx if len(pa) > 1 else _xx[:, 0]
        ksi = np.abs(r.normal(size=(2 * n2, d)))
        delta = np.arange(1, p + 1)
        np.random.seed(seed)  # gaussian_supremum draws multivariate_normal from global RNG
        res = V.gaussian_supremum(mock, delta, c, b, 0.10, 2 * n2, d, ksi, 0.05, n2, n2)
        oNGS, oUNGS, oNV, oUG = (c @ res[:, k] for k in range(4))
        tol = 1e-9
        assert oNV <= oUG + tol, f"NV should be <= UG (seed {seed}): {oNV} vs {oUG}"
        assert oUG <= oNGS + tol, f"UG should be <= NGS (seed {seed}): {oUG} vs {oNGS}"
        assert oUG <= oUNGS + tol, f"UG should be <= UNGS (seed {seed}): {oUG} vs {oUNGS}"

@test
def test_gaussian_supremum_deterministic():
    d, p, n2, b = 3, 5, 200, 3.0
    r = np.random.default_rng(0); xx = r.random((d, p))
    def mock(para, cc, bb, data, al, dd, n, _xx=xx):
        pa = np.atleast_1d(para); return _xx if len(pa) > 1 else _xx[:, 0]
    ksi = np.abs(r.normal(size=(2 * n2, d))); delta = np.arange(1, p + 1); c = -np.ones(d)
    np.random.seed(123); a = V.gaussian_supremum(mock, delta, c, b, 0.10, 2 * n2, d, ksi, 0.05, n2, n2)
    np.random.seed(123); b2 = V.gaussian_supremum(mock, delta, c, b, 0.10, 2 * n2, d, ksi, 0.05, n2, n2)
    assert np.array_equal(a, b2), "validator must be reproducible under a fixed seed"


def main():
    failed = 0
    for fn in TESTS:
        try:
            fn(); print(f"PASS {fn.__name__}")
        except Exception as e:
            failed += 1; print(f"FAIL {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(TESTS) - failed}/{len(TESTS)} passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()

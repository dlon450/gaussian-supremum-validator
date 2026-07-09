"""End-to-end smoke test of the non-solver infrastructure (no Gurobi needed).

Composes the whole pipeline with a **mock solver** so it runs anywhere with just
numpy/scipy:

    config -> per-replication RNG stream -> DGP sampling -> solution path (mock)
    -> validators (UG/NGS/UNGS/NV) -> true feasibility oracle -> metrics summary.

This is the plan's `smoke` entry point; the real runner swaps `mock_solve_path`
for the actual `gsv.formulations` solvers (which require Gurobi/cvxpy). Run:

    /tmp/gsv_venv/bin/python runners/smoke.py

Exits non-zero on any sanity-check failure.
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gsv import config, rng, dgp as D, oracle as O, metrics as M, validators as V


def mock_solve_path(para, c, b, data, alpha, d, n):
    """A stand-in solution path with a realistic conservativeness ordering.

    Returns candidates x*(s_j) that go from aggressive (x = 1, least feasible) to
    fully conservative (x = 0, always feasible) as j increases. Same signature as
    the real CCP_* solvers, so validators can call it unchanged.
    """
    pa = np.atleast_1d(para)
    p = len(pa)
    scales = np.linspace(1.0, 0.0, p)              # decreasing conservativeness parameter -> smaller x
    path = np.clip(np.outer(np.ones(d), scales), 0.0, 1.0)   # (d, p)
    return path if p > 1 else path[:, 0]


def run_smoke(cfg_name="paper_ro_ellipsoid", n=200, n_reps=40, p=25):
    cfg = config.get_config(cfg_name)
    d, alpha, beta = cfg.d, cfg.alpha, cfg.beta
    b = cfg.b_factor * d
    c = cfg.c_value * np.ones(d)
    n1 = n // 2
    n2 = n - n1
    delta = np.arange(1, p + 1)                     # mock parameter grid
    target = 1.0 - beta

    rows = {k: {"feas": [], "obj": [], "gap": [], "rel": []} for k in ["UG", "NGS", "UNGS", "NV"]}
    failures = 0
    for rep in range(n_reps):
        r = rng.make_stream(cfg.base_seed, cfg_name, rep)
        ksi = D.sample(cfg.dgp, n, d, r)
        try:
            # NV, UG, NGS, UNGS selections (validators draw the GS quantile from global RNG;
            # seed it per-rep for reproducibility of this smoke run).
            np.random.seed(rng.stable_hash(f"{cfg_name}:{rep}", bits=31))
            res = V.gaussian_supremum(mock_solve_path, delta, c, b, alpha, n, d, ksi, beta, n1, n2)
        except Exception as e:  # never silently drop
            failures += 1
            continue
        sel = {"NGS": res[:, 0], "UNGS": res[:, 1], "NV": res[:, 2], "UG": res[:, 3]}

        # Ground-truth feasibility + oracle on the SAME mock path (independent of the split).
        full_path = mock_solve_path(delta, c, b, ksi, alpha, d, n1)   # (d, p)
        feas_path = O.true_feasibility(full_path, cfg.dgp, d, b)      # closed-form (Gaussian)
        obj_path = c @ full_path
        po = O.path_oracle(delta, obj_path, feas_path, target=1.0 - alpha)

        for k, x in sel.items():
            feas = O.true_feasibility(x, cfg.dgp, d, b)               # scalar
            rows[k]["feas"].append(1.0 if feas >= 1.0 - alpha else 0.0)
            obj = float(c @ x)
            rows[k]["obj"].append(obj)
            if po["any_feasible"]:
                rows[k]["gap"].append(obj - po["best_obj"])
                denom = abs(po["best_obj"]) if abs(po["best_obj"]) > 1e-9 else 1.0
                rows[k]["rel"].append((obj - po["best_obj"]) / denom)

    print(f"smoke: config={cfg_name}  n={n}  reps={n_reps}  target(1-beta)={target:.2f}  "
          f"tolerance(1-alpha)={1-alpha:.2f}  failures={failures}")
    print(f"{'method':>6} {'coverage':>9} {'95% CI':>16} {'meets':>6} {'mean_obj':>9} {'oracle_gap':>10}")
    summaries = {}
    for k in ["NV", "UG", "NGS", "UNGS"]:
        s = M.summarize_method(k, feasible_flags=rows[k]["feas"], objectives=rows[k]["obj"],
                               target=target, oracle_gaps=rows[k]["gap"], rel_oracle_gaps=rows[k]["rel"])
        summaries[k] = s
        print(f"{k:>6} {s.coverage:>9.3f} [{s.coverage_lo:>6.3f},{s.coverage_hi:>6.3f}] "
              f"{str(s.meets_target):>6} {s.mean_obj:>9.3f} {s.mean_oracle_gap:>10.3f}")

    # sanity checks
    ok = True
    for k, s in summaries.items():
        if not (0.0 <= s.coverage <= 1.0 and np.isfinite(s.mean_obj)):
            print(f"  BAD summary for {k}"); ok = False
    # paper ordering: UG objective no worse than the supremum validators (nested margins)
    if not (summaries["NV"].mean_obj <= summaries["UG"].mean_obj + 1e-6
            <= summaries["NGS"].mean_obj + 1e-6):
        print("  WARN mean-objective ordering NV<=UG<=NGS not observed (can happen with fallback)")
    pd = M.paired_diff(rows["NGS"]["obj"], rows["UG"]["obj"])
    print(f"paired obj diff NGS-UG: mean={pd['mean_diff']:.3f}  95% CI [{pd['lo']:.3f},{pd['hi']:.3f}]  "
          f"(>=0 expected: UG less conservative)")
    return ok


if __name__ == "__main__":
    sys.exit(0 if run_smoke() else 1)

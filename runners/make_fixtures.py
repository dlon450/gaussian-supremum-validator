"""Generate / check regression fixtures for the validator statistics.

These fixtures pin the *validator* behavior (the numerically delicate part:
margins + selection) on fixed deterministic inputs, WITHOUT a solver — the
solution path and validation data are supplied directly. This is the regression
gate that guards against accidental changes to the focal Univariate Gaussian
validator (and NGS/UNGS/NV) during refactoring.

    /tmp/gsv_venv/bin/python runners/make_fixtures.py --generate   # write fixtures
    /tmp/gsv_venv/bin/python runners/make_fixtures.py --check      # compare (CI gate)

Solver-level fixtures (deterministic x*(s_j) for a fixed dataset per formulation)
require Gurobi/cvxpy and are generated in the solver environment by the analogous
`--solver` path (documented; not run here).
"""
import sys, os, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gsv import validators as V

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "fixtures", "validator_selections.json")


def _cases():
    """Metadata for deterministic cases spanning easy/hard feasibility regimes.

    Inputs (path xx, validation data ksi) are reconstructed from the seed in
    :func:`_run_case`, so generation and checking are guaranteed symmetric.
    """
    return [{"seed": s, "n2": n2, "p": p, "d": d}
            for s, n2, p, d in [(0, 300, 8, 4), (1, 120, 6, 5), (2, 500, 12, 3)]]


def _run_case(case):
    d, p, n2, seed = case["d"], case["p"], case["n2"], case["seed"]
    r = np.random.default_rng(seed)
    xx = r.random((d, p)); xx[:, 0] = 0.0                # always-feasible anchor (avoids fallback)
    ksi = np.abs(r.normal(size=(2 * n2, d)))
    c = -np.ones(d); b = 3.0; alpha, beta = 0.10, 0.05

    def solve(para, cc, bb, data, al, dd, n, _xx=xx):
        pa = np.atleast_1d(para); return _xx if len(pa) > 1 else _xx[:, 0]

    np.random.seed(1000 + seed)                          # GS quantile draw (global RNG) fixed
    res = V.gaussian_supremum(solve, np.arange(1, p + 1), c, b, alpha, 2 * n2, d, ksi, beta, n2, n2)
    # record objective per method (order NGS,UNGS,NV,UG) — robust to tiny vertex ties
    return [round(float(c @ res[:, k]), 10) for k in range(4)]


def generate():
    os.makedirs(os.path.dirname(FIXTURE), exist_ok=True)
    data = {"columns": ["NGS", "UNGS", "NV", "UG"],
            "cases": [{**{k: c[k] for k in ("seed", "n2", "p", "d")}, "expected_obj": _run_case(c)}
                      for c in _cases()]}
    with open(os.path.abspath(FIXTURE), "w") as f:
        json.dump(data, f, indent=2)
    print(f"wrote {len(data['cases'])} fixtures -> {os.path.relpath(FIXTURE)}")


def check():
    with open(os.path.abspath(FIXTURE)) as f:
        data = json.load(f)
    ok = True
    for case in data["cases"]:
        got = _run_case(case)
        exp = case["expected_obj"]
        match = all(abs(a - b) <= 1e-8 for a, b in zip(got, exp))
        print(f"case seed={case['seed']}: {'OK' if match else 'MISMATCH'}  got={got} exp={exp}")
        ok = ok and match
    print("REGRESSION OK" if ok else "REGRESSION FAILED")
    return ok


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "--check"
    if mode == "--generate":
        generate()
    else:
        sys.exit(0 if check() else 1)

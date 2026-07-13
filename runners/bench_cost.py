"""Wall-clock cost of each validation scheme at a common data budget.

Quantifies the reviewer response's computational-cost argument: the two-phase
Gaussian-supremum/univariate validators solve the path ONCE on Phase-1 data and
then only evaluate cheap sample statistics on Phase-2, whereas K-fold CV and
B-bootstrap resolve the optimization ~K (resp. B) times over the parameter grid.

    /tmp/gsv_venv/bin/python runners/bench_cost.py [--reps R]

Reports, per (formulation, n): seconds for two-phase (path+select), CV, BS,
sectioning, and the ratios CV/two-phase etc. Run when the machine is otherwise
idle for clean numbers. Single BLAS/solver thread.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, time, json, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gsv import config as C, rng as RNG, dgp as D, paths as P, select as S, validators as V
from gsv import experiment as E

# Sizes chosen so the (deliberately expensive) CV/BS calls still terminate for a
# clean ratio; the CV/two_phase ratio is the reported quantity, not absolute n.
CASES = [("paper_so", 500), ("paper_saa", 200), ("paper_dro_wasserstein", 100)]


def time_call(fn, reps):
    ts = []
    for _ in range(reps):
        t0 = time.time(); fn(); ts.append(time.time() - t0)
    return float(np.mean(ts)), float(np.std(ts))


def bench(cfg_name, n, reps):
    cfg = C.get_config(cfg_name)
    alpha, beta, d = cfg.alpha, cfg.beta, cfg.d
    b = cfg.b_factor * d; c = cfg.c_value * np.ones(d)
    split = 0.7; n1 = max(2, round(split * n)); n2 = n - n1
    r = RNG.make_stream(cfg.base_seed, f"bench:{cfg_name}:{n}", 0)
    data = D.sample(cfg.dgp, n, d, r)
    perm = r.permutation(n); phase1, phase2 = data[perm[:n1]], data[perm[n1:]]

    def two_phase():
        xx, s_values = P.build_path(cfg.formulation, phase1, c, b, alpha, d, n1, cfg.mesh, rng=r)
        objectives = c @ xx
        P_hat, sigma_hat, Sigma_hat = S.phase2_stats(xx, phase2, b)
        S.select_all(P_hat, sigma_hat, Sigma_hat, objectives, alpha=alpha, beta=beta, n2=n2, rng=r)

    cv_grid, bs_grid, solve = E._existing_grid(cfg.formulation, alpha)
    K = 10
    delta = P.WASSERSTEIN_GRID if cfg.formulation == "dro_wasserstein" else np.arange(1, n + 1)

    def cv():
        np.random.seed(0)
        V.cross_validation(solve, delta, c, b, alpha, n, d, data, beta, K, param_grid=cv_grid)

    def bs():
        np.random.seed(0)
        V.bootstrapping(solve, delta, c, b, alpha, n, d, data, beta, K, param_grid=bs_grid)

    def sec():
        np.random.seed(0)
        V.sectioning(solve, (cv_grid(n1) if cv_grid else delta), c, b, alpha, n, d, data, beta, n1, n2, K)

    out = {"config": cfg_name, "formulation": cfg.formulation, "n": n, "K": K}
    for name, fn in [("two_phase", two_phase), ("CV", cv), ("BS", bs), ("Sectioning", sec)]:
        m, sd = time_call(fn, reps)
        out[name + "_s"] = round(m, 3); out[name + "_sd"] = round(sd, 3)
    tp = out["two_phase_s"] or 1e-9
    out["CV/two_phase"] = round(out["CV_s"] / tp, 1)
    out["BS/two_phase"] = round(out["BS_s"] / tp, 1)
    out["Sectioning/two_phase"] = round(out["Sectioning_s"] / tp, 1)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--reps", type=int, default=3); a = ap.parse_args()
    results = []
    for cfg_name, n in CASES:
        o = bench(cfg_name, n, a.reps)
        results.append(o)
        print(f"{o['config']:22s} n={n:4d} K={o['K']}  two_phase={o['two_phase_s']:7.3f}s  "
              f"CV={o['CV_s']:8.3f}s ({o['CV/two_phase']}x)  BS={o['BS_s']:8.3f}s ({o['BS/two_phase']}x)  "
              f"Sec={o['Sectioning_s']:8.3f}s ({o['Sectioning/two_phase']}x)", flush=True)
    outdir = os.path.join(os.path.dirname(__file__), "..", "results", "analysis")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "compute_cost.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("wrote results/analysis/compute_cost.json")


if __name__ == "__main__":
    main()

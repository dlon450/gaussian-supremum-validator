"""Two-phase validation on the CVaR-constrained portfolio (the paper's second
problem class / a general stochastic constraint).

    /tmp/gsv_venv/bin/python runners/portfolio_experiment.py <matrix> [--reps R] [--workers W]

Matrices:
  ndgrid  (N x D) grid, coverage & return of NV/UG/NGS/UNGS + SAA benchmark (s=0)
  budget  Phase-1 fraction sweep
  gamma   risk-threshold sweep
  tails   heavy-tailed (multivariate-t) returns

Feasibility = true CVaR_alpha(-xi'x) <= gamma (closed form under Gaussian, large
sample under t). Objective = c'x = -expected return (lower is better = higher
return). Each replication is deterministic in (config, rep) via gsv.rng streams.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, argparse, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gsv import rng as RNG, portfolio as PF, metrics as M

ALPHA = 0.10          # CVaR tail level
BETA = 0.05           # target confidence 1-beta = 0.95
GAMMA = 0.15          # CVaR risk threshold (calibrated so s=0 overfits / is infeasible)
BASE_SEED = 2
RADII = PF.WASSERSTEIN_RADII
METHODS = ["NV", "UG", "NGS", "UNGS", "benchmark"]


def moments(d):
    """Return DGP moments: increasing mean & vol with asset index, moderate corr."""
    idx = np.arange(1, d + 1) / d
    mu = 0.08 + 0.06 * idx
    sig = 0.15 + 0.15 * idx
    Sigma = np.outer(sig, sig) * (0.2 + 0.8 * np.eye(d))
    return mu, Sigma


def sample(kind, mu, Sigma, n, d, r, df=6):
    if kind == "gaussian":
        return r.multivariate_normal(mu, Sigma, n)
    scale = Sigma * (df - 2.0) / df                          # matched covariance
    z = r.multivariate_normal(np.zeros(d), scale, n)
    w = r.chisquare(df, n)
    return mu + z / np.sqrt(w / df)[:, None]


def run_replication(n, split, d, gamma, dgp, rep):
    mu, Sigma = moments(d); c = -mu
    n1 = max(5, round(split * n)); n2 = n - n1
    r = RNG.make_stream(BASE_SEED, f"portfolio:{dgp}:n{n}:s{split}:d{d}:g{gamma}", rep)
    data = sample(dgp, mu, Sigma, n, d, r)
    eval_sample = None if dgp == "gaussian" else sample(dgp, mu, Sigma, 200_000, d, r)
    perm = r.permutation(n); p1, p2 = data[perm[:n1]], data[perm[n1:]]

    xx = PF.dro_cvar_path(RADII, c, gamma, p1, ALPHA, d, n1)      # (d, p)
    objectives = c @ xx
    cvar_hat, se, IF = PF.cvar_phase2_stats(xx, p2, ALPHA)
    picks = PF.select_cvar(cvar_hat, se, IF, objectives, gamma=gamma, beta=BETA, n2=n2, rng=r)

    def true_cvar(x):
        if dgp == "gaussian":
            return PF.true_cvar_gaussian(x, mu, Sigma, ALPHA)
        losses = -(eval_sample @ x); var = np.quantile(losses, 1 - ALPHA)
        return float(var + np.maximum(losses - var, 0).mean() / ALPHA)

    recs = []
    for m in ["NV", "UG", "NGS", "UNGS"]:
        idx = picks[m][0]; x = xx[:, idx]
        tc = true_cvar(x)
        recs.append({"method": m, "n": n, "n1": n1, "n2": n2, "d": d, "split": split, "gamma": gamma,
                     "dgp": dgp, "rep": rep, "feasible": 1.0 if tc <= gamma else 0.0,
                     "true_cvar": tc, "objective": float(c @ x), "return": float(mu @ x)})
    # benchmark = SAA (radius 0, aggressive) — what you'd do without validation
    xb = xx[:, 0]; tcb = true_cvar(xb)
    recs.append({"method": "benchmark", "n": n, "d": d, "split": split, "gamma": gamma, "dgp": dgp,
                 "rep": rep, "feasible": 1.0 if tcb <= gamma else 0.0, "true_cvar": tcb,
                 "objective": float(c @ xb), "return": float(mu @ xb)})
    return recs


def run_cell(n, split, d, gamma, dgp, reps, workers):
    rows = []
    if workers > 1:
        import concurrent.futures as cf
        with cf.ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(run_replication, n, split, d, gamma, dgp, rep) for rep in range(reps)]
            for f in cf.as_completed(futs):
                try: rows.extend(f.result())
                except Exception as e: rows.append({"method": "__error__", "error": str(e)[:200]})
    else:
        for rep in range(reps):
            rows.extend(run_replication(n, split, d, gamma, dgp, rep))
    summ = {}
    for m in METHODS:
        mr = [r for r in rows if r.get("method") == m]
        if not mr: continue
        s = M.summarize_method(m, feasible_flags=[r["feasible"] for r in mr],
                               objectives=[r["objective"] for r in mr], target=1 - BETA)
        d_ = s.to_dict(); d_["mean_return"] = float(np.mean([r["return"] for r in mr]))
        d_["mean_true_cvar"] = float(np.mean([r["true_cvar"] for r in mr]))
        summ[m] = d_
    summ["_n_replication_errors"] = sum(1 for r in rows if r.get("method") == "__error__")
    return rows, summ


def matrix(name):
    cells = []
    if name == "ndgrid":
        for d in (2, 5, 10, 20):
            for n in (100, 200, 500, 1000):
                cells.append((f"port_n{n}_d{d}", dict(n=n, split=0.5, d=d, gamma=GAMMA, dgp="gaussian")))
    elif name == "budget":
        for s in [round(i/10, 1) for i in range(1, 10)]:
            cells.append((f"port_n500_s{s}_d10", dict(n=500, split=s, d=10, gamma=GAMMA, dgp="gaussian")))
    elif name == "gamma":
        for g in (0.12, 0.15, 0.18, 0.22, 0.28):
            cells.append((f"port_n500_d10_g{g}", dict(n=500, split=0.5, d=10, gamma=g, dgp="gaussian")))
    elif name == "tails":
        for df in (3, 4, 6, 10):
            cells.append((f"port_n500_d10_t{df}", dict(n=500, split=0.5, d=10, gamma=GAMMA, dgp="multivariate_t")))
    else:
        raise SystemExit(f"unknown matrix {name!r}")
    return cells


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("matrix"); ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--workers", type=int, default=48); a = ap.parse_args()
    outdir = os.path.join(os.path.dirname(__file__), "..", "results", "experiments", f"portfolio_{a.matrix}")
    os.makedirs(outdir, exist_ok=True)
    cells = matrix(a.matrix)
    print(f"portfolio matrix={a.matrix} cells={len(cells)} reps={a.reps} workers={a.workers}", flush=True)
    for i, (key, kw) in enumerate(cells):
        sp = os.path.join(outdir, f"{key}_summary.json")
        if os.path.exists(sp): print(f"[skip] {key}", flush=True); continue
        t0 = time.time()
        rows, summ = run_cell(reps=a.reps, workers=a.workers, **kw)
        dt = time.time() - t0
        with open(os.path.join(outdir, f"{key}_raw.jsonl"), "w") as f:
            for r in rows: f.write(json.dumps(r) + "\n")
        with open(sp, "w") as f:
            json.dump({"key": key, **kw, "reps": a.reps, "elapsed_s": dt, "summaries": summ}, f, indent=2)
        cov = {m: round(summ[m]["coverage"], 3) for m in METHODS if m in summ}
        ret = {m: round(summ[m]["mean_return"], 3) for m in METHODS if m in summ}
        print(f"[done] ({i+1}/{len(cells)}) {key} {dt:.1f}s err={summ['_n_replication_errors']} cov={cov} ret={ret}", flush=True)
    print(f"portfolio {a.matrix} complete.", flush=True)


if __name__ == "__main__":
    main()

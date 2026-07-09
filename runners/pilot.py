"""Pilot experiment runner with REAL solvers (Gurobi) + oracle + metrics.

Runs the revised experiment for one formulation at pilot scale: samples the
(Gaussian) DGP, builds the Phase-1 solution path, validates on Phase-2 with the
focal Univariate Gaussian validator + NGS/UNGS/NV and the formulation benchmark,
then scores each selected solution against the closed-form Gaussian oracle
(feasibility coverage, objective, path-oracle gap, excess conservativeness).

Usage:
    /tmp/gsv_venv/bin/python runners/pilot.py <config> [--n N] [--reps R] [--split F] [--mesh-p P]

Example:
    /tmp/gsv_venv/bin/python runners/pilot.py paper_ro_ellipsoid --n 200 --reps 50 --split 0.5
"""
import sys, os, json, argparse, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dataclasses import replace
from gsv import config, rng as RNG, dgp as D, oracle as O, metrics as M, paths as P, select as S

METHODS = ["benchmark", "NV", "UG", "NGS", "UNGS"]


def run(cfg_name, n, reps, split, mesh_p=None):
    cfg = config.get_config(cfg_name)
    if mesh_p is not None:
        cfg = replace(cfg, mesh=replace(cfg.mesh, p=mesh_p))
    d, alpha, beta = cfg.d, cfg.alpha, cfg.beta
    b = cfg.b_factor * d
    c = cfg.c_value * np.ones(d)
    target = 1.0 - beta
    mu_true, sigma_true = D.moments(cfg.dgp, d)
    n1 = max(2, round(split * n)); n2 = n - n1
    stream_key = f"{cfg_name}:n{n}:s{split}"

    acc = {m: {"feas": [], "obj": [], "gap": [], "rel": [], "exc": []} for m in METHODS}
    failures = 0
    t0 = time.time()
    for rep in range(reps):
        r = RNG.make_stream(cfg.base_seed, stream_key, rep)
        data = D.sample(cfg.dgp, n, d, r)
        perm = r.permutation(n)
        phase1, phase2 = data[perm[:n1]], data[perm[n1:]]
        try:
            xx, s_values = P.build_path(cfg.formulation, phase1, c, b, alpha, d, n1, cfg.mesh, rng=r)
            objectives = c @ xx
            P_hat, sigma_hat, Sigma_hat = S.phase2_stats(xx, phase2, b)
            picks = S.select_all(P_hat, sigma_hat, Sigma_hat, objectives,
                                 alpha=alpha, beta=beta, n2=n2, rng=r)
            bench = P.benchmark_solution(cfg.benchmark, c=c, b=b, alpha=alpha, d=d,
                                         mu_true=mu_true, sigma_true=sigma_true,
                                         full_data=data, n1=n1) if cfg.benchmark else None
        except Exception as e:
            failures += 1
            print(f"  rep {rep}: FAILURE {type(e).__name__}: {str(e)[:80]}")
            continue

        # closed-form oracle feasibility of the whole path -> path oracle
        feas_path = O.gaussian_feasibility(xx, mu_true, sigma_true, b)
        po = O.path_oracle(s_values, objectives, feas_path, target=1 - alpha)

        for m in METHODS:
            if m == "benchmark":
                if bench is None or not np.all(np.isfinite(bench)):
                    continue
                x = np.asarray(bench, float); s_sel = np.nan
            else:
                idx = picks[m][0]; x = xx[:, idx]; s_sel = float(s_values[idx])
            feas = float(O.gaussian_feasibility(x, mu_true, sigma_true, b))
            obj = float(c @ x)
            acc[m]["feas"].append(1.0 if feas >= 1 - alpha else 0.0)
            acc[m]["obj"].append(obj)
            if po["any_feasible"]:
                acc[m]["gap"].append(obj - po["best_obj"])
                denom = abs(po["best_obj"]) or 1.0
                acc[m]["rel"].append((obj - po["best_obj"]) / denom)
                if np.isfinite(s_sel):
                    acc[m]["exc"].append(s_sel - po["min_feasible_s"])

    elapsed = time.time() - t0
    print(f"\n=== {cfg_name}  formulation={cfg.formulation}  n={n} (n1={n1},n2={n2})  reps={reps}  "
          f"benchmark={cfg.benchmark}  1-alpha={1-alpha:.2f}  target 1-beta={target:.2f} ===")
    print(f"    elapsed {elapsed:.1f}s  failures={failures}")
    hdr = f"{'method':>10} {'coverage':>9} {'95% CI':>16} {'meets':>6} {'mean_obj':>10} {'obj_SE':>7} {'oracle_gap':>10} {'excess_s':>9}"
    print(hdr)
    summaries = {}
    for m in METHODS:
        if not acc[m]["feas"]:
            continue
        s = M.summarize_method(m, feasible_flags=acc[m]["feas"], objectives=acc[m]["obj"],
                               target=target, oracle_gaps=acc[m]["gap"], rel_oracle_gaps=acc[m]["rel"],
                               excess_s=acc[m]["exc"] or None)
        summaries[m] = s
        print(f"{m:>10} {s.coverage:>9.3f} [{s.coverage_lo:>6.3f},{s.coverage_hi:>6.3f}] "
              f"{str(s.meets_target):>6} {s.mean_obj:>10.3f} {s.obj_se:>7.3f} {s.mean_oracle_gap:>10.3f} {s.mean_excess_s:>9.3f}")

    if "UG" in summaries:
        for other in ["benchmark", "NGS", "UNGS"]:
            if other in summaries and acc[other]["obj"] and len(acc[other]["obj"]) == len(acc["UG"]["obj"]):
                pd = M.paired_diff(acc[other]["obj"], acc["UG"]["obj"])
                print(f"  paired obj {other}-UG: mean={pd['mean_diff']:.3f} 95%CI [{pd['lo']:.3f},{pd['hi']:.3f}]")

    # persist raw + summary
    outdir = os.path.join(os.path.dirname(__file__), "..", "results", "pilot")
    os.makedirs(outdir, exist_ok=True)
    stamp = f"{cfg_name}_n{n}_s{split}_r{reps}"
    with open(os.path.join(outdir, f"{stamp}_summary.json"), "w") as f:
        json.dump({"config": cfg_name, "formulation": cfg.formulation, "n": n, "n1": n1, "n2": n2,
                   "reps": reps, "alpha": alpha, "beta": beta, "elapsed_s": elapsed, "failures": failures,
                   "summaries": {m: s.to_dict() for m, s in summaries.items()}}, f, indent=2)
    print(f"  saved -> results/pilot/{stamp}_summary.json")
    return summaries


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--split", type=float, default=0.5)
    ap.add_argument("--mesh-p", type=int, default=None)
    a = ap.parse_args()
    run(a.config, a.n, a.reps, a.split, a.mesh_p)

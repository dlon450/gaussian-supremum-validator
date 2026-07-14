"""Saturation campaign: sweep every remaining experimental axis.

    /tmp/gsv_venv/bin/python runners/run_saturation.py <matrix> [--reps R] [--workers W]

Matrices:
  alpha      tolerance 1-alpha sweep (alpha in {0.01..0.30}), per formulation
  beta       target-confidence 1-beta sweep (calibration: empirical vs nominal)
  tails      heavy-tailed multivariate-t, df in {3,4,6,10,30}
  corr       correlated Gaussian, corr in {0.0..0.8} (with a d=20 slice)
  foldsdeep  CV/BS/Sec at K in {2,3,5,10,20} (SO) + larger N for SAA/Wasserstein
  budgetfull Wasserstein d=20 n1-fraction sweep (deferred) + high-D SO/RO n1 sweep
  bigd       extreme dimension d in {100,200} for RO/SO

Each cell -> results/experiments/<matrix>/<key>_{summary.json,raw.jsonl}; done cells
skipped on resume. Single BLAS/solver thread per worker.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys, json, argparse, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gsv import experiment as E
from gsv.util import dump_json, dump_jsonl

SO_MESH_P = 60
SPLITS = [round(i / 10, 1) for i in range(1, 10)]
# per-formulation fixed operating point for the 1-D sweeps (alpha/beta/tails/corr)
FIXED = {"paper_so": (500, 10), "paper_ro_ellipsoid": (500, 10),
         "paper_saa": (400, 10), "paper_dro_wasserstein": (200, 10)}
FORMS = list(FIXED)


def cellkw(cfg, n, split, d, **kw):
    base = dict(include_existing=False, folds=None, mesh_p=(SO_MESH_P if cfg == "paper_so" else None),
                alpha=None, beta=None, dgp_params=None)
    base.update(kw)
    return dict(cfg=cfg, n=n, split=split, d=d, **base)


def matrix(name):
    cells = []
    if name == "alpha":
        for cfg in FORMS:
            n, d = FIXED[cfg]
            for a in (0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30):
                cells.append((f"{cfg}_n{n}_d{d}_a{a}", cellkw(cfg, n, 0.5, d, alpha=a)))
    elif name == "beta":
        for cfg in ("paper_so", "paper_ro_ellipsoid", "paper_saa", "paper_dro_wasserstein"):
            n, d = FIXED[cfg]
            for be in (0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50):
                cells.append((f"{cfg}_n{n}_d{d}_be{be}", cellkw(cfg, n, 0.5, d, beta=be)))
    elif name == "tails":
        for cfg in FORMS:
            n, d = FIXED[cfg]
            for df in (3, 4, 6, 10, 30):
                cells.append((f"{cfg}_n{n}_d{d}_t{df}",
                              cellkw(cfg, n, 0.5, d, dgp_params={"kind": "multivariate_t", "df": df})))
    elif name == "corr":
        for cfg in FORMS:
            n, d = FIXED[cfg]
            for co in (0.0, 0.2, 0.4, 0.6, 0.8):
                cells.append((f"{cfg}_n{n}_d{d}_c{co}", cellkw(cfg, n, 0.5, d, dgp_params={"corr": co})))
        for cfg in ("paper_so", "paper_ro_ellipsoid"):     # a higher-dim slice
            for co in (0.0, 0.4, 0.8):
                cells.append((f"{cfg}_n500_d20_c{co}", cellkw(cfg, 500, 0.5, 20, dgp_params={"corr": co})))
    elif name == "foldsdeep":
        for n in (200, 300, 500, 1000):
            cells.append((f"paper_so_n{n}_deepfolds",
                          cellkw("paper_so", n, 0.7, 10, include_existing=True, folds=(2, 3, 5, 10, 20))))
        for n in (100, 200, 300):
            cells.append((f"paper_saa_n{n}_folds", cellkw("paper_saa", n, 0.7, 10, include_existing=True, folds=(3, 5, 10))))
        for n in (100, 200):
            cells.append((f"paper_dro_wasserstein_n{n}_folds",
                          cellkw("paper_dro_wasserstein", n, 0.7, 10, include_existing=True, folds=(3, 5, 10))))
    elif name == "budgetfull":
        for d in (20,):                                    # deferred Wasserstein high-d n1 sweep
            for n in (100, 200):
                for s in SPLITS:
                    cells.append((f"paper_dro_wasserstein_n{n}_s{s}_d{d}", cellkw("paper_dro_wasserstein", n, s, d)))
        for cfg in ("paper_so", "paper_ro_ellipsoid"):     # high-D best-n1
            for n in (500, 1000):
                for s in SPLITS:
                    cells.append((f"{cfg}_n{n}_s{s}_d100", cellkw(cfg, n, s, 100)))
    elif name == "bigd":
        for d in (100, 200):
            for n in (200, 500, 1000, 2000):
                cells.append((f"paper_ro_ellipsoid_n{n}_d{d}", cellkw("paper_ro_ellipsoid", n, 0.5, d)))
            for n in (500, 1000, 2000):
                cells.append((f"paper_so_n{n}_d{d}", cellkw("paper_so", n, 0.5, d)))
    elif name == "highd":
        # high-dimension coverage for the MILP + SDP formulations, to match the paper's
        # d=50 RO/FAST tables and push moment-DRO beyond its d=10 table.
        for n in (200, 500):
            cells.append((f"paper_saa_n{n}_d50", cellkw("paper_saa", n, 0.5, 50)))
        for n in (100, 200):
            cells.append((f"paper_dro_wasserstein_n{n}_d50", cellkw("paper_dro_wasserstein", n, 0.5, 50)))
        for d in (20, 30):
            for n in (200, 500):
                cells.append((f"paper_dro_moment_n{n}_d{d}", cellkw("paper_dro_moment", n, 0.5, d)))
    elif name == "moment":
        # moment-based DRO (SDP) — the one paper formulation not yet in the pipeline.
        # Paper reports only d=10; we sweep N x D plus n1 / tolerance / tails slices.
        for d in (2, 5, 10, 15):
            for n in (100, 200, 500):
                cells.append((f"paper_dro_moment_n{n}_d{d}", cellkw("paper_dro_moment", n, 0.5, d)))
        for s in SPLITS:                                  # best-n1 for moment-DRO
            cells.append((f"paper_dro_moment_n200_s{s}_d10", cellkw("paper_dro_moment", 200, s, 10)))
        for a in (0.05, 0.10, 0.20):                      # tolerance slice
            cells.append((f"paper_dro_moment_n200_d10_a{a}", cellkw("paper_dro_moment", 200, 0.5, 10, alpha=a)))
        for df in (3, 6, 30):                             # heavy-tail slice
            cells.append((f"paper_dro_moment_n200_d10_t{df}",
                          cellkw("paper_dro_moment", 200, 0.5, 10, dgp_params={"kind": "multivariate_t", "df": df})))
    elif name == "extremes":
        # push D and N to the edges to confirm the qualitative trends saturate.
        for cfg in ("paper_ro_ellipsoid", "paper_so"):
            for d in (200, 500):                        # extreme dimension (+ high N so SO can catch up)
                for n in (2000, 5000):
                    cells.append((f"{cfg}_n{n}_d{d}", cellkw(cfg, n, 0.5, d)))
            for n in (5000, 10000):                     # extreme sample size (convergence)
                cells.append((f"{cfg}_n{n}_d10", cellkw(cfg, n, 0.5, 10)))
            for s in (0.05, 0.1, 0.2, 0.3, 0.5, 0.7):   # best-n1 at large N: does it keep shrinking?
                cells.append((f"{cfg}_n5000_s{s}_d10", cellkw(cfg, 5000, s, 10)))
    else:
        raise SystemExit(f"unknown matrix {name!r}")
    return cells


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix")
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--workers", type=int, default=88)
    a = ap.parse_args()

    outdir = os.path.join(os.path.dirname(__file__), "..", "results", "experiments", a.matrix)
    os.makedirs(outdir, exist_ok=True)
    cells = matrix(a.matrix)
    print(f"matrix={a.matrix}  cells={len(cells)}  reps={a.reps}  workers={a.workers}", flush=True)

    for i, (key, kw) in enumerate(cells):
        summ_path = os.path.join(outdir, f"{key}_summary.json")
        if os.path.exists(summ_path):
            print(f"[skip] ({i+1}/{len(cells)}) {key}", flush=True); continue
        t0 = time.time()
        try:
            rows, summ = E.run_cell(kw["cfg"], kw["n"], kw["split"], kw["d"], a.reps,
                                    include_existing=kw["include_existing"], workers=a.workers,
                                    folds=kw["folds"], mesh_p=kw["mesh_p"],
                                    alpha=kw["alpha"], beta=kw["beta"], dgp_params=kw["dgp_params"])
        except Exception as e:
            print(f"[FAIL] ({i+1}/{len(cells)}) {key}  {type(e).__name__}: {str(e)[:120]}", flush=True)
            continue
        dt = time.time() - t0
        dump_jsonl(rows, os.path.join(outdir, f"{key}_raw.jsonl"))
        meta = {"key": key, "config": kw["cfg"], "n": kw["n"], "split": kw["split"], "d": kw["d"],
                "reps": a.reps, "alpha": kw["alpha"], "beta": kw["beta"], "dgp_params": kw["dgp_params"],
                "folds": kw["folds"], "mesh_p": kw["mesh_p"], "elapsed_s": dt, "summaries": summ}
        dump_json(meta, summ_path)
        cov = {m: round(summ[m]["coverage"], 3) for m in summ
               if isinstance(summ[m], dict) and "coverage" in summ[m]}
        print(f"[done] ({i+1}/{len(cells)}) {key}  {dt:6.1f}s  errors={summ.get('_n_replication_errors',0)}  cov={cov}", flush=True)
    print(f"matrix {a.matrix} complete.", flush=True)


if __name__ == "__main__":
    main()

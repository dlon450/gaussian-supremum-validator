"""Comprehensive experiment matrices for the revision (best n1; N x D scaling;
CV/BS at multiple fold counts) across all formulations.

    /tmp/gsv_venv/bin/python runners/run_comprehensive.py <matrix> [--reps R] [--workers W]

Matrices:
  budget  Phase-1 fraction n1/n in {0.1..0.9} across N and D, per formulation
          (two-phase validators only) -> "what is the best n1, and how does it move
          with N and D?"
  ndgrid  full (N, D) grid per formulation at split=0.5 (two-phase + benchmark)
          -> "how do feasibility/objective/gap change with N and D?"
  folds   CV / bootstrap / sectioning at K in {3,5,10} vs the proposed validators,
          at split=0.7, per formulation and N -> "how do CV/BS at different fold
          counts compare to our method?"

Per cell we save results/experiments/<matrix>/<key>_{summary.json,raw.jsonl};
done cells are skipped on resume. Single BLAS/solver thread per worker.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys, json, argparse, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gsv import experiment as E

SPLITS = [round(i / 10, 1) for i in range(1, 10)]     # 0.1 .. 0.9
FOLDS = (3, 5, 10)
SO_MESH_P = 60   # cap the SO two-phase path (n1 LP solves) so wide sweeps stay cheap


def cell(cfg, n, split, d, mesh_p=None, include_existing=False, folds=None):
    return dict(cfg=cfg, n=n, split=split, d=d, mesh_p=mesh_p,
                include_existing=include_existing, folds=folds)


def matrix(name):
    cells = []
    if name == "budget":
        # best n1 vs (N, D); two-phase only. Cheap convex formulations get a wider
        # (N,D) net than the MILP ones.
        for d in (2, 10, 50):
            for n in (200, 500, 1000):
                for s in SPLITS:
                    cells.append(cell("paper_so", n, s, d, mesh_p=SO_MESH_P))
                    cells.append(cell("paper_ro_ellipsoid", n, s, d))
        for d in (5, 10, 20):
            for n in (200, 400):
                for s in SPLITS:
                    cells.append(cell("paper_saa", n, s, d))
        # Wasserstein big-M MILP path is the most expensive; keep d in {5,10} for the n1 sweep.
        for d in (5, 10):
            for n in (100, 200):
                for s in SPLITS:
                    cells.append(cell("paper_dro_wasserstein", n, s, d))
    elif name == "ndgrid":
        for d in (2, 5, 10, 20, 50):
            for n in (100, 200, 500, 1000, 2000):
                cells.append(cell("paper_so", n, 0.5, d, mesh_p=SO_MESH_P))
        for d in (2, 5, 10, 20, 50, 100):
            for n in (100, 200, 500, 1000, 2000):
                cells.append(cell("paper_ro_ellipsoid", n, 0.5, d))
        for d in (2, 5, 10, 20):
            for n in (100, 200, 400):
                cells.append(cell("paper_saa", n, 0.5, d))
                cells.append(cell("paper_dro_wasserstein", n, 0.5, d))
        for d in (2, 5, 10, 20, 50):
            for n in (200, 500, 1000):
                cells.append(cell("paper_fast", n, 0.5, d))
    elif name == "folds":
        # CV/BS/Sec at K in {3,5,10} vs proposed. SO (LP) goes wide; MILP formulations
        # are limited to smaller N (CV on big-M MILPs is the intractable regime).
        for n in (100, 200, 300, 500):
            cells.append(cell("paper_so", n, 0.7, 10, mesh_p=SO_MESH_P, include_existing=True, folds=FOLDS))
        for n in (100, 200):
            cells.append(cell("paper_saa", n, 0.7, 10, include_existing=True, folds=FOLDS))
        for n in (100,):
            cells.append(cell("paper_dro_wasserstein", n, 0.7, 10, include_existing=True, folds=FOLDS))
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

    for i, cc in enumerate(cells):
        key = f"{cc['cfg']}_n{cc['n']}_s{cc['split']}_d{cc['d']}"
        summ_path = os.path.join(outdir, f"{key}_summary.json")
        if os.path.exists(summ_path):
            print(f"[skip] ({i+1}/{len(cells)}) {key}", flush=True); continue
        t0 = time.time()
        try:
            rows, summ = E.run_cell(cc["cfg"], cc["n"], cc["split"], cc["d"], a.reps,
                                    include_existing=cc["include_existing"], workers=a.workers,
                                    folds=cc["folds"], mesh_p=cc["mesh_p"])
        except Exception as e:
            print(f"[FAIL] ({i+1}/{len(cells)}) {key}  {type(e).__name__}: {str(e)[:120]}", flush=True)
            continue
        dt = time.time() - t0
        with open(os.path.join(outdir, f"{key}_raw.jsonl"), "w") as f:
            for r in rows:
                f.write(json.dumps(r) + "\n")
        with open(summ_path, "w") as f:
            json.dump({"config": cc["cfg"], "n": cc["n"], "split": cc["split"], "d": cc["d"], "reps": a.reps,
                       "mesh_p": cc["mesh_p"], "include_existing": cc["include_existing"], "folds": cc["folds"],
                       "elapsed_s": dt, "summaries": summ}, f, indent=2)
        cov = {m: round(summ[m]["coverage"], 3) for m in summ
               if isinstance(summ[m], dict) and "coverage" in summ[m]}
        errs = summ.get("_n_replication_errors", 0)
        print(f"[done] ({i+1}/{len(cells)}) {key}  {dt:6.1f}s  errors={errs}  coverage={cov}", flush=True)
    print(f"matrix {a.matrix} complete.", flush=True)


if __name__ == "__main__":
    main()

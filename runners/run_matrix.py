"""Run a named experiment matrix in parallel, with checkpointing.

    /tmp/gsv_venv/bin/python runners/run_matrix.py <matrix> [--reps R] [--workers W]

Matrices (each a list of cells (config, n, split, d, include_existing)):
  dim      RO dimension sweep d=2..100 (the dimension-free headline)
  robust   heavy-tailed multivariate-t (where GS can beat UG)
  nsweep   sample-size sweep for RO/SO/FAST (coverage->target, gap->0)
  split    n1/n2 split-ratio sweep (the ~70% budgeting finding)
  existing CV / bootstrap / sectioning comparison on SO/SAA/Wasserstein

Per cell we save results/experiments/<matrix>/<key>_{summary.json,raw.jsonl};
existing cells are skipped on resume. One solver + BLAS thread per worker.
"""
import os
# limit BLAS/OpenMP threads BEFORE numpy import; inherited by spawned workers
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys, json, argparse, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gsv import experiment as E

RO_DIMS = [2, 5, 10, 20, 50, 100]


def matrix(name, reps):
    cells = []
    if name == "dim":
        for d in RO_DIMS:
            cells.append(("paper_ro_ellipsoid", 500, 0.5, d, False))
    elif name == "robust":
        for n in [200, 500, 1000]:
            cells.append(("robust_ro_ellipsoid_t", n, 0.5, 10, False))
        for n in [100, 200]:
            cells.append(("robust_dro_wasserstein_t", n, 0.5, 10, False))
    elif name == "nsweep":
        for n in [100, 200, 300, 500, 1000]:
            cells.append(("paper_ro_ellipsoid", n, 0.5, 10, False))
        for n in [100, 200, 300, 500, 800, 1000]:
            cells.append(("paper_so", n, 0.5, 10, False))
        for n in [200, 500]:
            cells.append(("paper_fast", n, 0.5, 10, False))
    elif name == "split":
        for s in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            cells.append(("paper_so", 500, s, 10, False))
            cells.append(("paper_saa", 400, s, 10, False))
    elif name == "existing":
        for n in [100, 200, 300, 500]:
            cells.append(("paper_so", n, 0.7, 10, True))
            cells.append(("paper_saa", n, 0.7, 10, True))
        for n in [100, 200, 300]:
            cells.append(("paper_dro_wasserstein", n, 0.7, 10, True))
    else:
        raise SystemExit(f"unknown matrix {name!r}")
    return cells


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix")
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()

    outdir = os.path.join(os.path.dirname(__file__), "..", "results", "experiments", a.matrix)
    os.makedirs(outdir, exist_ok=True)
    cells = matrix(a.matrix, a.reps)
    print(f"matrix={a.matrix}  cells={len(cells)}  reps={a.reps}  workers={a.workers}")

    for (cfg_name, n, split, d, incl) in cells:
        key = f"{cfg_name}_n{n}_s{split}_d{d}"
        summ_path = os.path.join(outdir, f"{key}_summary.json")
        if os.path.exists(summ_path):
            print(f"[skip] {key} (done)"); continue
        t0 = time.time()
        try:
            rows, summ = E.run_cell(cfg_name, n, split, d, a.reps, include_existing=incl, workers=a.workers)
        except Exception as e:
            print(f"[FAIL] {key}  {type(e).__name__}: {str(e)[:120]}  (continuing)")
            continue
        dt = time.time() - t0
        with open(os.path.join(outdir, f"{key}_raw.jsonl"), "w") as f:
            for r in rows:
                f.write(json.dumps(r) + "\n")
        with open(summ_path, "w") as f:
            json.dump({"config": cfg_name, "n": n, "split": split, "d": d, "reps": a.reps,
                       "include_existing": incl, "elapsed_s": dt, "summaries": summ}, f, indent=2)
        cov = {m: round(summ[m]["coverage"], 3) for m in summ
               if isinstance(summ[m], dict) and "coverage" in summ[m]}
        errs = summ.get("_n_replication_errors", 0)
        print(f"[done] {key}  {dt:5.1f}s  errors={errs}  coverage={cov}")
    print("matrix complete.")


if __name__ == "__main__":
    main()

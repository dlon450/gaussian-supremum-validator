"""Acceptance gate: serial and parallel runs are bit-identical per replication.

    /tmp/gsv_venv/bin/python tests/test_parallel.py
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gsv import experiment as E


def key(r):
    return (r["rep"], r["method"])


def main():
    rows_serial, _ = E.run_cell("paper_ro_ellipsoid", n=100, split=0.5, d=5, reps=8, workers=1)
    rows_par, _ = E.run_cell("paper_ro_ellipsoid", n=100, split=0.5, d=5, reps=8, workers=4)
    s = {key(r): (round(r["feasible"], 12), round(r["objective"], 10), round(r["true_feas"], 10)) for r in rows_serial}
    p = {key(r): (round(r["feasible"], 12), round(r["objective"], 10), round(r["true_feas"], 10)) for r in rows_par}
    assert set(s) == set(p), f"key sets differ: {set(s) ^ set(p)}"
    mismatches = [k for k in s if s[k] != p[k]]
    if mismatches:
        for k in mismatches[:5]:
            print("MISMATCH", k, "serial", s[k], "parallel", p[k])
        print(f"FAIL: {len(mismatches)}/{len(s)} records differ")
        sys.exit(1)
    print(f"PASS: serial == parallel across {len(s)} (rep,method) records (workers 1 vs 4)")


if __name__ == "__main__":
    main()

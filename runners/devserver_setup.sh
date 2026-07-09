#!/usr/bin/env bash
# Set up the venv and launch the full experiment program on a devserver.
# Run this ON the devserver, from the synced repo directory (e.g. ~/gsv).
#
#   bash runners/devserver_setup.sh            # setup + launch full run under tmux
#   bash runners/devserver_setup.sh setup      # just build the venv
#
# Transfer the code first, FROM YOUR LAPTOP (code only; never the paper/ or results/):
#   rsync -av --relative \
#     gsv runners tests sca.py fast.py \
#     devvm7564.hil0:~/gsv/
#
# NOTE: do NOT `git push` to origin — it is the public repo and paper/ is unpublished.

set -euo pipefail
cd "$(dirname "$0")/.."
VENV="${VENV:-$HOME/gsv_venv}"
WORKERS="${WORKERS:-$(( $(nproc) - 2 ))}"
PIP_TRUST="--trusted-host pypi.org --trusted-host files.pythonhosted.org"

echo "== building venv at $VENV =="
python3 -m venv "$VENV"
"$VENV/bin/python" -m pip install --quiet --upgrade pip $PIP_TRUST || true
"$VENV/bin/python" -m pip install --quiet $PIP_TRUST numpy scipy pandas gurobipy cvxpy matplotlib
"$VENV/bin/python" - <<'PY'
import numpy, scipy, pandas, gurobipy, cvxpy, matplotlib
print("deps OK; gurobipy", gurobipy.gurobi.version(), "| cvxpy", cvxpy.__version__,
      "| solvers", [s for s in ['CLARABEL','HIGHS','SCS'] if s in cvxpy.installed_solvers()])
PY

# smoke + acceptance gates
"$VENV/bin/python" tests/test_infra.py
"$VENV/bin/python" tests/test_parallel.py

if [ "${1:-run}" = "setup" ]; then
  echo "setup complete. Launch with: bash runners/devserver_setup.sh"
  exit 0
fi

# single-threaded BLAS/solver per worker
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# full program at paper scale (N=1000); checkpointed + resumable, so safe to re-run
LOG="results/experiments/run.log"; mkdir -p results/experiments
CMD="
set -e
for m in dim nsweep robust split; do
  nice -n 10 $VENV/bin/python runners/run_matrix.py \$m --reps 1000 --workers $WORKERS
done
nice -n 10 $VENV/bin/python runners/run_matrix.py existing --reps 300 --workers $WORKERS
$VENV/bin/python runners/make_figures.py
echo ALL_DONE
"
if command -v tmux >/dev/null; then
  tmux new-session -d -s gsv "bash -lc '$CMD' 2>&1 | tee $LOG"
  echo "launched under tmux session 'gsv' (workers=$WORKERS). Watch: tmux attach -t gsv   or   tail -f $LOG"
else
  nohup bash -lc "$CMD" > "$LOG" 2>&1 &
  echo "launched with nohup (workers=$WORKERS). Watch: tail -f $LOG"
fi

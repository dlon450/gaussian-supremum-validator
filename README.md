# Gaussian-Supremum Validation for Data-Driven Chance-Constrained Programs

Code accompanying the paper [arXiv:1909.06477](https://arxiv.org/pdf/1909.06477).

The experiments study how to **choose the conservativeness parameter** of a
data-driven solution to a chance-constrained program

```
min  c^T x
s.t. P( xi^T x <= b ) >= 1 - alpha
```

Each base formulation (`CCP_*`) produces a family of candidate solutions indexed
by a parameter `delta` (an uncertainty-set size, a scenario count, …). A
**validation method** then uses held-out data to pick the parameter that is
cheap while still (statistically) feasible. The paper's proposal is the
**Gaussian-supremum validator**; the code compares it against several baselines.

## Repository layout

```
.
├── gsv/                       # shared library (import as `gsv.*`; no install needed)
│   ├── formulations.py        # base solvers
│   └── validators.py          # validation methods + run_experiment
├── dro_moment.py  ro_ellipsoid.py  dro_wasserstein.py  so.py  saa.py   # experiment drivers
├── fast.py  sca.py  knapsack.py  portfolio_dro.py                      # standalone baselines / demos
├── main.sh                    # SLURM submit helper (run from repo root)
├── requirements.txt
├── analysis/
│   └── result.ipynb           # figure/table generation (run from this directory)
├── results/
│   ├── final/                 # canonical, paper-ready CSVs
│   │   └── figures/           # canonical figures (*_all, *_comparison, …)
│   └── raw/                   # per-experiment intermediate dumps
│       ├── dro_wasserstein/  saa/  so/
├── paper/                     # LaTeX source + reviewer response + summary tables
│   └── experiment/
├── legacy/                    # earlier MATLAB / prototype code + .mat data
│   ├── code_from_Jierong/  code_huajie/  final_code_Linyun/  matlab_data/
└── neuripsGenAI/              # separate, unrelated to this repository tidy-up
```

### Shared library

The importable `gsv/` package holds the shared code. The drivers are run from the
repository root (`python dro_moment.py`), so `import gsv` resolves without any
installation step.

| Module | Contents |
|------|----------|
| `gsv/formulations.py` | Base optimization formulations (solvers), all sharing the interface `solve(para, c, b, data, alpha, d, n)`. |
| `gsv/validators.py`   | Data-driven parameter-selection methods and the Monte-Carlo driver `run_experiment`. |

Solvers in `gsv/formulations.py`:

- `CCP_DRO_moment` — moment-based distributionally robust CC (cvxpy).
- `CCP_RO_ellipsoid` — robust ellipsoidal approximation (Gurobi).
- `CCP_SO` — scenario optimization (Gurobi).
- `CCP_SAA` — sample-average approximation with a big-M relaxation (Gurobi).
- `CCP_DRO_wasserstein` — Wasserstein-ball DRO CC, big-M MILP (Gurobi).
- `mean_risk_portfolio_dro` — mean-CVaR portfolio under a Wasserstein set (rsome).

Validation methods in `gsv/validators.py`:

- `cross_validation` — K-fold cross validation.
- `bootstrapping` — bootstrap resampling with out-of-bag validation.
- `sectioning` — a single train/validate split scored over `K` sections.
- `gaussian_supremum` — returns four candidate solutions, one per rule:
  - **NGS** — Normalized Gaussian Supremum (this paper),
  - **UNGS** — Unnormalized Gaussian Supremum (this paper),
  - **NV** — Naive validation, no statistical margin (baseline),
  - **UG** — Univariate-Gaussian margin (baseline).

### Experiment drivers

Each script seeds NumPy, builds the parameter grid, and calls `run_experiment`
with the matching solver, writing results to CSV:

| Script | Solver | Methods compared |
|--------|--------|------------------|
| `dro_moment.py`      | `CCP_DRO_moment`      | sectioning + Gaussian supremum |
| `ro_ellipsoid.py`    | `CCP_RO_ellipsoid`    | sectioning + Gaussian supremum |
| `dro_wasserstein.py` | `CCP_DRO_wasserstein` | Gaussian supremum |
| `so.py`              | `CCP_SO`              | cross validation + bootstrapping |
| `saa.py`             | `CCP_SAA`             | cross validation + bootstrapping |

`so.py` and `saa.py` use a *discrete* parameter (a scenario / constraint count),
so their candidate grid is regenerated inside the validators via the
`cv_grid` / `bs_grid` hooks.

### Standalone baselines / demos

- `fast.py` — FAST scenario-optimization method.
- `sca.py` — safe convex approximation (Bonferroni/Nemirovski surrogate).
- `knapsack.py` — DRO chance-constrained knapsack (Gurobi lazy constraints).
- `portfolio_dro.py` — demo of the mean-risk portfolio DRO solver.

### Results, analysis, and archives

- `results/final/` — the canonical, paper-ready CSVs and (under `figures/`) the
  figures used in the write-up. `analysis/result.ipynb` reads from and writes to
  this directory.
- `results/raw/<experiment>/` — the per-experiment intermediate dumps
  (`objectives_*.csv`, `feas_levels_*.csv`, per-fold / per-split sub-runs).
- `analysis/result.ipynb` — regenerates the comparison figures/tables. Run it
  from the `analysis/` directory so its `../results/...` paths resolve.
- `paper/` — LaTeX source (`paper/experiment/`), the reviewer response, and the
  summary spreadsheets.
- `legacy/` — earlier MATLAB / prototype code and the original `.mat` result
  files, kept for reference.
- `neuripsGenAI/` — separate, unrelated to this repository.

## Running

```bash
pip install -r requirements.txt        # Gurobi additionally needs a license
python dro_moment.py                    # or so.py, ro_ellipsoid.py, saa.py, dro_wasserstein.py
```

Each driver writes `experiment_results_<name>.csv` (a summary row per
configuration) plus per-replication `objectives_<n>.csv` / `feas_levels_<n>.csv`
to the repository root. These are treated as transient scratch output (they are
git-ignored); the curated copies used for the paper live in `results/final/`.

## Reproducibility

Every driver seeds the global NumPy RNG (`np.random.seed(2)`, or `seed(1)` for
`portfolio_dro.py`). The validation routines draw all randomness through that
global RNG, so the sequence of draws — and hence the results — is deterministic
given the seed. A few drivers intentionally generate an unused
`data = np.abs(np.random.normal(...))` right after seeding; that line is kept
(and commented) because it advances the RNG stream and removing it would change
the published numbers.

Optimizer outputs (Gurobi / cvxpy / rsome) may differ in the last digits across
solver versions; the NumPy-driven data generation and validation are exactly
reproducible.

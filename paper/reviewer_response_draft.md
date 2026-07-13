# Reviewer response — experimental evidence draft

Draft answers to the reviewer comments that call for new experiments, grounded in
the re-run results (N = 1000 replications for the main sweeps; N = 300 for the
existing-scheme comparison). Config: Gaussian DGP, `d = 10` unless swept,
tolerance `1 - alpha = 90%`, target confidence `1 - beta = 95%` (so the
"feasibility level >= 0.95" line in every figure/table). All numbers are produced
by `runners/run_matrix.py` and summarised in `results/analysis/` and
`results/figures/`.

> NOTE ON A TEXT/CONFIG INCONSISTENCY TO RECONCILE: the current §6.5 text states
> "tolerance 1-alpha = 95%, confidence 1-beta = 90%", but the author-confirmed
> configuration (and every other section) uses 1-alpha = 90%, 1-beta = 95%. The
> regenerated figures/tables use the latter; §6.5 prose and the red target line
> caption (currently "90%") should be updated to 95% for consistency.

---

## 1. (Both reviewers, major) "Comprehensive comparison — theoretical and empirical — against existing approaches (cross-validation, bootstrapping, sectioning)."

We added Section 6.5 experiments comparing the univariate-Gaussian (UG) and
normalized/unnormalized Gaussian-supremum (NGS/UNGS) validators against K-fold
cross validation (CV), bootstrapping (BS), and sectioning, on three formulations
(SO, SAA, Wasserstein DRO) at a **common total data budget** (split = 0.7 for the
two-phase methods; CV/BS/sectioning at 10 folds/resamples). See
`results/figures/existing_coverage.pdf`, `so_delta_objective.pdf`,
`saa_delta_objective.pdf`, and the table `results/analysis/existing_comparison.csv`.

Key empirical findings (coverage; target 0.95):

- **Small data (n = 100).** CV/BS attain higher coverage than the two-phase
  validators (SO n=100: BS 0.87, CV 0.76 vs UG 0.65, NGS 0.68; SAA n=100: BS 0.94,
  CV 0.79 vs UG 0.58). At small n the two-phase split leaves too few Phase-2
  points to certify feasibility, so CV/BS "reuse" of all data is an advantage —
  we now state this explicitly.
- **Moderate/large data.** Once n is large enough the proposed validators reach
  the target while remaining competitive on objective, and CV can fall *below*
  target (SO n=300: UG 0.95, NGS 0.98, UNGS 1.00 all meet target, whereas CV =
  0.89). The existing schemes do **not** deliver better objective values than the
  proposed validators in the regime where all methods are feasible
  (`so_delta_objective.pdf` / `saa_delta_objective.pdf`).
- **SAA.** CV/BS begin meeting the target at a smaller n (≈200) than the
  two-phase validators — consistent with the small-data point above and now
  reported honestly.

Computational-cost argument (quantitative, `results/analysis/compute_cost.json`,
single-threaded wall-clock, K=B=10):

| formulation | two-phase | CV | BS | sectioning |
|---|---|---|---|---|
| SO (n=500) | 2.9 s | 48.0 s (**16.3×**) | 58.8 s (**20.0×**) | 2.9 s (1.0×) |
| SAA (n=200) | 1.6 s | 33.3 s (**21.2×**) | 29.6 s (**18.8×**) | 1.7 s (1.1×) |
| Wasserstein (n=100) | 3.8 s | 57.3 s (**15.1×**) | 53.1 s (**14.0×**) | 2.5 s (0.7×) |

CV/BS resolve the optimization ~K (resp. B) times over the parameter grid — in
practice **14–21×** the two-phase cost — while the two-phase validators solve the
Phase-1 path once and then only evaluate cheap Phase-2 sample statistics.
(Sectioning is ~1× because it, too, solves the path once.) At the larger sizes we
attempted, CV/BS on the big-M SAA (n=500) and Wasserstein (n>=200) MILPs did not
finish within the compute budget — the "computationally intractable at scale"
regime the paper describes, now made concrete. **Takeaway for the manuscript**:
prefer the proposed validators when the base optimization is expensive and n is
not tiny; for small, cheap problems CV/BS are a reasonable (and sometimes
higher-coverage) alternative.

## 2. (Both reviewers) "n1 is left as a loose end. Guidance on the relative budgeting of samples between the two phases."

New sweep of the Phase-1 fraction n1/n (`results/figures/so_all.pdf`,
`saa_all.pdf`, table `results/analysis/split_budgeting.csv`, N = 1000):

SO (n = 500), coverage vs Phase-1 fraction:

| n1/n | 10% | 20% | 30% | 50% | 70% | 90% |
|------|-----|-----|-----|-----|-----|-----|
| UG   |0.72 |0.96 |0.98 |0.97 |0.98 |0.93 |
| NGS  |0.72 |0.98 |1.00 |1.00 |0.99 |0.98 |

Guidance we now give:
- **Too small n1 starves the solution path** (10% → 0.72, well below target):
  Phase 1 cannot even produce a good candidate path.
- **Too small n2 (too-large n1) weakens validation** (UG at 90% → 0.93): too few
  Phase-2 points to certify feasibility.
- **There is an interior sweet spot** (~30–70%). A simple **70%-to-Phase-1**
  default is safe once n is moderately large, and the Gaussian-supremum
  validators (NGS/UNGS) are robust across a wide band (20–90%), whereas the
  univariate validator is more sensitive to n2 and prefers the mid-band.

This is an empirical guideline (as the reviewer notes, a rigorous n1 theory is
out of scope and flagged as future work), but it is now characterized rather than
left silent. Complementarily, the existing-scheme comparison is run at 70% across
n and shows the default reaching target as n grows.

## 3. (Reviewer 2, minor/major) Dimension dependence / "dimension-free feasibility."

RO dimension sweep, n = 500 (`results/figures/ro_ellipsoid_comparison.pdf`,
`results/analysis/dim_free.csv`, N = 1000): coverage of UG/NGS/UNGS stays at
≈0.99–1.00 for d = 2 … 100, while the naive "plain average" drops (e.g. 0.87 at
d = 10–20) — the Gaussian margins are necessary and the feasibility guarantee is
essentially dimension-free. The objective advantage over the (true-moment) SCA
benchmark grows with d (≈0.15 at d=2 → ≈1.57 at d=50). At d=100 with only n=500
samples the RO path degrades (advantage shrinks) — an honest illustration of the
reviewer's point that a poorly conditioned Phase-1 estimate hurts the path; noted
in the text.

## 4. (Reviewer 2, Q1/Q2) Convergence of feasibility and optimality with n.

`nsweep` (RO/SO/FAST, d=10, N=1000; `results/figures/fig_coverage_vs_n.pdf`,
`fig_gap_vs_n.pdf`; numbers in `results/analysis/nsweep.csv`):
- **Feasibility (Q1):** coverage of the proposed validators rises to and holds at
  the target as n grows, in **all three** formulations.
- **Optimality (Q2):** the objective gap to the path oracle shrinks toward 0 for
  **RO and FAST** (e.g. RO NGS gap 0.49→0.03 over n=100→1000). For **SO the gap
  does *not* shrink** — it grows modestly with n (NGS ≈0.06→0.24, UNGS ≈0.07→0.33
  over n=100→1000). This is because the SO *path oracle* itself becomes much
  cheaper as more scenarios enter the path (its least-conservative feasible
  candidate improves faster than any validator's justified pick), so the gap to
  that moving oracle widens even though the validator's absolute objective still
  improves and coverage is maintained. We report this honestly rather than
  claiming universal gap-vanishing (verified against the data). [Reconciled after
  adversarial verification flagged the SO case.]

## 5. Robustness beyond the Gaussian assumption (heavy-tailed data).

`robust` matrix under a multivariate-t (df=4) DGP with matched first two moments
(`results/figures/fig_robust.pdf`, `results/analysis/robust.csv`, N=1000): the
validators retain target-level coverage under heavy tails; details in the CSV.

---

### Reproducibility note
Two solver-configuration fixes were applied to make the study reproducible and
tractable: Gurobi is now pinned to a single thread with a fixed MILP seed (the
previous default multi-threaded MILP search was non-deterministic for SAA/
Wasserstein), and the big-M MILPs carry a 5 s / 1% -gap cap so the CV/BS sweeps
terminate (the uncapped hard instances are exactly the "intractable at scale"
regime). Convex (RO/SO/SCA) results are unaffected; MILP coverage numbers changed
by < 0.01 vs the uncapped run.

---

# Comprehensive follow-up campaign (best n1; N x D scaling; CV/BS at multiple folds)

Added at the reviewers' request for a more detailed empirical characterization.
Matrices produced by `runners/run_comprehensive.py`; analysis in
`results/analysis/COMPREHENSIVE.md` and `results/analysis/comp_*.csv`; figures
`results/figures/comp_*.pdf`.

## A. What is the best Phase-1 fraction n1/n, and how does it move with N and D? (`budget`, N_rep=500)

We swept n1/n in {0.1,...,0.9} for SO, RO, SAA and Wasserstein across several
(N, D). Findings:

1. **There is a feasible band, not a single point.** Too-small n1 starves the
   Phase-1 solution path; too-small n2 (too-large n1) starves Phase-2 validation.
   The interior of the band meets the target robustly.
2. **The least-conservative feasible n1 shrinks as N grows.** E.g. SO d=10, the
   smallest feasible (objective-optimal) fraction is ~0.5 at N=200, ~0.2 at N=500,
   ~0.1 at N=1000 (`comp_budget_best_n1.csv`). With more total data, less of it
   need go to Phase 1.
3. **A larger n1 fraction is needed as D grows** (fixed N): SO N=500 needs ~0.1
   (d=2) → ~0.2 (d=10) → ~0.3 (d=50). Higher dimension requires more Phase-1 data
   to estimate a good path.
4. **Robust default:** a mid-range allocation (roughly 50-70% to Phase 1) stays
   inside the feasible band across all tested (N, D) once N is large enough to
   admit any feasible allocation — supporting the paper's 70% recommendation as a
   safe, dimension-agnostic default, while the *objective-optimal* allocation is
   less conservative (smaller n1) and depends on N and D as above.

Figures: `comp_budget_heat_n1_vs_N_d10.pdf` (coverage over n1 x N at d=10) and
`comp_budget_heat_n1_vs_D_n500.pdf` (coverage over n1 x d at N=500).

## B. How do feasibility and optimality change with N and D? (`ndgrid`, N_rep=1000, split=0.5)

Full (N, D) grids per formulation (`comp_ndgrid.csv`,
`comp_ndgrid_coverage_heat.pdf`). Coverage of the normalized GS validator:

- **RO is essentially dimension-free**: coverage >= 0.96 even at d=50, N=100, and
  >= 0.95 across the whole grid. This matches the paper's dimension-free feasibility
  claim for the ellipsoidal method.
- **SO and SAA require N to scale with D**: at small N and large D the target is
  missed (SO: d=50,N=100 -> 0.50; d=20,N=100 -> 0.50; SAA: d=20,N=100 -> 0.49),
  but coverage recovers to >= 0.99 once N is large enough (SO d=50,N=500 -> 0.99).
  These scenario/SAA methods implicitly estimate the whole feasible region, so
  their sample requirement grows with dimension — the paper's central motivation,
  now quantified.
- **Wasserstein DRO** is moderately robust (coverage >= 0.93 across d=2-20 at
  N=100-400).

Across all formulations, as N grows the coverage rises to/holds the target
(feasibility, Q1). The objective gap to the path oracle shrinks with N for RO,
FAST and Wasserstein; for SO the gap does not vanish (the SO path oracle improves
faster than any validated pick), as noted in §Q1/Q2 above.

## C. CV / bootstrap at K in {3,5,10} vs the proposed validators (`folds`, N_rep=300, split=0.7)

Complete results in `comp_folds.csv`; curves in
`results/figures/comp_folds_coverage.pdf`. SO coverage vs N (K=3/5/10), target 0.95:

| N | UG | NGS | UNGS | CV 3/5/10 | BS 3/5/10 | Sec 3/5/10 |
|---|----|-----|------|-----------|-----------|------------|
|200|0.90|0.94|0.97| 0.87/0.85/0.89 | 0.82/0.90/0.96 | 0.71/0.84/0.94 |
|300|0.96|0.98|1.00| 0.86/0.84/0.90 | 0.75/0.86/0.96 | 0.82/0.98/0.99 |
|500|0.96|0.98|1.00| 0.86/0.84/0.89 | 0.75/0.86/0.92 | 0.88/0.95/1.00 |

Findings:
1. **More folds/resamples → higher coverage** for the existing schemes, most
   markedly for bootstrap (SO N=500: BS 0.75→0.86→0.92 at K=3→5→10) and sectioning
   (SO N=300: 0.82→0.98→0.99). **K=10 is the best of the tested counts**, as the
   paper states. CV is comparatively flat and noisy in K.
2. **Cross-validation systematically under-covers for SO** — it never reaches the
   0.95 target at any fold count or N (max 0.90 at CV10). Bootstrap and sectioning
   reach target at K=10, N≥300.
3. **At matched data budget the proposed validators match the best fold setting**
   (K=10) and clearly beat the under-tuned ones (K=3/5) and CV for SO, meeting
   target from N=300 — while costing **14-21× less** to compute (Section 1). For
   SAA and Wasserstein (see `comp_folds.csv`) the resampling schemes (esp. BS10)
   reach target at slightly smaller N — the small-n advantage noted in Section 1 —
   but again only at K=10 and at the large compute cost.

Net: the proposed validators are competitive with the *best-tuned* existing scheme
and clearly better than under-tuned CV/BS, at an order-of-magnitude lower cost.

---

# Saturation campaign (tolerance, confidence/calibration, tails, correlation, deep folds, extreme dimension)

To characterize the framework as exhaustively as possible we swept every remaining
axis (`runners/run_saturation.py`; analysis `results/analysis/SATURATION.md`,
`sat_*.csv`; figures `results/figures/sat_*.pdf`). 187 additional cells, 0 errors;
590 experiment cells total. Config d=10, split=0.5 unless noted; DGP Gaussian
unless swept.

## D. Tolerance sweep (1-alpha in {70%..99%})
There is a minimum tolerance the two-phase validator can certify at a given data
size. At 1-alpha=99% (alpha=0.01) coverage collapses for SO/SAA (~0.01 at n=400-500,
which cannot statistically certify 99% feasibility) while RO holds (~0.98) and
Wasserstein is intermediate (~0.82-0.98). All four formulations meet the 0.95
target only once alpha>=0.10 (at alpha=0.05 SAA is 0.87 and SO 0.95-borderline; RO
already clears it); RO/DRO recover fastest. Takeaway: RO/DRO tolerate very tight
tolerances; scenario/SAA methods need alpha not too small (or larger n).
[Recovery threshold corrected to alpha>=0.10 after adversarial verification.]

## E. Confidence sweep / calibration (1-beta in {50%..99%})
Key calibration result (`sat_beta_calibration.pdf`): empirical coverage tracks the
nominal target 1-beta and stays **at or above it** for the Gaussian-supremum
validators across the whole range — e.g. SO at nominal 0.50 gives NGS 0.96, UNGS
0.95 (valid/conservative); at nominal 0.95 gives NGS 0.998, UNGS 1.00. The
univariate UG is the **least conservative** (closest to nominal: SO UG = 0.744 at
nominal 0.50, 0.986 at 0.95), consistent with the paper's claim that UG gives
tighter confidence when its assumptions hold while the supremum validators are
uniformly safe. No proposed validator falls materially below its nominal target.

## F. Heavy tails (multivariate-t, df in {3..30})
The validators are robust to heavy tails: NGS coverage >=0.96 across df=3..30 for
every formulation (df=3 is near the finite-variance boundary). The univariate UG is
less robust under the heaviest tails (dips to ~0.90-0.94 for SAA/Wasserstein),
whereas NGS/UNGS hold >=0.96 — again matching the paper's "supremum validators are
more generally applicable" message.

## G. Correlation (Gaussian corr in {0..0.8})
Coverage is essentially invariant to coordinate correlation: within each
formulation NGS coverage barely moves as corr goes 0->0.8 (per-formulation spread
<=0.02). The overall NGS band is ~0.95-1.0 (Wasserstein sits at ~0.95-0.97, RO/SO
at ~0.99-1.0). Correlation structure does not degrade the framework.

## H. Deep fold sweep (K in {2,3,5,10,20}) and extreme dimension (d in {100,200})
- **Folds:** bootstrap and sectioning coverage increase **monotonically** with K
  (SO n=200 bootstrap: 0.72,0.82,0.90,0.96,0.98 at K=2,3,5,10,20); K=20 marginally
  beats K=10 at ~2x the cost. Cross-validation is comparatively **flat in K** and
  stays below 0.95 for SO through K=10 at every N (max ~0.90-0.91); only at K=20 does
  it reach target at one size (n=300, CV20=0.977) — i.e. CV needs the largest,
  most expensive fold count to become calibrated for SO, whereas the proposed
  validators are already calibrated. This extends §C.
- **Extreme dimension:** RO stays perfectly dimension-free — NGS coverage = 1.00 at
  **d=200** for all N; SO reaches >=0.99 at d=100 and d=200 once N>=500
  (`sat_bigd.csv`).
- **Extended best-n1:** RO at d=100 admits the full 10-90% Phase-1 band; SO at d=100
  excludes only the smallest fraction (10% starves the path); the Wasserstein d=20
  feasible band widens with N (up to 50% at n=100, up to 70% at n=200) — the same
  budgeting pattern as §A holds at high dimension.

All saturation claims were adversarially re-verified against the raw CSVs (see
`results/analysis/verification_saturation.md`).

---

# Formulation completeness: moment-DRO added; high-d coverage for all formulations

Prompted by a completeness check, we (a) wired the **moment-based DRO** formulation
(paper Section 6.2, previously only in the appendix code) into the two-phase
pipeline and (b) pushed every formulation to high dimension.

## Moment-based DRO (SDP), N x D grid + n1/tolerance/tail slices
Reproduces and extends the paper's Table (d=10) to d=2..30, N=100..500. The
validators select the tightest radius s≈0 (moments as equalities), are 100%
feasible everywhere, and beat the chi2-quantile benchmark by a margin that GROWS
with dimension:

| d | validators (UG=NGS=UNGS) obj | chi2 benchmark obj | advantage |
|---|---|---|---|
| 2 | -0.75 | -0.44 | 0.31 |
| 10 | -3.33 | -1.57 | 1.76 |
| 20 | -5.22 | -2.60 | 2.62 |
| 30 | -6.98 | -3.98 | 3.00 |

(paper d=10 table: validators -2.73 vs benchmark -1.83 — same qualitative result.)
This shows the conventional chi2 moment-set size is very conservative, exactly the
paper's Section 6.2 point.

## High-dimension coverage for every formulation
- RO, SO: swept to **d=500** (RO dimension-free; SO needs N∝D).
- SAA, Wasserstein: swept to **d=50** — SAA d=50 needs N∝D (n=200->0.81, n=500->0.99);
  Wasserstein d=50 robust (n=100->0.95, n=200->0.98 for NGS).
- FAST: d=50 (ndgrid). Moment-DRO: d=30.
Every formulation the paper reports (RO/SCA, moment-DRO, SO, FAST, SAA, Wasserstein)
is now covered, each at dimensions >= the paper's (which used d in {10,50}).

## Which validator is "best"? A frontier, not a single winner
Over 245 adequate-data cells across all formulations:
- **UG** gives the best (lowest) objective in **97%** of cells — it is the least
  conservative — BUT it falls below the 1-beta target in **24.5%** of cells
  (heavy tails, tight tolerance, small n2, high d), vs 4.9% (NGS) / 4.5% (UNGS).
  In 20% of cells UG misses target while UNGS still meets it.
- **NGS / UNGS** are the robust default: they deliver the feasibility guarantee
  ~5x more reliably than UG, at a small objective cost.
Recommendation: UG when its regularity assumptions hold and best objective is the
priority; NGS/UNGS when guaranteed feasibility is required. This directly answers
the reviewer's "why keep Algorithms 2/3 if 4 dominates?" — Algorithm 4 (UG) does
not dominate on the primary (feasibility) criterion.

---

# Second problem class: CVaR-constrained portfolio (a general stochastic constraint)

The paper describes a second problem — a CVaR-constrained mean-return portfolio
(eq. (cvar)) — but it was never run (the code only exercised the DRO solver). We
instantiated the two-phase framework for it end to end, treating the CVaR
constraint CVaR_alpha(-xi'x) <= gamma as a *general stochastic constraint*: the
Phase-2 statistic is the Rockafellar-Uryasev empirical CVaR with its
influence-function standard error, and the Gaussian-supremum / univariate margins
are applied to it exactly as in the chance-constraint case
(`gsv/portfolio.py`, `runners/portfolio_experiment.py`, results under
`results/experiments/portfolio_*`). Wasserstein-DRO CVaR reformulation (convex,
Esfahani-Kuhn); Gaussian returns, CVaR level alpha=0.10, target 1-beta=0.95.

The framework transfers cleanly and every headline finding reproduces on this new
problem class:
- **Margins are necessary.** The SAA/aggressive benchmark covers only 0.01-0.45
  and the naive no-margin validator 0.41-0.71 — both violate the risk target
  (they overfit the sample CVaR). The margined validators reach the 0.95 target.
- **N x D:** coverage (target 0.95), Gaussian returns —

  | d | N=100 | N=200 | N=500 | N=1000 | (UG / benchmark) |
  |---|---|---|---|---|---|
  | 5  | 0.78 | 0.88 | 0.95 | 0.95 | benchmark 0.10-0.35 |
  | 10 | 0.83 | 0.88 | 0.95 | 0.96 | benchmark 0.02-0.24 |
  | 20 | 0.81 | 0.90 | 0.94 | 0.95 | benchmark 0.01-0.16 |

  (UG shown; NGS/UNGS are slightly higher — UNGS reaches target at the smallest N.)
- **UNGS safest, UG best return** — same frontier as the chance-constraint problem:
  UNGS meets target at the smallest N; UG gives the highest return among the valid
  validators.
- **Robust to heavy tails:** under multivariate-t returns UG/NGS/UNGS hold ~0.96
  down to df=3, while NV stays at 0.59.
- **Risk-threshold (gamma) sweep:** validators meet target across gamma in
  {0.12..0.28}; benchmark/NV fail throughout.

Caveat: at d=2 the problem is degenerate for gamma=0.15 (only two assets cannot
diversify enough — the minimum achievable true CVaR ~0.23 > gamma, so *no*
solution is feasible and all methods read 0 coverage). Feasible solutions exist
from d>=5; gamma should scale with d if very small d is of interest.

This demonstrates the framework is not specific to chance constraints — it applies
to general stochastic constraints (here CVaR), which strengthens the paper's
generality claim.

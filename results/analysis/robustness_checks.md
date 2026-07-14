# Robustness / validity checks
Target coverage 1-beta = 0.95.
## 1. Independent-stream confirmation of recommended Phase-1 splits
In-sample = same seed family used to pick the split; Out-of-sample = a disjoint RNG stream (seed_offset). If coverage holds out-of-sample, the recommendation is not an artifact of selection optimism.
- paper_so s=0.5 UG: in=0.972  out=0.952 (out CI_lo=0.93)
- paper_so s=0.5 NGS: in=0.994  out=0.99 (out CI_lo=0.977)
- paper_so s=0.5 UNGS: in=0.998  out=1.0 (out CI_lo=0.992)
- paper_ro_ellipsoid s=0.3 UG: in=0.996  out=0.992 (out CI_lo=0.98)
- paper_ro_ellipsoid s=0.3 NGS: in=1.0  out=0.998 (out CI_lo=0.989)
- paper_ro_ellipsoid s=0.3 UNGS: in=0.996  out=0.992 (out CI_lo=0.98)
- paper_saa s=0.5 UG: in=0.926  out=0.914 (out CI_lo=0.886)
- paper_saa s=0.5 NGS: in=0.988  out=0.988 (out CI_lo=0.974)
- paper_saa s=0.5 UNGS: in=0.988  out=0.992 (out CI_lo=0.98)

## 2. Solution-path mesh-size (p) sensitivity — RO d=10 n=500
- p=25: UG cov=0.99 NGS cov=0.993 UG obj=-7.448
- p=50: UG cov=0.983 NGS cov=0.993 UG obj=-7.582
- p=100: UG cov=0.963 NGS cov=0.997 UG obj=-7.649

## 3. Gaussian-supremum simulation-count (sim_num) sensitivity
- sim_num=1000: qhat mean=2.1694 std=0.0501
- sim_num=2000: qhat mean=2.1683 std=0.0498
- sim_num=5000: qhat mean=2.1682 std=0.0294
- sim_num=10000: qhat mean=2.1767 std=0.0197
(2000 draws already give a stable qhat; std shrinks ~1/sqrt(sim_num).)

## 4. Generality: negatively-correlated data (corr=-0.08 vs 0.0)
- paper_ro_ellipsoid corr=-0.08: UG=0.993 NGS=1.0 UNGS=0.993
- paper_ro_ellipsoid corr=0.0: UG=0.99 NGS=0.997 UNGS=0.99
- paper_so corr=-0.08: UG=0.97 NGS=0.997 UNGS=1.0
- paper_so corr=0.0: UG=0.98 NGS=0.997 UNGS=1.0
- paper_saa corr=-0.08: UG=0.923 NGS=0.987 UNGS=0.993
- paper_saa corr=0.0: UG=0.94 NGS=0.993 UNGS=0.99

## 5. MILP time-cap sensitivity (SAA two-phase, n=300 d=10)
- cap=5s : UG 0.957, NGS 0.980, UNGS 0.983 | UG obj -7.501
- cap=30s: UG 0.957, NGS 0.980, UNGS 0.983 | UG obj -7.501
Identical — the 5s cap does not bind at these sizes (solves finish to the 1% gap
well under 5s), so it does not affect the reported coverage/objective.

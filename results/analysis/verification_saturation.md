# Adversarial verification of saturation-campaign claims

7 claims (tolerance, calibration, tails, correlation, deep folds, extreme-d, extended best-n1)
each independently fact-checked against sat_*.csv by a separate agent instructed to refute.
Outcome: **4 SUPPORTED, 3 PARTIAL, 0 REFUTED**.

Corrections applied to the response draft:
1. alpha: all formulations meet target at alpha>=0.10 (not 0.05; at 0.05 SAA=0.87, SO=0.946).
2. corr: NGS coverage floor is ~0.95 (Wasserstein 0.954-0.974), not ~0.98.
3. foldsdeep/CV: CV under-covers for SO through K=10 but reaches target at K=20,n=300 (0.977);
   scope the 'CV never meets' statement to K<=10.

Raw workflow output:

```
{
  "summary": "Adversarially verify the saturation-campaign claims against sat_*.csv",
  "agentCount": 8,
  "logs": [],
  "result": {
    "verdicts": [
      {
        "claim_id": "alpha-tight-collapse",
        "verdict": "PARTIAL",
        "evidence": "Read results/analysis/sat_alpha.csv, method NGS. The direction of the collapse claim holds and Part-1 values match within tolerance: at alpha=0.01, coverage collapses for paper_so (0.006) and paper_saa (0.006) — exactly the quoted ~0.006 — while paper_ro_ellipsoid stays high at 0.984 (matches ~0.98). However Part 2 (\"by alpha>=0.05 all formulations meet the 0.95 target\") is refuted by the data: at alpha=0.05 paper_ro_ellipsoid=0.984 (meets=True) but paper_saa=0.87 (meets=False, 0.08 below target, beyond 0.03) and paper_so=0.946 (meets=False). All three formulations only clear 0.95 at alpha>=0.10 (paper_saa=0.992, paper_so=0.998, paper_ro_ellipsoid=0.996). So the recovery threshold should be alpha>=0.10, not alpha>=0.05 — the collapse narrative is correct but the \"meets by 0.05\" number is off.",
        "correction": "The tight-tolerance collapse is accurate (SO and SAA ~0.006, RO ~0.984 at alpha=0.01), but all formulations do NOT meet the 0.95 target by alpha=0.05 — at alpha=0.05 SAA is only 0.87 and SO is 0.946 (both meets=False). All three formulations meet 0.95 only at alpha>=0.10."
      },
      {
        "claim_id": "beta-calibration",
        "verdict": "SUPPORTED",
        "evidence": "From results/analysis/sat_beta.csv (config paper_so, nominal=1-beta), coverage per beta [nominal | NGS, UNGS, UG]: 0.01 [0.99 | 0.998, 1.000, 0.994]; 0.02 [0.98 | 0.994, 1.000, 0.990]; 0.05 [0.95 | 0.998, 1.000, 0.986]; 0.10 [0.90 | 0.990, 0.998, 0.964]; 0.15 [0.85 | 0.986, 0.998, 0.954]; 0.20 [0.80 | 0.982, 0.990, 0.930]; 0.30 [0.70 | 0.964, 0.982, 0.888]; 0.40 [0.60 | 0.966, 0.964, 0.838]; 0.50 [0.50 | 0.956, 0.954, 0.744]. NGS and UNGS are at or above nominal for every beta (valid/conservative). UG has the smallest gap to nominal at every beta (least conservative), with gap rising from +0.004 (beta=0.01) to +0.244 (beta=0.50), and UG dips lowest at large beta (0.744 at beta=0.50 vs NGS 0.956, UNGS 0.954). Quoted example matches: beta=0.50 NGS 0.956~0.96, UNGS 0.954~0.95, UG 0.744~0.74, nominal 0.50 (all within ~0.004)."
      },
      {
        "claim_id": "tails-robust",
        "verdict": "SUPPORTED",
        "evidence": "Source: /home/dereklong/local/scratch/gaussian-supremum-validator/results/analysis/sat_tails.csv (90 rows; methods NGS,UG,NV,UNGS,benchmark; df values {3,4,6,10,30} spanning the claimed {3..30}; 4 formulations: paper_dro_wasserstein, paper_ro_ellipsoid, paper_saa, paper_so).\n\nPart 1 (NGS >=0.96 across all df and formulations, robust near df=3): CONFIRMED. All 20 NGS rows have coverage in [0.966, 1.000], min=0.966. Every value is >=0.96 (min exceeds threshold by 0.006). At the finite-variance boundary df=3: wasserstein=0.966, ro_ellipsoid=0.998, saa=0.982, so=1.000 — all >=0.96. Per-formulation NGS minima: wasserstein 0.966, ro_ellipsoid 0.990, saa 0.978, so 0.994.\n\nPart 2 (UG less robust, dips ~0.90-0.94 for SAA/Wasserstein): CONFIRMED direction, quoted values match within tolerance. UG wasserstein = [0.908, 0.918, 0.930, 0.932, 0.936] (min 0.908, max 0.936 — squarely in 0.90-0.94). UG saa = [0.894, 0.912, 0.912, 0.928, 0.942] (min 0.894 at df=30, max 0.942). By contrast UG for ro_ellipsoid (0.982-0.998) and so (0.976-0.990) stay high, so the weakness is specific to SAA/Wasserstein, and UG is uniformly below NGS on those two.\n\nMinor note: the SAA UG minimum is 0.894 (df=30), just below the quoted \"0.90\" floor by 0.006 — well within the ~0.03 tolerance — so the direction and quoted range both hold."
      },
      {
        "claim_id": "corr-invariant",
        "verdict": "PARTIAL",
        "evidence": "From results/analysis/sat_corr.csv (method=NGS, d=10), coverage by corr (0.0,0.2,0.4,0.6,0.8): paper_dro_wasserstein = 0.968, 0.972, 0.954, 0.974, 0.960 (range 0.954-0.974, spread 0.020); paper_ro_ellipsoid = 0.996, 1.0, 1.0, 0.998, 0.994 (range 0.994-1.0); paper_saa = 0.992, 0.980, 0.984, 0.986, 0.978 (range 0.978-0.992); paper_so = 0.998, 0.996, 0.990, 0.994, 0.996 (range 0.990-0.998). The INVARIANCE direction clearly holds: within each formulation coverage barely changes as corr goes 0.0->0.8 (per-config spread only 0.006-0.020). However the quoted band \"~0.98-1.0\" is off for one whole formulation: all five paper_dro_wasserstein points (0.954-0.974) fall below 0.98. Overall NGS d=10 min=0.954, max=1.000; the min is ~0.026 below the quoted 0.98 lower bound (borderline within ~0.03 tolerance, but it is an entire formulation rather than a single stray point).",
        "correction": "Coverage IS essentially invariant to correlation (correct direction), but the true NGS d=10 range across formulations is ~0.95-1.0, not ~0.98-1.0: paper_dro_wasserstein sits at 0.954-0.974 for all corr values, below the stated 0.98 lower bound."
      },
      {
        "claim_id": "foldsdeep-bs-monotone",
        "verdict": "PARTIAL",
        "evidence": "From results/analysis/sat_foldsdeep.csv, config=paper_so, n=200 (K parsed from method name). Bootstrap (BS) coverage is strictly monotone increasing in K and matches the quoted values within ~0.003: BS2/K=2=0.7167 (~0.72), BS3/K=3=0.8167 (~0.82), BS5/K=5=0.8967 (~0.90), BS10/K=10=0.9633 (~0.96), BS20/K=20=0.9833 (~0.98). Cross-validation (CV) at n=200 is comparatively flat and stays below 0.95: CV2=0.8900, CV3=0.8700, CV5=0.8500, CV10=0.8933, CV20=0.8767 (max 0.893). However, the absolute clause \"CV never reaches 0.95 for SO at any K\" is false for the paper_so config beyond n=200: at n=300, CV20=0.9767 (>=0.95); n=500 CV20=0.9300, n=1000 CV2=0.9367 stay below. The paper_so CV max coverage across all n is 0.9767.",
        "correction": "The BS-monotone assertion and all five quoted n=200 values are correct. Scope the CV statement to n=200: \"CV stays below 0.95 at n=200 (max 0.893).\" As written, \"CV never reaches 0.95 for SO at any K\" is contradicted within the same config at n=300, where CV20 = 0.9767."
      },
      {
        "claim_id": "bigd-RO-dimfree-200",
        "verdict": "SUPPORTED",
        "evidence": "Source: /home/dereklong/local/scratch/gaussian-supremum-validator/results/analysis/sat_bigd.csv (method=NGS). Tested d={100,200}, n={200,500,1000,2000}.\n\nPart 1 (NGS coverage = 1.0 at d=200 for all tested N; config paper_ro_ellipsoid): d=200,n=200 -> 1.0; n=500 -> 1.0; n=1000 -> 1.0; n=2000 -> 1.0. All exactly 1.0. Confirmed.\n\nPart 2 (SO reaches >=0.99 at d=100 and d=200 once N>=500; config paper_so): d=100 -> n=500:0.996, n=1000:0.998, n=2000:1.0; d=200 -> n=500:0.992, n=1000:0.996, n=2000:1.0. All six values with N>=500 are >=0.99 (minimum 0.992 at d=200,n=500). Note paper_so has no n=200 row, so N>=500 covers every SO datapoint. Confirmed.\n\nBoth sub-claims match exactly (well within ~0.03 tolerance); direction and quoted values hold."
      },
      {
        "claim_id": "budgetfull-bands",
        "verdict": "SUPPORTED",
        "evidence": "Source: /home/dereklong/local/scratch/gaussian-supremum-validator/results/analysis/sat_budgetfull.csv, method=NGS, target=0.95, columns split (Phase-1 fraction) and coverage.\n\n(1) RO full band — config paper_ro_ellipsoid d=100: for BOTH n=500 and n=1000, splits 0.1 through 0.9 all have coverage=1.0 (meets=True). Full feasible band 10%-90%. SUPPORTED.\n\n(2) SO excludes smallest fraction — config paper_so d=100: n=500 split0.1 coverage=0.586 (fails), splits 0.2-0.9 coverage 0.954/0.98/0.992/0.996/0.992/0.99/0.982/0.962 (all meet). n=1000 split0.1 coverage=0.928 (fails), splits 0.2-0.9 all >=0.99 (meet). 10% fails, 20%-90% meet in both cases. SUPPORTED.\n\n(3) Wasserstein band widens with N — config paper_dro_wasserstein d=20: n=100 feasible up to split0.5 (0.5 coverage=0.95 meets=True; 0.6 coverage=0.948 fails) -> up to 50%; n=200 feasible up to split0.7 (0.7 coverage=0.956 meets=True; 0.8 coverage=0.92 fails) -> up to 70%. Band widens 50%->70% as n grows 100->200. SUPPORTED.\n\nAll quoted directions and boundary values match the CSV exactly (within ~0.03)."
      }
    ],
    "report": "# Saturation-Campaign Faithfulness Report\n\n## Verdict Counts (7 claims)\n- SUPPORTED: 4\n- PARTIAL: 3\n- REFUTED: 0\n\n## Claims Needing Correction (3)\n\n**1. alpha-tight-collapse [PARTIAL]**\n- What holds: The tight-tolerance collapse is real. At alpha=0.01, SO and SAA collapse to ~0.006 while RO_ellipsoid stays at 0.984 — direction and Part-1 numbers match.\n- Overstatement: The \"by alpha>=0.05 all formulations meet the 0.95 target\" claim is wrong. At alpha=0.05, SAA=0.87 (0.08 below target) and SO=0.946 both fail; only RO meets.\n- FIX: Change the recovery threshold from alpha>=0.05 to **alpha>=0.10**. All three formulations clear 0.95 only at alpha>=0.10 (SAA=0.992, SO=0.998, RO=0.996).\n\n**2. corr-invariant [PARTIAL]**\n- What holds: Invariance to correlation is clear — within each formulation, coverage spread across corr 0.0->0.8 is only 0.006-0.020.\n- Overstatement: The quoted band \"~0.98-1.0\" is wrong for an entire formulation. All five paper_dro_wasserstein points sit at 0.954-0.974, below 0.98.\n- FIX: State the true NGS d=10 range as **~0.95-1.0**, not ~0.98-1.0.\n\n**3. foldsdeep-bs-monotone [PARTIAL]**\n- What holds: BS coverage is strictly monotone in K, and all five quoted n=200 values match (0.72/0.82/0.90/0.96/0.98).\n- Overstatement: \"CV never reaches 0.95 for SO at any K\" is false — it is contradicted within the same config at n=300, where CV20=0.9767.\n- FIX: Scope the CV claim to n=200: **\"CV stays below 0.95 at n=200 (max 0.893).\"**\n\n## Fully Supported (4)\n- **beta-calibration** — NGS/UNGS at or above nominal for every beta; UG least conservative and dips lowest at large beta. Quoted beta=0.50 values match.\n- **tails-robust** — NGS >=0.96 across all df/formulations (min 0.966, robust at df=3); UG weakness isolated to SAA/Wasserstein (~0.90-0.94). Only nit: SAA UG floor is 0.894 vs quoted 0.90, well within tolerance.\n- **bigd-RO-dimfree-200** — RO=1.0 at d=200 for all N; SO >=0.99 at d=100/200 once N>=500 (min 0.992). Exact match.\n- **budgetfull-bands** — RO full 10%-90% band, SO excludes only the 10% split, Wasserstein band widens 50%->70% as n grows 100->200. All boundary values match.\n\n## Overall Faithfulness Verdict\nDirectionally faithful, with three quantitative overstatements. Every core narrative (collapse, calibration, tail robustness, correlation invariance, bootstrap monotonicity, dimension-free RO, budget bands) survives verification, and no claim was refuted outright. But the three PARTIALs share a common failure mode: **stated recovery thresholds and coverage bands are rosier than the data.** Specifically, formulations recover at alpha>=0.10 not 0.05; the correlation-invariance band bottoms at ~0.95 not ~0.98; and the absolute \"CV never reaches 0.95\" claim breaks at n=300. These are precision/scope errors, not fabrications — fix the three numbers above and the campaign is faithful."
  },
  "workflowProgress": [
    {
      "type": "workflow_phase",
      "index": 1,
      "title": "Verify"
    },
    {
      "type": "workflow_phase",
      "index": 2,
      "title": "Synthesize"
    },
    {
      "type": "workflow_agent",
      "index": 1,
      "label": "verify:alpha-tight-collapse",
      "phaseIndex": 1,
      "phaseTitle": "Verify",
      "agentId": "af1e013f915a23f4d",
      "agentType": "general-purpose",
      "model": "claude-opus-4-8[1m]",
      "state": "done",
      "startedAt": 1783946638957,
      "queuedAt": 1783946638942,
      "attempt": 1,
      "lastToolName": "StructuredOutput",
      "lastToolSummary": "alpha-tight-collapse",
      "promptPreview": "Adversarially fact-check a claim against saturation-experiment CSV output.\nRepo root: /home/dereklong/local/scratch/gaussian-supremum-validator. CLAIM (alpha-tight-collapse): At very tight tolerance (alpha=0.01, i.e. 1-alpha=99%) coverage collapses for SO and SAA (~0.006) while RO stays high (~0.98); by alpha>=0.05 all formulations meet the 0.95 target (NGS).\nRead: results/analysis/sat_alpha.csv (…",
      "lastProgressAt": 1783946692428,
      "tokens": 20984,
      "toolCalls": 6,
      "durationMs": 53471,
      "resultPreview": "{\"claim_id\":\"alpha-tight-collapse\",\"verdict\":\"PARTIAL\",\"evidence\":\"Read results/analysis/sat_alpha.csv, method NGS. The direction of the collapse claim holds and Part-1 values match within tolerance: at alpha=0.01, coverage collapses for paper_so (0.006) and paper_saa (0.006) — exactly the quoted ~0.006 — while paper_ro_ellipsoid stays high at 0.984 (matches ~0.98). However Part 2 (\\\"by alpha>=0.0…"
    },
    {
      "type": "workflow_agent",
      "index": 2,
      "label": "verify:beta-calibration",
      "phaseIndex": 1,
      "phaseTitle": "Verify",
      "agentId": "ade1655c2468bea78",
      "agentType": "general-purpose",
      "model": "claude-opus-4-8[1m]",
      "state": "done",
      "startedAt": 1783946638958,
      "queuedAt": 1783946638943,
      "attempt": 1,
      "lastToolName": "StructuredOutput",
      "lastToolSummary": "beta-calibration",
      "promptPreview": "Adversarially fact-check a claim against saturation-experiment CSV output.\nRepo root: /home/dereklong/local/scratch/gaussian-supremum-validator. CLAIM (beta-calibration): Calibration: as nominal target 1-beta decreases (beta rises 0.01->0.50), empirical coverage of NGS/UNGS stays at or above the nominal target (valid/conservative), while UG is the least conservative (closest to nominal, dipping lo…",
      "lastProgressAt": 1783946686072,
      "tokens": 27482,
      "toolCalls": 5,
      "durationMs": 47114,
      "resultPreview": "{\"claim_id\":\"beta-calibration\",\"verdict\":\"SUPPORTED\",\"evidence\":\"From results/analysis/sat_beta.csv (config paper_so, nominal=1-beta), coverage per beta [nominal | NGS, UNGS, UG]: 0.01 [0.99 | 0.998, 1.000, 0.994]; 0.02 [0.98 | 0.994, 1.000, 0.990]; 0.05 [0.95 | 0.998, 1.000, 0.986]; 0.10 [0.90 | 0.990, 0.998, 0.964]; 0.15 [0.85 | 0.986, 0.998, 0.954]; 0.20 [0.80 | 0.982, 0.990, 0.930]; 0.30 [0.70…"
    },
    {
      "type": "workflow_agent",
      "index": 3,
      "label": "verify:tails-robust",
      "phaseIndex": 1,
      "phaseTitle": "Verify",
      "agentId": "ab0e0a471b88d16d6",
      "agentType": "general-purpose",
      "model": "claude-opus-4-8[1m]",
      "state": "done",
      "startedAt": 1783946638958,
      "queuedAt": 1783946638943,
      "attempt": 1,
      "lastToolName": "StructuredOutput",
      "lastToolSummary": "tails-robust",
      "promptPreview": "Adversarially fact-check a claim against saturation-experiment CSV output.\nRepo root: /home/dereklong/local/scratch/gaussian-supremum-validator. CLAIM (tails-robust): Under heavy tails (multivariate-t), NGS coverage stays >=0.96 across df in {3..30} for all formulations, i.e. robust even near the finite-variance boundary df=3; the univariate UG is less robust (dips to ~0.90-0.94 for SAA/Wasserstei…",
      "lastProgressAt": 1783946706049,
      "tokens": 25106,
      "toolCalls": 7,
      "durationMs": 67091,
      "resultPreview": "{\"claim_id\":\"tails-robust\",\"verdict\":\"SUPPORTED\",\"evidence\":\"Source: /home/dereklong/local/scratch/gaussian-supremum-validator/results/analysis/sat_tails.csv (90 rows; methods NGS,UG,NV,UNGS,benchmark; df values {3,4,6,10,30} spanning the claimed {3..30}; 4 formulations: paper_dro_wasserstein, paper_ro_ellipsoid, paper_saa, paper_so).\\n\\nPart 1 (NGS >=0.96 across all df and formulations, robust ne…"
    },
    {
      "type": "workflow_agent",
      "index": 4,
      "label": "verify:corr-invariant",
      "phaseIndex": 1,
      "phaseTitle": "Verify",
      "agentId": "a58a930e74110efd9",
      "agentType": "general-purpose",
      "model": "claude-opus-4-8[1m]",
      "state": "done",
      "startedAt": 1783946638959,
      "queuedAt": 1783946638943,
      "attempt": 1,
      "lastToolName": "StructuredOutput",
      "lastToolSummary": "corr-invariant",
      "promptPreview": "Adversarially fact-check a claim against saturation-experiment CSV output.\nRepo root: /home/dereklong/local/scratch/gaussian-supremum-validator. CLAIM (corr-invariant): Coverage is essentially invariant to Gaussian coordinate correlation: NGS coverage stays ~0.98-1.0 as corr goes 0.0->0.8 across formulations (d=10).\nRead: results/analysis/sat_corr.csv (method NGS, d=10; column corr, coverage). Use…",
      "lastProgressAt": 1783946706006,
      "tokens": 21668,
      "toolCalls": 6,
      "durationMs": 67047,
      "resultPreview": "{\"claim_id\":\"corr-invariant\",\"verdict\":\"PARTIAL\",\"evidence\":\"From results/analysis/sat_corr.csv (method=NGS, d=10), coverage by corr (0.0,0.2,0.4,0.6,0.8): paper_dro_wasserstein = 0.968, 0.972, 0.954, 0.974, 0.960 (range 0.954-0.974, spread 0.020); paper_ro_ellipsoid = 0.996, 1.0, 1.0, 0.998, 0.994 (range 0.994-1.0); paper_saa = 0.992, 0.980, 0.984, 0.986, 0.978 (range 0.978-0.992); paper_so = 0.9…"
    },
    {
      "type": "workflow_agent",
      "index": 5,
      "label": "verify:foldsdeep-bs-monotone",
      "phaseIndex": 1,
      "phaseTitle": "Verify",
      "agentId": "a33dfd88a3fd8ed98",
      "agentType": "general-purpose",
      "model": "claude-opus-4-8[1m]",
      "state": "done",
      "startedAt": 1783946638959,
      "queuedAt": 1783946638943,
      "attempt": 1,
      "lastToolName": "StructuredOutput",
      "lastToolSummary": "foldsdeep-bs-monotone",
      "promptPreview": "Adversarially fact-check a claim against saturation-experiment CSV output.\nRepo root: /home/dereklong/local/scratch/gaussian-supremum-validator. CLAIM (foldsdeep-bs-monotone): Bootstrap coverage increases monotonically with the resample count K for SO (n=200: ~0.72,0.82,0.90,0.96,0.98 at K=2,3,5,10,20); cross-validation is comparatively flat and never reaches 0.95 for SO at any K.\nRead: results/an…",
      "lastProgressAt": 1783946741669,
      "tokens": 31308,
      "toolCalls": 6,
      "durationMs": 102710,
      "resultPreview": "{\"claim_id\":\"foldsdeep-bs-monotone\",\"verdict\":\"PARTIAL\",\"evidence\":\"From results/analysis/sat_foldsdeep.csv, config=paper_so, n=200 (K parsed from method name). Bootstrap (BS) coverage is strictly monotone increasing in K and matches the quoted values within ~0.003: BS2/K=2=0.7167 (~0.72), BS3/K=3=0.8167 (~0.82), BS5/K=5=0.8967 (~0.90), BS10/K=10=0.9633 (~0.96), BS20/K=20=0.9833 (~0.98). Cross-val…"
    },
    {
      "type": "workflow_agent",
      "index": 6,
      "label": "verify:bigd-RO-dimfree-200",
      "phaseIndex": 1,
      "phaseTitle": "Verify",
      "agentId": "ab8b9b73e48e8c256",
      "agentType": "general-purpose",
      "model": "claude-opus-4-8[1m]",
      "state": "done",
      "startedAt": 1783946638959,
      "queuedAt": 1783946638943,
      "attempt": 1,
      "lastToolName": "StructuredOutput",
      "lastToolSummary": "bigd-RO-dimfree-200",
      "promptPreview": "Adversarially fact-check a claim against saturation-experiment CSV output.\nRepo root: /home/dereklong/local/scratch/gaussian-supremum-validator. CLAIM (bigd-RO-dimfree-200): At extreme dimension the RO validator is still dimension-free: NGS coverage = 1.0 at d=200 for all tested N; SO reaches >=0.99 at d=100 and d=200 once N>=500.\nRead: results/analysis/sat_bigd.csv (method NGS; configs paper_ro_e…",
      "lastProgressAt": 1783946706028,
      "tokens": 23037,
      "toolCalls": 9,
      "durationMs": 67069,
      "resultPreview": "{\"claim_id\":\"bigd-RO-dimfree-200\",\"verdict\":\"SUPPORTED\",\"evidence\":\"Source: /home/dereklong/local/scratch/gaussian-supremum-validator/results/analysis/sat_bigd.csv (method=NGS). Tested d={100,200}, n={200,500,1000,2000}.\\n\\nPart 1 (NGS coverage = 1.0 at d=200 for all tested N; config paper_ro_ellipsoid): d=200,n=200 -> 1.0; n=500 -> 1.0; n=1000 -> 1.0; n=2000 -> 1.0. All exactly 1.0. Confirmed.\\n\\…"
    },
    {
      "type": "workflow_agent",
      "index": 7,
      "label": "verify:budgetfull-bands",
      "phaseIndex": 1,
      "phaseTitle": "Verify",
      "agentId": "aba51f5311fd20457",
      "agentType": "general-purpose",
      "model": "claude-opus-4-8[1m]",
      "state": "done",
      "startedAt": 1783946638960,
      "queuedAt": 1783946638943,
      "attempt": 1,
      "lastToolName": "StructuredOutput",
      "lastToolSummary": "budgetfull-bands",
      "promptPreview": "Adversarially fact-check a claim against saturation-experiment CSV output.\nRepo root: /home/dereklong/local/scratch/gaussian-supremum-validator. CLAIM (budgetfull-bands): Extended best-n1: RO at d=100 has a full feasible Phase-1 band (10%-90% all meet target); SO at d=100 excludes the smallest fraction (10% fails, 20%-90% meet); the Wasserstein d=20 feasible band widens as N grows (n=100 -> up to …",
      "lastProgressAt": 1783946707545,
      "tokens": 24816,
      "toolCalls": 10,
      "durationMs": 68585,
      "resultPreview": "{\"claim_id\":\"budgetfull-bands\",\"verdict\":\"SUPPORTED\",\"evidence\":\"Source: /home/dereklong/local/scratch/gaussian-supremum-validator/results/analysis/sat_budgetfull.csv, method=NGS, target=0.95, columns split (Phase-1 fraction) and coverage.\\n\\n(1) RO full band — config paper_ro_ellipsoid d=100: for BOTH n=500 and n=1000, splits 0.1 through 0.9 all have coverage=1.0 (meets=True). Full feasible band …"
    },
    {
      "type": "workflow_agent",
      "index": 8,
      "label": "synthesize",
      "phaseIndex": 2,
      "phaseTitle": "Synthesize",
      "agentId": "aead50b4e9b8f872b",
      "model": "claude-opus-4-8[1m]",
      "state": "done",
      "startedAt": 1783946741671,
      "queuedAt": 1783946741669,
      "attempt": 1,
      "promptPreview": "Adversarial verification verdicts on the saturation-campaign claims:\n- [PARTIAL] alpha-tight-collapse: Read results/analysis/sat_alpha.csv, method NGS. The direction of the collapse claim holds and Part-1 values match within tolerance: at alpha=0.01, coverage collapses for paper_so (0.006) and paper_saa (0.006) — exactly the quoted ~0.006 — while paper_ro_ellipsoid stays high at 0.984 (matches ~0.…",
      "lastProgressAt": 1783946761408,
      "tokens": 20832,
      "toolCalls": 0,
      "durationMs": 19737,
      "resultPreview": "# Saturation-Campaign Faithfulness Report\n\n## Verdict Counts (7 claims)\n- SUPPORTED: 4\n- PARTIAL: 3\n- REFUTED: 0\n\n## Claims Needing Correction (3)\n\n**1. alpha-tight-collapse [PARTIAL]**\n- What holds: The tight-tolerance collapse is real. At alpha=0.01, SO and SAA collapse to ~0.006 while RO_ellipsoid stays at 0.984 — direction and Part-1 numbers match.\n- Overstatement: The \"by alpha>=0.05 all form…"
    }
  ],
  "totalTokens": 195233,
  "totalToolCalls": 49
}
```

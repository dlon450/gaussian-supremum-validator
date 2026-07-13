export const meta = {
  name: 'verify-saturation-findings',
  description: 'Adversarially verify the saturation-campaign claims against sat_*.csv',
  phases: [{ title: 'Verify' }, { title: 'Synthesize' }],
}

const REPO = '/home/dereklong/local/scratch/gaussian-supremum-validator'
const A = 'results/analysis'
const CLAIMS = [
  { id: 'alpha-tight-collapse', claim: 'At very tight tolerance (alpha=0.01, i.e. 1-alpha=99%) coverage collapses for SO and SAA (~0.006) while RO stays high (~0.98); by alpha>=0.05 all formulations meet the 0.95 target (NGS).',
    files: `${A}/sat_alpha.csv (method NGS; configs paper_so, paper_saa, paper_ro_ellipsoid; column alpha, coverage)` },
  { id: 'beta-calibration', claim: 'Calibration: as nominal target 1-beta decreases (beta rises 0.01->0.50), empirical coverage of NGS/UNGS stays at or above the nominal target (valid/conservative), while UG is the least conservative (closest to nominal, dipping lowest at large beta). E.g. SO beta=0.50: NGS~0.96, UNGS~0.95, UG~0.74 (nominal 0.50).',
    files: `${A}/sat_beta.csv (config paper_so; methods UG,NGS,UNGS; columns beta, coverage; nominal target = 1-beta)` },
  { id: 'tails-robust', claim: 'Under heavy tails (multivariate-t), NGS coverage stays >=0.96 across df in {3..30} for all formulations, i.e. robust even near the finite-variance boundary df=3; the univariate UG is less robust (dips to ~0.90-0.94 for SAA/Wasserstein).',
    files: `${A}/sat_tails.csv (methods NGS and UG; column df, coverage)` },
  { id: 'corr-invariant', claim: 'Coverage is essentially invariant to Gaussian coordinate correlation: NGS coverage stays ~0.98-1.0 as corr goes 0.0->0.8 across formulations (d=10).',
    files: `${A}/sat_corr.csv (method NGS, d=10; column corr, coverage)` },
  { id: 'foldsdeep-bs-monotone', claim: 'Bootstrap coverage increases monotonically with the resample count K for SO (n=200: ~0.72,0.82,0.90,0.96,0.98 at K=2,3,5,10,20); cross-validation is comparatively flat and never reaches 0.95 for SO at any K.',
    files: `${A}/sat_foldsdeep.csv (config paper_so; methods BS2..BS20 and CV2..CV20; parse K from method name)` },
  { id: 'bigd-RO-dimfree-200', claim: 'At extreme dimension the RO validator is still dimension-free: NGS coverage = 1.0 at d=200 for all tested N; SO reaches >=0.99 at d=100 and d=200 once N>=500.',
    files: `${A}/sat_bigd.csv (method NGS; configs paper_ro_ellipsoid and paper_so; columns d, n, coverage)` },
  { id: 'budgetfull-bands', claim: 'Extended best-n1: RO at d=100 has a full feasible Phase-1 band (10%-90% all meet target); SO at d=100 excludes the smallest fraction (10% fails, 20%-90% meet); the Wasserstein d=20 feasible band widens as N grows (n=100 -> up to 50%, n=200 -> up to 70%). (method NGS)',
    files: `${A}/sat_budgetfull.csv (method NGS; configs paper_ro_ellipsoid/paper_so d=100, paper_dro_wasserstein d=20; columns split, coverage)` },
]

const SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['claim_id', 'verdict', 'evidence'],
  properties: {
    claim_id: { type: 'string' },
    verdict: { type: 'string', enum: ['SUPPORTED', 'PARTIAL', 'REFUTED'] },
    evidence: { type: 'string' },
    correction: { type: 'string' },
  },
}

phase('Verify')
const results = await parallel(CLAIMS.map(c => () =>
  agent(
    `Adversarially fact-check a claim against saturation-experiment CSV output.\n` +
    `Repo root: ${REPO}. CLAIM (${c.id}): ${c.claim}\n` +
    `Read: ${c.files}. Use Bash (python/pandas) to inspect the numbers. Try to REFUTE. ` +
    `SUPPORTED only if the data clearly backs the direction and quoted values match within ~0.03 (coverage). ` +
    `PARTIAL if direction holds but a number is off or holds only partially; REFUTED if wrong/missing. ` +
    `Quote the actual numbers you read in 'evidence'.`,
    { label: `verify:${c.id}`, phase: 'Verify', schema: SCHEMA, agentType: 'general-purpose' }
  ).catch(() => ({ claim_id: c.id, verdict: 'REFUTED', evidence: 'verifier error', correction: 'could not verify' }))
))

phase('Synthesize')
const table = results.filter(Boolean).map(r => `- [${r.verdict}] ${r.claim_id}: ${r.evidence}${r.correction ? ' | FIX: ' + r.correction : ''}`).join('\n')
const report = await agent(
  `Adversarial verification verdicts on the saturation-campaign claims:\n${table}\n\n` +
  `Write a short report: counts SUPPORTED/PARTIAL/REFUTED, list any needing correction with the fix, ` +
  `and an overall faithfulness verdict. Be blunt about any overstatement.`,
  { label: 'synthesize', phase: 'Synthesize' })
return { verdicts: results.filter(Boolean), report }

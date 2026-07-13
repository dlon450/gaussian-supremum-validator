export const meta = {
  name: 'verify-comprehensive-findings',
  description: 'Adversarially verify the comprehensive-campaign claims against comp_*.csv',
  phases: [{ title: 'Verify' }, { title: 'Synthesize' }],
}

const REPO = '/home/dereklong/local/scratch/gaussian-supremum-validator'
const A = 'results/analysis'
const CLAIMS = [
  { id: 'budget-n1-decreases-with-N', claim: 'The recommended (best feasible) Phase-1 fraction DECREASES as N grows: for SO d=10, ~0.5 at N=200, ~0.2 at N=500, ~0.1 at N=1000 (method NGS).',
    files: `${A}/comp_budget_best_n1.csv (formulation SO, method NGS, d=10)` },
  { id: 'budget-n1-increases-with-D', claim: 'The recommended Phase-1 fraction INCREASES with dimension d at fixed N: for SO at N=500, ~0.1 (d=2), ~0.2 (d=10), ~0.3 (d=50) (method NGS).',
    files: `${A}/comp_budget_best_n1.csv (formulation SO, method NGS, n=500)` },
  { id: 'budget-band-exists', claim: 'There is a feasible band of n1: for SO n=500 d=10, a too-small fraction (0.1) is below target while a mid band (~0.2-0.7) meets 0.95 (any of NGS/UNGS).',
    files: `${A}/comp_budget.csv (config paper_so, n=500, d=10) — coverage by split` },
  { id: 'ndgrid-RO-dimension-free', claim: 'RO is essentially dimension-free: NGS coverage >= 0.95 across the whole (N,D) grid, including d=50 at N=100 (>=0.96).',
    files: `${A}/comp_ndgrid.csv (config paper_ro_ellipsoid, method NGS)` },
  { id: 'ndgrid-SO-SAA-need-N-with-D', claim: 'SO and SAA need N to scale with D: at small N large D coverage is far below target (SO d=50 N=100 ~0.50; SAA d=20 N=100 ~0.49) but recovers to >=0.98 at larger N (SO d=50 N=500 ~0.99).',
    files: `${A}/comp_ndgrid.csv (configs paper_so and paper_saa, method NGS)` },
  { id: 'folds-more-folds-higher-coverage', claim: 'For the existing schemes, more folds/resamples give higher coverage; bootstrap on SO at N=500 rises 0.75 -> 0.86 -> 0.92 for K=3 -> 5 -> 10, and K=10 is the best of the tested counts.',
    files: `${A}/comp_folds.csv (config paper_so, methods BS3,BS5,BS10)` },
  { id: 'folds-CV-undercovers-SO', claim: 'Cross-validation systematically under-covers for SO: CV never reaches 0.95 at any fold count or N (max ~0.90 at CV10).',
    files: `${A}/comp_folds.csv (config paper_so, methods CV3,CV5,CV10)` },
  { id: 'folds-proposed-match-best', claim: 'At matched data budget the proposed validators meet target for SO from N=300 (UG~0.96, NGS~0.98, UNGS~1.0), comparable to the best fold setting BS10/Sec10 and strictly better than CV10.',
    files: `${A}/comp_folds.csv (config paper_so, N=300 and 500)` },
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
    `Adversarially fact-check a claim against comprehensive-experiment CSV output.\n` +
    `Repo root: ${REPO}. CLAIM (${c.id}): ${c.claim}\n` +
    `Read: ${c.files}. Columns include formulation/config, method, n, d, split, coverage, recommended_n1.\n` +
    `Use Bash (python/pandas or awk/grep) to read the numbers. Try to REFUTE. SUPPORTED only if the data ` +
    `clearly backs the direction and quoted values match within ~0.03 (coverage) or one grid-step (n1). ` +
    `PARTIAL if direction holds but a number is off; REFUTED if wrong/missing. Quote the actual numbers in 'evidence'.`,
    { label: `verify:${c.id}`, phase: 'Verify', schema: SCHEMA, agentType: 'general-purpose' }
  ).catch(() => ({ claim_id: c.id, verdict: 'REFUTED', evidence: 'verifier error', correction: 'could not verify' }))
))

phase('Synthesize')
const table = results.filter(Boolean).map(r => `- [${r.verdict}] ${r.claim_id}: ${r.evidence}${r.correction ? ' | FIX: ' + r.correction : ''}`).join('\n')
const report = await agent(
  `Adversarial verification verdicts on the comprehensive-campaign claims:\n${table}\n\n` +
  `Write a short report: counts of SUPPORTED/PARTIAL/REFUTED, list any needing correction with the fix, ` +
  `and an overall faithfulness verdict. Be blunt about overstatements.`,
  { label: 'synthesize', phase: 'Synthesize' })
return { verdicts: results.filter(Boolean), report }

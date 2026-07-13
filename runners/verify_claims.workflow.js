export const meta = {
  name: 'verify-reviewer-evidence',
  description: 'Adversarially verify each reviewer-response claim against the raw result files',
  phases: [
    { title: 'Verify' },
    { title: 'Synthesize' },
  ],
}

// Each claim names the exact numeric assertion in paper/reviewer_response_draft.md
// and the result file(s) that should back it. A verifier agent must READ the file(s)
// and try to REFUTE the claim; default to REFUTED if the data does not clearly support it.
const REPO = '/home/dereklong/local/scratch/gaussian-supremum-validator'
const CLAIMS = [
  { id: 'so-small-n', claim: 'For SO at n=100 (split=0.7), bootstrap coverage (~0.87) and CV (~0.76) EXCEED the proposed UG (~0.65) and NGS (~0.68).',
    files: 'results/analysis/existing_comparison.csv (config paper_so, n=100) or results/experiments/existing/paper_so_n100_s0.7_d10_summary.json' },
  { id: 'so-n300-cv-below', claim: 'For SO at n=300, UG/NGS/UNGS all meet the 0.95 target while CV is BELOW target (~0.89).',
    files: 'results/experiments/existing/paper_so_n300_s0.7_d10_summary.json' },
  { id: 'so-n500-cv-below', claim: 'For SO at n=500, the proposed validators meet the 0.95 target while CV is BELOW target (~0.86).',
    files: 'results/experiments/existing/paper_so_n500_s0.7_d10_summary.json' },
  { id: 'saa-cv-earlier', claim: 'For SAA, CV and bootstrap reach the 0.95 target at a SMALLER n than the proposed validators (e.g. at n=200 CV~0.96, BS~0.99 vs UG~0.82).',
    files: 'results/experiments/existing/paper_saa_n200_s0.7_d10_summary.json and paper_saa_n300_s0.7_d10_summary.json' },
  { id: 'n1-too-small-fails', claim: 'For SO n=500, a too-small Phase-1 fraction (10%) fails the target (coverage ~0.72), while a broad band (~20-90%) meets it for NGS.',
    files: 'results/analysis/split_budgeting.csv (config paper_so)' },
  { id: 'n1-too-large-fails', claim: 'For SO n=500, the univariate Gaussian (UG) validator degrades at a too-large Phase-1 fraction (90% -> ~0.93, below target), showing a two-sided sweet spot.',
    files: 'results/analysis/split_budgeting.csv (config paper_so)' },
  { id: 'dimension-free', claim: 'RO coverage of UG/NGS/UNGS stays >= ~0.99 across dimensions d=2..100, while the naive plain-average (NV) drops below target (e.g. ~0.87 at d=10-20).',
    files: 'results/analysis/dim_free.csv' },
  { id: 'sca-advantage-grows', claim: 'The UG objective advantage over the SCA benchmark GROWS with dimension d (roughly 0.15 at d=2 up to ~1.5 at d=50).',
    files: 'results/analysis/dim_free.csv (mean_obj of UG vs benchmark)' },
  { id: 'convergence-vs-n', claim: 'In the nsweep matrix, proposed-validator coverage rises to/holds the target and the mean oracle gap shrinks toward 0 as n grows.',
    files: 'results/analysis/nsweep.csv' },
  { id: 'robust-tails', claim: 'Under heavy-tailed multivariate-t data (robust matrix), the proposed validators retain approximately target-level coverage.',
    files: 'results/analysis/robust.csv' },
]

const SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['claim_id', 'verdict', 'evidence'],
  properties: {
    claim_id: { type: 'string' },
    verdict: { type: 'string', enum: ['SUPPORTED', 'PARTIAL', 'REFUTED'] },
    evidence: { type: 'string', description: 'Actual numbers read from the file(s) that justify the verdict' },
    correction: { type: 'string', description: 'If PARTIAL/REFUTED, the corrected statement; else empty' },
  },
}

phase('Verify')
const results = await parallel(CLAIMS.map(c => () =>
  agent(
    `You are adversarially fact-checking a claim in a paper's reviewer response against raw experiment output.\n` +
    `Repo root: ${REPO}\n` +
    `CLAIM (${c.id}): ${c.claim}\n` +
    `Read this/these file(s) and inspect the actual numbers: ${c.files}\n` +
    `Use Bash (cat/python) or Read. The result JSONs have summaries.<METHOD>.coverage and .mean_obj. ` +
    `CSVs have columns method,coverage,mean_obj,n,split,d.\n` +
    `Try to REFUTE the claim. Return SUPPORTED only if the data clearly backs the stated direction AND the quoted numbers are within ~0.03 (coverage) / ~15% (objective). ` +
    `Return PARTIAL if the direction holds but a quoted number is off; REFUTED if the direction is wrong or data is missing. ` +
    `Put the actual numbers you read in 'evidence'.`,
    { label: `verify:${c.id}`, phase: 'Verify', schema: SCHEMA, agentType: 'general-purpose' }
  ).catch(() => ({ claim_id: c.id, verdict: 'REFUTED', evidence: 'verifier error', correction: 'could not verify' }))
))

phase('Synthesize')
const table = results.filter(Boolean).map(r => `- [${r.verdict}] ${r.claim_id}: ${r.evidence}${r.correction ? ' | FIX: ' + r.correction : ''}`).join('\n')
const summary = await agent(
  `These are adversarial verification verdicts on the reviewer-response claims for the Gaussian-supremum-validator paper.\n${table}\n\n` +
  `Write a concise verification report: (1) how many SUPPORTED / PARTIAL / REFUTED; (2) list any PARTIAL/REFUTED with the correction needed; ` +
  `(3) an overall statement on whether the reviewer-response draft is faithful to the data. Be blunt about any overstatement.`,
  { label: 'synthesize', phase: 'Synthesize' }
)
return { verdicts: results.filter(Boolean), report: summary }

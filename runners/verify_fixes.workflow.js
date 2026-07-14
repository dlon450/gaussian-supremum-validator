export const meta = {
  name: 'verify-audit-fixes',
  description: 'Adversarially verify the audit fixes against the actual corrected files/data',
  phases: [{ title: 'Verify' }, { title: 'Synthesize' }],
}

const REPO = '/home/dereklong/local/scratch/gaussian-supremum-validator'
const CLAIMS = [
  { id: 'portfolio-tails-fixed', claim: 'Portfolio heavy-tail is fixed: the four results/experiments/portfolio_tails/*_raw.jsonl files are NO LONGER content-identical, and coverage now varies with df (df=3 lower than df=10 for the proposed validators).',
    how: `Compare the 4 raw files (sort each and diff, or compare summary coverage) under ${REPO}/results/experiments/portfolio_tails/. Summaries have summaries.{UG,NGS,UNGS}.coverage and a top-level "df".` },
  { id: 'moment-radius-fixed', claim: 'Moment-DRO radius now includes 1/sqrt(n): gsv/paths.py _moment_chi2_radius returns sqrt(chi2.ppf(0.95,q)/n) (matches legacy DRO2.m), the path includes s=0, and moment-DRO is a degeneracy (all validators incl. NV are ~1.0 coverage).',
    how: `Read ${REPO}/gsv/paths.py (_moment_chi2_radius and build_path dro_moment). Check moment summaries under ${REPO}/results/experiments/moment/ have NV coverage ~1.0 too (degeneracy).` },
  { id: 'failures-first-class', claim: 'Solver failures are first-class: gsv/experiment.py records a benchmark failure row when the benchmark is non-finite, masks NaN candidate columns, and failure_rate = failures/(successes+failures) with n_attempted reported.',
    how: `Read ${REPO}/gsv/experiment.py (run_replication benchmark else-branch + NaN masking; run_cell failure_rate). Confirm the code paths exist.` },
  { id: 'json-valid', claim: 'All result files are strict JSON: none of results/experiments/**/*.json or *.jsonl contain the tokens NaN or Infinity.',
    how: `grep -rl "NaN\\|Infinity" ${REPO}/results/experiments/ (should be empty). Also confirm a sampled summary parses with a strict JSON parser.` },
  { id: 'provenance-present', claim: 'Newly regenerated result files carry a _provenance block with a git_commit and package versions.',
    how: `Read a moment or portfolio summary JSON under ${REPO}/results/experiments/ and check for _provenance.git_commit and _provenance.versions.` },
  { id: 'pytest-passes', claim: 'pytest collects and passes the infra tests (the decorator/gurobipy collection failure is fixed).',
    how: `Run: cd ${REPO} && OMP_NUM_THREADS=1 /tmp/gsv_venv/bin/python -m pytest -q tests/test_infra.py  — expect "10 passed".` },
  { id: 'independent-confirmation', claim: 'Recommended splits hold out-of-sample: in results/analysis/robustness_checks.json, the NGS/UNGS out-of-sample coverage (out_cov) is >= ~0.95 for SO/RO/SAA (not selection optimism).',
    how: `Read ${REPO}/results/analysis/robustness_checks.json ("confirmation" list: in_cov vs out_cov per config/method).` },
  { id: 'sensitivity-stable', claim: 'Coverage is stable to mesh size p, GS sim_num, and the 5s MILP cap (the cap gives identical results at 5s vs 30s).',
    how: `Read ${REPO}/results/analysis/robustness_checks.md (mesh_p, sim_num sections) and ${REPO}/results/analysis/milp_cap_sensitivity.txt (5s vs 30s identical).` },
]

const SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['claim_id', 'verdict', 'evidence'],
  properties: {
    claim_id: { type: 'string' }, verdict: { type: 'string', enum: ['SUPPORTED', 'PARTIAL', 'REFUTED'] },
    evidence: { type: 'string' }, correction: { type: 'string' },
  },
}

phase('Verify')
const results = await parallel(CLAIMS.map(c => () =>
  agent(
    `Adversarially verify an audit-fix claim against the ACTUAL repo files/data. Try to REFUTE it.\n` +
    `CLAIM (${c.id}): ${c.claim}\nHOW TO CHECK: ${c.how}\n` +
    `Use Bash (grep/python/pytest/cat) and Read. SUPPORTED only if the files/data clearly confirm it; ` +
    `PARTIAL if partially; REFUTED if not. Quote the concrete evidence (numbers, code lines, command output).`,
    { label: `verify:${c.id}`, phase: 'Verify', schema: SCHEMA, agentType: 'general-purpose' }
  ).catch(() => ({ claim_id: c.id, verdict: 'REFUTED', evidence: 'verifier error', correction: 'could not verify' }))
))

phase('Synthesize')
const table = results.filter(Boolean).map(r => `- [${r.verdict}] ${r.claim_id}: ${r.evidence}${r.correction ? ' | FIX: ' + r.correction : ''}`).join('\n')
const report = await agent(
  `Adversarial verification of the audit fixes:\n${table}\n\nWrite a short report: counts SUPPORTED/PARTIAL/REFUTED, ` +
  `list anything still not fully fixed, and an overall verdict on whether the two invalid families are corrected and the ` +
  `traceability gaps closed. Be blunt.`, { label: 'synthesize', phase: 'Synthesize' })
return { verdicts: results.filter(Boolean), report }

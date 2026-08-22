export const meta = {
  name: 'verify-phase1-reports',
  description: 'Adversarially verify each Phase 1 scE2G characterisation report before Phase 2 consumes it',
  whenToUse: 'After Phase 1 characterisation agents have written phase1_reports/*.md, to separate honest "data cannot constrain this" findings from overreach, hallucinated citations, and unevaluable cost functions.',
  phases: [
    { title: 'Verify', detail: 'one skeptic per report: citations, overreach, function evaluability' },
    { title: 'Adjudicate', detail: 'second independent skeptic on reports flagged FAIL or SUSPECT' },
    { title: 'Synthesise', detail: 'cross-report consistency and a Phase 2 readiness verdict' },
  ],
}

const ROOT = '/oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G_optimize'
const ANALYSIS = '/scratch/users/kaybrand/scE2G_optimize_results/analysis'

// args: array of report filenames relative to phase1_reports/, e.g. ["07-make_kendall_pairs.md", ...]
const reports = Array.isArray(args) ? args : []
if (!reports.length) {
  log('No reports passed in args — nothing to verify.')
  return { error: 'empty report list' }
}
log(`Verifying ${reports.length} Phase 1 reports.`)

const VERDICT_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['report', 'verdict', 'honest_unconstrained', 'class_derivation_quality',
             'citation_errors', 'overreach_findings', 'class_authority_problems',
             'function_problems', 'variable_problems', 'benchmark_mismatches', 'notes'],
  properties: {
    report: { type: 'string' },
    verdict: {
      type: 'string',
      enum: ['SOUND', 'SOUND_WITH_CAVEATS', 'SUSPECT', 'FAIL'],
      description: 'SOUND = usable by Phase 2 as-is. SOUND_WITH_CAVEATS = usable, minor corrections listed. SUSPECT = a central claim could not be verified. FAIL = a load-bearing claim is wrong, or the cost functions are unusable.',
    },
    honest_unconstrained: {
      type: 'boolean',
      description: 'True if the report CORRECTLY states that the benchmarks cannot constrain its exponent. This is a virtue, not a defect — do not penalise it.',
    },
    citation_errors: {
      type: 'array', items: { type: 'string' },
      description: 'Each file:line cite that does not say what the report claims. Quote the actual line.',
    },
    overreach_findings: {
      type: 'array', items: { type: 'string' },
      description: 'Claims of empirical support that the data cannot bear — e.g. an exponent fitted to a variable spanning <1.2x, or a claim to have separated two correlated variables from 6 points.',
    },
    class_authority_problems: {
      type: 'array', items: { type: 'string' },
      description: 'The INVERTED error: places where the report let the benchmark timings determine or bound the complexity class instead of deriving it from the code. Includes a t_low/t_high pair whose EXPONENTS differ without a stated ambiguity in the code to justify it; a linear t_low chosen because the measurements look flat when the code shows a sort or nested loop; treating non-monotonic timings as evidence against a code-derived class rather than as a disagreement to explain; or omitting the code-derived class entirely in favour of a fitted curve.',
    },
    class_derivation_quality: {
      type: 'string',
      enum: ['FROM_CODE', 'PARTLY_FROM_CODE', 'FROM_BENCHMARKS', 'ABSENT'],
      description: 'Where the reported complexity class actually came from. FROM_CODE requires a cited sort / hash join / nested loop / interval tree that pins the exponent. FROM_BENCHMARKS means the class was inferred from timings and is not acceptable.',
    },
    function_problems: {
      type: 'array', items: { type: 'string' },
      description: 'Anything that stops Phase 2 evaluating t_est/t_low/t_high mechanically: unstated constants, undefined variables, wrong units, t_low > t_est, arithmetic that does not reproduce the reported values.',
    },
    variable_problems: {
      type: 'array', items: { type: 'string' },
      description: 'Named size variables that are not real size_manifest columns; use of the artifact columns n_frag_bytes / n_umi_bytes / n_cells_bytes; multiome values used for scATAC or vice versa.',
    },
    benchmark_mismatches: {
      type: 'array', items: { type: 'string' },
      description: 'Quoted measured numbers that disagree with the actual benchmark rows.',
    },
    notes: { type: 'string', description: 'Anything Phase 2 must know, in a few sentences.' },
  },
}

const sharedContext = `
Repo root: ${ROOT}
Analysis tables: ${ANALYSIS}/benchmarks_{multiome,scatac}.tsv and size_manifest_{multiome,scatac}.tsv

Ground truth you must read before judging:
  ${ROOT}/phase0_measurements/PHASE1_BRIEFING.md
  ${ROOT}/phase0_measurements/PHASE1_ADDENDUM.md

THE CENTRAL JUDGEMENT YOU ARE HERE TO MAKE. Across the six IGVF10 clusters only n_frag (15.1x),
n_umi (12.0x) and n_cells (6.0x) vary meaningfully, and those three are mutually correlated, so
even they cannot be separated from six points. Everything downstream of candidate-region calling
spans only 1.01x-1.17x. Therefore:

  * A report that says "class derived from code; the benchmarks cannot confirm it because this
    variable spans only 1.01x" is CORRECT and should be marked SOUND. Set honest_unconstrained
    true. Do NOT treat this as a failure to deliver.
  * A report that presents a fitted exponent, an R-squared, or a confident scaling claim for a
    near-constant variable is OVERREACH, even if it sounds rigorous. That is the specific failure
    mode this verify pass exists to catch.
  * Distinguish both of those from an agent that simply did not do the work.

THE SECOND JUDGEMENT, AND IT IS EQUALLY IMPORTANT — THE INVERTED ERROR.

The CODE is the authority on complexity class. The benchmarks are background context on OVERHEAD
(startup, I/O volume, memory, realised threading, and the numeric constant), and are NOT evidence
about the exponent in EITHER direction. O(n log n) versus O(n^2) is settled by finding the sort,
the hash join, the nested loop or the interval tree in the source — and by nothing else.

So as well as catching reports that fit an exponent to noise, catch reports that do the reverse:

  * A report that uses flat or non-monotonic timings to ARGUE DOWN an exponent it found in the
    code. If the code contains a nested loop over j_cand_elems, the class is O(j^2) even if the
    six points look linear. The report should state the code-derived class, then report the
    disagreement AS a disagreement with a hypothesised mechanism (a dominant constant at this
    scale, caching, I/O dominance, an inner loop over a bounded subset rather than all of j).
  * A t_low / t_high pair whose EXPONENTS differ, with no stated ambiguity in the CODE to justify
    the spread. t_low and t_high bound the class you derived from the code — they are not the
    envelope of a scatter plot. If the code pins the class, t_low and t_high should differ in
    their CONSTANT, not their exponent. Differing exponents are legitimate ONLY where the code
    itself leaves the class genuinely ambiguous (e.g. a pandas merge whose cost depends on key
    uniqueness that cannot be determined statically) — and the report must say which code
    ambiguity drives it.
  * A report whose stated class has no cited code construct behind it at all, only a curve.

Record these in class_authority_problems and set class_derivation_quality accordingly. Note the
asymmetry when you judge severity: a report that honestly says "code says O(n log n), the six
points cannot confirm it" is SOUND. A report that says "the timings look linear so I model it as
linear" when the code sorts is NOT sound, however well-calibrated its numbers are, because Phase 2
will extrapolate it to other cluster sizes.

Useful fact for calibrating this: over the full range of the most variable quantity in the dataset,
n_frag = 26.2M to 394.8M (15.1x), the factor ln(n_frag) itself spans only 1.16x. The log term is
nearly constant across every measurement available, so no dataset of this shape can separate O(n)
from O(n log n). That separation must come from the code, where it is available for free.

Also check:
  * Every file:line citation. Open the file at that line. If it does not say what the report
    claims, that is a citation_error — quote the actual line. Hallucinated cites are the classic
    failure of this kind of report.
  * That t_est / t_low / t_high are mechanically evaluable: every constant stated numerically,
    every variable a real size_manifest column name, units in seconds, and t_low <= t_est <=
    t_high over the range of interest. Recompute the report's own mean/min/max evaluations and
    confirm the arithmetic. Phase 2 will evaluate these functions automatically, so a function
    that cannot be evaluated is a FAIL regardless of how good the prose is.
  * That quoted measured numbers match the real benchmark rows:
      awk -F'\\t' -v r=RULE 'NR==1||$2==r' ${ANALYSIS}/benchmarks_multiome.tsv
  * n_frag_bytes / n_umi_bytes / n_cells_bytes are ARTIFACTS (the digit count of a count file),
    not data volumes. Using them as sizes is a variable_problem.
  * io_in / io_out claims. Verified against snakemake/benchmark.py in the run_snakemake9 env:
    io_in/io_out are summed over the main process and its STILL-ALIVE children each poll and then
    plain-ASSIGNED to the record (:404-405), unlike max_rss which takes a max (:399). So the value
    is a LAST-POLL SNAPSHOT and any child that already exited contributes nothing. The underlying
    read_bytes/write_bytes are BLOCK-DEVICE counters, so page-cache hits do not count either.
    Consequences you must check for:
      - A report treating io_in/io_out as total logical read/write volume, or comparing it directly
        against an on-disk file size as if they were the same quantity -> variable_problem.
      - Specifically: any report claiming io_out measures the UNCOMPRESSED byte volume written into
        a compression stream. That is FALSE — pipe traffic and in-process gzip never reach the block
        device. If a report asserts this, flag it. (Unit 22 asserted it; treat as unverified.)
      - Conversely, do NOT mark a report down for reporting a large io_out as real: write_bytes
        counts genuine device writes, so a rule that writes large uncompressed temp files really is
        doing that I/O (unit 9's per-chromosome BED splits are the likely genuine case).
  * n_called_all / n_called_distal / n_feat_cols differ substantially between modalities
    (see the addendum). Using the multiome figures for scATAC is a variable_problem.

You are a skeptic. Default to doubting the report and let it earn SOUND. But be fair: do not
manufacture problems, and do not mark a report down for honestly reporting that the measurement
is uninformative. Report edits are NOT your job — you only judge.
`

const verified = await pipeline(
  reports,
  (rep) => agent(
    `${sharedContext}

YOUR TASK: verify the Phase 1 report at ${ROOT}/phase1_reports/${rep}

Read the report, then independently check it against the code and the benchmark tables. Verify
every file:line citation by opening the file. Recompute the cost-function evaluations yourself.
Then return your structured verdict.`,
    { label: `verify:${rep}`, phase: 'Verify', schema: VERDICT_SCHEMA }
  ),
  // Second, independent skeptic only where the first flagged trouble. Cheap reports pass straight through.
  (v, rep) => {
    if (!v) return null
    if (v.verdict === 'SOUND' || v.verdict === 'SOUND_WITH_CAVEATS') return v
    const problems = JSON.stringify({
      citation_errors: v.citation_errors,
      overreach_findings: v.overreach_findings,
      class_authority_problems: v.class_authority_problems,
      function_problems: v.function_problems,
      variable_problems: v.variable_problems,
      benchmark_mismatches: v.benchmark_mismatches,
    }, null, 2)

    // Framing is deliberately asymmetric by severity.
    //  SUSPECT -> refute-framing. The risk is a reviewer over-flagging an honest
    //             "data cannot constrain this", so pressure-test the complaint.
    //  FAIL    -> NEUTRAL framing. Refute-framing here would let a real, load-bearing
    //             defect get talked out of existence by an adjudicator whose job was
    //             defined as finding the first reviewer wrong.
    const task = v.verdict === 'SUSPECT'
      ? `YOUR TASK: independently re-examine, and try to REFUTE the first reviewer. Check whether
each claimed problem is real by going to the code and the data yourself. Reviewers sometimes mark
a report down for correctly saying "the data cannot constrain this" — which is the right answer
here — or for misreading a citation. If the complaints do not survive contact with the source,
say so and raise the verdict. Return your own independent verdict, not a summary of theirs.`
      : `YOUR TASK: independently assess whether these load-bearing claims are actually wrong.
Do NOT set out to vindicate or to overturn the first reviewer — determine the truth from the code
and the data directly. A FAIL means something Phase 2 would propagate as a wrong number, so both
error directions are costly: confirming a defect that is not real, and dismissing one that is.
For each claimed problem, go to the cited file and line, or recompute the arithmetic, and state
what you actually found. Then give your own independent verdict, which may agree or disagree.`

    return agent(
      `${sharedContext}

A first reviewer marked ${ROOT}/phase1_reports/${rep} as ${v.verdict}. Their stated problems:

${problems}

${task}`,
      { label: `adjudicate:${rep}`, phase: 'Adjudicate', schema: VERDICT_SCHEMA }
    ).then(second => ({ ...second, first_pass_verdict: v.verdict, adjudicated: true }))
  }
)

// pipeline() returns results positionally aligned with `reports`, with null where an agent
// crashed, timed out, hit a session/rate limit, or was skipped. Those must NOT be silently
// dropped: an unverified report is an unknown, not a pass. Zip against the input to recover
// exactly which ones produced no verdict, and carry that list through to synthesis.
const results = []
const unverified = []
reports.forEach((rep, i) => {
  const v = verified[i]
  if (v && v.verdict) results.push(v)
  else unverified.push(rep)
})

const byVerdict = (v) => results.filter(r => r.verdict === v).map(r => r.report)

log(`SOUND ${byVerdict('SOUND').length} · CAVEATS ${byVerdict('SOUND_WITH_CAVEATS').length} · SUSPECT ${byVerdict('SUSPECT').length} · FAIL ${byVerdict('FAIL').length} · UNVERIFIED ${unverified.length}`)
if (unverified.length) {
  log(`WARNING — ${unverified.length} report(s) produced NO verdict and are unverified: ${unverified.join(', ')}`)
}

phase('Synthesise')
const synthesis = await agent(
  `${sharedContext}

${results.length} of ${reports.length} Phase 1 reports were verified. Here are the structured verdicts:

${JSON.stringify(results, null, 2)}

${unverified.length ? `CRITICAL — ${unverified.length} report(s) produced NO verdict at all, because the
verifying agent crashed, timed out, or hit a session limit:

${unverified.map(r => '  - ' + r).join('\n')}

These are UNVERIFIED, which is not the same as sound and not the same as failed — they are
unknown. Do not let them disappear from your assessment. Give them their own clearly labelled
section, list them by name, and state plainly that Phase 2 must either re-verify them or treat
their cost functions as unvalidated. Count them separately from the verdict tallies, and do not
describe coverage as complete.` : 'Every report received a verdict; none were dropped.'}

YOUR TASK — cross-report consistency, which no single-report reviewer could see. Produce a Phase 2
readiness assessment covering:

1. CONTRADICTIONS: two reports making incompatible claims about the same shared quantity. Watch
   especially for units that share a variable (n_cand_pairs, j_cand_elems, k_genes, n_frag) or
   sit adjacent in the DAG, and for the four rules whose cost differs by modality
   (add_external_features, run_e2g_qnorm, gen_final_features, get_stats_per_model_per_cluster).
2. DOUBLE-COUNTING RISK: places where two units describe the same physical work, so Phase 2 would
   sum the same seconds twice. The grouped Tier C units overlap with single-script units by
   design — bedSplitSort.sh (unit 9) sits inside rule frag_to_tagAlign (unit 24), and
   make_biosample_feature_table (unit 23) feeds checkpoint features_required (unit 26). Say
   exactly which node Phase 2 should attribute each block of time to.
3. COVERAGE GAPS: rules that ran and have benchmark data but that no report characterises, and
   DAG nodes with no cost function at all.
4. THE UNCONSTRAINED SET: list the rules whose cost class is code-derived only, with no empirical
   support. Phase 2 must propagate wide bounds for these, so name them explicitly.
5. WHICH REPORTS PHASE 2 SHOULD NOT TRUST as written, and what specifically to re-derive.

Be concrete and name rules. Write your assessment to
${ROOT}/phase1_reports/00-VERIFICATION_SYNTHESIS.md and also return it.`,
  { label: 'synthesis', phase: 'Synthesise' }
)

return {
  counts: {
    reports_in: reports.length,
    verified: results.length,
    sound: byVerdict('SOUND').length,
    caveats: byVerdict('SOUND_WITH_CAVEATS').length,
    suspect: byVerdict('SUSPECT').length,
    fail: byVerdict('FAIL').length,
    unverified: unverified.length,
  },
  unverified,
  unconstrained_but_honest: results.filter(r => r.honest_unconstrained).map(r => r.report),
  needs_attention: results.filter(r => r.verdict === 'SUSPECT' || r.verdict === 'FAIL'),
  synthesis,
}

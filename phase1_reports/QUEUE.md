# Phase 1 launch log — ALL 26 UNITS LAUNCHED

Concurrency ceiling is **20 subagents** (`CLAUDE_CODE_MAX_CONCURRENT_SUBAGENTS`), so the roster went
out as a wave of 20 plus backfill as slots freed. All 26 are now dispatched.

Tier A on Opus (6): 01 compute_kendall · 02 run_e2g_cv · 03 create_predictions ·
04 create_neighborhoods · 05 make_candidate_regions · 06 arc_e2g

Tier B on Sonnet (17): 07 make_kendall_pairs · 08 generate_atac_matrix · 09 bedSplitSort ·
10 make_external_features_config · 11 element_and_gene_summaries · 12 get_stats_per_cluster ·
13 plot_stats · 14 hover_plots · 15 threshold_e2g_predictions · 16 process_model_output ·
17 gen_num_candidate_enh_gene · 18 gen_num_tss_enh_gene · 19 gen_num_sum_nearby_enhancers ·
20 activity_only_features · 21 add_external_features · 22 gen_final_features ·
23 make_biosample_feature_table

Tier C on Sonnet, grouped (3): 24 fragment_shell_pipelines · 25 abc_external_binaries ·
26 run_blocks_and_checkpoints

Note: the `Agent` tool exposes no reasoning-effort parameter, only `model`. Tier A units were asked
for maximum reasoning effort in-prompt instead of being set to xhigh directly.

---

## Cross-cutting findings from returning reports — Phase 2 must account for these

**Agents corrected the launch-time priors about which variable drives their rule.** This was the
intended behaviour (derive from code, not from the suggested variable), and it means the roster's
rule↔variable guesses should not be trusted over the reports:
- **Unit 7** (`make_kendall_pairs`): neither `j_cand_elems` nor `k_genes` is an input at all. It
  re-anchors already-gene-tagged pairs onto raw peak intervals via `GenomicRanges::findOverlaps`.
  Real drivers: `n_cand_pairs` and `n_peaks`.
- **Unit 17** (`gen_num_candidate_enh_gene`): never opens `EnhancerList.txt` or `GeneList.txt`.
  Sole driver is `n_cand_pairs`.
- **Unit 10** (`make_external_features_config`): O(1) — it `grepl()`s the input *path* and never
  opens the ARC file it declares as a dependency. Also caught that the briefing's ~15× row labelled
  `add_external_features` belongs to the ENCODE_rE2G rule (unit 21), not this one.
- **Unit 14** (`hover_plots`): the tab template instantiates per **model** (M=1), not per cluster.

**A hard limit on how tight any `t_est` can be — from the only replicated rule.** Unit 9 found the
n=2 modality replicates disagree by **46% (telohaec_crispri) and 21% (k562_crispri)** on
byte-identical input, attributed to Lustre/Slurm contention that grows with job size. Phase 2's
`t_low`/`t_high` bounds cannot honestly be tighter than this for any fragment-stage rule, and the
critical path's bound must carry it.

**`O(n)` vs `O(n log n)` is not distinguishable anywhere in this dataset, even where `n` varies
15×.** Unit 9's reason is general and worth reusing: over `n_frag` = 26.2M → 394.8M, `ln(n_frag)`
itself spans only **1.16×**. So the log factor is nearly constant across the entire measured range.
Any report claiming to have empirically confirmed a log factor is overreaching.

**Two rules are I/O-shaped in a way the row counts hide.** Unit 9: `bedSplitSort.sh` round-trips
~24.3 GB in / ~27.4 GB out through **Lustre** (`-T` resolves to `RESULTS_DIR/tmp`, confirmed via
`mount`, not `$L_SCRATCH`) for a ~4 GB compressed input. Unit 17: benchmark `io_out` (15–60 MB)
undercounts the real 433–488 MB output — so `io_out` is not reliable as an I/O measure.

**Dead code confirmed present on the executed path:** `bedops`/`sort-bed` fast path in
`bedSplitSort.sh` never fires (`bedops` absent from `workflow/envs/sc_e2g.yml`); the
`Pairs.Kendall.tsv.gz` branch in `format_external_features_config_sc.R:40-42` is unreachable because
the checkpoint never emits `"Kendall"` alone.

**Zero-variance variables — cost classes in these are code-derived only, unfalsifiable here:**
`k_genes` (20,531, exactly constant), `n_feat_cols` (24 multiome / 18 scATAC, constant within
modality), cluster count `C` (always 6), model count `M` (always 1). Units 13, 14 and 26 depend on
`C` or `M`; none can be empirically confirmed even directionally.

---

## CONVERGENT CORRECTION — three independent agents overturned the same prompt prior

**My launch prompts pointed the post-scoring QC rules at `n_called_all` (~60k rows). That is wrong
for at least two of them.** Units 11 and 12 independently traced the rule `input:` and verified on
disk that the file actually read is `scE2G_predictions.tsv.gz` — the **unthresholded** output of
`run_e2g_qnorm`, whose row count equals **`n_cand_pairs`** (~10.2M–11.4M), about **170× larger**
than `n_called_all`. Unit 11 (`element_and_gene_summaries`) never touches the thresholded file at
all. Unit 12 (`get_stats_per_model_per_cluster`) reads *both*, but the unthresholded read dominates
by three orders of magnitude (~750–890 MB vs 4.7–5.5 MB), so `n_called_all` sets neither the cost
nor even the direction of the modality difference.

Phase 2 must use `n_cand_pairs`, not `n_called_all`, for units 11 and 12. Check unit 15
(`threshold_e2g_predictions`) and unit 16 against the same trap when they land — 15 was correctly
briefed on the input side, 16 genuinely is downstream of thresholding.

**A second table's width, not in the size manifest.** Units 11, 12 and 16 each measured the
`scE2G_predictions.tsv.gz` header directly and all three agree: **27 columns in multiome, 20 in
scATAC.** This is NOT the manifest's `n_feat_cols` (24 / 18), which describes
`genomewide_features.tsv.gz` — a different table. Both unit 11's ~1.33× and unit 12's ~1.25×
modality gaps are explained by this 27-vs-20 width at identical row count, which is why scATAC is
*faster* despite calling ~2× more links. Phase 2 needs this as a size variable and the manifest
does not carry it.

## `io_in` / `io_out` ARE UNUSABLE — six independent confirmations

Do not use the benchmark I/O columns for anything. Confirmed unreliable in both directions by units
9, 11, 12, 14, 16 and 17:
- **Undercounts:** unit 17 reported `io_out` 15–60 MB against a real 433–488 MB output; unit 16
  reported `io_out` 0.00 MB for jobs that wrote 6.5–13.2 MB; unit 14 reported `io_out ≈ 0` for a
  5.35 MB HTML; unit 12 saw `io_in` 0.2–30.6 MB against a true ~750 MB read (page-cache masking).
- **Overcounts:** unit 11 saw writes overcounted 100–400×; unit 12 saw `io_out` 1.5–3.4 GB for a
  one-row output file; unit 7 saw 2.7–3.1 GB `io_out` against an 83 MB file.

Same root cause as `hover_plots`' zero `cpu_time`: psutil does not see child process trees, and what
it does see is confounded by the page cache. **Use the on-disk `*_bytes` manifest columns, or
`ls`-measured file sizes, instead.**

## Job-level graph tables built for Phase 2

Parsed both resolved DAGs, propagating cluster identity along edges (labels alone don't carry it):

```
analysis/dag_{multiome,scatac}_nodes.tsv    node_id, rule, cluster, n_pred, n_succ
analysis/dag_{multiome,scatac}_edges.tsv    from_id, from_rule, to_id, to_rule, cluster
```

Cluster assignment resolved cleanly and symmetrically — **25 nodes per cluster in multiome, 21 per
cluster in scATAC**, with exactly 5 genuinely shared nodes in both: `all`,
`save_reference_configs`, `abc_generate_chrom_sizes_bed_file`, `hover_plots`, `plot_stats`. The
25 − 21 = 4 difference per cluster is precisely the multiome-only Kendall→ARC chain. The two
fan-in nodes (`plot_stats`, `hover_plots`) and the shared `abc_generate_chrom_sizes_bed_file` are
the only places cluster subgraphs touch, which means the six clusters are otherwise **fully
independent parallel branches** — important for the critical-path arithmetic.

## Verify pass

Workflow authored at `phase1_reports/verify_workflow.js`. Launch with
`Workflow({scriptPath: ".../phase1_reports/verify_workflow.js", args: [<report filenames>]})`
once the reports are on disk. It runs one skeptic per report (citations opened and checked, cost
functions re-evaluated, variable names validated against real manifest columns), a second
independent adjudicator prompted to *refute* the first on anything marked SUSPECT/FAIL, then a
cross-report synthesis for contradictions, double-counting, coverage gaps and the unconstrained set.
It is explicitly told that "the data cannot constrain this" is a virtue to be marked SOUND.

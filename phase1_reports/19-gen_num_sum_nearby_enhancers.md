# Unit 19 — `gen_num_sum_nearby_enhancers.py`

## 1. Unit identity

- Script: `ENCODE_rE2G/workflow/scripts/feature_tables/gen_num_sum_nearby_enhancers.py` (131 lines).
- Rule: `generate_num_sum_enhancers`, `ENCODE_rE2G/workflow/rules/genomewide_features.smk:116-151`.
- Modalities: **both**. No `threads:` declared — only `resources: mem_mb =
  partial(ABC.determine_mem_mb, min_gb=16)` (`:126-127`) — so it runs at Snakemake's default of 1
  requested CPU, consistent with the briefing's "every ABC-adjacent rule declares no threads."
- **Extra wildcard confirmed and job count corrected/confirmed**: the rule has wildcards
  `{biosample}` and `{kb}` (`NumEnhancersEG{kb}kb.txt` / `SumEnhancersEG{kb}kb.txt`,
  `genomewide_features.smk:129-130`). Pulling every row for this rule from both benchmark tables:
  ```
  awk -F'\t' -v r=generate_num_sum_enhancers 'NR==1||$2==r' benchmarks_multiome.tsv
  awk -F'\t' -v r=generate_num_sum_enhancers 'NR==1||$2==r' benchmarks_scatac.tsv
  ```
  gives exactly 6 rows per modality, and **every `wildcards` value is `<cluster>~5`** — e.g.
  `jurkat~5`, `telohaec_crispri~5`. **The roster's "6 jobs" is 6 clusters × exactly one `kb`
  value (5), not 6 clusters spread across multiple `kb` values.** Phase 0 never exercised any
  other window size for this rule, in either modality — 12 measured jobs total, all at `kb=5`.
  I confirmed `kb=5` is also the only value ever *consumed* downstream: `get_numNearbyEnhancers_file`
  / `get_numSumEnhancersEG_file`-equivalents (`genomewide_features.smk:54-65`) hard-code
  `NumEnhancersEG5kb.txt` / `SumEnhancersEG5kb.txt`, so `kb=5` is the only value the rest of the
  pipeline ever requests in this config; the wildcard exists for generality but is unexercised at
  any other value in these runs.

## 2. What the code actually does

`main()` (`:97-127`) reads `enhancer_list` (`EnhancerList.txt`, `j_cand_elems` rows after the
header) with `pd.read_csv(..., usecols=["chr","start","end","name"])` (`:110-112`) and converts
`--distance_threshold_kb` to bp (`:114`). All real work is in `generate_num_sum_enhancers()`
(`:7-83`), which does **four shell-outs via `os.system`** plus pandas post-processing:

1. **`:21-28`** — compute `midpoint = (start+end)/2` then `.apply(np.int64)` (row-wise Python
   call, not a vectorized cast) over all `j_cand_elems` rows, and write a 4-column BED
   (`chr, midpoint, midpoint, name`) to `enh_midpoint`. O(`j_cand_elems`).
2. **`:31-35`** — `bedtools slop -b {kb*1000} -i enh_midpoint -g chr_sizes > enh_expanded`:
   expands every element's midpoint by `kb` kilobases on each side, one line in → one line out.
   O(`j_cand_elems`), independent of `kb`'s value (slop is a per-row arithmetic pass regardless of
   window width).
3. **`:38-42`** — `zcat {abc_predictions} | csvtk cut -t -f chr,start,end,name,activity_base |
   sed '1d' > pred_slim`: decompresses and column-selects the **entire** predictions file. I
   verified `wc -l` on the resulting `pred_slim` equals `n_cand_pairs` from the size manifest
   exactly for every one of the 6 multiome clusters (e.g. `thp1_1`: 10,173,077 lines both ways;
   `jurkat`: 11,036,283 both ways) — this step is a single **O(`n_cand_pairs`) linear pass** over
   the full predictions file with zero filtering, independent of `kb`.
4. **`:45-51`** — `bedtools intersect -a enh_expanded -b pred_slim -wa -wb | sort -u >
   enh_pred_int`. **This is the dominant, `kb`-sensitive step and the one the unit brief asked me
   to characterise specifically.** It is neither a sorted-sweep/cumulative-sum implementation nor
   a naive Python nested loop — it delegates to `bedtools intersect`, an off-the-shelf spatial
   join, piped directly into `sort -u`. Two sub-costs:
   - `bedtools intersect` is called **without `-sorted`**, so per bedtools' own algorithm it
     builds an in-memory searchable interval structure over the **entire** `-b` file
     (`pred_slim`, `n_cand_pairs` rows) rather than doing a coordinate sweep — O(`n_cand_pairs`
     log `n_cand_pairs`) to build — then answers each of the `j_cand_elems` window queries in
     roughly O(log `n_cand_pairs` + hits-for-that-window).
   - The **hit volume itself is the `kb`-dependent term.** I measured it directly (this number
     exists nowhere in the benchmark table — the raw pre-dedup stream is piped straight into
     `sort -u` and never touches disk in production) by re-running the intersect half of the
     command on the smallest cluster's retained intermediates:
     ```
     bedtools intersect -a enhancerMidpoint_exp5kb.txt -b EnhancerPredictionsAllPutative_5kb.slim.txt -wa -wb | wc -l
     # thp1_1, kb=5: 22,395,901 raw hit rows (88 s wall for this diagnostic alone)
     ```
     `n_cand_pairs` for `thp1_1` is 10,173,077, so raw hits ≈ **2.20×** `n_cand_pairs` at `kb=5`.
     This is because `pred_slim` rows are (element, gene) pairs — every element repeats once per
     candidate gene it's paired to (~64–70× per element on average, `n_cand_pairs`/`j_cand_elems`)
     — and each element's ±5 kb window catches itself plus its immediate neighbours, each also
     repeated ~64–70× in `pred_slim`. Widening `kb` only increases how many *distinct elements*
     fall inside the window; it does not change the per-element repeat multiplicity. So raw hit
     volume ≈ `elements_per_window(kb) × (n_cand_pairs / j_cand_elems)`, and
     `elements_per_window(kb)` grows with `kb` (more genomic space per window → more neighbouring
     candidate elements captured). **`kb` therefore directly controls the size of the stream that
     `sort -u` must then sort** — `sort -u` over `H_raw(kb)` lines is O(`H_raw(kb)` log
     `H_raw(kb)`), and `H_raw(kb)` scales with `kb`.
   - `sort -u` then **collapses the duplication back down**: because a `pred_slim` row's first
     five columns (`chr,start,end,name,activity_base`) are identical for the same element across
     every gene it pairs with, exact-duplicate lines from different genes collapse to one. The
     post-dedup file (`enh_pred_int`) has only 322,411–360,643 lines across the 6 multiome
     clusters at `kb=5` — i.e. it is *not* `n_cand_pairs`-scaled at all, it's ≈2.03–2.29×
     `j_cand_elems` (self + a small number of genuinely-distinct nearby elements at this window
     width). **The 12–140× blow-up from `n_cand_pairs`-scale to `j_cand_elems`-scale happens
     entirely inside the intersect|sort pipe and leaves no trace in the benchmark table or on
     disk** — only the wall-clock/CPU cost of doing it is observable.
5. **`:54-58`** — `pd.read_csv(enh_pred_int, header=None)`: O(deduped hits), ≈2.0–2.3×
   `j_cand_elems` at `kb=5`.
6. **`:61-69`** — filter `int_df[3] != int_df[7]` (drop self-matches), then
   `.groupby([3]).size()` (count) and `.groupby([3])[8].sum()` (sum `activity_base`). Both
   groupbys are hash-based, O(deduped hits); no nested loop.
7. **`:71-83`** — write `out_num`/`out_sum`, O(`j_cand_elems`).

## 3. Size variable(s)

Named exactly as in `size_manifest_multiome.tsv` / `size_manifest_scatac.tsv`:

- **`n_cand_pairs`** (10,173,077–11,447,564, Tier 2, 1.13×) — governs step 3 (linear
  decompress+cut, `:38-42`) and, via the per-element repeat multiplicity `n_cand_pairs /
  j_cand_elems` (~64–70×), the raw hit volume of step 4's intersect (before `sort -u` collapses
  it).
- **`j_cand_elems`** (156,420–158,548, Tier 3, 1.01×) — governs steps 1, 2, 7 directly (O(j) each)
  and is the number of intersect *queries* (`-a` side) in step 4, and (after dedup) the
  approximate scale of steps 5–6.
- **`kb`** (rule wildcard; only value ever run is **5**, i.e. zero range, not even the narrow
  ranges the briefing warns about for Tier 3 variables) — controls how many distinct neighbouring
  elements fall in each window, and therefore the size of the raw stream `sort -u` must sort in
  step 4. This is the parameter the brief specifically asked me to surface: **the cost of this
  rule is parameterised by a config value, `distance_threshold_kb`, not purely by dataset size.**
- **`k_genes`** (20,531, constant) — **does not appear anywhere in this script's control flow.**
  It never reads `GeneList.txt` or loops over genes. It enters only indirectly, as the thing that
  sets the roughly-constant per-element multiplicity `n_cand_pairs/j_cand_elems` in the *upstream*
  ABC predictions file this script consumes — i.e. `k_genes` shapes the input this rule reads, but
  is not a variable this rule's own code touches. I'm flagging this the same way unit #17's report
  flagged it for `gen_num_candidate_enh_gene.py`, per that report's note to units #18/#19.

No literal `O(n_cand_pairs · k_genes)` cross-term exists in the code; the closest real cross-term
is `O(elements_per_window(kb) · (n_cand_pairs / j_cand_elems))` for the raw intersect volume,
which is not a manifest-named product because `elements_per_window(kb)` is not itself a manifest
column (see §7 for how it's estimated).

## 4. Asymptotic class derived from code

Per-step classes (all cited above): O(j) [steps 1,2,7] + O(n) [step 3] + **O(n log n) [interval
build] + O(H_raw(kb) log H_raw(kb)) [step 4's `sort -u`, dominant and `kb`-dependent]** + O(H_dedup)
[steps 5-6, H_dedup ≈ 2.0–2.3× j at kb=5].

The rule's overall class is:

**O(n_cand_pairs) + O(n_cand_pairs · log n_cand_pairs) + O(H_raw(kb) · log H_raw(kb))**,
where **H_raw(kb) ≈ (n_cand_pairs / j_cand_elems) · elements_per_window(kb)**, and
`elements_per_window(kb)` is non-decreasing in `kb` (a wider window cannot capture fewer
neighbouring candidate elements). This is neither of the two hypotheses the brief flagged as the
two extremes to distinguish (sorted-sweep/cumsum, or naive per-element neighbour rescanning) —
it's a general-purpose interval-join tool (`bedtools intersect`) whose cost structure sits between
them: it avoids literal O(j²) all-pairs comparison via an interval-tree/tree-query design, but
because it is invoked without `-sorted` and its `-b` side carries `n_cand_pairs`-scale duplication
(not `j_cand_elems`-scale), its *practical* cost is driven by `n_cand_pairs`, not by `j_cand_elems`
alone, and its **output volume before `sort -u`** is a real, measured `kb`-dependent quantity
(§2 point 4) that is completely invisible in the benchmark table.

Bounding the class (see §7 for the numeric forms):
- **Lower bound**: if the intersect+sort behaved as a single linear coordinate-sweep pass with
  negligible dependence on `kb` in this regime (e.g. if hit volume stayed of the same order as
  `n_cand_pairs` regardless of window width) — O(`n_cand_pairs`).
- **Upper bound**: the absolute ceiling on intersect output, regardless of `kb`, is one row per
  (window, predictions-row) pair — O(`j_cand_elems` · `n_cand_pairs`) — reached only in the
  pathological case where `kb` is wide enough that every window contains every predictions row.
  This is a hard structural ceiling from the code (a `-wa -wb` join can never emit more than
  `|a| × |b|` rows), not a fitted extrapolation.

## 5. Benchmark cross-check

All 12 measured jobs (all at `kb=5`):

| Modality | Cluster | `n_cand_pairs` | `j_cand_elems` | wall s | cpu_time s | io_in MB | io_out MB |
|---|---|---:|---:|---:|---:|---:|---:|
| multiome | jurkat | 11,036,283 | 158,397 | 167.81 | 83.62 | 136.78 | 830.49 |
| multiome | jurkat_pma_cd3_4hr | 10,854,698 | 157,395 | 100.83 | 65.39 | 43.61 | 737.86 |
| multiome | k562_crispri | 11,447,564 | 156,420 | 114.03 | 93.75 | 306.71 | 768.16 |
| multiome | telohaec_crispri | 11,106,687 | 157,551 | 180.36 | 93.13 | 148.18 | 759.89 |
| multiome | thp1_1 | 10,173,077 | 158,548 | 146.95 | 79.04 | 170.56 | 701.39 |
| multiome | thp1_2 | 10,278,325 | 157,498 | 157.08 | 70.33 | 442.00 | 706.95 |
| scATAC | jurkat | 11,036,283 | 158,397 | 104.69 | 66.72 | 59.76 | 749.41 |
| scATAC | jurkat_pma_cd3_4hr | 10,854,698 | 157,395 | 100.89 | 61.96 | 69.65 | 738.02 |
| scATAC | k562_crispri | 11,447,564 | 156,420 | 114.39 | 93.93 | 31.81 | 783.73 |
| scATAC | telohaec_crispri | 11,106,687 | 157,551 | 220.96 | 114.03 | 443.74 | 758.90 |
| scATAC | thp1_1 | 10,173,077 | 158,548 | 142.87 | 65.23 | 436.96 | 702.97 |
| scATAC | thp1_2 | 10,278,325 | 157,498 | 147.47 | 81.47 | 96.77 | 707.98 |

**`n_cand_pairs` spans only 1.13×; `j_cand_elems` spans only 1.01×; `kb` spans 1.00× (it never
varies).** With these ranges, the six/twelve points are a weak check at best. Sorting by
`n_cand_pairs`: the smallest-`n` cluster (`thp1_1`, 10.17M) is *not* the fastest (146.95 s / 142.87 s
— mid-pack), and the largest-`n` cluster (`k562_crispri`, 11.45M) is one of the *fastest* (114.03 s
/ 114.39 s). This is the same qualitative pattern unit #17 found for the same input file
(`EnhancerPredictionsAllPutative.tsv.gz`) — **not a coincidence, since this rule reads that same
file in step 3.** `cpu_time`/wall ratios (`mean_load` equivalent) are 0.22–0.64 (`k562_crispri`
scATAC: 93.93/114.39=0.82 is the high outlier; most sit near 0.4–0.6), meaning 20–80% of wall time
is not active CPU — some real I/O/subprocess-startup cost is present but it is not monotonic in
either size variable across these 12 points.

**Verdict: the benchmarks cannot confirm the O(n log n) class in `n_cand_pairs`, and they cannot
possibly say anything at all about the `kb` dependence, since `kb` never varied.** The one
`n_cand_pairs`-vs-cpu_time relationship that is directionally consistent (`k562_crispri` has both
the highest `n_cand_pairs` and among the highest `cpu_time` in both modalities: 93.75/93.93) is
suggestive but sits alongside `jurkat` (second-highest `n_cand_pairs`, mid-range cpu_time) and
`telohaec_crispri` (third-highest `n_cand_pairs`, highest wall time in both modalities — 180.36 s
multiome, 220.96 s scATAC, the two slowest jobs in the whole set) breaking any clean monotonic
story. I do not fit an exponent to this; per the brief's rule, I report the code-derived class and
flag it as unconfirmed.

## 6. Concrete overhead sources

- **The `kb`-dependent raw hit stream never touches disk and is invisible to `io_in`/`io_out`.**
  I directly measured it for `thp1_1` (§2 point 4): 22,395,901 raw lines piped from
  `bedtools intersect` into `sort -u`, collapsing to 322,411 lines on disk. Any cost model built
  only from the benchmark table's I/O columns would miss this entirely — it is pure intermediate
  pipe volume, real CPU/memory/sort cost with zero footprint in `io_in`/`io_out` or on-disk file
  size.
- **`io_in` substantially undercounts real read volume.** For `thp1_1` multiome, the on-disk
  bytes this job actually reads are: `enhancer_list` (33.3 MB) + compressed `abc_predictions`
  (281.5 MB, `n_cand_pairs_bytes`, verified against `ls -la` on
  `Predictions/EnhancerPredictionsAllPutative.tsv.gz`) + a re-read of `pred_slim`
  (664.5 MB, verified on disk) by the intersect step + a re-read of `enh_expanded` (9.0 MB) by the
  same step + a re-read of `enh_pred_int` (39.3 MB) by pandas — roughly **1.03 GB** of logical
  reads. The benchmark table reports `io_in = 170.56 MB` for that exact job — a **~6×
  undercount**, and the same pattern holds for `telohaec_crispri` (computed ≈1.12 GB vs. reported
  148.18 MB, ~7.6×). This is consistent with psutil's Linux `io_counters` reporting actual
  block-device bytes, not logical `read()` bytes — the intermediates this rule writes
  (`pred_slim`, `enh_expanded`) are re-read moments later while still warm in the page cache, so
  the re-read never reaches the block device and isn't counted, even though it costs real wall
  time (memory copies, decompression-adjacent parsing).
- **`io_out` is much closer to the true write volume for this rule** (unlike unit #17's finding
  for a sibling script). Summing the on-disk sizes of this rule's actual outputs for `thp1_1`
  (`pred_slim` 664.5 MB + `enh_midpoint` 9.0 MB + `enh_expanded` 9.0 MB + `enh_pred_int` 39.3 MB +
  `NumEnhancersEG5kb.txt` 3.4 MB + `SumEnhancersEG5kb.txt` 4.1 MB) ≈ 729.3 MB, versus the
  benchmark's `io_out = 701.39 MB` — within ~4%. The difference from unit #17 is plausibly that
  these writes go through external processes (`bedtools`, `sort`, shell redirection) that flush to
  the OS before the wrapping Python process exits, rather than through a single buffered
  `pandas.to_csv()` call.
- **Redundant full re-materialisation of the predictions file.** Step 3 (`:38-42`) writes a
  5-column, decompressed copy of the *entire* `n_cand_pairs`-row predictions file to disk
  (`pred_slim`, 664.5–722.7 MB across clusters) purely so `bedtools intersect` can read it back in
  step 4. This is a full write-then-reread of a table roughly 2.3–2.6× the size of the original
  compressed input, solely to strip 24 of 29 columns — the same "gzip in, uncompressed textual
  intermediate out" pattern flagged for the sibling script in unit #17, but here compounded by
  being read a second time immediately afterward.
- **`sort -u` runs without `--parallel`.** GNU `sort`'s default behaviour can use more than one
  core for the merge phase depending on the node's visible CPU count, even though the rule
  declares no `threads:` (so Snakemake requests 1 CPU from Slurm). The observed `cpu_time`/wall
  ratios (mostly 0.2–0.6, not >1) give no evidence of effective multi-core speedup actually being
  realised inside the Slurm cgroup for this job — consistent with the briefing's general finding
  that declared/realised parallelism gaps are common in this pipeline, though here I cannot
  isolate `sort`'s specific contribution from the combined 4-subprocess pipeline.
- **Four separate `os.system()` shell-outs** (`:31-35`, `:38-42`, `:45-51`, implicitly a fourth for
  the pipe in step 4) each pay conda/shell subprocess-spawn overhead; with wall times of
  100–220 s this is a small fraction of total cost but is non-zero, repeated four times per job.

## 7. Three functions: `t_est`, `t_low`, `t_high`

Let `n` = `n_cand_pairs`, `j` = `j_cand_elems`, `kb` = the rule's window half-width in kilobases
(only `kb=5` was ever measured — every calibration below is anchored at that single point).

From §2's direct measurement on `thp1_1`, define the raw-hit-volume ratio at the one measured
point: `ρ(5) = 22,395,901 / 10,173,077 = 2.201`. I model `ρ(kb)` as linear through the origin-ish
point `ρ(0)=1` (a zero-width window degenerates to self-only matches) and the one measured point:
`ρ(kb) ≈ 1 + 0.2402·kb`. **This is a single-point-calibrated straight line, not a fit — there is
no second `(kb, ρ)` pair anywhere in Phase 0 to check linearity against.**

```
t_est(n, kb) = 2.598e-7 * n * log2(n) * (1 + 0.2402*kb)         [seconds]
t_low(n)     = 1.337e-5 * n                                     [seconds]
t_high(n, j) = 8.477e-11 * n * j                                [seconds]
```

- `t_est` combines the O(n log n) interval-build/sort term with the measured `kb`-scaling factor
  from step 4; it collapses the separate O(n) linear-scan term (step 3) into the same coefficient
  because, as in unit #17, `log2(n)` barely moves across the observed 1.13× range of `n`, so O(n)
  and O(n log n) are not separable at this data density and there's no value in pretending
  otherwise with a second free parameter.
- `t_low` is the class floor from §4 — a pure O(n) pass, `kb`-independent, representing the
  best case where intersect+sort cost is not appreciably worse than the single linear
  decompress-and-cut pass that also happens in this job.
- `t_high` is the class ceiling from §4 — the structural O(n·j) bound on `-wa -wb` join output
  size, which holds for **any** `kb` (it's the absolute worst case the join operator could produce,
  not a `kb`-specific extrapolation), so no `kb` term is needed or justified in this bound.
- All three coefficients are calibrated to coincide at the pooled multiome mean (`n =
  10,816,106`, `j = 157,635`, `kb=5`, mean wall = 144.51 s) **by construction**, per the brief's
  instruction not to blend point estimates: they are meant to diverge away from this point, not
  average out at it.

## 8. Evaluated at mean / min / max

All at `kb = 5` (the only value ever run):

| Point | cluster | `n_cand_pairs` | `j_cand_elems` | `t_low` | `t_est` | `t_high` |
|---|---|---:|---:|---:|---:|---:|
| mean (6 multiome clusters) | — | 10,816,106 | 157,635 | 144.6 s | 144.6 s | 144.5 s |
| min | `thp1_1` | 10,173,077 | 158,548 | 135.9 s | 135.5 s | 136.7 s |
| max | `telohaec_crispri` | 11,106,687 | 157,551 | 148.5 s | 148.6 s | 148.3 s |

For reference, the actual measured wall times: `thp1_1` 146.95 s (multiome) / 142.87 s (scATAC) —
both bracket the ~135–137 s model band loosely; `telohaec_crispri` 180.36 s (multiome) / 220.96 s
(scATAC) — **both well above the ~148 s model band**, `telohaec_crispri` scATAC being the single
slowest job in the whole 12-job set by a wide margin. Because `n` and `j` span so little range, all
three functions sit within ~10 s of each other at every evaluated point — the observed 100–220 s
spread across real jobs is **larger than the spread the model can produce from size alone**,
reinforcing §5's conclusion that non-size noise (scheduling, filesystem state) dominates the
variance at this scale, not that the model is well-calibrated.

## 9. Does the class differ between modalities?

**No — same class, and the constant looks the same modulo run-to-run noise.** `n_cand_pairs` and
`j_cand_elems` are numerically identical between the multiome and scATAC size manifests for every
cluster (verified directly from both `size_manifest_*.tsv` files, e.g. `telohaec_crispri`:
`n_cand_pairs`=11,106,687 and `j_cand_elems`=157,551 in both). This is **not** because the two
modalities share the ABC output on disk — `ABC_BIOSAMPLES_DIR` resolves independently per config
(`ENCODE_rE2G/workflow/rules/utils.smk:121-131`, defaulting to `RESULTS_DIR/{dataset}`, and
`results_dir` differs between `configs/igvf10_multiome_config.yaml` and
`configs/igvf10_scatac_config.yaml`) — it's because both configs run ABC candidate-calling
independently but deterministically from the **same `atac_frag_file`** for a given cluster
(`configs/tables/igvf10_cell_clusters_{multiome,scatac}.tsv` list identical `atac_frag_file` paths
per cluster, differing only in `model_dir`), so the two independent ABC runs converge on the same
candidate set and predictions row count. This rule sits entirely upstream of the ARC/Kendall
branch point (`checkpoint features_required`) that the briefing identifies as the real source of
modality-dependent cost for `add_external_features` et al. — this rule is not in that four-rule
list, and the benchmarks agree it doesn't belong there: same-cluster multiome-vs-scATAC wall times
differ by anywhere from −0.03 s (`jurkat_pma_cd3_4hr`, near-identical) to +40.6 s
(`telohaec_crispri`, scATAC slower) with no consistent sign or magnitude relationship to any size
variable — the signature of run-to-run noise on an identical workload, not a modality effect.

## 10. What the measurement cannot tell us

- **Nothing about the `kb` dependence at all.** Every one of the 12 measured jobs ran at
  `kb=5`. There is no second data point to check even the *sign* of the modeled `ρ(kb)` slope,
  let alone its linearity. §2's direct 22.4M-line measurement establishes that the raw
  pre-dedup hit volume is real and large at `kb=5`, but says nothing about how it would change at,
  say, `kb=1` or `kb=25` — the linear model in §7 is a code-motivated guess about a single
  unmeasured axis, not a calibrated relationship.
- **Cannot confirm O(n log n) in `n_cand_pairs`** — the 1.13× range and the non-monotonic
  wall/cpu pattern in §5 (same underlying issue as unit #17's report on the same input file) mean
  the six/twelve points cannot separate O(n) from O(n log n), or either from a constant dominated
  by I/O/scheduling noise.
- **Cannot separate the four `os.system()` subprocess stages' individual contributions to wall
  time.** The benchmark table reports one wall/cpu/I/O measurement for the whole job; there is no
  per-subprocess breakdown, so I cannot say what fraction of the 100–220 s is the step-3 decompress
  pass versus the step-4 intersect+sort versus pandas post-processing. My 88 s stand-alone
  re-measurement of the intersect-only half of step 4 on `thp1_1` (§2) is suggestive that it is a
  large fraction of that job's 146.95 s total, but it was run outside the Snakemake job context
  (different node, different cache state), so it is an approximate, not exact, decomposition.
- **Cannot rule out page-cache effects as the explanation for the `io_in` gap in §6.** I
  hypothesised warm-cache re-reads of just-written intermediates based on the ~6–7.6× gap between
  computed logical read volume and the benchmark's `io_in`, but I have no direct cache-hit
  instrumentation to confirm this versus some other undercounting mechanism in the psutil-based
  benchmark plugin.
- **These are single (n=1) measurements per cluster per modality** (12 points, no replicates) —
  the scatter in §5 and §9 is consistent with pure measurement/scheduling noise but cannot be
  proven to be only that.
- **Cannot say anything about behaviour outside the observed 10.17M–11.45M `n_cand_pairs` /
  156,420–158,548 `j_cand_elems` band**, or at any `kb` other than 5. The O(n log n) +
  `kb`-dependent-hit-volume class in §4 is a code-derived asymptotic statement; §7's numeric
  bounds are anchored to the single measured band and should not be trusted far outside it.

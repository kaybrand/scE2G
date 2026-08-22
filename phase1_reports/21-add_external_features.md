# Unit 21 — `add_external_features`

## 1. Unit identity

- Script: `ENCODE_rE2G/workflow/scripts/feature_tables/merge_external_features.R` (134 lines)
- Rule: `add_external_features`, `ENCODE_rE2G/workflow/rules/genomewide_features.smk:179-193`
- Modality: **both**, 6 jobs each (jurkat, jurkat_pma_cd3_4hr, k562_crispri, telohaec_crispri,
  thp1_1, thp1_2). 12 measured jobs total.
- **Disambiguation confirmed**: this is NOT `make_external_features_config`
  (`workflow/rules/add_external_features.smk:37-52`, unit #10, a trivial config-formatting rule
  with sub-15 s runtimes). I verified against the benchmark tables:

  ```
  awk -F'\t' -v r=add_external_features 'NR==1||$2==r' benchmarks_multiome.tsv
  ```
  gives 756.96–1275.25 s (12.6–21.3 min) across the 6 multiome clusters, and the scATAC table
  gives 52.53–75.87 s — an exact match to the range quoted in the assignment. This is my rule.

## 2. What the code actually does

Upstream mechanism (owned by unit #26, confirmed here by direct inspection): `checkpoint
features_required` (`workflow/rules/add_external_features.smk:1-22`) reads each cluster's
`feature_table.tsv` and writes `"ARC"` or `"Neither"` to `to_generate.txt`. `features_to_generate()`
(`add_external_features.smk:26-34`) resolves the ARC-branch path to
`.../ARC/EnhancerPredictionsAllPutative_ARC.tsv.gz`; `make_external_features_config` then either
emits a 5-row config (multiome) or a 0-row, header-only config (scATAC) — **confirmed on disk**:
`external_features_config.tsv` is 868–928 bytes / 5 data rows for all 6 multiome clusters and
exactly 60 bytes / 0 data rows for all 6 scATAC clusters (also independently confirmed by unit
#10). The 5 multiome rows all share one `source_file` (the ARC table) and list
`join_by = "overlap"` for all 5 (`format_external_features_config_sc.R:32-36`).

`merge_external_features.R` main script (`:100-134`):

1. `abc <- fread(predictions_extended)` (`:103`) — reads `ActivityOnly_features.tsv.gz`
   (**19 columns**, confirmed by `zcat … | head -1 | tr '\t' '\n' | wc -l` on both modalities —
   multiome 19, scATAC 18; row count = `n_cand_pairs`).
2. `features <- fread(feature_table_file)` (`:104`) and `ext_ft <- fread(external_features_config)`
   (`:105`) — both tiny (config is 60 B scATAC / ~900 B multiome).
3. `needed = setdiff(input_features, colnames(abc))` (`:108-110`) then `ext_ft = filter(ext_ft,
   input_col %in% needed)` (`:111`). I confirmed on disk that none of the 5 candidate columns
   (`mean_log_normalized_rna`, `RnaPseudobulkTPM`, `RnaDetectedPercent`, `Kendall`,
   `ARC.E2G.Score`) exist in `ActivityOnly_features.tsv.gz`'s 19 columns, so for multiome all 5
   survive the filter; for scATAC `ext_ft` was already 0 rows, so it stays 0 rows regardless.
4. `unique_sources = na.omit(unique(ext_ft$source_file))` (`:112`) — **1** for multiome (all 5 rows
   share one path), **0 (`character(0)`)** for scATAC.
5. `if (length(unique_sources>0)){ for (i in 1:length(unique_sources)) {…} }` (`:115-131`). Despite
   the odd-looking `unique_sources>0` (a character-vs-numeric comparison), `length()` of that
   comparison equals `length(unique_sources)`, so the conditional is functionally `if
   (length(unique_sources) > 0)`. **For scATAC this is `if (0)` — the entire loop body never
   executes.** No file is read, no join happens, no rows or columns change.
   For multiome the loop runs exactly **once** (one unique source):
   - `this_source <- fread(unique_sources[1])` (`:119`) — reads the **entire ARC table**,
     `EnhancerPredictionsAllPutative_ARC.tsv.gz`, 36 columns, same row count as `abc`
     (confirmed: `wc -l` on both `ActivityOnly_features.tsv.gz` and the ARC file for thp1_1 both
     give 10,173,078 lines).
   - `colnames(this_source)[...] <- this_ext$input_col` (`:121`) — cosmetic rename (source_col ==
     input_col here, so a no-op in practice, but still a full column-name scan).
   - `join_by=="overlap"` for all 5 rows → `overlap_feature_with_abc(abc, this_source,
     feature_score_cols = <5 names>, ...)` is called **once**, covering all 5 feature columns
     together (`:123-124`). The `merge_feature_by_gene` / `dplyr::left_join` branch (`:88-98`,
     used only when `join_by == "TargetGene"`) is **not exercised** in this project — every row in
     both configs uses `"overlap"`.
6. Inside `overlap_feature_with_abc` (`:14-68`):
   - Builds two `GRanges` objects, one per table, each ~10.17M–11.45M ranges, with
     `seqnames = paste0(chr, ":", TargetGene)` (`:22-27`) — i.e. overlaps are only ever computed
     within matching chr:gene groups.
   - Unifies `seqlevels` across both objects (`:29-32`) — a full-vector `unique()` and factor-level
     reassignment over ~2×n elements.
   - `findOverlaps(abc_gr, feat_gr)` (`:35`) — an interval-tree overlap search, sort/hash-based,
     not a nested pairwise loop.
   - `merged <- cbind(abc[queryHits(ovl)], feature[subjectHits(ovl), feature_score_cols])` (`:38`)
     — materializes a new ~n-row, 24-column table (a full copy).
   - `aggregate_features(merged, feature_score_cols, agg_cols, agg_fun)` (`:42-43`, defined
     `:71-85`): `agg_cols <- setdiff(colnames(merged), feature_score_cols)` — **all 19 non-feature
     columns of `merged`** become the group-by key. `mapply` over the 5 `feature_score_cols`
     (`:74-78`) runs **5 separate `data.table` `by=agg_cols` group-bys**, each over ~10-11M rows
     keyed on that 19-column composite key (cardinality ≈ `n_cand_pairs`, since each row is
     already essentially unique on `chr,start,end,TargetGene,CellType,...` — this is a wide-key
     "group-by" that mostly produces groups of size 1). The 5 resulting tables are then combined
     via `Reduce(function(df1,df2) merge(df1,df2,by=agg_cols), agg_list)` (`:81`) — **4 sequential
     merges**, again keyed on the same 19-column composite key.
   - `missing <- abc[setdiff(seq_len(nrow(abc)), queryHits(ovl)), ]` (`:50`) — a `setdiff` over up
     to ~11M integers; benign if `perc_overlap` is near 100% (consistent with the confirmed
     matching row counts, though the actual percentage isn't captured in the benchmark files).
   - `output <- output[order(chr, start, end, TargetGene), ]` (`:64`) — an explicit **full sort**
     of the merged ~n-row table by 4 columns, materializing yet another full copy.
7. `fwrite(abc, plus_external_features)` (`:134`) writes the result — for scATAC this is the
   **unmodified input table** (confirmed: 18 columns in, 18 columns out, byte size essentially
   unchanged, `189.8–223.5 MB` → `189.3–222.9 MB`); for multiome this is `abc` **plus the 5 new
   columns** (confirmed: 19 → 24 columns, `203–238 MB` → `494–605 MB`, verified directly with
   `zcat … | head -1 | tr '\t' '\n' | wc -l` on thp1_1 for both stages).

**Library-load note:** the script unconditionally loads `data.table`, `tidyverse`, and
`GenomicRanges` (`:7-9`) even on the scATAC path, where none of `GenomicRanges`'s functionality
(`GRanges`, `findOverlaps`) is ever called. This is a real, code-grounded fixed cost paid by
every job regardless of modality (see §6).

## 3. Size variables

Named exactly as in `size_manifest_{multiome,scatac}.tsv`:

- `n_cand_pairs` — rows of `ActivityOnly_features.tsv.gz` = rows of the ARC table = rows of the
  output. Governs the group-by/merge/sort cost in multiome; governs nothing in scATAC (loop body
  never runs).
- `arc_bytes` — compressed size of `EnhancerPredictionsAllPutative_ARC.tsv.gz`, read **only** in
  multiome (`this_source <- fread(unique_sources[1])`, `:119`). Absent entirely from scATAC's cost.
- `actonly_bytes` — compressed size of `ActivityOnly_features.tsv.gz`, read (`:103`) in **both**
  modalities.
- `feat_bytes` — compressed size of `genomewide_features.tsv.gz` (the *next* rule's output, not
  this rule's). I use it here only as a documented-manifest proxy for this rule's own output
  volume, because I measured it directly: at thp1_1, `feat_bytes` = 494,990,740 B vs. the
  actually-measured `ActivityOnly_plus_external_features.tsv.gz` size of 494,468,191 B (within
  0.1%) — close enough to use as a named stand-in, flagged as an approximation, not an identity.

Cross term used below: `O(n_cand_pairs · log(n_cand_pairs))` for the sort/join machinery, and
`O(arc_bytes + actonly_bytes + feat_bytes)` for I/O. These do **not** collapse into one variable —
`arc_bytes` exists as an I/O term only in multiome, and the `n_cand_pairs`-driven CPU term exists
only when `unique_sources` is non-empty (multiome only).

## 4. Asymptotic class derived from code

**scATAC: O(n_cand_pairs) [equivalently O(actonly_bytes)], with a zero-valued join/sort term.**
The `if (length(unique_sources) > 0)` gate at `:115` is `if (0)` on every scATAC job (0-row
config, confirmed on disk), so the entire loop body (`:116-130`) — the only place any join, sort,
or GRanges call exists in this script — never executes. What remains is exactly one `fread` and
one `fwrite` of the same ~n-row, 18-column table: a single linear pass, no comparison-based
operation at all.

**Multiome: O(arc_bytes + actonly_bytes + feat_bytes) I/O term + O(n_cand_pairs · log(n_cand_pairs))
CPU term.** Justification, cited above: `findOverlaps` (`:35`, interval-tree, sort/hash-based),
5× `data.table` `by=`-group aggregation (`:74-78`, radix-sort based), 4× `merge()` on the same
composite key (`:81`, dispatches to `merge.data.table`, also sort/hash-based), and 1× `order()`
full-table sort (`:64`). None of the operations actually invoked is a nested O(n²) loop over row
pairs — every one is sort- or hash-based — but there are **~10 such passes** (5 aggregates + 4
merges + 1 order) stacked sequentially, each over a wide **19-column composite key** whose
cardinality is ≈ `n_cand_pairs` itself (rows are already unique on that key, so each "group-by" is
mostly groups of size 1 — an expensive way to do what is, given the confirmed 1:1 row-count match
between `abc` and the ARC table, functionally a straight column-bind). This wide-key, 10-pass
design is the code-grounded reason the *constant* is unusually large, not evidence of a worse
exponent (data.table's radix-based grouping/merging remains O(n log n) regardless of key width).

## 5. Benchmark cross-check

```
multiome  jurkat               1275.25 s   max_rss 20234.2 MB   cpu_time 1221.23  (cpu/wall .957)
multiome  jurkat_pma_cd3_4hr    805.77 s   max_rss 19914.0 MB   cpu_time  796.81  (cpu/wall .988)
multiome  k562_crispri          851.95 s   max_rss 20971.1 MB   cpu_time  825.00  (cpu/wall .968)
multiome  telohaec_crispri      993.36 s   max_rss 20358.1 MB   cpu_time  953.66  (cpu/wall .960)
multiome  thp1_1                756.96 s   max_rss 18705.4 MB   cpu_time  731.54  (cpu/wall .966)
multiome  thp1_2                774.01 s   max_rss 18891.9 MB   cpu_time  760.40  (cpu/wall .983)

scatac    jurkat                 52.53 s   max_rss  2880.8 MB   cpu_time   44.83  (cpu/wall .853)
scatac    jurkat_pma_cd3_4hr     55.91 s   max_rss  1497.4 MB   cpu_time   41.75  (cpu/wall .747)
scatac    k562_crispri           55.77 s   max_rss  2496.8 MB   cpu_time   44.35  (cpu/wall .795)
scatac    telohaec_crispri       75.87 s   max_rss  1521.7 MB   cpu_time   51.95  (cpu/wall .685)
scatac    thp1_1                 56.06 s   max_rss  1425.6 MB   cpu_time   37.01  (cpu/wall .660)
scatac    thp1_2                 56.65 s   max_rss  1436.3 MB   cpu_time   36.97  (cpu/wall .653)
```

`n_cand_pairs` spans only 10,173,077–11,447,564 (**1.13×**) and `arc_bytes` spans
630,109,119–778,223,611 (**1.23×**) — both Tier 2, narrow. With 6 collinear points I cannot fit an
exponent, and the data show it: **jurkat is the slowest multiome job (1275.25 s) despite NOT
having the largest `arc_bytes` (732.2 MB, 5th of 6), `n_cand_pairs` (11,036,283, 4th of 6), or
`max_rss` (20,234 MB, 3rd of 6)** — k562_crispri has the largest `arc_bytes` (778.2 MB) and
`n_cand_pairs` (11,447,564) and *max_rss* (20,971 MB) but is only 851.95 s, well below jurkat.
This is a genuine code-vs-data disagreement for the ranking (not just noise in the exponent): I
cannot explain jurkat's outlier status from any manifest variable. Most likely hypothesis (not
verifiable with n=1 per cluster): Lustre I/O contention, GC-timing variance in a script that
builds ~10 full-table copies (see §6), or node co-tenancy — flagged, not resolved.

What the data *do* confirm cleanly, independent of the narrow-range problem:
- **Multiome is CPU-bound** (cpu/wall 0.957–0.988), consistent with the sort/join/aggregate
  machinery in §4, not I/O wait.
- **scATAC is markedly less CPU-bound** (cpu/wall 0.653–0.853) and its `cpu_time` (37.0–51.9 s) is
  large relative to a job whose only work is one `fread`+`fwrite` of a ~200 MB gzip table — this
  is consistent with the unconditional `library(GenomicRanges)`/`library(tidyverse)` load noted in
  §2 dominating scATAC's runtime (see §6).
- **Bytes moved badly under-explain the ratio.** Using `io_in`+`io_out` as a byte-moved proxy:
  scATAC moves `io_in`+`io_out` ≈ 215.6+1803.1 = 2018.7 MB (mean); multiome moves
  488.75+5913.39 = 6402.1 MB (mean) — only **3.2×** more bytes than scATAC. But multiome takes
  909.55/58.80 = **15.5×** longer (mean wall). A 3.2× I/O increase cannot explain a 15.5× time
  increase — the residual ~4.8× must come from the CPU-bound join/aggregate/sort pipeline that
  scATAC's zero-length loop skips entirely. This is the strongest available support for "class
  differs, not just constant" (see §9).
- `io_out` (5,913 MB mean, multiome; 1,803 MB mean, scATAC) is **~10–11×** larger than the
  actual compressed output file (~551 MB / ~208 MB respectively) in *both* modalities — this ratio
  is roughly consistent across modality, so it is not a multiome-specific effect; I cannot explain
  the mechanism (possibly uncompressed intermediate buffering inside `fwrite`'s gzip connection, or
  the benchmark sampler counting page-cache writeback separately from the final file) — flagged in
  §10.
- `io_in` (488.75 MB mean, multiome) is *less* than half the ~929 MB (`arc_bytes`+`actonly_bytes`)
  this rule is declared to read, while scATAC's `io_in` (215.6 MB mean) closely matches its single
  ~206 MB input. Plausible explanation: the ARC table was just written by the immediately upstream
  `arc_e2g` rule and is likely still warm in page cache when this rule reads it, so much of that
  read never reaches the block device counter that `io_in` tracks — a hypothesis, not confirmed.

## 6. Concrete overhead sources

- **An entire extra ~700 MB table read, present only in multiome.** `arc_bytes` (630–778 MB
  compressed) is read in full (`:119`) only when `unique_sources` is non-empty; scATAC never opens
  this file at all. This alone is a categorical (not incremental) I/O difference between
  modalities.
- **Ten full-table sort/hash passes over a 19-column composite key** (§4) where a single pass
  would functionally suffice, given the confirmed 1:1 row correspondence between `abc` and the ARC
  table (equal row counts, `perc_overlap` message at `:47` not captured in the benchmark but
  consistent with near-100%). This is architecturally wasteful but not a worse asymptotic class —
  it inflates the constant, and directly explains why `n_cand_pairs`'s narrow 1.13× range still
  costs 750–1275 s rather than a few seconds.
- **Multiple simultaneous full-size in-memory copies plausibly explain the ~18.7–21.0 GB peak
  RSS.** Rough R object-size accounting (not a profiler trace — an estimate): `abc` (19 cols,
  ~10-11M rows) ≈ 1.9 GB; `this_source`/ARC table (36 cols) ≈ 3.3 GB; two `GRanges` objects
  ≈ 0.6–1.0 GB combined; the `cbind` `merged` result (24 cols) ≈ 2.4 GB; and — critically —
  `mapply` at `:74` builds `agg_list`, a list of **all 5** per-feature aggregate tables
  (~2.0 GB each ≈ 10 GB total) **before** the `Reduce(merge, …)` step at `:81` can consume them,
  so all 5 are resident simultaneously by construction of `mapply`. Summing what is plausibly
  concurrently live at that moment (1.9 + 3.3 + ~1.0 + 2.4 + 10 ≈ 18.6 GB) lands within the
  observed 18.7–21.0 GB band. This is consistent with — not proof of — "the whole ARC table plus
  the whole feature table plus (several) copies are resident simultaneously"; I have no
  object-level memory trace to confirm the exact composition, only that the order of magnitude
  matches and the code structure (a `mapply` that must materialize its whole result list before
  the caller can reduce it) makes multiple simultaneous full-size copies unavoidable, not just
  possible.
- **Unconditional heavy library loads** (`library(data.table)`, `library(tidyverse)`,
  `library(GenomicRanges)` at `:7-9`) regardless of modality. For scATAC, where the only real work
  is one `fread`+`fwrite` of a ~200 MB gzip table, `cpu_time` of 37–52 s is disproportionate to
  that data volume and is consistent with `GenomicRanges`/`tidyverse` load time dominating —
  `GenomicRanges` is never actually used on the scATAC path (its functions are called only inside
  `overlap_feature_with_abc`, which is never reached, `:115`).
- **No `threads:` declared** (`genomewide_features.smk:179-193`) — runs single-threaded by
  Snakemake's default. cpu/wall ≈ 0.96–0.99 for multiome confirms the work is genuinely
  single-core CPU-bound (not multi-core, not I/O-wait); there is no declared-vs-realized threading
  gap here because none is declared, consistent with the project-wide "every ABC-adjacent rule
  declares no `threads:`" pattern.
- **Resource request is config-driven, not size-driven**: `mem_mb=partial(ABC.determine_mem_mb,
  min_gb=min_mem)` where `min_mem` is 32 GB if `config["final_score_col"] ==
  "E2G.Score.qnorm"` else 8 GB (`genomewide_features.smk:174-178`) — the floor depends on which
  scoring mode the whole pipeline run uses, not on any per-cluster size variable, so it cannot
  adapt to the ~1.24× `arc_bytes` spread actually observed.
- **Redundant re-derivation of `colnames` matching** (`:121`) is a no-op here (`source_col ==
  input_col` for all 5 ARC rows, confirmed from the printed config) but still executes a full
  column-name vector scan every job.

## 7. Three functions

I/O rate calibration: from scATAC's pure linear pass (no join term at all), total bytes moved
(read `actonly_bytes` + write ≈ same size) ≈ 2 × 206 MB = 400 MB mean, in ≈ (58.80 − 35) = 23.8 s
of size-dependent time (35 s taken as the fixed library/conda/Slurm floor, close to the observed
`cpu_time` floor of 37.0 s at the smallest cluster) → **R_io ≈ 1.68×10⁷ B/s**. This is a rough,
single-pipeline-stage calibration, not a general filesystem benchmark — stated explicitly as a
guess.

CPU-term calibration: solving the multiome mean point (mean wall 909.55 s; mean `n_cand_pairs`
10,816,106; mean `arc_bytes`+`actonly_bytes`+`feat_bytes` = 1,481,034,927 B) for the one free
constant `c_sj` in `c_sj · n_cand_pairs · log2(n_cand_pairs)`, after subtracting the I/O term
(88.16 s at `R_io`), gives **c_sj ≈ 3.249×10⁻⁶ s per (row·log2(row)) unit**.

```
t_est_multiome(n_cand_pairs, arc_bytes, actonly_bytes, feat_bytes) =
    (arc_bytes + actonly_bytes + feat_bytes) / 1.68e7
  + 3.249e-6 * n_cand_pairs * log2(n_cand_pairs)

t_est_scatac(actonly_bytes) =
    35.0
  + (2 * actonly_bytes) / 1.68e7
    [the "2×" reflects one read + one write of the same table; no join/sort term — the
     loop body at merge_external_features.R:115-131 has a zero-valued coefficient in scATAC,
     confirmed by the 0-row config and the unchanged 18-column output]
```

Class bounds (multiome only — scATAC's class is not in doubt, see §4):

```
t_low_multiome(n_cand_pairs, arc_bytes, actonly_bytes, feat_bytes) =
    (arc_bytes + actonly_bytes + feat_bytes) / 1.68e7
  + 7.593e-5 * n_cand_pairs
    [floor: treats the ~10 sort/join/aggregate passes as if collapsed into one efficient
     linear pass — the best a correctly-implemented single join could do, given the
     confirmed 1:1 row correspondence between abc and the ARC table]

t_high_multiome(n_cand_pairs, arc_bytes, actonly_bytes, feat_bytes) =
    (arc_bytes + actonly_bytes + feat_bytes) / 1.68e7
  + 10 * 3.249e-6 * n_cand_pairs * log2(n_cand_pairs)
    [ceiling: charges the full n·log(n) rate once for each of the ~10 sort/hash passes
     actually present in the code (5 mapply aggregates at merge_external_features.R:74-78,
     4 Reduce merges at :81, 1 order() at :64), i.e. assumes none of the redundant passes
     benefit from any shared work or caching]
```

Both `t_low` and `t_high` are O(n log n)-or-better by construction — nothing in the code (§4)
behaves worse than sort/hash-based, so I am not proposing an O(n²) ceiling; the uncertainty here
is in the *constant* (how much the 19-column composite key and the ~10 redundant passes actually
cost), not in the exponent.

## 8. Evaluated at mean / min (`thp1_1`) / max (`telohaec_crispri`)

| Point | cluster | `n_cand_pairs` | `arc_bytes` (MB) | `actonly_bytes` (MB) | `feat_bytes` (MB) | t_low | t_est | t_high | measured |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| mean | (6-cluster avg) | 10,816,106 | 709.3 | 220.1 | 551.7 | 909.6 s* | 909.6 s* | 5,689 s | mean 909.55 s |
| min | thp1_1 | 10,173,077 | 630.1 | 203.0 | 495.0 | 851.2 s | 848.9 s | 5,275 s | 756.96 s |
| max | telohaec_crispri | 11,106,687 | 750.5 | 226.6 | 579.1 | 933.6 s | 937.3 s | 5,779 s | 993.36 s |

\* `t_low` and `t_est` are both calibrated to hit the mean exactly by construction (one free
constant each, fit to one point) — the min/max columns are the actual cross-validation, and both
land within ~12% of measured. `t_high` is a deliberate ceiling (§7) and is *not* expected to match
— it bounds the class/constant uncertainty, not the point estimate; it lands ~5–7× above measured
at all three points, which is consistent with the "~10 redundant wide-key passes" charge being
pessimistic relative to what data.table actually achieves in practice.

scATAC, same three points (`t_est_scatac` only depends on `actonly_bytes`):

| Point | cluster | `actonly_bytes` (MB) | t_est_scatac | measured |
|---|---|---:|---:|---:|
| mean | (6-cluster avg) | 206.0 | 59.5 s | mean 58.80 s |
| min | thp1_1 | 189.8 | 57.6 s | 56.06 s |
| max | telohaec_crispri | 212.5 | 60.3 s | 75.87 s (outlier — see §5) |

## 9. Does the class differ between modalities, or only the constant?

**The class differs — this is one of the genuine cases, not an overreach.** For scATAC, the
join/sort machinery's contribution is not "small," it is **exactly zero**: the 0-row
`external_features_config.tsv` (confirmed 60 bytes/header-only on all 6 scATAC clusters) makes
`unique_sources` empty, so the `if (length(unique_sources) > 0)` gate at
`merge_external_features.R:115` is never entered — no `GRanges` object is built, no
`findOverlaps` call happens, no `data.table` group-by or `merge()` runs, no `order()` sort
happens, and the ~700 MB ARC file is never opened. The scATAC path is a single O(n) read+write
pass, full stop. For multiome, all of that machinery runs, once, over ~10-11M rows on a
19-column composite key (§4). This is corroborated quantitatively in §5: bytes moved
(`io_in`+`io_out`) increase only 3.2× from scATAC to multiome, but wall time increases 15.5× — a
gap far too large to attribute to the constant on a shared linear I/O term; it requires an
additional term (the sort/join machinery) that literally does not exist in the scATAC code path.
Additionally, multiome pays a categorical extra I/O cost (`arc_bytes`, ~700 MB) that has no
scATAC counterpart at all — not a "same file, different size" scaling, but a whole extra file
that is read in one modality and never referenced in the other.

## 10. What the measurement cannot tell us

- **Cannot separate the I/O-rate constant from the sort/join constant.** `arc_bytes` and
  `n_cand_pairs` co-vary across the 6 multiome clusters (both grow with cluster depth), so the
  split I made in §7 (`R_io` from scATAC extrapolation, `c_sj` from the multiome residual) is a
  modeling choice calibrated to fit the mean, not an independently-verified decomposition. A
  different split would fit the mean equally well.
- **Cannot confirm an n log n vs. a plain linear exponent from these 6 points.** `n_cand_pairs`
  spans only 1.13× and `arc_bytes` only 1.23× across multiome clusters — Tier 2, narrow — so the
  exponent claim in §4 rests on reading the code (GRanges/data.table/merge/order are sort-and-hash
  based), not on curve-fitting. `t_low` and `t_high` in §7 bound this uncertainty explicitly;
  they are not distinguishable from the benchmark data alone.
- **Cannot explain the jurkat multiome outlier (1275.25 s, the single highest value despite not
  having the largest `arc_bytes`, `n_cand_pairs`, or `max_rss` — see §5).** With n=1 per cluster
  there is no way to tell whether this is a genuine property of that cluster's data or
  node/scheduler/cache noise on the day it ran.
- **Cannot verify the ~18.6 GB memory-composition estimate in §6** beyond order-of-magnitude
  plausibility — I have no object-level memory profiler trace, only R object-size arithmetic from
  column counts and row counts. The true peak composition (which objects are alive at the exact
  moment of peak RSS, and whether R's garbage collector delayed freeing any of them) is not
  something this benchmark instrumentation (wall time, `max_rss`, `io_in`/`io_out`, `cpu_time`)
  can resolve.
- **Cannot explain the `io_out` ≈ 10–11× final-file-size mismatch**, in either modality. It is
  present in both (not multiome-specific), so it is likely a general artifact of how `fwrite`'s
  gzip connection or the benchmark sampler counts bytes — but I cannot pin down the mechanism from
  this data.
- **Cannot confirm `perc_overlap` (the fraction of ARC rows found by `findOverlaps`) is actually
  ~100%** — the `message()` at `:47` reporting it is not captured in the benchmark files. I infer
  it is high from the confirmed equal row counts between `abc`, the ARC table, and the output, but
  a near-100% match with a small "missing" fallback (`:50-58`) executing on a handful of rows is
  also consistent with that observation and would not be visible in these numbers.

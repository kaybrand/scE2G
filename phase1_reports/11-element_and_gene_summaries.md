# Unit 11 — `generate_element_gene_lists.R`

## 1. Unit identity

- Script: `workflow/scripts/prediction_qc/generate_element_gene_lists.R` (67 lines).
- Rule: `element_and_gene_summaries`, `workflow/rules/sc_predictions.smk:136-155`.
- Modalities: **both**, 6 jobs each = 12 measured jobs total (`{cluster} × {model_name}`, one
  `model_name` per modality: `multiome_powerlaw_v3` / `scATAC_powerlaw_v3`). No `threads:`
  declared — runs at 1, consistent with the briefing's "every ABC rule declares no threads"
  cross-cutting finding. `resources: mem_mb = encode_e2g.ABC.determine_mem_mb`.
- Inputs (`sc_predictions.smk:137-141`): `Neighborhoods/GeneList.txt` (`abc_gene_list`),
  `Neighborhoods/EnhancerList.txt` (`abc_element_list`), `{model_name}/scE2G_predictions.tsv.gz`
  (`prediction_file`), and a conditional `gene_expr_file` (`Kendall/gene_expression_metrics.tsv.gz`
  if the ARC/Kendall branch ran, else `RESULTS_DIR` — a directory path, never actually read because
  the R script only reads it inside the multiome-only `if` branch, see §2).
- Outputs: `{model_name}/scE2G_gene_list.tsv.gz`, `{model_name}/scE2G_element_list.tsv.gz`.

## 2. What the code actually does — and a correction to this unit's briefed expectation

**Correction up front:** this unit's brief hypothesised `n_called_all` (the *thresholded*
predictions row count) as the governing variable, and asked me to check whether the code filters
`class == "promoter"` to decide between `n_called_all`/`n_called_distal`. Neither applies. I
verified directly:

- The rule's `prediction_file` input (`sc_predictions.smk:140`) is
  `{cluster}/{model_name}/scE2G_predictions.tsv.gz` — the **unthresholded** output of
  `run_e2g_qnorm`, not the thresholded output of `filter_sc_e2g_predictions`
  (`scE2G_predictions_threshold{threshold}.tsv.gz`, a different rule/file entirely,
  `sc_predictions.smk:87` region).
- I counted rows directly on disk (`thp1_1`, multiome): `scE2G_predictions.tsv.gz` has
  **10,173,078** lines (header + 10,173,077 data rows) — exactly `n_cand_pairs` for that cluster
  (manifest: 10,173,077). The thresholded file `scE2G_predictions_threshold0.177.tsv.gz` has only
  60,436 lines, matching `n_called_all` (60,435) — but that file is **not** an input to this rule.
- So the governing variable for the one large table this script reads is **`n_cand_pairs`**
  (10.17M–11.45M rows, Tier 2, range 1.13×), roughly **170×** larger than `n_called_all`. The
  script never touches `n_called_all` or `n_called_distal`, and there is no `class == "promoter"`
  filter anywhere in the 67 lines. This is exactly the kind of "pairs vs links" trap the briefing
  warns about, just one level removed from where the unit brief expected it.

Control flow, in order, all citations `generate_element_gene_lists.R:line`:

1. `:12` `abc_gene_list <- fread(GeneList.txt)` — `k_genes` rows (20,531, constant).
2. `:13` `abc_element_list <- fread(EnhancerList.txt)` — `j_cand_elems` rows (156,420–158,548).
3. `:14` `pred <- fread(prediction_file)` — **`n_cand_pairs`** rows (10.17M–11.45M), 27 columns in
   multiome / 20 columns in scATAC (verified on disk, see §9). This is a single linear parse pass,
   O(`n_cand_pairs`).
4. `:31-33` — gene side: `dplyr::select(pred, TargetGene, TargetGeneEnsembl_ID, any_of(gene_columns))
   %>% distinct() %>% rename(...)`. This is a full pass over all `n_cand_pairs` rows of `pred`,
   reduced by `distinct()` down to one row per `TargetGene` — i.e. `k_genes` output rows. `distinct()`
   here is `dplyr::distinct()` (package loaded at `:6`), which is **hash-based** (`vctrs::vec_unique_loc()`
   under the hood) — a single-pass, expected-linear dedup, **not** a `data.table`-keyed
   `unique()`/`setkey()` operation and **not** a `group_by()`/`split()` (neither is called anywhere
   in this script). O(`n_cand_pairs`).
5. `:35-39` — `abc_gene_list %>% select(...) %>% mutate(%in%) %>% left_join(pred_genes, by=...) %>%
   mutate(...)`. All of these operate on `abc_gene_list` (`k_genes` rows) joined against
   `pred_genes` (also `k_genes`-scale after the dedup in step 4) — O(`k_genes`), negligible next to
   `n_cand_pairs`.
6. `:41-52` — **branch on column presence**: `if ("RNA_pseudobulkTPM" %in% colnames(pred))`. True
   only when `pred` carries RNA columns, i.e. multiome. In that branch (`:42-49`): a second small
   `fread` of `gene_expr_file` (`gene_expression_metrics.tsv.gz`, `k_genes` rows, 390 KB on disk),
   a `left_join` on `name` (`k_genes`-scale), and two `mutate()` calls — all O(`k_genes`). In the
   `else` branch (`:51`, scATAC), only a single `mutate()` runs and `gene_expr_file` is never read
   at all despite being declared as a rule input.
7. `:55` `fwrite(gene_list, genes_out, ...)` — writes `k_genes` rows (~20,532 incl. header).
8. `:61-62` — element side: `dplyr::select(pred, chr, start, end, any_of(element_columns)) %>%
   distinct()`. Same shape as step 4 but keyed on `(chr, start, end, ...)`, reducing all
   `n_cand_pairs` rows to `j_cand_elems` distinct rows. Second full O(`n_cand_pairs`) hash-dedup
   pass over `pred`.
9. `:64` `left_join(abc_element_list, pred_elements, by = c("chr","start","end"))` — join of two
   `j_cand_elems`-scale tables, O(`j_cand_elems`).
10. `:65` `fwrite(element_list, elements_out, ...)` — writes `j_cand_elems` rows (~158,549 incl.
    header).

**No sort, no nested loop, no cross join against a large table anywhere.** The entire cost is
concentrated in one `fread` and two `select()+distinct()` passes over the full `n_cand_pairs`-row
`pred` table (steps 3, 4, 8); everything else operates on the small `k_genes`/`j_cand_elems`-scale
tables and is asymptotically irrelevant at this data's scale (`k_genes`/`n_cand_pairs` ≈ 0.002,
`j_cand_elems`/`n_cand_pairs` ≈ 0.015).

## 3. Size variable(s)

**Primary: `n_cand_pairs`** — 10,173,077–11,447,564 rows (Tier 2, 1.13× range), read once
(`:14`) and scanned twice more for dedup (`:31`, `:61`). Confirmed **identical between multiome
and scATAC for the same cluster** (e.g. `thp1_1`: 10,173,077 in both `size_manifest_multiome.tsv`
and `size_manifest_scatac.tsv`) — this rule's dominant cost driver does not change with modality,
even though the rule appears in the launcher's "runs in both, materially different cost" family of
concerns. See §9.

**Secondary, output-scale only: `k_genes`** (20,531, constant) and **`j_cand_elems`**
(156,420–158,548, Tier 3) — these bound the sizes of `abc_gene_list`/`abc_element_list` and the
post-`distinct()` tables, and are the row counts of the two output files. They do not drive the
dominant cost term; the two `left_join()`s and the `fwrite()`s that touch them are O(`k_genes`)
and O(`j_cand_elems`) respectively, both ≪ `n_cand_pairs`.

**No cross term.** `n_cand_pairs` is never multiplied against `k_genes` or `j_cand_elems` in this
script — the joins in steps 5 and 9 are keyed joins between two already-small tables, not a scan
of `pred` against `abc_gene_list`/`abc_element_list`.

## 4. Asymptotic class derived from code

**O(`n_cand_pairs`)** — a small, fixed number of linear passes over the full prediction table:
one `fread` parse (`:14`) and two `dplyr::select()+distinct()` hash-dedup passes (`:31-33`,
`:61-62`). `dplyr::distinct()` is hash-based (`vctrs::vec_unique_loc`), which is expected-linear,
not `O(n log n)` sort-based — there is no `sort()`, `order()`, `arrange()`, `setkey()`, or
`data.table`-style keyed `unique()` anywhere in the file. The remaining operations (`%in%`/`match()`
at `:37`, the two `left_join()`s at `:38` and `:64`, `fwrite()` at `:55`/`:65`, the small `fread` at
`:43`) are all O(`k_genes`) or O(`j_cand_elems`), strictly smaller terms that do not change the
class.

Full model: `T = c1·n_cand_pairs (parse) + c2·n_cand_pairs (gene distinct) + c3·n_cand_pairs
(element distinct) + O(k_genes) + O(j_cand_elems)`, which collapses to **O(`n_cand_pairs`)** since
the `k_genes`/`j_cand_elems` terms are 2–3 orders of magnitude smaller.

## 5. Benchmark cross-check

All 12 measured jobs (`benchmarks_{multiome,scatac}.tsv`, filtered to `element_and_gene_summaries`):

| Modality | Cluster | `n_cand_pairs` | wall s | cpu_time s | mean_load | max_rss MB | io_in MB | io_out MB |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| multiome | jurkat | 11,036,283 | 38.86 | 6.31 | 16.21 | 178.15 | 85.38 | 1239.78 |
| multiome | jurkat_pma_cd3_4hr | 10,854,698 | 31.41 | 13.67 | 43.51 | 197.01 | 118.59 | 3013.62 |
| multiome | k562_crispri | 11,447,564 | 34.90 | 11.27 | 32.27 | 187.30 | 281.51 | 2365.17 |
| multiome | telohaec_crispri | 11,106,687 | 36.83 | 9.13 | 24.73 | 196.94 | 322.99 | 1840.60 |
| multiome | thp1_1 | 10,173,077 | 28.10 | 14.84 | 52.32 | 582.28 | 112.31 | 3078.26 |
| multiome | thp1_2 | 10,278,325 | 28.33 | 14.83 | 51.74 | 589.42 | 276.89 | 3112.40 |
| scATAC | jurkat | 11,036,283 | 21.17 | 14.92 | 69.71 | 2563.97 | 97.54 | 2100.99 |
| scATAC | jurkat_pma_cd3_4hr | 10,854,698 | 21.33 | 15.04 | 69.72 | 2517.68 | 245.32 | 2191.61 |
| scATAC | k562_crispri | 11,447,564 | 22.54 | 14.35 | 63.31 | 2109.62 | 37.07 | 2245.74 |
| scATAC | telohaec_crispri | 11,106,687 | 33.46 | 4.74 | 14.17 | 177.87 | 330.79 | 886.92 |
| scATAC | thp1_1 | 10,173,077 | 25.50 | 8.95 | 35.09 | 227.51 | 310.46 | 1939.17 |
| scATAC | thp1_2 | 10,278,325 | 25.42 | 9.18 | 36.10 | 228.86 | 288.32 | 1960.15 |

**`n_cand_pairs` spans only 1.13× — too narrow to fit or confirm an exponent, and the data doesn't
even show clean monotonicity.** Sorted by `n_cand_pairs`: `k562_crispri` has the *largest*
`n_cand_pairs` (11,447,564) of all 6 clusters but is **not** the slowest job in either modality
(34.90 s multiome, ranked 3rd of 6; 22.54 s scATAC, ranked 4th of 6) — `jurkat` (11,036,283,
smaller `n`) is the slowest multiome job (38.86 s), and `telohaec_crispri` (11,106,687) is the
slowest scATAC job (33.46 s). This is the same pattern the briefing predicts for near-constant/
mildly-variable Tier 2/3 variables: **the wall-clock ranking does not track the size-variable
ranking**, so the six-per-modality points can be used to sanity-check the class is *not* wildly
superlinear, but cannot confirm O(`n_cand_pairs`) over, say, O(`n_cand_pairs`·log(`n_cand_pairs`)).

**A real, code-consistent modality signal, separate from `n_cand_pairs`:** mean wall time is
33.07 s (multiome) vs 24.90 s (scATAC) — a **1.33× ratio**. This lines up well with the column-count
difference confirmed on disk (§9): `pred`'s `fread` parses **27 columns in multiome vs 20 in
scATAC** (ratio 1.35×, almost exactly the observed wall-time ratio), not with the (absent) `n`
difference — `n_cand_pairs` is identical between modalities for the same cluster. This is a clean
agreement between code (more columns tokenized per row in the multiome branch) and data (a
consistent ~1.33× constant-factor gap), even though the size-variable-vs-wall-time relationship
*within* a modality is noisy.

**`cpu_time`/wall ratios are low and erratic** (0.16–0.71, no relation to `n_cand_pairs`; matches
`mean_load`/100 almost exactly in every row, e.g. `thp1_1` multiome 14.84/28.10=0.528 ≈
`mean_load`=52.32/100), meaning 30–85% of wall time in every job is not active CPU — consistent
with gzip-decompression I/O wait and Slurm scheduling latency dominating the noise budget at this
job length (20–40 s), on top of whatever the ~O(`n_cand_pairs`) compute itself costs.

**Disagreement flagged, not smoothed:** `scATAC`/`telohaec_crispri` (33.46 s) is an outlier relative
to its four scATAC siblings with near-identical `n_cand_pairs` (21.2–25.5 s) despite not having a
distinguishably larger `n_cand_pairs` (11,106,687, 3rd of 6) — its `mean_load` (14.17) and
`cpu_time` (4.74 s) are also the lowest of the 12 jobs, meaning this job spent the *most* wall time
waiting per unit of CPU work of any of the 12 — most plausibly a slow/contended read of the input
predictions file on that particular run (cold Lustre cache, competing I/O on the node) rather than
anything algorithmic; I cannot confirm this from a single measurement (no replicates).

## 6. Concrete overhead sources

- **Two full-table hash dedups whose output is 0.02–1.5% of their input.** The `distinct()` calls
  at `:32` and `:62` each scan all `n_cand_pairs` rows (10.17M–11.45M) to produce `k_genes`
  (20,531) or `j_cand_elems` (≤158,548) output rows — 99.8% and 98.5% of the scanned rows are
  discarded. This is a real, code-confirmed inefficiency (two full linear passes purely to recover
  small per-gene/per-element tables that are logically derivable from `abc_gene_list`/
  `abc_element_list` plus a smaller feature source) — noted per scope as a mechanism, not a fix
  proposal.
- **Column-count asymmetry drives the modality constant.** `pred`'s on-disk header has 27 columns
  in multiome (includes `RNA_meanLogNorm`, `RNA_pseudobulkTPM`, `RNA_percentCellsDetected`,
  `Kendall`, `ARC.E2G.Score`, `normalizedATAC_enh`, `E2G.Score.qnorm.ignoreTPM`) vs 20 in scATAC —
  verified directly (`zcat .../scE2G_predictions.tsv.gz | head -1`). `fread` must tokenize every
  column of every row regardless of what's later selected, so this ~1.35× column ratio is the most
  plausible driver of the ~1.33× mean wall-time ratio between modalities (§5, §9) — a much cleaner
  match than the compressed-byte-size ratio (see next point).
- **Compressed input volume does not track the wall-time ratio as well as column count does.**
  On-disk `scE2G_predictions.tsv.gz` sizes: multiome 742.7–889.2 MB (mean 819.3 MB across the 6
  clusters) vs scATAC 414.0–476.5 MB (mean 445.8 MB) — a **1.84× byte ratio**, larger than the
  1.33× wall-time ratio and the 1.35× column ratio. This suggests gzip-decompression bandwidth is
  *not* the bottleneck (if it were, wall time should scale closer to bytes); per-row/per-column
  tokenization and hashing cost is the more consistent explanation.
- **Benchmark `io_in`/`io_out` do not match real I/O volume for this rule — use on-disk file sizes
  instead.** Real read volume ≈ `pred` (742.7–889.2 MB multiome, 414.0–476.5 MB scATAC) +
  `EnhancerList.txt` (`enh_list_bytes`, 32.9–33.3 MB) + `GeneList.txt` (`k_genes_bytes`, 5.57–5.62
  MB) + (multiome only) `gene_expression_metrics.tsv.gz` (~0.39 MB) ≈ 781–928 MB (multiome) /
  453–515 MB (scATAC). The benchmark's `io_in` column reports only 37–385 MB (undercounts the true
  read volume by roughly 2–20×). Real write volume ≈ `scE2G_element_list.tsv.gz` +
  `scE2G_gene_list.tsv.gz` ≈ 7.0–8.6 MB (multiome) / 7.9–8.1 MB (scATAC, verified `thp1_1`: 7.08 MB
  + 1.56 MB = 8.64 MB multiome; 6.81 MB + 1.12 MB = 7.93 MB scATAC), but the benchmark's `io_out`
  column reports 887–3112 MB — overcounting the true write volume by roughly **100–400×**. Neither
  direction is trustworthy for this rule; the on-disk file sizes are.
- **No threading.** No `threads:` on the rule (`sc_predictions.smk:136-155`); `data.table::fread()`
  can use multiple threads internally by default, but with no `threads:` declared Snakemake
  requests 1 CPU from Slurm, and the observed `cpu_time`/wall ratios (0.16–0.71, mostly well under
  1) show no evidence of effective multi-core work — consistent with (though not proof of) the
  briefing's "declared threads are largely fictional" pattern extending here even where no threads
  were ever declared.
- **The `gene_expr_file` rule input is provided but never read on the scATAC path.** Snakemake
  still resolves and depends on `get_gex_file()` (`sc_predictions.smk:132-135`) for both
  modalities, but the R script (`:41`) only opens it inside the multiome-only branch — on scATAC
  it is a DAG dependency with no corresponding read, a harmless but real mismatch between declared
  and used inputs.

## 7. Three functions

Class from §4 is **O(`n_cand_pairs`)**, confirmed by the absence of any sort/nested-loop/cross-join
in the code. The six-per-modality benchmark points (§5) span only 1.13× in `n_cand_pairs` and do
not even rank monotonically with wall time, so they cannot confirm the exponent — `t_low`/`t_high`
below bound the **class**, not a point estimate: `t_low` is a pure-linear floor (the fastest
observed per-row rate), `t_high` is an `n·log(n)` ceiling as a safety margin against undetected
superlinear hashing/rehash costs that a 1.13×-range benchmark could not expose even if present
(there is no sort in the code, but I cannot rule out hash-table resize costs behaving worse than
linear at scales well beyond what was measured).

```
t_est(n_cand_pairs; multiome) = 3.06e-6 * n_cand_pairs                       [seconds]
t_est(n_cand_pairs; scATAC)   = 2.30e-6 * n_cand_pairs                       [seconds]
t_low(n_cand_pairs)           = 1.90e-6 * n_cand_pairs                      [seconds]
t_high(n_cand_pairs)          = 1.505e-7 * n_cand_pairs * log2(n_cand_pairs) [seconds]
```

Calibration:
- `t_est`'s two coefficients are each the modality's mean observed wall time divided by the mean
  `n_cand_pairs` across the 6 clusters (10,816,106): multiome 33.07 s / 10,816,106 = 3.06e-6 s/row;
  scATAC 24.90 s / 10,816,106 = 2.30e-6 s/row. The 1.33× ratio between them matches the 27-vs-20
  column ratio (1.35×, §6), which is why they are kept as two separate constants rather than one
  pooled value — the class is shared (§9) but the constant genuinely differs by modality.
- `t_low`'s coefficient is the minimum observed wall-time/`n_cand_pairs` rate across all 12 jobs
  (scATAC `jurkat`: 21.17 s / 11,036,283 = 1.918e-6 s/row, rounded down slightly to 1.90e-6 as a
  floor).
- `t_high`'s coefficient is fit through the single slowest observed job (multiome `jurkat`: 38.86 s
  at `n_cand_pairs` = 11,036,283, log2(n) = 23.398): 38.86 / (11,036,283 × 23.398) = 1.505e-7.

These are **not** an exponent fit to the 12 points (§5 shows that would fit noise) — `t_low`/
`t_high` are calibrated only so the linear floor and n·log(n) ceiling pass through the observed
extremes, with the class itself (linear, from code) coming from §4.

## 8. Evaluated at mean / min / max

| Point | cluster | `n_cand_pairs` | `t_low` | `t_est` (multiome) | `t_est` (scATAC) | `t_high` |
|---|---|---:|---:|---:|---:|---:|
| mean (across 6 clusters) | — | 10,816,106 | 20.6 s | 33.1 s | 24.9 s | 38.1 s |
| min | `thp1_1` | 10,173,077 | 19.3 s | 31.1 s | 23.4 s | 35.6 s |
| max | `telohaec_crispri` | 11,106,687 | 21.1 s | 34.0 s | 25.6 s | 39.2 s |

(For reference, the actual measured wall times at these named clusters: `thp1_1` 28.10 s multiome
/ 25.50 s scATAC — the multiome value is below `t_est` (31.1 s) and above `t_low`; `telohaec_crispri`
36.83 s multiome / 33.46 s scATAC — both close to or (scATAC) above the modality-blind `t_est`
scATAC value, reflecting the outlier flagged in §5. All 12 real points fall inside the
[`t_low`, `t_high`] envelope except none exceed `t_high` and none fall below `t_low` — the envelope
holds, but it is wide relative to the tight `t_est` numbers, which is the expected shape given how
little the benchmarks actually constrain the exponent.)

## 9. Does the class differ between modalities?

**No — only the constant differs, and I can point to exactly why.** `n_cand_pairs` is numerically
identical between multiome and scATAC for every cluster (verified against both
`size_manifest_{multiome,scatac}.tsv`), because the input `pred` file's row count is fixed
upstream by the shared ABC candidate-calling stack, independent of the `checkpoint
features_required` ARC/Kendall branch. The code *does* take a genuinely different path per
modality (§2 step 6: the `if ("RNA_pseudobulkTPM" %in% colnames(pred))` branch at `:41-52` runs
only for multiome), but that branch operates entirely at `k_genes` scale (one extra small `fread`,
one join, two `mutate()`s) — 3–4 orders of magnitude below `n_cand_pairs`, so it cannot explain a
population-level wall-time difference on its own.

What *does* explain the observed ~1.33× multiome/scATAC wall-time ratio (§5) is the column-count
difference in the shared, dominant `fread`/`distinct()` passes: 27 columns (multiome) vs 20
columns (scATAC) in the same `n_cand_pairs`-row `pred` table (verified on disk, §6) — a 1.35×
ratio that lines up with the 1.33× time ratio much better than the underlying compressed-byte-size
ratio (1.84×) does. So: **same class (O(`n_cand_pairs`)), same dominant operations
(parse + two hash-dedups), different constant driven by how many columns those operations touch
per row** — this rule does not belong to the launcher's four-rule "genuinely different cost by
modality" list (`add_external_features`, `run_e2g_qnorm`, `gen_final_features`,
`get_stats_per_model_per_cluster`), and the modality gap here (~1.3×) is much smaller than that
list's largest entries (up to ~15×).

## 10. What the measurement cannot tell us

- **Cannot confirm O(`n_cand_pairs`) over any nearby exponent (e.g. `n_cand_pairs`·log(`n_cand_pairs`))
  from these 12 points.** The range is 1.13× and the wall-time ranking is not even monotonic in
  `n_cand_pairs` within a modality (§5: `k562_crispri` has the largest `n` but is not the slowest
  job in either modality). The O(`n_cand_pairs`) call in §4 is a code-reading conclusion, not an
  empirical fit; `t_low`/`t_high` in §7 bound the class as a hedge against that.
- **Cannot distinguish "column count causes the modality gap" from "modality gap is coincidental
  noise"** with certainty — the 1.35× column ratio matching the 1.33× time ratio (§6, §9) is
  suggestive, not proof, since I only have 6 jobs per modality and no replicates.
- **Cannot explain the `telohaec_crispri`/scATAC outlier (33.46 s vs 21–25 s for its four
  size-matched siblings)** from the benchmark table alone — its low `mean_load`/`cpu_time` point to
  I/O wait, but I cannot attribute that to cold cache, Slurm node contention, or genuine
  algorithmic behavior without a replicate run.
- **Cannot trust the benchmark's `io_in`/`io_out` columns for this rule at all** (§6) — they
  undercount real reads by ~2–20× and overcount real writes by ~100–400× relative to the on-disk
  file sizes I measured directly. Any cost model built from this rule's benchmark I/O columns
  should be discarded in favor of the file-size-based figures in §6.
- **These are single (n=1) measurements per cluster per modality**, not replicates (per the
  briefing's general caveat) — so I cannot separate "this job hit a slow node/cold cache" from
  "this cluster is systematically more expensive" for any of the scatter described in §5.
- **Cannot say anything about behavior at `n_cand_pairs` far outside the observed 10.17M–11.45M
  range** (e.g. an order of magnitude more candidate pairs) — whether `dplyr::distinct()`'s
  hash-table growth stays linear at that scale is a code-level plausibility argument (§4, §7), not
  something these clusters can validate.

# Unit 22 — `gen_final_features` (Tier B)

## 1. Unit identity

- Scripts:
  - `ENCODE_rE2G/workflow/scripts/feature_tables/gen_final_features.R` (54 lines)
  - `ENCODE_rE2G/workflow/scripts/feature_tables/get_fill_values.R` (33 lines)
- Rule: `gen_final_features`, `ENCODE_rE2G/workflow/rules/genomewide_features.smk:196-208`. No
  `threads:` declared (runs at 1, consistent with every other ABC/feature rule already noted
  in the briefing).
- Modality: **both**. 6 jobs each (one per cluster), 12 jobs total.
- File-identity note (per instructions): there are two copies of `get_fill_values.R` on disk —
  `ENCODE_rE2G/workflow/scripts/feature_tables/get_fill_values.R` (mine) and
  `workflow/scripts/model_application/get_fill_values.R` (another unit's). `diff` between them
  is empty — **byte-identical**. I characterise only the ENCODE_rE2G copy below.

## 2. What the code actually does

`gen_final_features.R`:

1. `gen_final_features.R:13` — `fread` the tiny `feature_table.tsv` config (12 rows in
   multiome, 6 in scATAC for every measured cluster — see §5).
2. `gen_final_features.R:15` — `fread` the actual input, `ActivityOnly_plus_external_features.tsv.gz`
   (`n_cand_pairs` rows × `n_feat_cols` columns). This is the dominant read.
3. `gen_final_features.R:18-26` — presence check of `config$input_col`/`second_input` against
   `colnames(df)`; operates on the column-name vector only (`O(n_feat_cols)`), never touches
   row data.
4. `gen_final_features.R:29-32` — interaction-term loop: for each row of `intx` (config rows
   with a non-blank `second_input`), computes `df[[a]] * df[[b]]` — one full `n_cand_pairs`-length
   vector multiply per interaction term. **Empirically dormant**: in all 12 measured jobs
   `second_input` is blank for every config row (verified below), so `nrow(intx) == 0` and this
   loop runs zero iterations in every measured job.
5. `gen_final_features.R:35-38` — rename loop over `single` (the non-interaction config rows,
   12 multiome / 6 scATAC): reassigns `names(df)[...]`, i.e. touches only the character vector
   of column names, `O(k_config)` where `k_config` is tiny; no row-data cost.
6. `gen_final_features.R:41` → `get_fill_values.R` (see below).
7. `gen_final_features.R:42` — `tidyr::replace_na(df, replace = fill_values)`: a full pass over
   every column named in `fill_values` (12 of 24 columns in multiome, 6 of 18 in scATAC —
   roughly half of `n_feat_cols` either way), scanning all `n_cand_pairs` rows of each to
   substitute `NA`. This is the second full-table pass.
8. `gen_final_features.R:44-51` — builds `output <- select(df, all_of(core_cols), all_of(config$feature))`.
   **This variable is never used again.** Line 54 writes `df`, not `output`. This `select()`
   call is a dead full-table copy: it re-materialises all `n_cand_pairs × n_feat_cols` cells for
   no purpose. It is real, measured wall/CPU time and real transient memory, wasted on an
   unused object — a genuine (if modest) overhead source, not a class-changing one, since it's
   just one more `O(n_cand_pairs · n_feat_cols)` pass alongside the others.
9. `gen_final_features.R:54` — `fwrite(df, ..., sep="\t", na="NA", quote=FALSE)` to a `.gz`
   path: full serialization of `df` to text plus gzip compression. Confirmed by direct
   measurement (§6) that the compression step's real I/O cost tracks the **uncompressed** text
   volume, not the final compressed file size.

`get_fill_values.R` (called once, at `gen_final_features.R:41`):

- **What it computes**: per-feature imputation ("fill") values used to replace `NA`s before the
  final table is written. Two sources, both metadata-only unless the "mean" branch fires:
  - `get_fill_values.R:4-8` — pulls the literal `fill_value` from the config table for every
    feature present in `colnames(features)`; `filter`/`select`/`distinct`/`deframe` operate on
    the config table only (`k_config` rows, e.g. 12/6) — **no pass over the big table**.
  - `get_fill_values.R:11-12` — flags which features have `fill_value == "mean"`.
  - `get_fill_values.R:15-27` — **only if any feature uses `"mean"`**: subsets `df` to those `m`
    columns (`features[, ..mean_features]`, `get_fill_values.R:18`, a copy of `n_cand_pairs × m`
    cells) and calls `apply(mean_features, MARGIN = 2, FUN = ...)` (`get_fill_values.R:21-23`) to
    compute one mean per column, excluding `Inf`/`-Inf`/`NA`. This is a genuine full pass over
    the `m` "mean" columns — `O(n_cand_pairs · m)` — but it is a **per-column summary statistic**
    (one scalar output per column), not a join, sort, or anything superlinear.
  - `get_fill_values.R:30` — reshapes the (tiny) fill-value vector into a named list.
- **Empirically dormant in all 12 measured jobs**: every `feature_table.tsv` for every cluster
  in both modalities has `fill_value == "0"` for every feature (checked directly — see §5), so
  `mean_fill_values` is `all FALSE` and `m = 0` in every measured run. The full-table branch
  (`get_fill_values.R:15-27`) never executes on Phase 0's data.
- **Class**: as written, `get_fill_values.R` is `O(k_config)` (metadata only) when no feature
  uses `"mean"`, and `O(k_config + n_cand_pairs · m)` in general, where `m` ≤ number of config
  rows with `fill_value == "mean"`. It is never `O(n_cand_pairs · n_feat_cols)` — even in the
  worst case it only touches the subset of columns that request a mean, not the whole table.
  In the runs actually measured, its cost is effectively zero (`k_config` ≈ 6–12, negligible).

## 3. Size variables (named as in `size_manifest_*.tsv`)

- `n_cand_pairs` — rows of the input/output tables (Tier 2, 10.17M–11.45M, 1.13× range).
  **Identical per cluster across modalities** (verified directly: e.g. `thp1_1` has
  `n_cand_pairs = 10173077` in both `size_manifest_multiome.tsv` and `size_manifest_scatac.tsv`)
  — the ABC candidate set doesn't change; only the feature columns differ.
- `n_feat_cols` — columns of `genomewide_features.tsv.gz`. Constant **within** a modality:
  **24** (multiome) vs **18** (scATAC) — confirmed by reading the actual output headers on disk
  (`igvf10_multiome/thp1_1/genomewide_features.tsv.gz` has 24 columns; the scATAC equivalent has
  18). The 6-column difference is exactly `normalizedATAC_enh, RNA_meanLogNorm,
  RNA_pseudobulkTPM, RNA_percentCellsDetected, Kendall, ARC.E2G.Score` — the RNA/ARC-derived
  features that don't exist on the scATAC path.
- `feat_bytes` — compressed bytes of `genomewide_features.tsv.gz` (495–605 MB multiome,
  189–223 MB scATAC). Used below as an I/O cross-check, not as the primary size variable (the
  code's control flow is expressed in rows and columns, not bytes).

No cross term with `k_genes` or any other Tier 1/3 variable — this script never touches genes,
cell counts, or fragments; it only reshapes/imputes/writes the feature table it's handed.

## 4. Asymptotic class, derived from code

**`O(n_cand_pairs · n_feat_cols)`**, i.e. linear in the number of table cells. Justification,
line by line:

- `fread` at `gen_final_features.R:15` — must deserialize every cell: `Θ(n_cand_pairs · n_feat_cols)`.
- `replace_na` at `gen_final_features.R:42` — full pass over ~half of `n_feat_cols` columns,
  each of length `n_cand_pairs`: `Θ(n_cand_pairs · n_feat_cols)` up to the constant ½.
- The dead `select()` copy at `gen_final_features.R:44-51` — re-copies all `n_cand_pairs ×
  n_feat_cols` cells (confirmed the selected column set has the same width as `df`, see §5):
  another `Θ(n_cand_pairs · n_feat_cols)` pass, wasted.
- `fwrite`/gzip at `gen_final_features.R:54` — serializes and compresses all cells:
  `Θ(n_cand_pairs · n_feat_cols)`.
- The interaction-term loop (`gen_final_features.R:29-32`) and the mean-fill branch
  (`get_fill_values.R:15-27`) are *also* linear in `n_cand_pairs` when active (each is a single
  vectorized pass per affected column, no nesting), so even in the worst case where every one
  of the config's features required an interaction term and a mean fill, the class would still
  be `O(n_cand_pairs · n_feat_cols)` — just with a larger constant. **There is no join, sort, or
  nested loop over rows anywhere in either script.** Both scripts only ever index columns by
  name or apply vectorized per-column operations; nothing indexes rows against other rows.

This is one of the more code-certain classes among Phase 1 units: not just "which exponent" but
"which class" is settled by inspection — there is no candidate superlinear mechanism in the code
at all, so I am not choosing between plausible classes here (contrast with rules whose class is
genuinely ambiguous between, say, O(n) and O(n log n)).

## 5. Verification of the "dead code" and "dormant branch" claims

Verified directly rather than inferred:

- **Dead `select()`**: `ActivityOnly_plus_external_features.tsv.gz` (rule input) and
  `genomewide_features.tsv.gz` (rule output) have the **same column count** for every cluster
  (24/24 multiome, 18/18 scATAC — checked via header line counts on disk), and their on-disk
  gzip sizes differ by <0.2% per cluster (e.g. `thp1_1` multiome: input 494,468,191 B, output
  494,990,740 B). The rename step changes names, not shape; `output` (the reordered/selected
  copy) is structurally redundant with `df` and is simply never referenced after line 51.
- **Dormant interaction terms**: `cut -f3 feature_table.tsv | tail -n +2` (the `second_input`
  column) is blank for **every row, every cluster, every modality** (12/12 jobs checked).
  `nrow(intx)` is 0 in every measured run.
- **Dormant mean-fill**: `cut -f5 feature_table.tsv` (`fill_value`) is `"0"` for every row,
  every cluster, every modality (12/12 jobs checked) — never `"mean"`.
- **Config size**: `feature_table.tsv` has 12 feature rows in multiome, 6 in scATAC, for every
  cluster (`k_config`), against `n_feat_cols` of 24/18 — i.e. roughly half of the output columns
  are "core" pass-through columns (`chr, start, end, name, class, TargetGene, ...`) untouched by
  `replace_na`, and half are config-driven features that are.

## 6. Benchmark cross-check

Pulled from `/scratch/users/kaybrand/scE2G_optimize_results/analysis/benchmarks_{multiome,scatac}.tsv`
and `size_manifest_{multiome,scatac}.tsv`.

| cluster | modality | wall s | cpu_time s | cpu/wall | io_in MB | io_out MB | n_cand_pairs | feat_bytes MB |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| jurkat | multiome | 129.22 | 93.57 | 0.72 | 298.7 | 3104.7 | 11,036,283 | 568 |
| jurkat_pma_cd3_4hr | multiome | 126.01 | 93.12 | 0.74 | 298.7 | 3186.1 | 10,854,698 | 558 |
| k562_crispri | multiome | 125.86 | 101.71 | 0.81 | 245.8 | 3319.7 | 11,447,564 | 605 |
| telohaec_crispri | multiome | 114.00 | 105.03 | 0.92 | 240.0 | 3338.8 | 11,106,687 | 579 |
| thp1_1 | multiome | 107.96 | 102.68 | 0.95 | 244.9 | 2941.7 | 10,173,077 | 495 |
| thp1_2 | multiome | 104.13 | 75.09 | 0.72 | 240.0 | 2823.3 | 10,278,325 | 505 |
| jurkat | scATAC | 48.91 | 45.21 | 0.92 | 3.8 | 1811.0 | 11,036,283 | 208 |
| jurkat_pma_cd3_4hr | scATAC | 48.89 | 45.11 | 0.92 | 3.9 | 1908.8 | 10,854,698 | 208 |
| k562_crispri | scATAC | 51.09 | 45.26 | 0.89 | 212.6 | 1940.8 | 11,447,564 | 223 |
| telohaec_crispri | scATAC | 49.29 | 44.99 | 0.91 | 202.2 | 1927.8 | 11,106,687 | 212 |
| thp1_1 | scATAC | 48.36 | 43.13 | 0.89 | 213.4 | 1677.8 | 10,173,077 | 189 |
| thp1_2 | scATAC | 48.43 | 43.54 | 0.90 | 183.9 | 1695.3 | 10,278,325 | 192 |

**Raw ranking against `n_cand_pairs` is not monotonic within multiome**: sorted by
`n_cand_pairs`, `thp1_1` (10.17M, 107.96 s) < `thp1_2` (10.28M, 104.13 s) — time *drops* as `n`
rises — and `telohaec_crispri` (11.11M, 114.00 s) < `jurkat_pma_cd3_4hr` (10.85M, 126.01 s) —
time drops again despite higher `n`. Taken at face value this looks like disagreement with a
linear-in-`n` model. But `n_cand_pairs` only spans **1.13×** across these clusters while raw
wall time spans **1.24×** (104.13–129.22 s) — a comparable range — so a small amount of
per-cluster noise (Slurm scheduling variance, filesystem/cache state, gzip time depending on
data entropy) is enough to scramble the rank order without contradicting a linear model.

**Normalizing removes most of the apparent disagreement.** Computing `wall_s / (n_cand_pairs ×
n_feat_cols)` — the implied constant of a pure `O(n·w)` model — gives:

- multiome: 4.22×10⁻⁷ – 4.88×10⁻⁷ s/cell (spread **1.16×**, i.e. tighter than the 1.13×–1.24×
  spread in the raw inputs)
- scATAC: 2.46×10⁻⁷ – 2.64×10⁻⁷ s/cell (spread **1.07×**, even tighter)

A near-constant per-cell rate at this precision is consistent with the code-derived linear class
and is the strongest evidence available from six points — but six points spanning ~1.13× with a
noise band of similar size (1.07–1.16×) **cannot distinguish O(n) from, e.g., O(n log n)** at
this `n` range; log-linear terms would look identical to noise here. This is exactly the
"near-constant variable, don't fit an exponent" situation the briefing warns about — I am
reporting a class from the code (§4), and separately noting the benchmark can only confirm that
the data are *not inconsistent* with it, not that it *proves* linearity over any alternative
that is also ≈linear at this scale.

**Cross-modality comparison (the useful n=2-per-cluster check) supports "I/O/serialization-
dominated," as the briefing's addendum hypothesized, but with a nuance**: `n_cand_pairs` is
literally identical between modalities per cluster, so the multiome/scATAC ratio isolates the
effect of `n_feat_cols` (24 vs 18, a 1.33× ratio) and column composition. The observed wall-time
ratio per cluster is 2.15×–2.64×, and the `feat_bytes` ratio is 2.62×–2.74× — much closer to each
other than either is to the raw column-count ratio (1.33×). **This says the constant scales with
bytes moved, not with column count alone**: the 6 multiome-only columns (`normalizedATAC_enh,
RNA_meanLogNorm, RNA_pseudobulkTPM, RNA_percentCellsDetected, Kendall, ARC.E2G.Score`) are
floating-point values with long decimal representations, so they cost disproportionately more
per cell to parse/format/compress than the 18 shared columns (mostly short strings, ints, or
booleans). Directly measured average bytes per (row×col) cell: multiome ≈10.6 B/cell, scATAC
≈8.6 B/cell — the extra columns' implied average width (≈16.7 B/cell) is nearly double the
shared-column average. **The class is unchanged (still `O(n_cand_pairs · n_feat_cols)`), but the
constant is column-composition-dependent, not just column-count-dependent** — this is a
constant-level effect, not a different exponent, and I flag it as such rather than reify it into
an `O(n·w·bytes_per_cell)` "class."

**CPU/wall ratios differ systematically by modality**: multiome averages ≈0.81 (range
0.72–0.95), scATAC averages ≈0.91 (range 0.89–0.92) — multiome jobs spend proportionally more
wall time not running CPU, consistent with more bytes moved causing more I/O/gzip wait relative
to compute. This is a modality-level pattern (n=6 pairs), not a single-job artifact.

## 7. Concrete overhead sources

1. **`io_in` is measured far below the actual compressed input size and should not be read as
   "bytes read from disk."** Multiome `io_in` is 240–299 MB, but the actual on-disk
   `ActivityOnly_plus_external_features.tsv.gz` is 494–605 MB (checked via `stat -c%s`) — i.e.
   `io_in` is roughly half the file size. The most likely explanation is page-cache locality:
   this file was written by the immediately preceding rule (`add_external_features`) moments
   before, so a large fraction of the read is served from cache rather than the block device,
   and psutil's I/O counters (which the benchmark reports) only count actual block-device reads.
   **Don't use `io_in` here as a proxy for data volume; it understates it.**
2. **`io_out` is ~5.5–8.8× the final compressed file size, and matches the *uncompressed* text
   volume almost exactly.** Directly measured for `thp1_1`: `genomewide_features.tsv.gz`
   decompresses to 2,597,006,553 bytes (multiome) / 1,578,123,187 bytes (scATAC), against a
   compressed size of 494,990,740 / 189,044,508 bytes — implied compression ratios of 5.25× and
   8.35×, which line up closely with the `io_out`-to-`feat_bytes` ratios of 5.5–5.9× (multiome)
   and 8.7–9.2× (scATAC) measured independently in the benchmark table. **`fwrite`'s gzip
   compression cost (and the I/O it generates) tracks the uncompressed byte volume, not the
   `feat_bytes` compressed size** — using `feat_bytes` alone to estimate this rule's write-side
   cost understates it by 5–9×, and understates it *more* for scATAC than multiome (because
   scATAC's data compresses better, being narrower and less floating-point-heavy).
3. **A full-table copy is computed and discarded** (`gen_final_features.R:44-51`, `output <-
   select(...)`, never used — line 54 writes `df`). One wasted `O(n_cand_pairs · n_feat_cols)`
   pass and a transient ~500 MB (multiome) / ~190 MB (scATAC) in-memory copy, every job, in both
   modalities.
4. **Three to four effectively-full-size passes per job**: `fread` (read+parse), `replace_na`
   (scan ~half the columns), the dead `select()` copy (all columns), `fwrite`+gzip (all
   columns). None of these are algorithmically necessary to run back-to-back — they are separate
   language-level operations that each re-touch the whole table.
5. **No parallelism requested or used.** The rule declares no `threads:` (runs at 1, per the
   already-established cross-cutting finding). `cpu_time` never exceeds wall time in any of the
   12 jobs, consistent with single-threaded execution throughout (fread/fwrite/gzip in `R`'s
   `data.table` do not multi-thread compression by default here).
6. **Two dormant code branches add zero measured cost but are real latent overhead sources**:
   the interaction-term loop (`gen_final_features.R:29-32`) and the mean-fill computation
   (`get_fill_values.R:15-27`). Neither executes in any of Phase 0's 12 configs (verified in
   §5), so none of the measured wall time includes them — if a future feature-table config used
   `fill_value == "mean"` or a non-blank `second_input`, this rule's cost would increase by
   additional `O(n_cand_pairs)`-per-affected-column passes, still within the same overall class.

## 8. Formulas

Let `n = n_cand_pairs`, `w = n_feat_cols` (constant per modality: 24 multiome, 18 scATAC).
Fitted per-cell rate `k = wall_s / (n · w)`, averaged/bounded across the 6 measured jobs per
modality (§6):

**Multiome** (`w = 24`): `k_est = 4.5365e-7`, `k_low = 4.220e-7`, `k_high = 4.879e-7` s/cell

```
t_est_multiome(n)  = 1.0888e-5 * n     [seconds]
t_low_multiome(n)  = 1.0128e-5 * n
t_high_multiome(n) = 1.1710e-5 * n
```

**scATAC** (`w = 18`): `k_est = 2.528e-7`, `k_low = 2.463e-7`, `k_high = 2.641e-7` s/cell

```
t_est_scatac(n)  = 4.5504e-6 * n
t_low_scatac(n)  = 4.4334e-6 * n
t_high_scatac(n) = 4.7538e-6 * n
```

General cross-modality form, same class, constant absorbs `w` and the modality-specific average
byte-width per cell (§6): `t(n, w) ≈ k(w) · n · w`, where `k(w)` is **not** a single
modality-independent constant — it differs ~1.8× between modalities at equal `n · w` because the
multiome-only columns are wider on average (§6) — so I give the two modality-specific linear
forms above rather than one combined formula; forcing a single `k` would hide a real,
measured, constant-level effect.

**These bounds cover only measurement noise around the branches actually exercised in Phase 0**
(no interaction terms, no mean-fill). They do **not** bound the dormant-branch worst case. If
exercised, add (order-of-magnitude, unmeasured, using the same per-cell rate order as the
established passes since both are single vectorized column operations):

```
+ p_intx  * k_est_modality * n     # interaction-term loop, gen_final_features.R:29-32
                                    # p_intx = number of config rows with non-blank second_input
+ m_mean  * k_est_modality * n     # mean-fill full-column pass, get_fill_values.R:21-23
                                    # m_mean = number of config rows with fill_value == "mean"
```

both currently `p_intx = m_mean = 0` in every measured job (§5).

## 9. Evaluated at mean/min/max cluster size

`n_cand_pairs`: mean of the 6 clusters = 10,816,106; `thp1_1` (designated min) = 10,173,077;
`telohaec_crispri` (designated max) = 11,106,687. (Note: `telohaec_crispri` is not the single
largest cluster by `n_cand_pairs` — `k562_crispri` is slightly larger, 11,447,564 — but I use the
designated min/max clusters per the task instructions rather than re-picking by this specific
variable.)

| cluster (n) | modality | t_low (s) | t_est (s) | t_high (s) | actual observed (s) |
|---|---|---:|---:|---:|---:|
| mean (10,816,106) | multiome | 109.55 | 117.77 | 126.66 | — (not a real cluster) |
| thp1_1 min (10,173,077) | multiome | 103.05 | 110.76 | 119.08 | 107.96 |
| telohaec_crispri max (11,106,687) | multiome | 112.50 | 120.93 | 130.06 | 114.00 |
| mean (10,816,106) | scATAC | 47.97 | 49.22 | 51.42 | — (not a real cluster) |
| thp1_1 min (10,173,077) | scATAC | 45.09 | 46.29 | 48.35 | 48.36 |
| telohaec_crispri max (11,106,687) | scATAC | 49.25 | 50.55 | 52.82 | 49.29 |

Both real-cluster actuals fall inside `[t_low, t_high]` in both modalities, and are close to
`t_est`, as expected — these constants were fit from the same 6 points, so this is a consistency
check on the linear form, not an independent validation.

## 10. Does the class differ between modalities, or only the constant?

**Only the constant.** Both scripts run the identical control flow regardless of modality — the
same `fread` → rename → `get_fill_values` → `replace_na` → dead `select()` → `fwrite` sequence,
with no modality-conditional branching in either script. What differs:

- `n_feat_cols`: 24 (multiome) vs 18 (scATAC) — a structural constant of the class, `w`.
- The average byte-width per cell, because the 6 multiome-only columns are floating-point with
  long decimal representations, inflating the effective per-cell constant ~1.8× beyond what the
  1.33× column-count ratio alone would predict (§6).
- `k_config` (12 vs 6 feature-table rows) — trivial, doesn't touch `n_cand_pairs`-scaled work.

No new loop, join, sort, or algorithmic step appears or disappears between modalities. The
~2×–2.6× time ratio the briefing flagged for this rule (§ briefing table: 1.7–2.2 min multiome
vs 0.8–0.9 min scATAC) is explained almost entirely by this constant-level width/byte effect at
essentially identical `n_cand_pairs`, not by a different complexity class.

## 11. What the measurement cannot tell us

- **Cannot confirm the exponent on `n_cand_pairs` independently of noise.** The variable spans
  only 1.13× across the 6 clusters, and the per-cell-rate noise band (1.07×–1.16×) is comparable
  in size. A linear model is consistent with the data and is what the code supports, but six
  points this close together cannot rule out a mild superlinear or sublinear correction at this
  scale — only that if one exists, it's small enough to hide inside measurement noise here.
- **The interaction-term loop and mean-fill branch are entirely unmeasured.** All 12 Phase 0
  configs have blank `second_input` and `fill_value == "0"` for every feature (§5), so neither
  `gen_final_features.R:29-32` nor `get_fill_values.R:15-27` executes even once in this dataset.
  Their cost estimates in §8 are order-of-magnitude extrapolations from the *other* per-column
  passes in the same scripts, not measurements.
- **`io_in` cannot be used to infer this rule's true read volume** — it's suppressed by page-
  cache locality from the immediately preceding rule's write, and the benchmark table gives no
  way to separate "read from cache" from "read from disk" bytes.
- **`io_out` cannot be attributed to a specific line without deeper profiling** — it reflects the
  combined effect of `fwrite`'s text serialization and gzip compression on the uncompressed byte
  stream; the benchmark has no per-line breakdown, so the individual contribution of the dead
  `select()` copy at `gen_final_features.R:44-51` vs. the `fwrite` call at line 54 cannot be
  separated from the single wall-clock total. I inferred the dead copy exists and roughly what
  it costs (§7.3) from reading the code, not from isolating it in the timing data.
- **n=1 per cluster** (no replicates, per the briefing) — the ~1.24× multiome wall-time spread
  that doesn't rank-order with `n_cand_pairs` (§6) could reflect genuine per-cluster data
  differences (NA density, float-value entropy affecting gzip time) rather than pure scheduling
  noise, and this dataset cannot distinguish the two explanations.
- **`feat_bytes` alone does not fully explain the multiome/scATAC time ratio** — it tracks it
  more closely than raw column count does, but the six-point comparison isn't precise enough to
  rule out some additional non-byte-proportional constant (e.g. gzip's compression effort scaling
  nonlinearly with entropy rather than pure byte count).

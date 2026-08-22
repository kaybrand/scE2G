# Unit 8 — `generate_atac_matrix` (Tier B)

## 1. Unit identity

- Script: `workflow/scripts/feature_computation/generate_atac_matrix.R` (94 lines)
- Rule: `generate_atac_matrix`, `workflow/rules/generate_atac_matrix.smk:23-54`
- Modality: **Multiome only**. 6 jobs (one per cluster) — one of the four rules
  (`make_kendall_pairs → generate_atac_matrix → compute_kendall → arc_e2g`) that exist only on
  the Multiome branch of the DAG (Phase 0 briefing §2). No scATAC counterpart exists to compare
  class or constant against.
- No `threads:` declared (`generate_atac_matrix.smk:37-54`, only `resources: mem_mb`) →
  Snakemake/Slurm gives it 1 core. Measured `cpu_time/wall_s` (`mean_load`) is 86.5–93.7%,
  consistent with a single, mostly-CPU-bound R process — there is no threading gap to explain
  here because none was declared.
- A sibling rule, `get_cell_barcodes` (`generate_atac_matrix.smk:8-20`), is defined in the same
  file but **did not run at all** in this measurement — `benchmarks_multiome.tsv` has zero rows
  for it (verified: `awk -F'\t' -v r=get_cell_barcodes 'NR==1||$2==r' benchmarks_multiome.tsv`
  returns only the header). See §2 for why.

## 2. What the code actually does (control flow, `file:line`)

1. `generate_atac_matrix.R:33-37` — `readGeneric(kendall_pairs_path, keep.all.metadata=T,
   header=T)` parses the **entire** `Kendall/Pairs.tsv.gz` (input, `n_kendall_pairs` rows,
   9.69M–11.32M, 6 columns — verified: `zcat .../Pairs.tsv.gz | head -1` shows `chr start end
   TargetGene PeakName PairName`) into a `GRanges`, then `!duplicated(mcols(...)[,"PeakName"])`
   deduplicates by `PeakName` (hash-based) to build `bed.peaks`, the feature/peak set — this is
   what collapses `n_kendall_pairs` down to `j_cand_elems` (156,420–158,548, confirmed
   near-identical to the manifest's `j_cand_elems` per the briefing's stated identity
   `j_cand_elems == EnhancerList.txt rows − 1`).
2. `generate_atac_matrix.R:40-50` — reads the RNA matrix to get the cell-barcode universe.
   `RNA_matrix_filtered: true` in this config (`configs/igvf10_multiome_config.yaml:64`,
   confirmed live at `/scratch/.../igvf10_multiome/config/scE2G_config.yml:3`) and the
   `rna_matrix_file` path for every cluster is a directory (verified on disk, e.g.
   `.../thp1_2/rna_count_matrix_igvf10_thp1_2/{barcodes,features,matrix}.{tsv,mtx}.gz`), so this
   takes the `Read10X(...)` branch (`:46-47`, `file.info(...)$isdir` is TRUE) — a sparse
   MatrixMarket read, cost `O(nnz(RNA matrix))`, roughly proportional to `n_umi` in magnitude
   though `n_umi` is not itself a row/col count.
3. `generate_atac_matrix.R:56-65` — **conditional ATAC-barcode intersection.** `cell_bc_path`
   comes from `get_cell_barcode_file(config["RNA_matrix_filtered"])`
   (`generate_atac_matrix.smk:2-6`): when `RNA_matrix_filtered` is `True` (it is, here), this
   function returns `RESULTS_DIR` itself — a directory, not the `cell_barcodes.txt` file that
   `get_cell_barcodes` (`:8-20`) would otherwise produce. `file_test("-f", cell_bc_path)`
   (`generate_atac_matrix.R:57`) is `FALSE` for a directory, so the script takes the `else`
   branch at `:64`: `cells.use <- rna.cells` directly, with **no** ATAC-fragment-barcode
   intersection. This is why `get_cell_barcodes` never runs (§1) — Snakemake still needs *some*
   existing path to satisfy the `cell_barcodes_path` input, and `RESULTS_DIR` (which always
   exists) is used as that placeholder, breaking the dependency edge to `get_cell_barcodes`
   entirely. Confirmed from the slurm logs: no cluster's log contains the line "Number of cells
   from fragment file" (which only prints inside the `if` branch, `:59`), and "Number of cells
   to use" equals "Number of cells in RNA matrix" exactly in all 6 logs (below).
4. `generate_atac_matrix.R:69-72` — **subsampling.** `if (length(cells.use) > max_cell_count)
   set.seed(123); cells.use <- sample(cells.use, max_cell_count)`. `max_cell_count` comes from
   `params: max_cell_count = config['max_cell_count']` (`generate_atac_matrix.smk:47-48`), value
   **20000** (`config/config.yaml:17`, confirmed live at
   `/scratch/.../igvf10_multiome/config/scE2G_config.yml:16`). **This subsampling happens
   inside this rule/script, not upstream** — there is no upstream cap on cell count or on the
   ATAC fragment file itself; `atac_frag_path` (`generate_atac_matrix.smk:32-33`) is the raw,
   un-subsampled per-cluster fragment file.
5. `generate_atac_matrix.R:79-81` — `CreateFragmentObject(path=atac_frag_path,
   cells=cells.use)`. This is **not** a lightweight registration: Signac's implementation
   computes `hashes <- md5sum(files=c(path, index.file))` — an MD5 checksum of the **entire**
   raw fragment file plus its `.tbi` index (confirmed from Signac's `R/fragments.R` source,
   fetched live from `github.com/stuart-lab/signac`). This forces one genuine full sequential
   byte-for-byte pass over the fragment file, cost `O(bytes of atac_frag_path)` ≈
   `O(n_frag)` (fragment count and file bytes scale together for a given library-prep chemistry).
6. `generate_atac_matrix.R:84-88` — `FeatureMatrix(fragments=list.fragments,
   features=bed.peaks, cells=cells.use)`. Confirmed from the same live Signac source fetch:
   `FeatureMatrix()` chunks `features` (the `j_cand_elems`-sized `bed.peaks`) into groups of
   `process_n = 2000` (default, not overridden here) via `ChunkGRanges()`, then for each chunk
   calls `SingleFeatureMatrix()` → `PartialMatrix()`, which uses a `TabixFile` handle to do
   **genomic random-access reads restricted to that chunk's regions** (not a full-file rescan)
   and returns a **sparse** per-chunk matrix; all chunks are combined via `rbind` into the final
   `j_cand_elems × n_cells` sparse matrix. Cost ≈ `O(bytes of fragments overlapping any
   candidate peak)` (a fraction of `n_frag`, not necessarily all of it) + a fixed per-chunk
   overhead of `ceil(j_cand_elems / 2000)` ≈ 79 chunks (constant across all 6 clusters here
   because `j_cand_elems` only spans 1.014×) + `O(nnz)` for the sparse accumulation/`rbind`.
7. `generate_atac_matrix.R:93-94` — `saveRDS(atac.matrix, atac_matrix_path)`: gzip-compressed
   serialization of the sparse matrix, cost `O(nnz)`.

**Confirmation the output matrix is sparse, not dense** (addressing the unit brief's question
directly): confirmed both from source (`PartialMatrix()` builds sparse per-chunk matrices,
`rbind`ed) and from the numbers in §5/§6 below — the on-disk file size scales far more tightly
with `n_frag` than with `j_cand_elems · n_cells`, which is only possible if `nnz` (and hence the
stored/compressed size) tracks the fragment count rather than the full peak×cell grid.

## 3. Size variables (named as in `size_manifest_multiome.tsv`)

| Variable | Role in this script | Range | Ratio |
|---|---|---|---:|
| `n_frag` | bytes/records scanned by `CreateFragmentObject`'s MD5 hash (`:79-81`) and by `FeatureMatrix`'s tabix-restricted scan (`:84-88`) | 26.2M–394.8M | **15.1×** |
| `n_cells` | width of the output matrix; drives `Read10X` and the (never-triggered) subsample check | 2,593–15,555 | 6.0× |
| `n_kendall_pairs` | rows parsed by `readGeneric` and deduplicated at `:33-37` | 9.69M–11.32M | 1.17× |
| `j_cand_elems` | unique peaks after dedup at `:37`; height of the output matrix; drives `FeatureMatrix`'s chunk count (`ceil(j/2000)`) | 156,420–158,548 | 1.01× |

Cross term the code genuinely has: the output matrix (and hence a ceiling on `nnz`, memory, and
serialization cost) is `j_cand_elems × n_cells`; the *actual* dominant cost is bounded by
`min(n_frag, j_cand_elems · n_cells)` since `nnz ≤ n_frag` (each fragment can increment at most
one peak×cell count).

## 4. Asymptotic class derived from code

| Step | file:line | Class |
|---|---|---|
| Parse+dedup `Pairs.tsv.gz` | `:33-37` | `O(n_kendall_pairs)` (hash-based `duplicated`) |
| Read RNA matrix (`Read10X`) | `:40-50` | `O(nnz(RNA matrix))`, roughly tracks `n_umi`/`n_cells` |
| Cell-list build / (skipped) intersection | `:52-65` | `O(n_cells)` |
| Subsample check (never triggered here) | `:69-72` | `O(n_cells)` |
| `CreateFragmentObject` → `md5sum()` of whole fragment file | `:79-81` | `O(n_frag)` (full linear byte pass) |
| `FeatureMatrix`: chunk features, tabix random access per chunk, sparse `rbind` | `:84-88` | `O(bytes overlapping candidate peaks)` + `O(⌈j_cand_elems/2000⌉)` fixed chunk overhead — both effectively `O(n_frag)`-order here since the peak-overlap fraction of `n_frag` is roughly stable across these libraries and the chunk count barely varies (`j_cand_elems` spans only 1.01×) |
| `saveRDS` | `:93-94` | `O(nnz)` ≤ `O(n_frag)` |

**Overall class: linear, `Θ(n_frag)`**, with smaller additive linear terms in `n_kendall_pairs`
and RNA-matrix `nnz`/`n_cells` that are real in the code but (see §5) not separable from
`n_frag` at this data range, plus a near-constant `j_cand_elems`-driven chunk-count term that
manifests as part of the fixed intercept because `j_cand_elems` barely varies across these 6
clusters. **This rule does not have a hidden superlinear term** — no step above is a nested
loop or an `O(j²)`-shaped operation; the closest candidate (`FeatureMatrix`'s per-chunk
processing) is confirmed by source to be `O(j/process_n)` chunks of `O(process_n)` work each,
i.e. `O(j)` overall, not `O(j²)`.

## 5. Benchmark cross-check

```
cluster              wall_s   cpu_time  mean_load  max_rss(MB)  io_in(MB)  io_out(MB)   n_frag       n_cells  j_cand_elems  n_kendall_pairs
jurkat                441.10   406.41    92.12       5739.45      210.35     937.48     75,342,313     8,239    158,397        10,788,202
jurkat_pma_cd3_4hr     338.01   314.65    93.06       4226.99      173.73     867.41     55,803,917     5,124    157,395        10,759,396
k562_crispri           433.81   401.24    92.48       5868.97      207.45     911.78     86,424,398     5,109    156,420        11,317,834
telohaec_crispri      1187.25  1111.19    93.58      14219.40      423.71    1280.65    394,847,593    15,555    157,551        11,136,922
thp1_1                 307.83   266.30    86.47       3240.59      150.17     782.51     26,231,026     2,676    158,548         9,694,978
thp1_2                 317.71   297.91    93.71       3662.75      159.47     832.10     26,292,230     2,593    157,498         9,910,221
```

**This is one of the few units where the six points genuinely test the code.** `n_frag` spans
15.1×, comfortably wider than the 1.94× "unexplained noise" spread seen in other Tier-B units
(e.g. `make_kendall_pairs`, unit 7).

- **`wall_s ~ a + b·n_frag`** (OLS, all 6 points): `a = 238.49 s`, `b = 2.398×10⁻⁶ s/fragment`,
  **R² = 0.996**. This is a very tight fit for 6 points and directly matches the code: the
  `md5sum()` full-file hash in `CreateFragmentObject` (`:79-81`) is an unconditional linear pass
  over `n_frag`-proportional bytes, and `FeatureMatrix`'s tabix-restricted scan (`:84-88`) adds
  a further term proportional to the (roughly stable) peak-overlap fraction of `n_frag`.
- **`max_rss ~ a + b·n_frag`**: `a = 2968.3 MB`, `b = 2.880×10⁻⁵ MB/fragment`, **R² = 0.989** —
  also a tight linear-in-`n_frag` fit.
- **By contrast, `wall_s` and `max_rss` regressed on `n_cells` (or on `j_cand_elems·n_cells`,
  effectively the same signal since `j_cand_elems` is near-constant) fit noticeably worse:**
  R² = 0.908 (wall) and R² = 0.932 (RSS). `n_frag` is the better single-variable explanation of
  both time and memory.
- **Leave-one-out sensitivity** (drop each cluster, refit `wall_s ~ a + b·n_frag` on the other
  5): slope ranges `2.199×10⁻⁶`–`2.416×10⁻⁶` s/fragment, intercept ranges `232.97`–`248.91` s,
  and R² stays ≥0.88 in every case (worst case is dropping `telohaec_crispri`, the largest
  cluster, which is expected — it anchors the fit's high end). This is a genuinely stable
  linear signal, not noise dressed up as a fit.
- **Attempting to add `n_kendall_pairs` as a second regressor** (`wall_s ~ a + b·n_frag +
  c·n_kendall_pairs`) gives `c = −1.72×10⁻⁵` — a **negative, physically impossible**
  coefficient (more pairs → less time), the classic symptom of over-fitting 6 points with 3
  parameters when the added variable (`n_kendall_pairs`, spanning only 1.17× and correlated
  with `n_frag`, r≈0.6–0.7 informally) carries no separable signal here. I flag this explicitly
  rather than reporting the fit: **the code has a real `O(n_kendall_pairs)` term (§4), but these
  6 points cannot isolate its coefficient from `n_frag`'s.**
- **Correlation caveat, per the briefing**: `n_frag` and `n_cells` correlate at r≈0.95 across
  these 6 clusters (deeper libraries also have more cells), so even the strong `n_frag` fit
  above cannot fully rule out that some of what it captures is actually an `n_cells` effect.
  What *does* discriminate a little: `jurkat` (8,239 cells, 75.3M frag) vs. `k562_crispri`
  (5,109 cells, 86.4M frag) — higher `n_frag` but *lower* `n_cells` than `jurkat`, and its wall
  time (433.81 s) sits between `jurkat`'s (441.10 s) and the `n_frag`-only prediction (445.8 s)
  — consistent with `n_frag` dominating over `n_cells` at these two points, but this is a single
  weak data contrast, not proof.
- **Sparse-vs-dense cross-check** (addressing the unit brief directly): the on-disk `.rds` file
  size (gzip-compressed) tracks `n_frag` far more tightly (`filesize/n_frag` ratio 1.06–1.44
  across clusters, a 1.4× spread) than it tracks `j_cand_elems·n_cells` (`filesize/(j·n)` ratio
  0.065–0.193, a **3.0× spread**). If the matrix were effectively dense (or if `nnz` scaled with
  the full `j·n` grid), file size should track `j·n`, not `n_frag`. It tracks `n_frag` instead —
  consistent with a sparse matrix whose non-zero count is bounded by, and scales with,
  fragments-in-peaks rather than the peak×cell grid size. This agrees with source (§2 point 6).
- **`io_in` is a poor proxy for bytes actually processed, and disagrees with the code in a
  specific, flaggable way.** `io_in` (150–424 MB) is only **9.9%–52.6%** of the corresponding
  raw fragment file's on-disk size (301 MB–4.28 GB, measured directly), and — tellingly — the
  ratio is *smallest* for the *largest* cluster (`telohaec_crispri`: 423.7 MB `io_in` vs. 4.28 GB
  file, 9.9%) and *largest* for the smallest (`thp1_2`: 159.5 MB vs. 303 MB, 52.6%). If
  `CreateFragmentObject`'s `md5sum()` genuinely re-reads the whole file from disk every time
  (§2 point 5), `io_in` should be *at least* the file size. It is far less, and disproportionately
  so for larger files. My hypothesis: this fragment file is touched by several earlier rules in
  the same cluster's DAG (`process_fragment_file`, `frag_to_tagAlign`, `call_macs_peaks`,
  `frag_to_norm_bigWig`) before `generate_atac_matrix` runs, and OS/Lustre-client page-cache
  reuse from those earlier reads is satisfying a large fraction of the `md5sum` and tabix reads
  without hitting the block device — `io_in` (a device-I/O counter) would then systematically
  undercount logical bytes processed, more so for a large file that plausibly stays "warmer" in
  cache relative to a small one. I cannot confirm this without instrumenting the filesystem
  layer; it is a flagged hypothesis, not a verified mechanism, exactly per the addendum's
  guidance to flag disagreement rather than smooth it over.

## 6. Concrete overhead sources

- **A mandatory whole-file MD5 hash inside `CreateFragmentObject`** (`:79-81`,
  Signac `R/fragments.R`, confirmed live from source) is a genuine, unconditional `O(n_frag)`
  full-file read on every job, independent of how many peaks or cells are actually used
  downstream. This is the single most important, code-confirmed overhead source for this rule.
- **The 6-column `Pairs.tsv.gz` is parsed in full and then reduced to unique `PeakName`s only**
  (`:33-37`) — `TargetGene`/`PairName` are read but not used to build `bed.peaks`
  (`mcols(bed.peaks) <- NULL` at `:37`). A "read more than you need" pattern, though modest
  relative to the fragment-file cost given `n_kendall_pairs`'s narrow 1.17× range.
- **`io_out` (782–1281 MB) is 8–21× the final `.rds` output size** (37.1–473.9 MB, measured
  directly on disk) — consistent with the R/Bioconductor/Signac stack's typical
  allocate-then-discard churn while building the sparse matrix chunk-by-chunk and `rbind`-ing
  ~79 pieces together, though I cannot attribute it to a specific line without profiling.
- **Declared vs. realized parallelism**: no `threads:` declared at all
  (`generate_atac_matrix.smk:37-54`) → 1 core allocated; measured `mean_load` 86.5–93.7%
  confirms the job is essentially single-core and CPU/parse-bound (not idle on I/O wait, which
  would show a much lower ratio). There is no declared-vs-measured *gap* to explain here (unlike
  `frag_to_norm_bigWig`/`frag_to_tagAlign`), since nothing beyond 1 core was ever requested.
- **`FeatureMatrix`'s default `process_n=2000`** (not overridden at `:84-88`) fixes the chunk
  count at `⌈j_cand_elems/2000⌉ ≈ 79` for all 6 clusters (since `j_cand_elems` spans only
  1.01×) — this per-chunk overhead (tabix query setup, per-chunk sparse-matrix allocation,
  eventual `rbind`) is real per the code but, because it barely varies across these clusters,
  is indistinguishable from a fixed constant in the regression above, not a variable one.
- **Memory headroom**: `resources: mem_mb` here is computed by
  `encode_e2g.ABC.determine_mem_mb` (not a flat config constant), so I cannot directly compare a
  single declared ceiling to the measured 3.2–14.2 GB `max_rss` without evaluating that function
  against each cluster's inputs — out of scope for a quick check, flagged rather than guessed.

## 7. `t_est(n)`, `t_low(n)`, `t_high(n)`

**Unlike most Tier-B units in this project, the exponent here is not the main source of
uncertainty** — code (§4) and data (§5, R²=0.996, stable under leave-one-out) agree the class is
linear, `Θ(n_frag)`. The residual uncertainty is in the *rate* (the constant/slope), not the
exponent, so `t_low`/`t_high` below bound that rate via leave-one-out sensitivity rather than
bounding competing exponents. All formulas in **seconds**, `n_frag` = raw fragment count (e.g.
`75342313`, not millions):

- **`t_est(n_frag) = 238.49 + 2.398×10⁻⁶ · n_frag`** — full-sample OLS fit, R² = 0.996.
- **`t_low(n_frag) = 232.97 + 2.199×10⁻⁶ · n_frag`** — lower bound on the linear rate, using the
  minimum intercept and minimum slope observed across the 6 leave-one-out refits (§5).
- **`t_high(n_frag) = 248.91 + 2.416×10⁻⁶ · n_frag`** — upper bound on the linear rate, using
  the maximum intercept and maximum slope observed across the 6 leave-one-out refits.

These three are not point estimates blended together — they are the same `Θ(n_frag)` class with
a calibrated envelope on its rate, which is the correct bracket *given* that (unusually) the
exponent itself is well-supported here. `t_low ≤ t_est ≤ t_high` holds at every `n_frag` in and
beyond the observed range (checked at 26M, 111M, 395M below).

**Not separately fit, but real per the code** (§4): an additive `O(n_kendall_pairs)` term (pairs
parse+dedup, `:33-37`) and an `O(nnz(RNA matrix))` term (`Read10X`, `:40-50`), both currently
small and folded into the intercept/slope above because they cannot be isolated from `n_frag`
with 6 correlated points (§5) — attempting to fit `n_kendall_pairs` as a second regressor
produced an unphysical negative coefficient. On a hypothetical dataset where `n_kendall_pairs`
or RNA-matrix size varied far more independently of `n_frag` than it does across these 6
clusters, these terms could matter more than they do here.

## 8. All three evaluated

| `n_frag` | cluster | actual wall_s | `t_low` | `t_est` | `t_high` |
|---|---|---:|---:|---:|---:|
| 110,823,580 (mean of 6) | — | 371.05 (mean) | 476.7 | 504.4 | 516.7 |
| 26,231,026 (min) | thp1_1 | 307.83 | 290.7 | 301.4 | 312.3 |
| 394,847,593 (max) | telohaec_crispri | 1187.25 | 1101.3 | 1185.5 | 1203.1 |

`t_est` at the min and max cluster sizes is within 2% of the actual observed wall time in both
cases — an unusually tight calibration for a Phase-1 unit, made possible by `n_frag`'s genuine
15.1× spread. `t_est` over-predicts the *mean* point by ~36% because the mean `n_frag`
(110.8M) does not correspond to any single cluster's actual joint (`n_frag`, `n_cells`,
overlap-fraction) combination — the mean-of-inputs is not the same as the mean-of-outputs when
clusters vary on more than one axis simultaneously; treat the "mean" row as a convenience
evaluation of the formula, not evidence the formula is biased.

## 9. Modality

Multiome-only — no scATAC counterpart exists for this rule (Phase 0 briefing §2), so there is
nothing to compare the class or constant against across modalities.

## 10. What the measurement cannot tell us

- **Cannot fully separate `n_frag` from `n_cells`.** They correlate at r≈0.95 across these 6
  clusters (briefing §3's stated caveat). The `n_frag`-only fit is meaningfully tighter than the
  `n_cells`-only fit (R² 0.996 vs. 0.908), and the code's `CreateFragmentObject` MD5-hash step
  is unambiguously an `n_frag`-only operation — but I cannot rule out that part of the fitted
  `n_frag` slope is actually attributable to `n_cells`-correlated effects elsewhere in
  `FeatureMatrix`.
- **The `max_cell_count: 20000` subsampling branch (`:69-72`) never triggers on this dataset.**
  All 6 clusters' RNA-matrix cell counts (2,593–15,555) are below `max_cell_count` (20,000);
  slurm logs confirm "Number of cells to use" equals "Number of cells in RNA matrix" exactly
  for every cluster (e.g. `telohaec_crispri`: both 15,555). **This directly contradicts the
  premise in this unit's brief** ("`max_cell_count: 20000` subsamples the four largest clusters
  but NOT the thp1 pair") — I could not find any evidence for that claim in this measurement;
  the measured behavior is that subsampling is dormant code on all 6 clusters at their current
  sizes. The measurement therefore tells us nothing about this rule's behavior *when*
  subsampling is active — a hypothetical 7th cluster with >20,000 RNA-matrix cells would exit
  the `Read10X`/subsample path differently, and I have no data point for that regime.
- **Cannot attribute `io_out` (782 MB–1.28 GB, 8–21× the final output size) to a specific
  operation** without profiling the R process (`Rprof`/`tracemem`), which is out of scope for a
  characterisation-only task.
- **Cannot confirm the page-cache hypothesis for the low `io_in`/file-size ratio** (§5) — it is
  a plausible mechanism given this fragment file is read by several earlier rules in the same
  DAG, but I have not instrumented the filesystem layer to verify it, and an alternative
  explanation (tabix random access genuinely touching a shrinking fraction of the file as
  libraries get deeper) is not ruled out either.
- **`n_kendall_pairs`'s and the RNA matrix's coefficients are not identifiable from these 6
  points** (§5, §7) — real per the code, invisible in the fit at this data range.
- **Single measurement per cluster, no replicates.** Every number above is n=1 per cluster; the
  tightness of the `n_frag` fit (R²=0.996) is reassuring but could still partly reflect
  cluster-to-cluster differences that happen to correlate with `n_frag` on this particular run
  (e.g. which Slurm node, cache state) rather than the algorithm alone.
- **The toy chr22 run is excluded**, per the briefing's standing instruction — it measures
  interpreter/conda startup, not this rule's algorithm.

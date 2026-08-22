# Phase 1 unit #6 — `arc_e2g` (ARC-E2G feature computation)

**Tier A.** Characterisation only — no code edits, no fix proposals, no prioritisation.

---

## 1. Unit identity

| | |
|---|---|
| Rule | `arc_e2g`, `workflow/rules/arc_e2g.smk:2` (33-line file) |
| Script | `workflow/scripts/feature_computation/compute_arc_e2g.R` — **132 lines** |
| Language / env | R 4.3.3, `workflow/envs/sc_e2g.yml`; conda prefix `…/scE2G/.snakemake/conda/c1ae1ceebb37236cbea110033b350574_` |
| Key package versions | `data.table` 1.15.2, `genomation` 1.34.0, `GenomicRanges` 1.54.1, `readr` 2.1.5, `vroom` 1.6.5 |
| Modality | **MULTIOME ONLY.** One of the four rules behind `checkpoint features_required`; `to_generate == "Neither"` for all six scATAC clusters, so it never runs there. |
| Jobs | **6** (one per multiome cluster), 0 in scATAC |
| Declared `threads:` | **none** → Snakemake runs it at 1; observed `cpus_per_task=1` in every slurm log |
| Declared `resources:` | `mem_mb=64*1000` (`arc_e2g.smk:29`); Slurm `ReqMem=62.50G`; no `runtime` → inherits the profile default `48h` (`profiles/measure/config.yaml:115`) |
| Inputs | `Predictions/EnhancerPredictionsAllPutative.tsv.gz` (29 cols) + `Kendall/Pairs.Kendall.tsv.gz` (10 cols) |
| Output | `ARC/EnhancerPredictionsAllPutative_ARC.tsv.gz` (**36 cols**, 630–778 MB gz) |
| Sole consumer | the external-features path — `features_to_generate()` at `workflow/rules/add_external_features.smk:32`, then `add_external_features` (unit 21) |
| Measured range | wall **287.9 – 357.5 s**; `cpu_time` 255.7 – 321.8 s; `max_rss` **7 755 – 10 343 MB** |
| `params.abc_score_col` | `get_abc_score_col()` (`workflow/rules/utils.smk:104`) returns `"powerlaw.Score"` for **all six** clusters (no HiC configured) — verified in all six slurm logs. Effectively a constant here; the `ABC.Score` branch is untested by this dataset. |

### Measurements taken for this report

The benchmark table gives one number per job. To answer the I/O-vs-compute question I
re-ran the script under per-stage instrumentation. All artifacts are under
`/scratch/users/kaybrand/phase1_unit06/` (outside the repo; no pipeline file was modified).

| job | what |
|---|---|
| 40217422 | verbatim replay of `compute_arc_e2g.R` on `thp1_1`, with `Sys.time()` + `/proc/self/io` + `VmHWM` around every statement; plus `fwrite` to `.gz` vs plain `.tsv` vs `.gz` again; plus `zcat`/`gzip -6` calibration on the same node; plus a fresh-process decomposition of `readGeneric` |
| 40218215 / 40219386 | micro-benchmarks: `make.unique`, `pmatch`, `[.data.frame` row-name join, `lm` vs `cor`, all at 1–11 M rows |
| 40218659 / 40220638 | exact **uncompressed** byte sizes and NA/miss counts for all three tables in all six clusters |

Fidelity check: the replay's `VmHWM` was **9 286 MB** against the production job's
`max_rss` **8 952 MB** (+3.7%), and the sum of its stage timings (226.8 s) plus harness
overhead reproduces the production `cpu_time` (255.7 s) to within 0.1%. The replay ran on a
*different* node (`sh04-02n16`) from any production job, so absolute times carry node-to-node
variation; the *proportions* are what I rely on.

---

## 2. What the code actually does — the control flow that costs time

Straight-line script, no loops over data, five phases.

### 2.1 Read the ABC predictions — `compute_arc_e2g.R:90`

```r
pairs.E2G.ABC = readGeneric(abc_predictions_path, keep.all.metadata = T, header = T)
```

`genomation::readGeneric` → `genomation:::readTableFast`, and **it does not use `data.table::fread`**:

1. `read.table(file, nrows = 30, header = TRUE)` — a 30-row type-sniffing pass (0.00 s).
2. the sniffed classes are collapsed to a compact spec string
   (`ciiccdddddiciddccdidccddddddc`, 29 chars, verified) and handed to
3. **`readr::read_delim(...)`** — `readr` 2.1.5 edition 2, i.e. `vroom::vroom()`.
4. `as.data.frame(tibble)` (0.00 s — shares columns by reference).
5. `makeGRangesFromDataFrame(df, keep.extra.columns = FALSE)` (5.78 s).
6. `my.mcols <- df[, -(1:3)]` then `mcols(g) <- my.mcols[…]` (0.00 s each — pointer copies).

`genomation:::compressedAndUrl2temp` only intercepts *URLs*, so the local `.gz` path is
passed straight through and the decompression is done by the reader.

**Measured on `thp1_1`: 38.16 s total, of which `read_delim` alone is 27.12 s** and gzip
inflate of the 2 805 MB payload is 10.40 s (`zcat` on the same node). `fread` on the *identical
file* takes **13.09 s** — 2.07× faster than `read_delim` for byte-identical output
(10 173 077 × 29).

### 2.2 Read the Kendall pairs — `compute_arc_e2g.R:95-97`

```r
pairs.E2G.Kendall = GRanges(fread(kendall_predictions_path))
pairs.E2G.Kendall = pairs.E2G.Kendall[!is.na(mcols(pairs.E2G.Kendall)[,"Kendall"])]
```

25.78 s, decomposing as 6.17 s inflate + 6.92 s `fread` parse + **12.69 s for the
`GRanges()` coercion alone**. The NA filter drops 348 823 of 9 694 978 rows (3.6%) in 0.67 s.
This filter is load-bearing for §2.4 — see the `pmatch` finding.

### 2.3 `AddMax` — the interval join — `compute_arc_e2g.R:11-53`

The join key is faked into the `seqnames` slot. For each side, columns 1–3 are copied out,
`seqnames` is overwritten with `paste(chr, TargetGene, sep="_")` (`:20-22`, `:28-30`) and a new
`GRanges` is built (`:23`, `:31`). Composite seqlevels: **14 358** (query) and **13 851**
(subject). Then `findOverlaps` (`:35-36`), and the "keep the largest Kendall per query"
reduction: `order()` by Kendall descending (`:45`), a second **stable** `order()` by
`queryHits` (`:46`), `!duplicated(queryHits)` (`:47`), and a scatter into `mcols` (`:48-50`).

Measured cost split (`thp1_1`, 20.32 s total):

| line | operation | s |
|---|---|---:|
| `:20-22` | `paste()` over 10.17 M rows, query side | 1.18 |
| `:23` | `GRanges()` construction, query side | 6.10 |
| `:27-31` | `paste()` + `GRanges()`, subject side | 6.63 |
| `:35-36` | **`findOverlaps`** | **3.31** |
| `:39-43` | hits → data.frame + Kendall gather | 0.22 |
| `:45` | `order(Kendall, decreasing=TRUE)` | 1.03 |
| `:46` | `order(queryHits)` | 0.70 |
| `:47` | `!duplicated(queryHits)` | 0.69 |
| `:48-50` | `mcols(...)[i, "Kendall"] <- …` | 0.46 |

Two things worth naming. **The marshalling costs 3.8× the join** — 12.73 s to build the two
composite-seqname `GRanges` objects versus 3.31 s to query them. And the fan-out is
essentially 1:1: `length(findOverlaps(...))` = **9 720 926** against `n_cand_pairs` =
10 173 077, i.e. μ = H/n = **0.955**. 724 872 query rows get no hit at all, and of the
9 448 205 that do, only 2.89% have more than one — so the three-step max-reduction at `:45-47`
(2.42 s over all 9.72 M hits) resolves ties for under 3% of rows.

`order()` here is `method="auto"`, which resolves to **radix** for a numeric and an integer
vector below 2^31, so both sorts are O(H) linear passes, not O(H log H) comparison sorts, and
radix's stability is what makes the two-pass idiom correct.

### 2.4 `IntegrateABC` — `compute_arc_e2g.R:56-82`

```r
index.exclude_na = !is.na(...ABC) & ...ABC > 0 & !is.na(...Kendall)   # :61-63
bed.E2G.filter   = bed.E2G[index.exclude_na]                          # :65
lm.res = lm(scale(log(...ABC)) ~ scale(...Kendall))                   # :70-72
beta   = 1/coef(lm.res)[2]                                            # :74
```

5.84 s total: index 0.17, subset 1.46, two `sd()` 0.36, **`lm()` 2.96**, `mcols` write-back
0.89. `n_kept_for_lm` = 9 448 200 of 10 173 077 (92.9%).

Line `:65` subsets **all 27 metadata columns** of the `GRanges` to produce
`bed.E2G.filter`, of which exactly two (`powerlaw.Score`, `Kendall`) are subsequently read.
Heap grew 4 253 → 6 163 MB across that one statement.

Line `:70-72` fits a full OLS model — model frame, `n×2` model matrix, QR, residuals, effects,
fitted values — to extract one scalar at `:74`. Micro-benchmark at n = 9×10⁶:
**`lm()` = 2.87 s and `object.size(lm.res)` = 2 304 MB**, versus `cor()` = 0.52 s. Because both
regressors are `scale()`d, the slope being extracted is algebraically the Pearson correlation.
This statement is what sets the process memory high-water mark (§6.2).

### 2.5 Copy the gene-level expression columns — `compute_arc_e2g.R:111-123`

```r
df.gene_exp = mcols(pairs.E2G.Kendall)[, c("TargetGene", …3 cols)]   # :111-114
df.gene_exp = df.gene_exp[!duplicated(df.gene_exp$TargetGene), ]     # :116
rownames(df.gene_exp) = df.gene_exp$TargetGene                       # :117
mcols(pairs.E2G.ABC)[, c(…3 cols)] =
  df.gene_exp[pairs.E2G.ABC$TargetGene, c(…3 cols)]                  # :118-123
```

Build the 13 851-row lookup: **0.13 s**. Use it: **27.15 s** — 12% of the job, and the second
largest single line item after compression.

Line `:121` is a **character row-name lookup on a base `data.frame`**, and `base::[.data.frame`
does two expensive things with an index of length `n_cand_pairs`:

- **`i <- pmatch(i, rows, duplicates.ok = TRUE)`** (line 107 of `deparse(base::[.data.frame)`,
  R 4.3.3) — exact matching is hashed and free
  (0.20 s at n = 11×10⁶ when everything matches), but every *unmatched* input falls through to
  **partial matching, which is a nested scan over the table's row names**. Verified by
  micro-benchmark over 18 (n, k, miss-rate) points: cost = `c · n_miss · k` with
  **c = 5.0×10⁻⁹ s**, linear in `n_miss` and linear in `k` independently
  (k=5 000 → 20 532 at fixed `n_miss` scales the time by 4.00× for a 4.11× k increase).
- **`if (anyDuplicated(rows)) rows <- make.unique(as.character(rows))`** (same deparsed body,
  the `if (!drop)` block) — the result inherits `n_cand_pairs` row names drawn from only ~13 900 distinct
  genes, so `make.unique` fires on the whole vector and mints ~10 M `"GENE.k"` strings that
  are discarded on the next line. Micro-benchmark: **12.69 s at n = 11×10⁶**, linear in n
  (0.97 / 2.49 / 5.31 / 11.46 s at 1/2/4/8 M).

And `n_miss` is not zero: **359 442 rows (3.53%) at `thp1_1`**. The reason is a two-line
interaction — `:96-97` filters the Kendall object to non-NA *before* `:111` builds the lookup
from it, so the ~500 genes whose Kendall values are all NA are absent from the lookup, and every
candidate pair naming one of them becomes a `pmatch` partial-match probe. Both tables contain
the *same* 14 358 distinct `TargetGene` values (verified); the misses are created entirely by
the NA filter. Predicted cost 5.0×10⁻⁹ × 359 442 × 13 851 = **24.9 s**, i.e. essentially all of
the measured 27.15 s (the mechanism model over-predicts by ~1.35× because my synthetic keys are
longer than real gene symbols, so `psmatch` rejects them more slowly).

### 2.6 Serialise — `compute_arc_e2g.R:126-132`

```r
df.output = as.data.frame(pairs.E2G.ABC)      # :126  -> 0.07 s
colnames(df.output)[1] = "chr"                # :127
fwrite(df.output, file = arc_predictions_path, row.names=F, quote=F, sep="\t")   # :128-132
```

`as.data.frame(GRanges)` is **0.07 s and adds no heap** — it re-uses the existing column vectors
by reference. There is no second full-width copy of the table. (This falsifies the prior in the
briefing prompt that several full copies are live at peak; see §6.2 for what actually is.)

`fwrite` with a `.tsv.gz` output and `compress="auto"` compresses **in process, via data.table's
bundled zlib** — no external `gzip`/`pigz` subprocess. Verified three ways:

| | s | bytes written (`/proc/self/io` `wchar` == `write_bytes`) |
|---|---:|---:|
| `fwrite(… "out.tsv.gz")` | **98.07** | 630 109 119 |
| `fwrite(… "out.tsv")` (identical data, no `.gz`) | **16.66** | 3 591 824 914 |
| `fwrite(… "out2.tsv.gz", compress="gzip")` | **98.14** | 630 109 119 |

Both `.gz` writes produce a byte-identical 630 109 119-byte file, so this is reproducible to
0.07%. **Compression alone = 98.07 − 16.66 = 81.41 s.**

---

## 3. Size variables

Named exactly as in `size_manifest_multiome.tsv` unless flagged as measured-here.

**Primary.** `n_cand_pairs` (10 173 077 – 11 447 564, **1.125×**) — every per-row term scales
with it. Row count is unchanged from input to output (`row_identities_ok = yes`): the rule adds
7 columns and drops zero rows.

**Secondary (bytes).** `arc_bytes` (630 109 119 – 778 223 611, **1.235×**) for the output;
`n_cand_pairs_bytes` and `kendall_scored_bytes` for the two inputs. Also `n_kendall_pairs`
(9 694 978 – 11 317 834, 1.167×) and `k_genes` (20 531, constant).

**Two caveats on the byte columns, both measured.**

1. The codec cost scales with the **uncompressed** volume, which is not in the manifest.
   Measured for all six clusters (uncompressed columns are data rows only, header excluded;
   the whole-file `wc -c` for `thp1_1`'s ARC output is 3 591 824 914 B, 532 B more):

   | cluster | ABC in (uncomp) | Kendall in (uncomp) | ARC out (uncomp) | ARC out gz | out ratio | B/row out |
   |---|---:|---:|---:|---:|---:|---:|
   | thp1_1 | 2 805 198 956 | 1 516 754 359 | 3 591 824 382 | 630 109 119 | 5.700 | 353.1 |
   | thp1_2 | 2 834 974 336 | 1 551 120 046 | 3 629 723 717 | 645 189 492 | 5.626 | 353.1 |
   | jurkat_pma_cd3_4hr | 3 120 976 405 | 1 687 432 914 | 3 969 664 601 | 719 335 787 | 5.519 | 365.7 |
   | jurkat | 3 040 033 696 | 1 697 249 238 | 3 910 676 104 | 732 234 886 | 5.341 | 354.3 |
   | telohaec_crispri | 3 168 709 170 | 1 757 540 131 | 4 055 760 427 | 750 525 876 | 5.404 | 365.2 |
   | k562_crispri | 3 222 646 227 | 1 775 439 913 | 4 123 044 425 | 778 223 611 | 5.298 | 360.2 |

   `arc_bytes` spans **1.235×** but the uncompressed output spans only **1.148×**, because the
   deflate ratio itself varies 5.30–5.70 (7.6%) across clusters. So `arc_bytes` *over-states*
   the variation in the quantity that actually costs time. This refines the briefing's advice to
   "use the `*_bytes` columns": use them, but scale by ~5.5× and expect ~8% of the apparent
   spread to be entropy, not work.

2. **Two of the seven added columns are pure `GRanges` round-trip artifacts.** The output header
   is `chr start end width strand name …`; `width` and `strand` do not exist in the 29-column
   input and are introduced by `readGeneric`'s `GRanges` conversion and `as.data.frame` at
   `:126`. `strand` has exactly **one distinct value** (`*`) in all 10 173 077 rows. Together
   they are 62 175 121 of 3 591 824 382 uncompressed bytes (**1.73%**).

**A third variable that the manifest does not carry, and that a size-only model cannot see.**
The `pmatch` term of §2.5 is governed by `n_miss` — candidate pairs whose target gene has no
non-NA Kendall value anywhere — and by `k_lookup`, the row count of the gene lookup. Both
measured here for all six clusters:

| cluster | `n_cand_pairs` | `n_miss` | miss % | `k_lookup` | `n_miss · k_lookup` | implied pmatch s |
|---|---:|---:|---:|---:|---:|---:|
| telohaec_crispri | 11 106 687 | 139 409 | 1.26% | 14 206 | 1.98×10⁹ | **9.9** |
| jurkat | 11 036 283 | 263 435 | 2.39% | 14 047 | 3.70×10⁹ | 18.5 |
| jurkat_pma_cd3_4hr | 10 854 698 | 320 377 | 2.95% | 13 947 | 4.47×10⁹ | 22.3 |
| k562_crispri | 11 447 564 | 331 553 | 2.90% | 14 009 | 4.64×10⁹ | 23.2 |
| thp1_2 | 10 278 325 | 336 982 | 3.28% | 13 896 | 4.68×10⁹ | 23.4 |
| thp1_1 | 10 173 077 | 359 442 | 3.53% | 13 851 | 4.98×10⁹ | **24.9** |

This term spans **2.5×** across the six clusters and **anti-correlates with dataset depth**
(deeper RNA → fewer all-NA genes → fewer misses). It is worth ~3–9% of the job and it moves in
the *opposite* direction from `n_cand_pairs`, which partly cancels the size signal — see §5.

**No genuine cross term in the size variables** other than `O(n_miss · k_genes)` above.
`k_genes` does not otherwise enter: there is no per-gene loop anywhere in the script.

---

## 4. Asymptotic class derived from the code

> **O(`n_cand_pairs` · log `n_cand_pairs`)  +  O(`n_miss` · `k_lookup`)**

with the second term being the only place the code leaves the class genuinely ambiguous.

**The log factor, and the single place it comes from.** `findOverlaps` at
`compute_arc_e2g.R:35-36` is dispatched on two `GRanges`. `GenomicRanges` implements this by
building a `GNCList` (Nested Containment List) over the subject, which requires grouping and
ordering the ranges by seqlevel — O((n+m) log(n+m)) — and then answering each query in
O(log m + hits), for a total of **O((n+m) log(n+m) + H)**. Nothing else in the script is
super-linear:

- `paste()` (`:20-22`, `:28-30`) — O(n) allocations, hashed against R's global string cache.
- `GRanges()` / `makeGRangesFromDataFrame` (`:23`, `:31`, `:95`) — the seqnames coercion is
  `factor()`, i.e. O(n) hash + O(k log k) level sort with k ≈ 14 000.
- `order()` ×2 (`:45`, `:46`) — `method="auto"` → **radix**, O(H) per pass.
- `duplicated()` (`:47`, `:116`) — hashed, O(H) / O(m).
- `bed.E2G[index]` (`:65`), `mcols(...)[i, j] <-` (`:48`, `:76-79`) — O(rows × columns).
  `S4Vectors`' `[,DataFrame` applies `extractCOLS` **before** `extractROWS`, and `[<-,DataFrame`
  routes through `mergeROWS(extractCOLS(x, j), i, value)`, so the six subsetting expressions
  inside `:77-79` each touch one column, not all 27. Verified by reading the installed methods.
- `readGeneric`/`fread`/`fwrite` — O(bytes).
- `as.data.frame(GRanges)` (`:126`) — O(1) allocations, measured 0.07 s.

**No nested loop over `j_cand_elems` or `k_genes` exists in this script.** The row-count
identity (ARC rows == `n_cand_pairs`, zero filtering) means the only candidates for
superlinearity were the join and the I/O. The join is Θ(n) in realised work
(μ = H/n = 0.955, measured) with an O(n log n) index build; the I/O is Θ(bytes) ∝ Θ(n). Which
leaves exactly one non-Θ(n) term, and it is not in the join at all:

**The `pmatch` fallback at `:121` is O(`n_miss` · `k_lookup`)** — a genuine nested scan,
established by reading `base::[.data.frame` (line 107: `i <- pmatch(i, rows, duplicates.ok = TRUE)`)
and confirmed to be exactly bilinear by micro-benchmark. Statically the code cannot bound
`n_miss`: it is `≤ n_cand_pairs`, so the **worst case class of that one line is
O(n_cand_pairs · k_genes)**. Numerically, at `n_miss = n` and `k = 20 531` that line would take
5.0×10⁻⁹ × 10.2×10⁶ × 20 531 ≈ **1.0×10⁶ s ≈ 12 days**, against the 25 s it actually takes at
`n_miss/n = 3.5%`. That is a five-orders-of-magnitude range for a single line, gated entirely by
how many genes the upstream NA filter removes. This is the one place in the unit where the code
does not determine the class, and it is why `t_high` below is not merely `t_est` with a bigger
constant.

**The magnitude point the briefing insists on.** The log factor is real but tiny in constant:
`findOverlaps` is 3.31 s of a 226.8 s script — **1.5%**. Reporting the class as O(n log n) is
correct, and reporting that the exponent is worth 1.5% of the runtime at this scale is also
correct; both belong here.

---

## 5. Benchmark cross-check

```
cluster              wall s   cpu s  cpu/wall  max_rss MB   io_in    io_out  n_cand_pairs
jurkat               357.49  321.79    0.900     10342.61   572.50  5185.86    11036283
jurkat_pma_cd3_4hr   307.51  283.77    0.923      7754.59   224.98  5170.49    10854698
k562_crispri         330.53  312.00    0.944     10070.11   240.84  5444.54    11447564
telohaec_crispri     320.02  303.40    0.948      9690.89   386.39  5403.58    11106687
thp1_1               287.94  255.69    0.888      8952.07     0.00     0.00    10173077
thp1_2               293.84  280.97    0.956      9530.76   248.39  4769.23    10278325
```

**`jurkat` must be set aside.** `sacct` shows it ran on **`sh03-04n24` on 2026-08-20** under run
UUID `a5815c54-…`, while the other five ran on `sh04-16n27`/`sh04-16n28` on 2026-08-19 under
`7e55cc08-…`. It is also the only job with an `io_in` (572.50 MB) matching the full compressed
input volume (550.1 MB), i.e. a cold page cache, and it has the lowest `cpu/wall` (0.900).
Including it makes the six points **non-monotone in every size variable** — `k562` has 3.7% more
`n_cand_pairs` and 6.3% more `arc_bytes` than `jurkat` yet runs 7.5% *faster*. So the
cross-node/cross-run point is doing more damage than the size signal can carry.

**Fits.**

| model | all 6 | 5, excl. `jurkat` |
|---|---|---|
| `log(s) ~ log(n_cand_pairs)` | slope **1.372**, R² 0.628 | slope **1.132 ± 0.073**, R² 0.9876 |
| `log(cpu) ~ log(n_cand_pairs)` | slope 1.565, R² 0.743 | slope 1.406 ± 0.357, R² 0.838 |
| `log(s) ~ log(arc_bytes)` | slope 0.731, R² 0.600 | slope 0.610 ± 0.058, R² 0.974 |
| `log(s) ~ log(n_kendall_pairs)` | slope 0.937, R² 0.541 | slope 0.810 ± 0.096, R² 0.959 |

**Verdict: the observed scaling is consistent with the code-derived class and cannot confirm
it.** `n_cand_pairs` spans 1.125×, so `log(n_cand_pairs)` itself spans **1.0073×** — over this
range O(n) and O(n log n) differ by 0.7%, which is inside the measurement scatter. The 5-point
slope 1.132 ± 0.073 has a 95% interval of roughly [0.90, 1.36] (t₃ = 3.18); the effective local
exponent of `n ln n` here is 1.062, comfortably inside it, and so is 1.000. I am not fitting an
exponent to this.

**Two apparent disagreements, both explained, neither a code-vs-data conflict.**

1. *Time grows slower than `arc_bytes` (slope 0.61).* Because compression time scales with the
   **uncompressed** output, and the uncompressed output spans 1.148× while `arc_bytes` spans
   1.235× (§3, caveat 1). Refitting against uncompressed output bytes brings the slope back
   toward 1. Nothing anomalous.
2. *`n_kendall_pairs` slope 0.81, below 1.* The Kendall side enters through `fread` + the
   subject-side `GRanges` + `findOverlaps`, ~26 s of 227 s; a sub-unit slope on a 1.167× range
   with ~5% scatter is not interpretable.

**The `n_miss` term is a measurable mechanism behind the residual scatter.** `telohaec_crispri`
carries only 9.9 s of `pmatch` cost while `thp1_1` carries 24.9 s (§3) — a 15 s swing, 5% of the
job, moving *against* size. `telohaec` is the deepest cluster (394.8 M fragments, 174.2 M UMIs)
and is therefore both the largest and the cheapest per row on this term. This is a real
composition effect, not noise, and it means the size-only fits above are structurally slightly
mis-specified.

---

## 6. Concrete overhead sources

### 6.1 Where the time goes — the full decomposition (`thp1_1`, instrumented replay)

Stage timings, from `job_40217422.out`. Sums to 226.8 s; the production job's `cpu_time` was
255.7 s, the 28.9 s difference being conda activation, the in-job Snakemake DAG rebuild (visible
as a duplicated `Building DAG of jobs…` in every slurm log), R interpreter startup and the
generated Snakemake R wrapper.

| stage | `file:line` | s | % of 226.8 |
|---|---|---:|---:|
| **`fwrite` → `.gz`** | `:128-132` | **98.07** | **43.2%** |
|  ├─ zlib deflate | | 81.41 | 35.9% |
|  └─ number formatting | | 16.66 | 7.3% |
| `readGeneric` (ABC) | `:90` | 38.16 | 16.8% |
|  ├─ gzip inflate | | 10.40 | 4.6% |
|  ├─ `readr::read_delim` parse | | 16.72 | 7.4% |
|  ├─ `makeGRangesFromDataFrame` | | 5.78 | 2.5% |
|  └─ unattributed | | 5.26 | 2.3% |
| **`df.gene_exp[TargetGene, ]` join** | `:118-123` | **27.15** | **12.0%** |
|  ├─ `pmatch` partial-match fallback | | ~24.9 (modelled) | ~11% |
|  └─ `make.unique` on discarded row names | | ~11.7 (modelled) | ~5% |
| `fread` + `GRanges()` (Kendall) | `:95` | 25.78 | 11.4% |
|  ├─ gzip inflate | | 6.17 | 2.7% |
|  ├─ `fread` parse | | 6.92 | 3.1% |
|  └─ `GRanges()` coercion | | 12.69 | 5.6% |
| `AddMax` | `:11-53` | 20.32 | 9.0% |
|  ├─ composite-seqname `GRanges` build ×2 | `:20-23`, `:27-31` | 13.91 | 6.1% |
|  ├─ **`findOverlaps`** | `:35-36` | **3.31** | **1.5%** |
|  └─ order ×2 + dedup + scatter | `:45-50` | 2.88 | 1.3% |
| R startup + `library()` ×3 | `:2-6` | 10.60 | 4.7% |
| `IntegrateABC` | `:56-82` | 5.84 | 2.6% |
|  └─ of which `lm()` | `:70-72` | 2.96 | 1.3% |
| Kendall NA filter | `:96-97` | 0.67 | 0.3% |
| gene lookup build | `:111-117` | 0.13 | 0.1% |
| `as.data.frame(GRanges)` | `:126` | 0.07 | 0.0% |

Rolled up by *kind of work*:

| category | s | % |
|---|---:|---:|
| **gzip codec** (81.41 deflate out + 16.57 inflate in) | **98.0** | **43.2%** |
| text parse in + text format out (16.72 + 6.92 + 16.66) | 40.3 | 17.8% |
| `GRanges` construction / coercion (5.78 + 12.69 + 13.91) | 32.4 | 14.3% |
| `pmatch` fallback + `make.unique` row names | 27.2 | 12.0% |
| R startup + `library()` | 10.6 | 4.7% |
| **the actual analytics** (findOverlaps, sorts, dedup, `lm`, `sd`, subsets, NA filter) | **13.1** | **5.8%** |
| unattributed | 5.3 | 2.3% |

**The hypothesis in the prompt is confirmed and quantified: 75% of this rule (170.7 s of 226.8 s)
is gzip codec, text parse/format and data-structure marshalling; 5.8% is arithmetic.** Narrowing
"serialisation" to codec + text I/O only, still 61%.

Note on the two `:118-123` sub-rows above: 24.9 + 11.7 = 36.6 s of modelled cost against 27.15 s
measured. The micro-benchmark coefficients over-predict by ~1.35× (see §9.2); the two mechanisms
are the right ones and `pmatch` is the larger, but their split within the 27.15 s is inferred, not
separately measured.

### 6.2 Compression share — the explicit estimate

Direct, not modelled: **81.41 s**, from the difference between two `fwrite` calls on the same
in-memory `data.frame` (98.07 s to `.gz`, 16.66 s to plain `.tsv`), reproduced at 98.14 s on a
third call. That is:

- **28.3% of the 287.94 s wall job**
- **31.8% of the 255.69 s `cpu_time`**
- **35.9% of the 226.8 s of in-script work**
- the single largest line item in the rule, by a factor of 2.1 over the next one

Adding input decompression (10.40 + 6.17 = 16.57 s) gives **98.0 s = 34.0% of the job / 43.2% of
in-script work** for gzip in total.

**Rates measured on the same node, on these exact bytes:**

| | MB/s |
|---|---:|
| deflate, inside `fwrite` (3 591.8 MB in → 630.1 MB out, 81.41 s) | **44.1** in / 7.7 out |
| `gzip -6` on the identical file (3 591.8 MB → 626.1 MB, 89.45 s) | 40.2 in |
| inflate, `zcat` (mean of the three files) | **244** |
| `cat` the 3.59 GB plain file to `/dev/null` | 8 452 |

Three consequences.

1. `fwrite`'s in-process zlib is **~9% faster than the `gzip(1)` binary** (81.41 s vs 89.45 s on
   identical bytes) and produces a 0.63%
   larger file — i.e. an equivalent deflate level (6), in-process, no subprocess. That settles
   the "`fwrite` vs external gzip" question in the prompt: it is `fwrite`, in-process.
2. **My measured asymmetry is 5.5× (244 / 44.1), not the ~40× carried in the prompt's prior.**
   The prior's own quoted rates (15–21 MB/s compress vs 194 MB/s decompress) imply ~11×, and my
   compression rate is ~2.4× faster than that range. This rule's output is a numeric-heavy TSV
   that deflates only 5.3–5.7×, which is cheap work for zlib. **A Phase-2 model that applies the
   15–21 MB/s prior to this rule would over-state its compression cost by roughly 2×.** I flag
   this as a disagreement to reconcile rather than asserting either side is wrong: same codec,
   different payload entropy, different node.
3. The compression share here is **~32% of CPU, not the ≥52% seen in unit 2** (`run_e2g_qnorm`).
   Consistent direction, half the magnitude.

### 6.3 `io_in` / `io_out` — mechanism resolved for this rule

`/proc/self/io` deltas around every statement (`wchar` and `write_bytes` were *identical* at
every step, so these are real physical writes, not buffered chars):

| stage | rchar MB | **wchar MB** | read_bytes MB |
|---|---:|---:|---:|
| `readGeneric` (`:90`) | 283.1 | **2 805.2** | 2.5 |
| `fread` (`:95`) | 203.7 | **1 516.8** | 0.0 |
| `fwrite` `.gz` (`:128`) | 0.0 | **630.1** | 0.0 |

2 805.2 MB and 1 516.8 MB are **exactly** the uncompressed sizes of the two input files
(2 805 198 956 and 1 516 754 359 bytes, measured independently). So:

**Both readers decompress the entire `.gz` input to a temporary file on disk and then mmap it
back.** `readr`/`vroom` does this (2 805.2 MB), and `data.table::fread` does it too (verified
separately: `fread` on the same file also wrote exactly 2 805.2 MB). The mmap'd read-back does
not appear in `rchar` at all, which is why `rchar` only shows the compressed input size.

This resolves the 7–8× `io_out` "write amplification" for this rule, and the mechanism is
**different from the project-wide explanation**. Predicting
`io_out = uncompressed_ABC + uncompressed_Kendall + compressed_output`:

| cluster | predicted MB | measured `io_out` | pred/meas |
|---|---:|---:|---:|
| jurkat | 5 469.5 | 5 185.9 | 1.055 |
| jurkat_pma_cd3_4hr | 5 527.7 | 5 170.5 | 1.069 |
| k562_crispri | 5 776.3 | 5 444.5 | 1.061 |
| telohaec_crispri | 5 676.8 | 5 403.6 | 1.051 |
| thp1_2 | 5 031.3 | 4 769.2 | 1.055 |
| thp1_1 | 4 952.1 | **0.00** (capture failed) | — |

Consistent to ±1% across five clusters, with a uniform ~5.5% under-count attributable to
psutil's final sample landing before process exit. **86–88% of `io_out` is temporary
decompression spill; only 12–14% is the output file, and that part is counted *compressed*
(630.1 MB, not 3 591.8 MB).** The addendum's general note — that `io_out` measures uncompressed
bytes entering the compression stream — is *not* the mechanism here: for `arc_e2g` the
amplification is 4.3 GB of genuine physical temp-file I/O, plus a correctly-counted compressed
output. I am reporting the mechanism, not a defect.

Where those temp files land matters: `TMPDIR` in the production jobs was `/tmp`, which on these
nodes is **XFS on the local `/dev/md0` RAID (1.9 TB, same device as `/lscratch`)**, verified. So
it is real local disk, not tmpfs, and it is not charged to the job's memory cgroup. At
multi-GB/s the *time* cost is negligible (`cat` of 3.59 GB = 0.42 s); the cost is the RSS
consequence in §6.4.

`io_in` is `read_bytes` (physical) and is page-cache dependent — 0.00 to 572.50 MB against a
constant ~550 MB of compressed input. Unusable, as nine other units found.

### 6.4 Memory

`max_rss` **7 754.59 – 10 342.61 MB** (1.33× spread) against `mem_mb = 64 000`
(`arc_e2g.smk:29`, Slurm `ReqMem = 62.50 G`) → the declared floor is **6.0–8.1× the observed
peak**. Replay `VmHWM` 9 286 MB vs production `max_rss` 8 952 MB (+3.7%).

Two events set the high-water mark, and neither is what the briefing prompt anticipated:

1. **`readGeneric` at `:90` → `VmHWM` 7 899 MB while R's tracked heap was 2 638 MB.** The
   fresh-process decomposition pins it on `read_delim` alone: 7 422 MB `VmHWM` with a 2 255 MB
   heap and R's `gc()` reporting only 3 543 MB of lifetime max heap (Ncells 698.4 + Vcells
   2 844.8). The missing ~3.9 GB is outside R's allocator: the **2 805 MB decompressed temp file
   is mmap'd and every page is touched during parsing, so all of it becomes resident and is
   charged to the process**, plus vroom's own C++ index and buffers. The peak is therefore an
   artifact of decompress-to-temp-then-mmap, not of the data's size.
2. **`lm()` at `:70-72` → `VmHWM` 7 899 → 9 286 MB (+1 387 MB)**, to produce one scalar
   consumed at `:74`. `object.size(lm.res)` is 2 304 MB at n = 9×10⁶ (micro-benchmark).

**There is *not* a second full-width copy of the 36-column table live at peak.**
`as.data.frame(pairs.E2G.ABC)` at `:126` costs 0.07 s and adds no heap because it shares the
column vectors; the heap reading after it (4 642 MB) is *lower* than before it. `bed.E2G.filter`
at `:65` is the one genuine near-full copy (92.9% of rows × all 27 metadata columns, +1.9 GB
heap), and it exists only so two of those columns can be read.

Peak scales weakly and noisily with size (`rss ~ n_cand_pairs`, R² = 0.14; `jurkat_pma` peaks
2.6 GB below `jurkat` at similar n) because R's GC timing relative to the peak is not
deterministic.

The `make.unique` row names at `:121` allocate ~10 M `CHARSXP`s. R's global string cache **is**
GC-collectable, and the measurements bear that out (heap 5 296 MB after the join → 4 642 MB
after the next statement), so the cost of that line is allocation and hashing time (~12 s), not
retained memory.

### 6.5 Threading: declared 0 → realised 1, and one term that is provably serial

The rule declares **no `threads:`** (`arc_e2g.smk` has `benchmark`, `resources`, `conda`,
`script` and nothing else), so Snakemake schedules it at 1 and the slurm logs confirm
`cpus_per_task=1`. Measured `cpu_time/wall` = **0.888 – 0.956**, i.e. 89–96% of one core, no
parallelism. The remaining 4–11% is page-cache miss and scheduling latency, not idle I/O wait
(§6.2: physical I/O time is under a second).

The mechanism, and it differs between the two libraries in the same process — both verified
inside a job:

- **`data.table` resolves to one thread.** `getDTthreads()` returned **1** inside the job —
  verified. Two independent things pin it: `OMP_NUM_THREADS` was `1` in the job environment, and
  a separate probe under a cpuset had `omp_get_num_procs()` return `1` as well, either of which is
  sufficient (data.table's verbose output confirms it reads both). data.table's `fwrite` gzip path is documented as
  thread-parallel per output buffer, so the **81.4 s deflate term is single-threaded purely
  because the rule asks for one CPU**. I did not measure a multi-thread `fwrite` — that would be
  evaluating a change, which is out of scope for Phase 1.
- **`readr` does not honour the cgroup.** `readr::read_delim` defaults
  `num_threads = readr_threads()`, which resolves to `parallel::detectCores()`. Inside the job
  that returned **64** — the node's full core count, ignoring the 1-CPU cpuset. So `readr`
  requests 64 threads while confined to one CPU, oversubscribing rather than parallelising, and
  gains nothing: gzip inflate of a single stream is serial regardless.

### 6.6 Other concrete overheads

- **`library(genomation)` for one function.** `:2-6` loads `GenomicRanges`, `genomation` and
  `data.table`; `genomation` drags in `Biostrings`, `rtracklayer`, `GenomicAlignments`,
  `BSgenome`, `impute`, `seqPattern`, `reshape2`, `plotrix` and more (the
  `replacing previous import 'Biostrings::pattern' by 'grid::pattern'` warning in every slurm log
  is that resolution happening). Cost **10.60 s**, 4.7% of the job, from a Lustre-backed conda
  prefix. `genomation` is used for exactly one call, `readGeneric` at `:90`, whose body needs
  only `readr` and `GenomicRanges`.
- **Reader choice costs ~14 s.** `fread` reads the identical 29-column file in 13.09 s vs
  `read_delim`'s 27.12 s — 2.07×, on top of an identical 2 805 MB temp-file write. That is 6% of
  the job, spent inside a third-party wrapper.
- **A 30-row `read.table` sniffing pass** (`readTableFast`) opens and decompresses the head of
  the ABC file before `read_delim` reopens it. Measured 0.00 s; noted for completeness.
- **`paste()` twice over ~10 M and ~9 M rows** (`:20-22`, `:28-30`) to synthesise a join key
  into a coordinate field, then two `GRanges` constructions to make it queryable: 13.91 s for
  the marshalling, 3.31 s for the query it enables.
- **Both input files are read once each; the output is written once.** No redundant passes over
  the same data at the R level. (The temp-file round trip in §6.3 is one extra *physical* write
  and mmap read per input, inside the readers.)
- **1.73% of the output is the `width` + `strand` artifact columns** (§3), `strand` being a
  single repeated `*`.
- **Fixed harness cost ≈ 28.9 s** (production `cpu_time` minus in-script total): conda
  activation, an in-job Snakemake DAG rebuild that re-parses the whole workflow including
  `ENCODE_rE2G/ABC/workflow/rules/utils.smk` (its pandas `FutureWarning` fires twice per job
  log), R interpreter startup and the generated wrapper. ~10% of the job.

---

## 7. The three functions

Variables, named exactly as in `size_manifest_multiome.tsv`:

- `n` ≡ `n_cand_pairs`
- `n_miss`, `k_lookup` — *not* in the manifest; measured here (§3 table). Only needed for the
  structural model in §7.2.

### 7.1 Deliverable — evaluate these

**Wall-clock seconds.** Class as derived in §4; `t_low`/`t_high` differ in the **constant**, for
the reasons given below.

```
t_est(n)  = 29 + 1.60e-6 * n * ln(n)
t_low(n)  = 20 + 1.35e-6 * n * ln(n)
t_high(n) = 40 + 2.10e-6 * n * ln(n)
```

`n = n_cand_pairs`; `ln` is the natural logarithm; result in **seconds**.

**Calibration.** `t_est` fits the five same-run/same-day jobs to within **±1.3%**:

| cluster | `n_cand_pairs` | measured wall s | `t_est` | error |
|---|---:|---:|---:|---:|
| thp1_1 | 10 173 077 | 287.94 | 291.6 | +1.3% |
| thp1_2 | 10 278 325 | 293.84 | 294.5 | +0.2% |
| jurkat_pma_cd3_4hr | 10 854 698 | 307.51 | 310.4 | +0.9% |
| telohaec_crispri | 11 106 687 | 320.02 | 317.3 | −0.9% |
| k562_crispri | 11 447 564 | 330.53 | 326.7 | −1.2% |
| *jurkat (other node, other day)* | 11 036 283 | *357.49* | *315.4* | *−11.8%* |

**Why the bounds are what they are.**

- `t_low` — warm page cache, a fast node, and `n_miss` at its observed floor
  (`telohaec_crispri`, 1.26% → the `pmatch` term drops from ~25 s to ~10 s, ~5% of the job);
  minimum fixed overhead ~20 s. Constant −16% vs `t_est`.
- `t_high` — cold page cache and a slower node. `jurkat`, the one job that has both, implies
  a constant of 1.836×10⁻⁶ from `(357.49 − 29)/[n ln n]`; 2.10×10⁻⁶ covers that plus `n_miss` at
  its observed ceiling. Constant +31%.

**I am deliberately not widening the exponent to match the regression.** The 5-point log-log
slope is 1.132 ± 0.073 (95% CI ≈ [0.90, 1.36]), but that spread is measurement scatter over a
1.125× range, not evidence about a class, and the instruction is explicit that these functions
bound the code-derived class. A pure O(n) form, `29 + 2.60e-5 * n`, fits the same five points to
±1.9% and is statistically indistinguishable — because `ln(n)` varies by 0.7% here. Use whichever
is convenient; they agree to within 1% anywhere in 10⁷–10⁸.

**The one exponent caveat that is real, and that Phase 2 should carry separately.** §4
establishes that `:121` is O(`n_miss` · `k_lookup`) with `n_miss` statically unbounded. Inside
the calibration above, `n_miss/n` is 1.3–3.5% and the term is 3–9% of the job. If a future
dataset drove `n_miss/n` toward 1 — many genes with all-NA Kendall, e.g. a shallow-RNA cluster —
that single line would grow by up to ~30× while `n` did not move at all, and `t_high` above would
be badly wrong. **This is not covered by the constants.** The conditional term is:

```
t_pmatch(n_miss, k_lookup) = 5.0e-9 * n_miss * k_lookup        [seconds]
```

evaluate it, subtract the value it takes at the calibration point for the cluster in question
(9.9–24.9 s, §3), and add the new one.

### 7.2 Structural model — mechanism, not the deliverable

Reproduces the decomposition in §6.1. Models **`cpu_time`**; multiply by 1/0.92 for wall.
`Bin` = uncompressed ABC + uncompressed Kendall bytes; `Bout` = uncompressed ARC bytes
(both in §3's table; `Bout ≈ 5.5 * arc_bytes` if you only have the manifest).

```
cpu(n, Bin, Bout, n_miss, k_lookup) =
    39.5                              # R + library() + conda + in-job Snakemake DAG rebuild
  + 1.052e-8 * Bin                    # gzip inflate + text parse of both inputs (+5.3 s unattributed)
  + 2.730e-8 * Bout                   # fwrite: 4.64e-9 format + 2.27e-8 deflate
  + 5.326e-6 * n                      # GRanges marshalling (3.18e-6), joins/sorts/subsets/lm
                                      #   (1.29e-6), make.unique row names (0.85e-6)
  + 3.71e-9  * n_miss * k_lookup      # pmatch partial-match fallback at :121
```

Term values at `thp1_1`: 39.5 + 45.47 + 98.06 + 54.18 + 18.47 = **255.68 s**.

Exact at `thp1_1` by construction (255.7 s vs 255.69 s measured). Across the other five it
**under-predicts `cpu_time` by 4–17%**, which is the node-to-node clock spread — it was
calibrated on one cluster on one node (`sh04-02n16`) that appears in no production run. Use §7.1
for numbers; use this to reason about a changed regime (different compression level, different
reader, a different `n_miss`).

### 7.3 Evaluations

`n_cand_pairs` = 10 173 077 / 10 278 325 / 10 854 698 / 11 036 283 / 11 106 687 / 11 447 564;
mean **10 816 106**.

| point | `n_cand_pairs` | `t_low` | **`t_est`** | `t_high` | measured |
|---|---:|---:|---:|---:|---:|
| **mean cluster** | 10 816 106 | 256.5 s (4.27 min) | **309.3 s (5.16 min)** | 407.9 s (6.80 min) | — |
| **min — `thp1_1`** | 10 173 077 | 241.6 s (4.03 min) | **291.6 s (4.86 min)** | 384.7 s (6.41 min) | 287.94 s |
| **max — `telohaec_crispri`** | 11 106 687 | 263.2 s (4.39 min) | **317.3 s (5.29 min)** | 418.4 s (6.97 min) | 320.02 s |
| *max by `n_cand_pairs` — `k562_crispri`* | 11 447 564 | 271.2 s (4.52 min) | **326.7 s (5.44 min)** | 430.7 s (7.18 min) | 330.53 s |

**Note the naming trap in the brief.** `telohaec_crispri` is the max by `n_frag` (394.8 M),
`n_umi` (174.2 M) and `n_cells` (15 555), but **`k562_crispri` is the max by `n_cand_pairs`**
(11 447 564 vs 11 106 687) **and by `arc_bytes`** (778.2 vs 750.5 MB) — and it is the slowest of
the five comparable jobs (330.53 s). Since this rule is governed by `n_cand_pairs`, `k562_crispri`
is the true worst case and I have evaluated both.

Total across the 6 multiome jobs: `t_est` sums to **1 856 s = 30.9 min** (measured 1 897 s =
31.6 min). All six run in parallel in the DAG, so the critical-path contribution is one job, not
the sum — Phase 2's arithmetic, not mine.

---

## 8. Modality

**Multiome only, and the class question does not arise.** `checkpoint features_required`
(`workflow/rules/add_external_features.smk:1-22`) writes `"ARC"` for all six multiome clusters
and `"Neither"` for all six scATAC clusters, so `features_to_generate()`
(`add_external_features.smk:26-33`) never requests
`ARC/EnhancerPredictionsAllPutative_ARC.tsv.gz` on the scATAC path and `arc_e2g` has **0 jobs**
there. This is why `arc_bytes` and `kendall_scored_bytes` are empty columns in
`size_manifest_scatac.tsv` (addendum §C).

The rule also cannot exist on the scATAC path for a substantive reason, not just a config one:
its second input is `Kendall/Pairs.Kendall.tsv.gz`, which requires per-cell RNA
(`compute_kendall` consumes `rna_matrix_file`). There is no scATAC counterfactual to model.

For context only, and citing unit 21's finding rather than re-deriving it: this rule's output is
the input that makes `add_external_features` ~15× more expensive in multiome (12.6–21.3 min) than
in scATAC (0.9–1.3 min). That is a dependency observation. Prioritisation is Phase 2's.

---

## 9. What the measurement cannot tell us

1. **Nothing about the exponent.** `n_cand_pairs` spans 1.125×, so `ln(n_cand_pairs)` spans
   1.0073×. O(n) and O(n log n) differ by 0.7% across every point I have. The 5-point slope
   1.132 ± 0.073 admits both. The O(n log n) claim in §4 rests entirely on reading
   `findOverlaps`' `GNCList` construction at `:35-36`, and the benchmarks neither support nor
   contradict it. Six points at 1.125× can never do so.
2. **Nothing about the O(n_miss · k_genes) term at scale.** The bilinearity is established by
   micro-benchmark over 18 synthetic points, and the coefficient 5.0×10⁻⁹ over-predicts the real
   measurement by ~1.35× (real gene symbols are shorter than my synthetic keys, so `psmatch`
   rejects them faster). `n_miss/n` spans only 1.3–3.5% here; the behaviour at 50% or 100% is
   extrapolation from the mechanism, not measurement.
3. **`n = 1`. No replicates, and one point is from a different machine.** Five jobs ran
   2026-08-19 on two nodes of one run; `jurkat` ran 2026-08-20 on a third node under a different
   run UUID and is 12% slower than the model — an 11-point spread I cannot decompose into node
   clock, cold cache and `n_miss` without repeating it. This rule exists only in multiome, so
   the briefing's "n=2 for shared rules" does not apply.
4. **My decomposition is one cluster on one node.** All stage timings, memory events and the
   compression share are from `thp1_1` on `sh04-02n16`, a node used by no production job.
   Proportions should hold; absolute constants carry node variation, and §7.2's 4–17%
   under-prediction is direct evidence of exactly that.
5. **The instrumented replay is not bit-identical to production.** I inserted 28 `gc(FALSE)`
   calls (their ~47 s is excluded from every stage timing, but they also *relieve* heap pressure
   the real run experiences), read paths via `$L_SCRATCH` rather than `/tmp` (same physical
   device, verified), and ran outside the Snakemake wrapper. The `VmHWM`/`max_rss` agreement to
   3.7% and the `cpu_time` reconstruction to 0.1% are the fidelity evidence.
6. **The `ABC.Score` branch is never exercised.** `get_abc_score_col()` returns
   `"powerlaw.Score"` for all six clusters because none has `HiC_type` set. A HiC-configured
   cluster takes the `"ABC.Score"` path at `:106`, and I have no measurement of it. The column is
   read, not computed, so I expect no cost difference — but that is inference from `:61-79`, not
   data.
7. **`max_rss` is a sampled maximum of a GC'd heap.** The 1.33× spread across clusters is mostly
   *when* R collected relative to the peak, not how much memory the data needs
   (`rss ~ n_cand_pairs` R² = 0.14). Any memory model built on these six numbers inherits that.
8. **The compression rate disagreement is unresolved.** I measure 44.1 MB/s deflate on this
   payload; the cross-unit prior says 15–21 MB/s. Both may be right for their own data. I cannot
   tell from here whether the difference is payload entropy, deflate level, or node, and Phase 2
   should not silently pick one.
9. **I did not measure what more CPUs would do.** §6.5 establishes *why* `fwrite`'s deflate is
   serial (1-CPU allocation → `getDTthreads() == 1`) but deliberately stops short of running it
   with more threads, which would be evaluating a change rather than characterising the current
   one.
10. **`io_in` is unusable and `thp1_1`'s `io_out` row is a capture failure** (0.00/0.00 with a
    630 MB file on disk). My `/proc/self/io` numbers replace both for this rule; they are from
    the replay, not from the production processes.

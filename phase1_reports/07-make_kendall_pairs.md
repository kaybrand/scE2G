# Unit 7 — `make_kendall_pairs` (Tier B)

## 1. Unit identity

- Script: `workflow/scripts/feature_computation/make_kendall_pairs.R` (75 lines)
- Rule: `make_kendall_pairs`, `workflow/rules/make_kendall_pairs.smk:3-31`
- Modality: **Multiome only**. 6 jobs (one per cluster).
- No `threads:` declared in the rule → Snakemake gives it the default of 1 core. `resources: mem_mb=32*1000` (`workflow/rules/make_kendall_pairs.smk:26-27`).

## 2. Correction to the premise in the unit brief

The unit brief hypothesizes the rule "generates ~10.7M pairs from `j_cand_elems` (~157k
elements) and `k_genes` (20,531 genes)" and asks whether the enumeration is `O(j·k)`,
`O(j · genes-per-window)`, or an interval-overlap join. Having read the code, **that premise
is wrong for this script**: neither `j_cand_elems` nor `k_genes` is an input here at all.

The two inputs declared at `workflow/rules/make_kendall_pairs.smk:5-16` are:

- `narrowPeak` = `Peaks/macs2_peaks.narrowPeak.sorted` → this is `n_peaks` (432,039–595,372
  raw MACS2 peaks, **pre-cap**, not the post-cap `j_cand_elems` candidate set).
- `allPutative` = `Predictions/EnhancerPredictionsAllPutative.tsv.gz` → this is `n_cand_pairs`
  (10.17M–11.45M rows), which is **already** a table of (enhancer region, `TargetGene`) pairs
  — the gene assignment (verified: 29 columns, includes `TargetGene`, `TargetGeneTSS`,
  `distance`, `gene_idx`, `ABC.Score`, etc.) happened upstream, in whatever ABC step produced
  `EnhancerPredictionsAllPutative.tsv.gz`, using its own distance-window logic — not in this
  script.

So the actual job of `make_kendall_pairs.R` is: **re-anchor the already gene-tagged candidate
pairs from `EnhancerPredictionsAllPutative.tsv.gz` onto the raw (non-extended) MACS2 peak
intervals**, via a genomic interval overlap, then deduplicate. The script's own header comment
says this explicitly: "Make enhancer-gene pairs for computing Kendall correlation using
non-extended peaks" (`make_kendall_pairs.R:1-2`).

**Answer to the specific question asked: this is an interval-overlap join** —
`GenomicRanges::findOverlaps()` — not `O(j·k)` and not a `data.table` non-equi/`foverlaps`
join. See §4.

## 3. What the code actually does (control flow, `file:line`)

1. `make_kendall_pairs.R:18` — `readGeneric(narrowPeak_path)`: parse `n_peaks` (432k–595k)
   BED-like rows into a `GRanges`.
2. `make_kendall_pairs.R:21-23` — `readGeneric(allPutative_path, keep.all.metadata=T,
   header=T)`: parse the full `EnhancerPredictionsAllPutative.tsv.gz` — `n_cand_pairs` rows
   (10.17M–11.45M) × 29 columns, gzip-compressed (`n_cand_pairs_bytes` 268.5–319.8 MB on disk)
   — into a `GRanges` with all 29 columns carried as `mcols`. This is by far the largest object
   built in the script and, I believe, the dominant cost (see §6); I cannot fully verify this
   from the pipeline's own code since `genomation::readGeneric`'s internals live outside this
   repo, but the I/O and RSS numbers below are consistent with it.
3. `make_kendall_pairs.R:26-27` — `findOverlaps(bed.narrowPeak, bed.allPutative)`: an
   interval-tree (NCList) overlap join between the `n_peaks`-sized query and the
   `n_cand_pairs`-sized subject. This is the step that actually produces the pairing; it is
   **not** a nested loop over peaks × genes.
4. `make_kendall_pairs.R:30-31` — index `bed.narrowPeak` by the hit list and attach
   `TargetGene` from the matched `allPutative` row. Output length = number of overlap hits,
   call it `H` (before dedup).
5. `make_kendall_pairs.R:34-44` — build `PeakName` and `PairName` by `paste()`-ing genomic
   coordinates and gene name for each of the `H` hits (vectorized string concatenation).
6. `make_kendall_pairs.R:47` — `order()` the `H` rows by `PairName` (sorts ~10–11M strings).
7. `make_kendall_pairs.R:48` — `duplicated()` on `PairName`, then subset to unique rows. This
   is the step that turns `H` (raw overlap hits) into the reported `n_kendall_pairs`
   (9.69M–11.32M, the row count of the final `Pairs.tsv.gz`).
8. `make_kendall_pairs.R:51-57` — `as.data.frame()` coerces the deduplicated `GRanges` (with 6
   selected columns) to a plain data frame.
9. `make_kendall_pairs.R:69-74` — `data.table::fwrite()` writes the final table
   (`n_kendall_pairs` rows × 6 cols) to `Kendall/Pairs.tsv.gz` (`n_kendall_pairs_bytes`,
   70.9–82.7 MB).

## 4. Size variables and asymptotic class

Named exactly as in `size_manifest_multiome.tsv`:

- `n_peaks` (432,039–595,372, ratio 1.38×) — query side of the overlap join; the MACS2 peak
  count, **pre-cap** (not `j_cand_elems`).
- `n_cand_pairs` (10.17M–11.45M, ratio 1.13×) — subject side of the join; also the dominant
  parse cost (29-column table).
- `n_kendall_pairs` (9.69M–11.32M, ratio 1.17×) — the output row count, i.e. the number of
  *unique* overlap hits after dedup; this drives the `paste`/`order`/`duplicated`/`fwrite`
  steps. `n_kendall_pairs` and `n_cand_pairs` are close in magnitude (ratio ~0.95–1.06 within
  a cluster) but **not identical** — this rule is explicitly excluded from the collapsed
  row-count identities in the briefing.

Class, by step:

| Step | file:line | Class |
|---|---|---|
| Parse `narrowPeak` | `:18` | `O(n_peaks)` |
| Parse `allPutative` (29 cols) | `:21-23` | `O(n_cand_pairs · 29)` — 29 is a fixed column count, not a manifest variable |
| `findOverlaps` (NCList interval tree) | `:26-27` | `O((n_peaks + n_cand_pairs) · log(n_cand_pairs))` — GenomicRanges builds an interval tree on one range set and queries it with the other; this is the standard NCList complexity, not `O(n_peaks · n_cand_pairs)` |
| Attach gene, build names | `:30-44` | `O(H)`, `H` = raw hit count, same order as `n_kendall_pairs` |
| `order()` by `PairName` | `:47` | `O(H log H)` worst case (comparison sort); R's default `order()` on character vectors uses a radix method that is closer to `O(H)` in practice, but I have not verified which method fires here |
| `duplicated()` + subset | `:48` | `O(H)` (hash-based) |
| `as.data.frame` | `:51-57` | `O(n_kendall_pairs)` |
| `fwrite` | `:69-74` | `O(n_kendall_pairs)`, fast C writer |

**Overall class**: dominated by the `n_cand_pairs`-sized parse and the
`(n_peaks + n_cand_pairs)·log(n_cand_pairs)` overlap join, with a secondary
`n_kendall_pairs·log(n_kendall_pairs)` sort/dedup pass. Since `n_cand_pairs` and
`n_kendall_pairs` are the same order of magnitude and `n_peaks` is 20–25× smaller than
`n_cand_pairs`, the whole script is well described as:

**`O(n_cand_pairs · log(n_cand_pairs))`**, with `n_peaks` entering additively as a much
smaller `O(n_peaks · log(n_cand_pairs))` query-side term that the six data points cannot
separate out (see §5).

This is the interval-overlap join class asked about in the brief — confirmed at
`make_kendall_pairs.R:26-27` — and it is a real, non-trivial improvement over a naive
`O(n_peaks · n_cand_pairs)` nested comparison (which would be ~4.8×10¹² operations at these
sizes and clearly is not what is happening, given the observed 3–6 minute runtimes).

## 5. Benchmark cross-check

```
cluster            wall_s  cpu_time  max_rss(MB) io_in(MB) io_out(MB)  n_peaks  n_cand_pairs  n_kendall_pairs
jurkat              307.19   250.32     8484.12     140.29    2899.20   520543     11036283        10788202
jurkat_pma_cd3_4hr   214.82   193.38     7368.10     259.35    2976.40   595372     10854698        10759396
k562_crispri         241.78   224.59     7625.26     277.79    3073.36   503068     11447564        11317834
telohaec_crispri     347.59   317.69     7468.78     383.89    3099.64   561437     11106687        11136922
thp1_1               196.12   145.53     7935.86     370.04    2675.25   432039     10173077         9694978
thp1_2               178.84   141.33     4735.17     371.94    2703.64   435637     10278325         9910221
```

- Wall time spans **178.8–347.6 s (1.94×)**, larger than the spread of *any* of the size
  variables at play: `n_cand_pairs` 1.13×, `n_kendall_pairs` 1.17×, `n_peaks` 1.38×. A
  variable that ranges 1.1–1.4× cannot mechanically produce a 1.94× output spread, so
  something other than problem size — Slurm node contention, page-cache state, GC pauses —
  is contributing comparably to, or more than, size itself at this data range.
- Direct evidence of this noise: `jurkat_pma_cd3_4hr` has the **largest** `n_peaks` (595,372)
  of all six clusters but a below-median wall time (214.82 s); `k562_crispri` has the largest
  `n_cand_pairs` (11,447,564) but is not the slowest job — `telohaec_crispri` is, despite
  `k562_crispri`'s `n_cand_pairs` being 3% larger.
- I fit a simple OLS line of `wall_s` on `n_cand_pairs` across the six points anyway, to be
  concrete about how weak the signal is: slope ≈ 86.8 s per million `n_cand_pairs`, intercept
  ≈ **−692 s** (physically impossible — negative time at n=0), R² ≈ 0.42. A negative intercept
  from a linear fit is the standard symptom of fitting noise over too narrow a range; I am
  reporting this fit only to show it should not be trusted, not as a calibration.
- `cpu_time/wall_s` ratios: jurkat 0.82, jurkat_pma_cd3_4hr 0.90, k562_crispri 0.93,
  telohaec_crispri 0.91, thp1_1 0.74, thp1_2 0.79. All comfortably below 1 and consistent with
  a single R process (no `threads:` declared → 1 core reserved) — no evidence of hidden
  multi-core parallelism from `fwrite` or `data.table`'s internal thread pool materially
  affecting wall time here.
- **Conclusion**: the code-derived class (`O(n_cand_pairs · log(n_cand_pairs))`, dominated by
  parse + interval join + sort) is plausible and not contradicted by the data, but the six
  clusters cannot confirm or refute it — the benchmark spread is dominated by something other
  than the input-size variables measured.

## 6. Concrete overhead sources

- **I/O asymmetry is large and worth flagging**: `io_in` is only 140–384 MB (matches reading
  a ~270–320 MB gzip `allPutative` file plus a ~32–44 MB `narrowPeak` file), but `io_out` is
  **2.7–3.1 GB** — 30–40× larger than the final `Pairs.tsv.gz` output (70.9–82.7 MB,
  `n_kendall_pairs_bytes`). The `psutil`-based `io_out` counter is almost certainly capturing
  R's internal churn — repeated allocation/GC of intermediate `GRanges`/`DataFrame` objects,
  temp files, and swap-like paging — not the declared output file. I cannot pin the exact
  source further without instrumenting the R process directly (out of scope here), but this is
  a real, measurable "hidden write volume" 30–40× the nominal output size.
- **Declared vs. realized threading**: no `threads:` line at `workflow/rules/make_kendall_pairs.smk:1-31`
  → Snakemake/Slurm allocates 1 core. `cpu_time/wall` of 0.74–0.93 confirms the job is
  essentially single-threaded and CPU/parse-bound for most of its wall time, not idle waiting
  on I/O (ratio would be much lower if I/O-wait dominated).
- **Memory is heavily overprovisioned**: `resources: mem_mb=32*1000` (32 GB,
  `workflow/rules/make_kendall_pairs.smk:26-27`) against a measured `max_rss` of 4.7–8.5 GB —
  roughly 4–7× headroom unused. This matches the pattern already documented for
  `compute_kendall`'s stale 63 GB floor.
- **Repeated full-table materialization**: the 29-column `allPutative` table is parsed in full
  (`:21-23`) even though only `TargetGene` and the interval coordinates are ultimately used
  (final output keeps `chr, start, end, TargetGene, PeakName, PairName`, `:51-57`) — the other
  ~25 columns are read, stored in `mcols`, and then discarded. This is a real "read more than
  you need" cost, though I have not measured what fraction of the 240–380 MB `io_in` or the
  4.7–8.5 GB RSS this specifically accounts for.
- **Two full passes over ~10–11M rows for sort+dedup** (`order()` at `:47`, `duplicated()` at
  `:48`) rather than a single grouped operation — plausible but modest overhead; `data.table`
  or `dplyr::distinct()` equivalent would likely fuse this, though that is a fix-shaped
  observation and out of scope per the characterisation-only mandate.
- **Subprocess/interpreter startup**: not separately measured for this rule; at 178–348 s wall
  time, R/Bioconductor package load (`genomation`, `GenomicRanges`, `data.table`) is a small
  fixed cost relative to the dominant parse+join+sort, unlike the sub-second/toy-run rules
  called out in the briefing.

## 7. `t_est(n)`, `t_low(n)`, `t_high(n)`

Primary variable: `n_cand_pairs` (drives both the dominant parse and the overlap-join
complexity; `n_peaks` is a real but much smaller additive term in the join and cannot be
separated from `n_cand_pairs` with six points — see §5 — so it is folded into the constant
below rather than fit as a separate coefficient).

All formulas in seconds, `n` = `n_cand_pairs` (raw count, e.g. 10816106, not millions).

- **`t_low(n) = 1.9285e-5 · n`** — linear-class lower bound (optimistic: assumes the
  `log(n)` factor in the join/sort is swamped by the per-row parse constant at this scale).
  Constant calibrated so `t_low` matches the smallest cluster exactly:
  `t_low(10173077) = 196.12` (thp1_1, actual wall time).

- **`t_high(n) = 1.3369e-6 · n · log2(n)`** — log-linear-class upper bound (pessimistic:
  the NCList interval-tree join and the comparison-based sort dominate). Constant calibrated
  so `t_high` matches the largest cluster exactly:
  `t_high(11106687) = 347.59` (telohaec_crispri, actual wall time).

- **`t_est(n) = 9.797e-7 · n · log2(n)`** — best-guess class, same functional form as
  `t_high` (I judge `O(n log n)` the more code-faithful class, given `findOverlaps`'s NCList
  algorithm and `order()`'s comparison/radix sort at `:26-27` and `:47`), but calibrated at the
  **mean** cluster rather than the max:
  `t_est(10816106) = 247.72` (mean of the six observed wall times).

These are **not** blended point estimates; `t_low` and `t_high` are different functional
classes (linear vs. log-linear) each anchored at one real observation, and `t_est` is the
log-linear class anchored at the mean. Do not average them.

## 8. All three evaluated

| n (`n_cand_pairs`) | cluster | actual wall_s | `t_low(n)` | `t_est(n)` | `t_high(n)` |
|---|---|---:|---:|---:|---:|
| 10,816,106 (mean of 6) | — | 247.72 (mean) | 208.6 | 247.7 | 338.0 |
| 10,173,077 (min) | thp1_1 | 196.12 | 196.1 | 231.9 | 316.4 |
| 11,106,687 (max) | telohaec_crispri | 347.59 | 214.2 | 254.7 | 347.6 |

By construction `t_low` matches thp1_1 exactly and `t_high` matches telohaec_crispri exactly;
`t_est` over-predicts thp1_1 by ~18% and under-predicts telohaec_crispri by ~27%. That
±20–30% miss, on a variable that only spans 1.13×, is exactly the kind of gap the Phase 0
briefing warns about — it reflects run-to-run noise (Slurm scheduling, I/O contention) at
least as much as it reflects `n_cand_pairs` itself.

## 9. Modality

Multiome-only rule — one of the four rules (`make_kendall_pairs → generate_atac_matrix →
compute_kendall → arc_e2g`) that exist only on the Multiome branch of the DAG (Phase 0
briefing §2). There is no scATAC counterpart to compare class or constant against.

## 10. What the measurement cannot tell us

- **Cannot confirm the `n_cand_pairs` exponent.** `n_cand_pairs` spans only 1.13× across the
  six clusters, while wall time spans 1.94×; the benchmarks cannot distinguish `O(n)` from
  `O(n log n)` from `O(n^1.5)` at this range. The `t_low`/`t_high` bracket in §7–8 is a bound
  on the class derived from the code (parse cost vs. interval-tree-join-plus-sort cost), not
  something the six points can verify or falsify.
- **Cannot attribute the ~2.7–3.1 GB `io_out` to a specific R operation.** I flagged it as
  likely GC/temp-object churn from constructing and discarding large `GRanges`/`DataFrame`
  objects, but confirming that would require profiling the R process directly (e.g.
  `Rprof`/`tracemem`), which is outside a characterisation-only, no-code-execution task.
  Whatever is generating it, it is 30–40× the actual output file size and is a real,
  measured phenomenon, not a guess about its magnitude — only its cause is a guess.
- **Cannot separate `n_peaks`'s contribution from `n_cand_pairs`'s.** Both vary across
  clusters, are not obviously correlated with each other in this data (`jurkat_pma_cd3_4hr`
  has max `n_peaks` but near-median `n_cand_pairs` and below-median wall time), and the join
  cost genuinely depends on both per `findOverlaps`'s NCList algorithm
  (`make_kendall_pairs.R:26-27`). Six points, two candidate drivers, one noisy outcome — this
  cannot be decomposed empirically.
- **Cannot verify `genomation::readGeneric`'s internal algorithm** from this repo's code —
  it is a Bioconductor dependency. I have inferred its cost is `O(n_cand_pairs · 29-columns)`
  from its documented contract (parse a flat file into a `GRanges` with all metadata columns
  attached) and from the observed `io_in`/RSS numbers, not from reading its source, since it
  is not part of this repository and I did not find its installed copy in a scoped,
  non-filesystem-scanning search.
- **Cannot rule out that the negative-intercept linear fit in §5 is masking a real but small
  size effect.** A weak, noisy trend is consistent with either "size barely matters at this
  range" or "size matters somewhat but node-to-node variance is larger" — the six points
  cannot distinguish these.
- **Single measurement per cluster, no replicates** (per the addendum's methodological rule
  6) — every wall/cpu/IO number above is n=1 per cluster; there is no way to know how much of
  the 1.94× wall-time spread would survive on a second run of the same clusters.

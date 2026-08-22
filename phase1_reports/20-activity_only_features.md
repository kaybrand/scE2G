# Unit 20 — `activity_only_features.R`

## 1. Unit identity

- Script: `ENCODE_rE2G/workflow/scripts/feature_tables/activity_only_features.R` (92 lines).
- Rule: `activity_only_features`, `ENCODE_rE2G/workflow/rules/genomewide_features.smk:153-171`.
- Modalities: **both**, 6 jobs each, 12 measured jobs total. No `threads:` declared
  (`resources: mem_mb = ABC.determine_mem_mb`), confirmed at runtime in the Slurm log
  (`cpus_per_task=1`) — consistent with the briefing's "every ABC rule declares no threads"
  cross-cutting finding.
- Output: `{biosample}/ActivityOnly_features.tsv.gz`. Confirmed on disk (`thp1_1`,
  `zcat | head -1`): **19 columns** — `chr, start, end, name, class, TargetGene, TargetGeneTSS,
  TargetGeneEnsembl_ID, isSelfPromoter, isSelfGenic, CellType, distance, normalized_atac_prom,
  ABC.Score, normalized_atac_enh, numCandidateEnhGene, numTSSEnhGene, numNearbyEnhancers,
  is_ubiquitous_uniform`. Row count is `n_cand_pairs` (tautology, per briefing — not re-derived
  here, only used).

## 2. What the code actually does

The script's dependencies mix two R data idioms with different copy semantics, and this mixture
is the central fact needed to explain the memory profile.

1. `:12-16` — `config <- fread(...)`, `abc <- fread(snakemake@input$abc)` (the
   `EnhancerPredictionsAllPutative.tsv.gz` file, `n_cand_pairs` rows, 29 columns), then
   `abc$ABC.Numerator = abc$ABC.Score.Numerator; abc$ABC.Denominator = abc$ABC.Score/abc$ABC.Numerator`.
   These are `fread`/`data.table` operations. The two `$<-` column replacements are **whole-column
   replacements**, which since R 3.1+ (bundled data.table 1.15.4's own
   `datatable-reference-semantics.Rmd:71-82`, verified in the conda env at
   `.../conda/49efddd2545177ca66d0194982a8c760_/lib/R/library/data.table/doc/`) trigger only a
   **shallow** copy of the column-pointer vector, not a deep copy of `abc`'s ~2.8 GB of
   underlying data (see §6 for that number). This part of the script is cheap.
2. `:19-20` — build `input_features`, the (deduplicated) list of feature names actually requested
   by `feature_table.tsv` for this biosample. This gates which of the merges below add real
   columns vs. contribute nothing (see §6).
3. `:23-31` — subset `abc` to `core_cols` via `select()` (dplyr). This is the **first dplyr
   verb**: from here on the script is in tibble/data.frame land, not data.table's reference
   semantics.
4. `:36-38` — `output <- abc %>% select(...) %>% left_join(output, ., by = c("name","TargetGene"))`.
   `abc` is used here for the only time after load. **No `rm(abc)` follows.** `abc` (29 columns,
   `n_cand_pairs` rows, loaded from a file whose *uncompressed* size is ~2.8 GB for `thp1_1` —
   measured directly, see §6) stays live in the R global environment for the rest of the script,
   long after its one useful column-subset has been extracted.
5. `:41-49`, `:52-60`, `:63-71`, `:74-82` — four `if (file_test("-f", ...))` blocks, each `fread`s
   one auxiliary feature file and `left_join`s it onto `output` by `name`+`TargetGene` or `name`
   alone. On disk for every one of the 12 jobs inspected, **all four files exist** (verified for
   `igvf10_multiome/thp1_1/new_features/`: `NumCandidateEnhGene.tsv` 433.6 MB,
   `NumTSSEnhGene.tsv` 425.3 MB, `NumEnhancersEG5kb.txt` 3.4 MB, `SumEnhancersEG5kb.txt` 4.1 MB),
   so all four branches execute for all 12 measured jobs — this is not a conditional cost, it is
   the actual cost. None of these four intermediate tables is ever `rm()`'d either.
   - **A concrete, verified waste**: for `thp1_1`'s `feature_table.tsv`, the requested-feature
     list (`input_features`) contains `numNearbyEnhancers` but **not** `sumNearbyEnhancers` — so
     the `SumEnhancersEG5kb.txt` fread+join at `:74-82` contributes zero new columns (confirmed:
     the 19-column output header has no `sumNearbyEnhancers` field). The file is small (4.1 MB)
     so this costs little I/O, but it is read and joined unconditionally based on *file
     existence*, not on whether `feature_table.tsv` actually requests the feature — a real,
     if minor, wasted pass.
6. `:85-89` — one more `fread` (small `geneClasses` reference file, 1.0 MB on disk) and
   `left_join` by `TargetGene`.
7. `:92` — `fwrite(output, file = snakemake@output[[1]], sep="\t", na="NA", quote=FALSE)` — the
   output path ends in `.tsv.gz`; `fwrite()` in this env's data.table (1.15.2/1.15.4, confirmed
   via `DESCRIPTION` in the three conda envs under `.snakemake/conda/`) compresses internally via
   bundled zlib (no external `gzip` subprocess — confirmed no such call in the R script, and
   `NEWS.md` documents zlib-based `fwrite()` gz support, with an OpenBSD/zlib bug entry as of
   line 595 that only makes sense if the compression happens in-process).

**No loop, no cross join, no non-unique-key merge anywhere in this script.** Every join key
(`name`+`TargetGene`, or `name` alone) is unique per side (all of `abc`, `NumCandidateEnhGene`,
`NumTSSEnhGene` are one row per candidate pair; `NumEnhancersEG5kb`/`SumEnhancersEG5kb` are one
row per element `name`), so each `left_join` is a 1:1 merge — no row-count blow-up, consistent
with the row-count tautology in the briefing.

## 3. Size variable(s)

**Primary: `n_cand_pairs`** (10.17M–11.45M rows, Tier 2, range 1.13×) — every `fread`/`left_join`/
`fwrite` in the script scales in row count with `n_cand_pairs` (the input, the four merges, and
the output are all the same row count, by the documented identity).

**Secondary/multiplicative: effective column width.** There is no single "`n_feat_cols`" term
here (that variable names the *later* `genomewide_features.tsv.gz` file, not this rule's output),
but the *number of columns simultaneously resident in memory* across `abc` (29 cols) +
`NumCandidateEnhGene` (3 cols) + `NumTSSEnhGene` (3 cols) + `NumEnhancersEG5kb`/`SumEnhancersEG5kb`
(2 cols each) + `output` growing to 19 cols is what turns a `n_cand_pairs`-row, ~200 MB-compressed
table into a multi-GB memory footprint. This isn't a manifest column, it's a property of the code
(sum of columns across live copies), cited at the specific lines in §2.

No cross term with `k_genes` or `j_cand_elems`: the script never opens `EnhancerList.txt` or
`GeneList.txt`, and `geneClasses` (read at `:85`) is a small, near-constant reference file
(1.0 MB) that does not scale with any manifest variable.

## 4. Asymptotic class derived from code

**O(`n_cand_pairs`)** — linear. Every operation in the script is one of:
- a single linear pass over an `n_cand_pairs`-row (or, for the two auxiliary files, `n_cand_pairs`-
  row) table (`fread`, `$<-` column replace, column `select`), or
- a **hash join** (dplyr's `left_join`, on unique keys on both sides — 1:1, not 1:many) which is
  O(rows) to build the hash table plus O(rows) to probe it, not O(n²) and not O(n log n) (no sort
  is performed by dplyr's default hash-join implementation; unlike unit #17's
  `gen_num_candidate_enh_gene.py`, which explicitly calls `sort_values` three times, this script
  contains zero `sort`/`order`/`arrange` calls — confirmed by reading all 92 lines), or
- a single `fwrite` pass.

There is no nested loop, no cross join, and no per-row R-level iteration (everything is
vectorized data.table/dplyr). **The class is unambiguously linear from the code — I am not
hedging this one the way Tier 3 rules must be hedged.** What is genuinely uncertain is not the
exponent but the **constant**: the number of full-table deep copies dplyr materializes (five
sequential `left_join` calls, each allocating an entirely new object) and the fact that none of
the large intermediates (`abc`, `NumCandidateEnhGene`, `NumTSSEnhGene`) is ever freed with `rm()`
or `gc()`, so they all remain simultaneously resident through to the final `fwrite()`. That
constant, not the class, is what the ~7 GB RSS is telling us (§6).

## 5. Benchmark cross-check

`n_cand_pairs` spans only 1.13× across the six clusters (10.17M–11.45M), and it is **identical
between multiome and scATAC for the same cluster** (both modalities read the same upstream ABC
`EnhancerPredictionsAllPutative.tsv.gz` — verified: `n_cand_pairs` values in
`size_manifest_multiome.tsv` and `size_manifest_scatac.tsv` match exactly per cluster, e.g.
`jurkat` = 11,036,283 in both). This is the same trap the briefing warns about for Tier 2
variables: not enough range to fit an exponent, and here the fit genuinely fails —

| cluster (multi) | n_cand_pairs | wall s | rate (µs/pair) |
|---|---:|---:|---:|
| thp1_1 | 10,173,077 | 91.81 | 9.02 |
| thp1_2 | 10,278,325 | 95.48 | 9.29 |
| jurkat_pma_cd3_4hr | 10,854,698 | 89.79 | 8.27 |
| jurkat | 11,036,283 | 106.80 | 9.68 |
| telohaec_crispri | 11,106,687 | 99.58 | 8.97 |
| k562_crispri | 11,447,564 | 93.22 | 8.14 |

| cluster (scATAC) | n_cand_pairs | wall s | rate (µs/pair) |
|---|---:|---:|---:|
| thp1_1 | 10,173,077 | 82.19 | 8.08 |
| thp1_2 | 10,278,325 | 82.30 | 8.01 |
| jurkat_pma_cd3_4hr | 10,854,698 | 94.56 | 8.71 |
| jurkat | 11,036,283 | 85.18 | 7.72 |
| telohaec_crispri | 11,106,687 | **138.42** | **12.46** |
| k562_crispri | 11,447,564 | 88.93 | 7.77 |

A least-squares fit of wall time against `n_cand_pairs` gives **R² = 0.081 (multiome)** and
**R² = 0.162 (scATAC)** — essentially no linear signal, exactly as the briefing predicts for a
1.13× range: `k562_crispri` has the *largest* `n_cand_pairs` in both modalities but is not the
slowest job in either. **This is expected, not a contradiction of the code-derived class**: a
linear rule over a variable that ranges only 1.13× produces a time range dominated by run-to-run
noise (Slurm scheduling, filesystem latency, GC timing), not by the size term itself.

**Flagged disagreement**: `telohaec_crispri`/scATAC at 138.42 s is a clear outlier — 1.5–1.7×
every other job in either modality, despite `telohaec_crispri` not having the largest
`n_cand_pairs` (`k562_crispri` does, and runs in 88.93 s on scATAC). I cannot explain this from
the code (nothing in the script branches on cluster identity or modality); the most plausible
hypotheses, in order of how much they're actually supported here, are (a) Slurm/Lustre I/O
contention specific to that job's wall-clock window (io_in for this job is 1198.54 MB, the
highest of the 12 — 26% above the median ~940 MB, which is at least directionally consistent with
an I/O-contention explanation) or (b) GC-timing noise in a single-run measurement (per the
briefing, these are not replicated). I am not able to distinguish these from one measurement, and
say so rather than picking one.

**Modality**: excluding the `telohaec_crispri`/scATAC outlier, the per-pair rate is **8.0–9.7
µs/pair pooled across both modalities** — the multiome and scATAC rates interleave rather than
separate into two clusters (e.g. `k562_crispri`/scatac at 7.77 µs/pair is faster than
`thp1_1`/multiome at 9.02 µs/pair). This is expected: at this point in the DAG, multiome and
scATAC read the *same* ABC outputs (`EnhancerPredictionsAllPutative.tsv.gz`, same
`NumCandidateEnhGene`/`NumTSSEnhGene` byte sizes per cluster — verified identical in both
manifests for `thp1_1`: 433,609,152 bytes on both paths) and run the *identical* R script; the
ARC-feature branch point (which the briefing identifies as the source of large multiome/scATAC
divergence for `add_external_features`, `run_e2g_qnorm`, `gen_final_features`) has not yet been
reached. **Neither the class nor the constant differs by modality for this rule.**

## 6. Concrete overhead sources

**Memory (~6.5–8.1 GB RSS against a 19-column, `n_cand_pairs`-row, 189–238 MB compressed output):**
Measured directly: `EnhancerPredictionsAllPutative.tsv.gz` for `thp1_1` is 281,522,213 bytes
compressed on disk, and decompresses (`zcat | wc -c`) to **2,805,198,956 bytes (2.8 GB) of raw
text** — a ~10:1 gzip ratio. That 2.8 GB text file is `fread`'d in full as `abc` (29 columns) and,
per §2 item 4, is never released (`rm(abc)`/`gc()` do not appear anywhere in the script) even
though it is used only once, at `:36-38`. Simultaneously and also never released:
`NumCandidateEnhGene.tsv` (433.6 MB of raw text, `fread`'d at `:42-45`) and `NumTSSEnhGene.tsv`
(425.3 MB, `fread`'d at `:53-56`). All three of these — the largest input by far, plus the two
next-largest — remain live in R's global environment through every subsequent `left_join` and
through the final `fwrite()`. On top of that, `output` is rebuilt from scratch **five times**
(once per `left_join` call at `:38`, `:48`, `:59`, `:70`, `:81`, `:89` — six total merges
including the initial `select`), and dplyr's `left_join`/`select` are not reference-semantic the
way data.table's `:=` is (per data.table's own `datatable-reference-semantics.Rmd:71-82`, only a
whole-column `$<-` replacement gets the R-3.1+ shallow-copy optimization; a `left_join` is a fresh
tibble every time). The measured 6.5–8.1 GB is consistent with: `abc` (~2.8 GB text, likely
somewhat less as a typed data.table due to R's string-interning of the low-cardinality `chr`,
`class`, `CellType` columns, but `name`/`TargetGene` have up to ~158K/20.5K distinct values so
interning savings are partial) + `NumCandidateEnhGene`/`NumTSSEnhGene` (~0.86 GB combined text) +
one or more live `output` copies (up to 19 columns × `n_cand_pairs` rows) all resident
simultaneously, with none of it freed until the process exits after `fwrite()`. **This is a
missing-cleanup problem, not an algorithmic one**: the class is linear (§4), but the script never
drops an object it no longer needs, so the memory high-water mark is the *sum* of every large
table it ever loaded, not the size of the table it is currently working on.

**I/O amplification — write side.** `io_out` (3,105–3,220 MB across the 12 jobs) against the
actual compressed output (`actonly_bytes`, 189.0–238.3 MB) gives a **write amplification of
13.5×–14.7× (multiome) and 14.4×–15.1× (scATAC)** — remarkably *consistent* across all 6+6 jobs
(not noisy, unlike the wall-time rates above), which argues for a structural cause rather than
measurement jitter:

| cluster | modality | io_out (MB) | actonly_bytes (MB) | ratio |
|---|---|---:|---:|---:|
| jurkat | multi | 3105.75 | 211.83 | 14.66 |
| jurkat_pma_cd3_4hr | multi | 3120.96 | 214.59 | 14.54 |
| k562_crispri | multi | 3210.04 | 227.28 | 14.12 |
| telohaec_crispri | multi | 3105.09 | 216.14 | 14.37 |
| thp1_1 | multi | 2783.35 | 193.55 | 14.38 |
| thp1_2 | multi | 2792.44 | 195.80 | 14.26 |
| jurkat | scatac | 3051.22 | 199.03 | 15.33 |
| jurkat_pma_cd3_4hr | scatac | 3068.16 | 199.48 | 15.38 |
| k562_crispri | scatac | 3217.72 | 213.15 | 15.10 |
| telohaec_crispri | scatac | 3214.77 | 202.68 | 15.86 |
| thp1_1 | scatac | 2825.75 | 181.00 | 15.61 |
| thp1_2 | scatac | 2856.53 | 183.42 | 15.57 |

I looked for an explicit uncompressed intermediate written to disk and found none: the rule
declares exactly one `output:` (`ActivityOnly_features.tsv.gz`), the R script contains exactly one
`fwrite()` call, and this env's `data.table` (1.15.2/1.15.4) implements gz output via bundled
zlib in-process (confirmed via `NEWS.md`'s OpenBSD/zlib entry for `fwrite()`'s gz path; there is
no shell-out to an external `gzip`/`pigz` binary in the script, though both are present in
`encode_re2g.yml`, presumably for other rules). **I cannot fully confirm the mechanism from
available evidence** — the two candidate explanations are (a) `fwrite`'s internal buffered
write-then-compress path genuinely issuing several times more bytes to the storage layer than the
final compressed size (e.g. via preallocated per-chunk buffers or a page-cache accounting effect
on Lustre), or (b) the process's own memory pressure (§'s 6.5–8.1 GB RSS) causing swap writes that
get attributed to this PID's `io_counters().write_bytes` alongside the genuine output write. I
lean toward (a) being at least partly responsible, because the ratio is *tight* (13.5–15.9× with
no relationship to cluster size) rather than variable the way swap-driven writes typically would
be if they depended on transient memory-pressure timing; but I have not instrumented the process
directly and am reporting this as a flagged, quantified, but mechanistically unresolved overhead
source, not a confirmed one. For comparison, `io_in` (880–1,199 MB) tracks real reads very well:
summing the on-disk sizes of `abc` (compressed, 281–335 MB) + `NumCandidateEnhGene.tsv`
(uncompressed, ~410–488 MB) + `NumTSSEnhGene.tsv` (uncompressed, ~400–480 MB) + the two small
`NumEnhancersEG5kb`/`SumEnhancersEG5kb` files (~7–8 MB combined) for `thp1_1` gives ≈1,146 MB,
matching the measured `io_in` of 879.72 MB–1,126.29 MB across clusters closely enough that
`io_counters()` is clearly measuring something real for reads — which makes the write-side
amplification more, not less, notable, since the same accounting mechanism does not show
comparable inflation on the input side.

**Wasted read+join** (minor, quantified in §2 item 5): the `SumEnhancersEG5kb.txt` (4.1 MB)
fread+join executes on all 12 jobs but contributes zero output columns for at least the `thp1_1`
config, because `feature_table.tsv` requests `numNearbyEnhancers` but not `sumNearbyEnhancers`.
Small in absolute bytes, but it is representative of the script's general pattern: branches key
off *file existence*, not off whether `feature_table.tsv` actually needs the feature.

**Declared-vs-realised resources**: the Slurm log for `thp1_2` (jobid `39917090`) shows
`mem_mb=35502` (34.7 GiB) requested via `ABC.determine_mem_mb`, against a measured peak of
6,460–8,069 MB across the 12 jobs — a **4.3×–5.5× over-allocation**. This doesn't cost wall time
(it's a reservation, not a throttle) but is a real resource-accounting overhead worth noting for
Phase 2's cluster-utilization arithmetic, not for this rule's own runtime.

## 7. Three functions: `t_est(n)`, `t_low(n)`, `t_high(n)`

Variable: `n` = `n_cand_pairs` (rows). Class is linear, confirmed unambiguously from the code
(§4) — unlike Tier 3 rules, the uncertainty here is in the **constant** (I/O + copy overhead per
row), not the exponent, so `t_low`/`t_high` bound the realised per-row rate rather than a
different asymptotic class. All three below pool multiome and scATAC, since §5 established the
constant does not differ by modality for this rule.

```
t_est(n)  = 8.9e-6  * n      # seconds; mean observed rate (95.7 s / 10.816M pairs, pooled 12-job mean)
t_low(n)  = 7.7e-6  * n      # seconds; fastest observed rate (jurkat, scATAC: 85.18 s / 11,036,283 pairs)
t_high(n) = 1.25e-5 * n      # seconds; slowest observed rate, INCLUDING the telohaec_crispri/scATAC
                              #   outlier (138.42 s / 11,106,687 pairs). Excluding that one outlier,
                              #   the realistic high bound is 9.7e-6 * n (jurkat, multiome: 106.80 s /
                              #   11,036,283 pairs) — reported here because §5 could not confirm the
                              #   outlier is size-driven rather than a one-off scheduling/I/O artifact.
```

No additive constant term is fit separately: with only a 1.13× range in `n` and R² ≈ 0.08–0.16
for a linear fit, splitting the observed time into a slope and an intercept would manufacture
false precision. The rate already includes fixed costs (conda activation, interpreter startup,
the small `geneClasses`/`NumEnhancersEG5kb`/`SumEnhancersEG5kb` reads) folded into the per-row
constant, because at this data range they cannot be separated from the size-dependent term.

## 8. All three evaluated

| | `n_cand_pairs` | `t_low` | `t_est` | `t_high` (incl. outlier) | `t_high` (excl. outlier) |
|---|---:|---:|---:|---:|---:|
| mean (6-cluster mean) | 10,816,106 | 83.3 s | 96.3 s | 135.2 s | 104.9 s |
| min (`thp1_1`) | 10,173,077 | 78.3 s | 90.5 s | 127.2 s | 98.7 s |
| max (`telohaec_crispri`) | 11,106,687 | 85.5 s | 98.8 s | 138.8 s | 107.7 s |

(For reference: `thp1_1` is also the min-size cluster per the briefing's convention;
`telohaec_crispri` is the max — same clusters, same convention as other units' reports.)

## 9. Modality: class or constant?

**Neither differs.** Both the asymptotic class (linear, from identical code) and the realised
constant (8.0–9.7 µs/pair pooled, excluding one outlier) are the same between multiome and
scATAC for this rule. This is expected given §2/§5: `activity_only_features` runs before the
ARC/Kendall branch point the briefing identifies as the source of multiome/scATAC divergence for
downstream rules (`add_external_features`, `run_e2g_qnorm`, `gen_final_features`); at this stage
of the DAG both modalities are reading byte-identical upstream ABC files and running the same
script.

## 10. What the measurement cannot tell us

- **Cannot confirm linear scaling empirically.** `n_cand_pairs` spans only 1.13× across the six
  clusters, and the six-point fit has R² = 0.08–0.16 — indistinguishable from noise at this
  range. The linear class in §4 rests entirely on reading the code (no sorts, only hash joins on
  unique keys, no nested loops); the benchmarks neither confirm nor contradict it, they are just
  too flat a range to say anything about the exponent at all.
- **Cannot attribute the `telohaec_crispri`/scATAC 138.42 s outlier to a specific cause.** I have
  one data point per (cluster, modality) — no replicates (per the briefing) — so I cannot tell
  scheduling/I/O contention apart from GC-timing noise apart from a genuine but unmeasured
  size effect. I flagged this rather than smoothing it into the fit.
- **Cannot fully confirm the mechanism behind the 13.5×–15.9× write amplification.** I ruled out
  an external `gzip` subprocess and a second declared output file, and confirmed `io_in` tracks
  real reads accurately for the same jobs (so `io_counters()` is not generically unreliable here)
  — but I could not directly instrument the `fwrite()` call itself (e.g. with `strace`) to
  distinguish "gzip's internal buffering genuinely writes more bytes to storage than the final
  file" from "swap activity under memory pressure gets attributed to the same PID." Both remain
  live hypotheses; I lean toward the former given how tight and size-independent the ratio is,
  but say so as a lean, not a conclusion.
- **Cannot decompose the ~7 GB RSS into per-object contributions.** I estimated `abc`'s
  in-memory size is somewhat less than its 2.8 GB raw-text equivalent due to R's string-interning
  of low-cardinality columns, but I did not instrument the R session (e.g. with `lobstr::obj_size`
  or `Rprofmem`) to get an exact breakdown across `abc`, `NumCandidateEnhGene`, `NumTSSEnhGene`,
  and the five successive `output` copies. The claim that all of these are simultaneously live
  and none is freed is fully supported by reading the code (no `rm()`/`gc()` anywhere in the
  script); the claim that this *specific set* of objects sums to exactly the measured 6.5–8.1 GB
  is a plausibility argument, not a verified accounting.
- **Cannot separate a fixed per-job overhead from the size-dependent term** in `t_est`/`t_low`/
  `t_high` — see §7. At `n_cand_pairs`'s observed range, any split between "intercept" and
  "slope" would be curve-fitting noise, not signal.

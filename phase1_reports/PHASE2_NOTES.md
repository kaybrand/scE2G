# Phase 2 working notes

Findings established during the Phase 1 launch that change how Phase 2 must be computed.

## 0. TWO CORRECTIONS THAT AFFECT EVERY UNIT — from unit 4

### 0a. Cross-cluster timing is confounded by CPU generation. This supersedes the "46% noise floor" note.

Unit 4 pulled `sacct` node assignments and stratified by CPU model. The result:

- **9 jobs on AMD Siena 8224P: CV 0.49% across a 3.30× `n_frag` range.** Perfectly linear.
- Both `telohaec_crispri` jobs landed on **AMD Rome 7502**; the fast `k562_crispri` outlier on
  **Sapphire Rapids**.
- **Same-hardware replicate agreement is ≤1.3%.** The 28% gap is *pure CPU generation*, not a
  size-dependent noise floor.

**This supersedes §"A hard limit on how tight any `t_est` can be" in `QUEUE.md`,** which took unit
9's 46%/21% replicate disagreement as an irreducible contention floor. It is largely *explainable*,
not irreducible — and it plausibly explains the `telohaec_crispri` outlier that units 9, 12, 15, 18,
20, 21 and 25 each independently flagged and could not account for. Phase 2 should treat unexplained
cross-cluster excess as a **hardware covariate first**, and only then as contention.

**Consequence for the deliverable:** the apparent `n^1.049`–`n^1.07` superlinearity reported by units
5 and 4 is a node-generation artifact. On fixed hardware both rules are **dead linear** in `n_frag`.
Do not propagate a fitted exponent above 1.0 for the fragment-stage rules.

### 0b. `max_rss` is ALSO unreliable — it undercounts 1.5–4.7×

Unit 4 compared Snakemake's `max_rss` against the true cgroup peak: on `telohaec_crispri`
`create_neighborhoods` really peaked at **1.33 GB, not the reported 290 MB**. Same
sampling/child-process blind spot as `cpu_time` and `io_*` (see §3a) — the polling monitor simply
misses peaks between samples and in short-lived children.

So the *only* benchmark columns that survive scrutiny are **`s` (wall) and `cpu_time`**, and even
`cpu_time` is zero for `hover_plots`. Treat `max_rss` as a lower bound, never as "the" peak. The
290 MB figure was also never telohaec-specific — `thp1_2`, the *smallest* cluster, reports 287 MB.

### 0c. `cpu_time/wall ≈ 0.98` does NOT prove the work is serial

The baton records `create_neighborhoods` as "CPU-bound on one core, not I/O-bound," which is right
as far as it goes. Unit 4 sharpens it: **`AllocCPUS=1` for all 12 jobs, so CPU/wall ≈ 1.0 is a
cpuset ceiling, not a measurement of serialisability.** Run the same `count_tagalign` pipe with 4
CPUs and it reaches **CPU/wall = 1.98** (36.2 s wall / 71.5 s CPU) versus 74.6 s wall = 74.5 s CPU
at 1 CPU. The pipeline stages *are* concurrent; the single-CPU allocation forces them to time-share
so their costs **add**, a 1.66× rule-wide penalty.

This generalises to every rule in the pipeline measured at `AllocCPUS=1`: a CPU/wall near 1.0 is
consistent with *either* genuinely serial work *or* concurrent work squeezed onto one core, and the
benchmark cannot distinguish them. Phase 2 must not infer "unparallelisable" from CPU/wall ≈ 1.

## 1. "min / median / max cluster size" is NOT well defined globally

The plan specifies recomputing the critical path at the min (`thp1_1`, 287 MiB ATAC) and max
(`telohaec_crispri`, 4.0 GiB ATAC) cluster. **That ranking holds only for the fragment stage.**
Computed argmin/argmax per size variable from both manifests:

| Variable group | argmin | argmax | ratio |
|---|---|---|---:|
| `n_frag`, `n_umi`, `n_cells` | `thp1_1` / `thp1_2` | **`telohaec_crispri`** | 15.1× / 12.0× / 6.0× |
| `n_peaks` | `thp1_1` | **`jurkat_pma_cd3_4hr`** | 1.38× |
| `n_cand_pairs`, `n_kendall_pairs`, and every `*_bytes` table | `thp1_1` | **`k562_crispri`** | 1.13–1.26× |
| `j_cand_elems`, `enh_list_bytes` | `k562_crispri` | **`thp1_1`** | 1.01× (**inverted**) |
| `n_called_all` (multiome) | `k562_crispri` | `thp1_2` | 1.03× |
| `n_called_all` (scATAC) | `thp1_1` | `k562_crispri` | 1.13× |

Three consequences:

1. **There are three different "biggest clusters"** depending on pipeline stage:
   `telohaec_crispri` for the fragment/MACS2 stack, `jurkat_pma_cd3_4hr` for peak sorting, and
   `k562_crispri` for everything from candidate pairs onward. Unit 2 caught this independently for
   `run_e2g_qnorm` (max is `k562_crispri`, not `telohaec_crispri`).
2. **`j_cand_elems` is inversely ranked**: `thp1_1`, the *smallest* cluster by fragments, has the
   *most* candidate elements (158,548 vs 156,420 for `k562_crispri`). This is the `nStrongestPeaks`
   cap working as designed, and it means a rule whose cost is driven by `j_cand_elems` gets very
   slightly *slower* on the smallest dataset.
3. So the deliverable's "critical path at avg / min / max cluster size" must be computed
   **per-cluster across all six**, and the min/max presented as *which cluster's path is longest*,
   not as an assumed size ordering. Presenting `telohaec_crispri` as "the max case" would be wrong
   for every rule downstream of `make_candidate_regions`.

## 1a. The plan's `max_cell_count` premise is wrong, and the conclusion drawn from it is too

`cozy-finding-candle.md` § Phase 2 states: *"`max_cell_count: 20000` subsamples the top four but
not the thp1 pair, so downstream Kendall/ATAC-matrix cost does **not** track fragment-file size
linearly."* Unit 8 checked both halves and **both are false**:

- **The subsampling never triggers.** All six clusters have 2,593–15,555 cells, every one below the
  20,000 threshold. Slurm logs confirm `"Number of cells to use" == "Number of cells in RNA matrix"`
  for all six. The branch (`generate_atac_matrix.R:69-72`) is dormant in this dataset — and it lives
  *inside* this rule, not upstream of it.
- **ATAC-matrix cost does track `n_frag`, almost perfectly linearly.** OLS gives **R² = 0.996**
  (`t = 238.49 + 2.398e-6·n_frag`), stable under leave-one-out (R² ≥ 0.88, slope 2.20–2.42e-6).

This is one of the very few units where the exponent is genuinely well-supported empirically,
because `n_frag` really does span 15.1×. Unit 8 also showed the code reason: Signac's
`CreateFragmentObject` runs an unconditional `md5sum()` over the entire raw fragment file *plus its
`.tbi` index* (verified against Signac's source), which is a mandatory full linear byte pass before
any matrix work begins.

So Phase 2 should not carry the plan's assumption that downstream cost is decoupled from fragment
volume. For `generate_atac_matrix` it is tightly coupled. Unit 8 also confirmed the matrix is
**sparse**, not dense — on-disk `.rds` size tracks `n_frag` (1.4× spread) far better than
`j_cand_elems · n_cells` (3.0× spread), so nnz scales with fragments rather than the peak×cell grid.

Caveat unit 8 stated honestly: adding `n_kendall_pairs` as a second regressor yields a physically
impossible negative coefficient, confirming that term is real in the code but **not separable** from
`n_frag` at six correlated points (`n_frag` vs `n_cells`, r ≈ 0.95).

## 2. Job-level graph tables

Parsed from the regenerated, checkpoint-resolved DAGs (see `PHASE1_ADDENDUM.md` §F):

```
analysis/dag_{multiome,scatac}_nodes.tsv    node_id, rule, cluster, n_pred, n_succ
analysis/dag_{multiome,scatac}_edges.tsv    from_id, from_rule, to_id, to_rule, cluster
```

Cluster identity was propagated along edges, since dot labels carry a wildcard only where it is
first introduced. Result is clean and symmetric: **25 nodes per cluster (multiome), 21 per cluster
(scATAC)**, with exactly **5 genuinely shared nodes** in both — `all`, `save_reference_configs`,
`abc_generate_chrom_sizes_bed_file`, `hover_plots`, `plot_stats`.

**The six clusters are otherwise fully independent parallel branches.** They touch only at the two
fan-in nodes (`plot_stats`, then `hover_plots`) and at the shared chrom-sizes node. So the critical
path is: the longest single-cluster chain, plus the two fan-in tails. Under the plan's
unconstrained-parallelism assumption, the six branches do not contend — which makes the *shared
resource* caveat (Lustre bandwidth during the fragment stage, when all six clusters are reading and
writing simultaneously) the main threat to that assumption, and worth flagging explicitly.

Note the checkpoints (`features_required`, `basic_features_required`) are absent from the resolved
DAGs but are real benchmarked barrier nodes — 6 jobs each per modality. They must be inserted
manually, per unit 26:

### 2a. Checkpoint barrier structure (from unit 26) — insert these into both graphs

```
make_biosample_feature_table  ──►  checkpoint features_required   (per cluster)
                                        ├──►  make_external_features_config
                                        │        └──►  [if ARC] make_kendall_pairs
                                        │                → generate_atac_matrix
                                        │                → compute_kendall → arc_e2g   (24 jobs, multiome only)
                                        ├──►  element_and_gene_summaries
                                        └──►  get_stats_per_model_per_cluster

                                   checkpoint basic_features_required  (per cluster)
                                        └──►  generate_num_candidate_enh_gene
                                              generate_num_tss_enh_gene
                                              generate_num_sum_enhancers     (18 jobs per modality)

both checkpoints ──►  add_external_features   (needs outputs from both branches)
```

Corrections this makes to my earlier framing:
- `features_required` has **three** out-edges, not one. Only the first transitively pulls in the
  ARC chain.
- `basic_features_required` is **always True in this dataset**, so the three `generate_num_*` rules
  run in *both* modalities — 18 jobs per modality, 36 total. They are **not** part of the 24-job
  multiome/scATAC delta, which is purely the ARC chain. (They were missing from the stale multiome
  DAG only because the checkpoint was unresolved, not because they don't run.)
- `save_reference_configs` has **no `input:` at all** — zero upstream edges — and nothing consumes
  its outputs, but `rule all` (`Snakefile:123-125`, via `Snakefile:89`) depends on
  `expanded_biosample_config.tsv` directly. So it is a leaf that still feeds `all`: schedulable at
  any time, full slack, cannot be on the critical path.

### 2b. Checkpoint bookkeeping latency is real critical-path time

Unit 26 pulled Slurm log timestamps rather than relying on the benchmark row, and found both
checkpoints execute in a **nested double pattern**: outer sbatch dispatch, a 2–6 s gap, then a
nested "storing output in storage" DAG-update step of 1–2 s. **Total user-visible latency is 3–7 s
per checkpoint job, of which under 0.5 s is the actual work.**

This matters because the benchmark row shows only the ~0.03–0.4 s of compute. With 6 jobs of each
checkpoint per modality sitting as hard barriers on every cluster's chain, Phase 2 must charge the
3–7 s dispatch latency, not the benchmarked compute — otherwise it will score these nodes as free
when they are not. Unit 26 flags its `t_high = 60 s` as an unverified placeholder for queue
contention beyond what this run's logs show.

## 2c. Attribution rule — no double counting, resolved

**Snakemake's `benchmark:` measures the whole `shell:` block, including every subprocess.** So all
`bedSplitSort.sh` and `sort` seconds are **already inside** `frag_to_tagAlign`'s and
`frag_to_norm_bigWig`'s measured totals. Phase 2 must **not** create a separate node for
`bedSplitSort.sh`, and must not add unit 9's numbers to unit 24's.

Unit 24 handled this correctly by adopting unit 9's already-fitted `O(n_frag log n_frag)`
calibration for `frag_to_tagAlign` rather than re-fitting, since both units are measuring the same
DAG node from the same 12 points. Treat unit 9 as the *mechanism* report for a node that unit 24
*owns*. The same logic applies to unit 23 (`make_biosample_feature_table`) versus unit 26's
checkpoint — those are genuinely distinct rules with distinct benchmark rows, so both are charged.

## 2d. `add_external_features` — the class genuinely differs by modality

The one place in the pipeline where modality changes the complexity class, not just the constant.
Unit 21 traced it to `merge_external_features.R:112-115`: `unique_sources` is empty for scATAC
(its external-features config is header-only, 60 bytes, verified on all 6 clusters), so the
`if (length(unique_sources) > 0)` guard means the **only** block containing any join or sort code
(`:116-130`) never executes on the scATAC path.

- **scATAC:** join/sort term is exactly **zero**. Column count stays at 18 before and after.
- **multiome:** one pass doing a full `fread` of the 630–778 MB ARC table, a
  `GRanges`/`findOverlaps` overlap join, 5× `data.table` group-bys on a 19-column composite key,
  4× `merge()`, and a full `order()` sort — roughly **10 sort/hash passes over ~10–11M rows**.
  Columns go 19 → 24.

The decisive arithmetic against a pure I/O explanation: bytes moved rise only **3.2×** from scATAC
to multiome, but wall time rises **15.5×**. So this is compute in an extra code path, not volume.
Peak RSS 18.7–21.0 GB is accounted for by simultaneous residency of both tables, two `GRanges`
objects, the cbind result, and the 5-element `agg_list` that `mapply` must fully materialise before
`Reduce`.

For Phase 2 this means `add_external_features` needs **two different cost functions**, not one
function with a modality constant.

## 3. Do not use `io_in` / `io_out`

Confirmed unreliable in both directions by units 2, 7, 9, 11, 12, 14, 15, 16 and 17 — undercounts
of 2–41× (page-cache hits) and overcounts of 100–400×, plus `0.00` for jobs that demonstrably wrote
megabytes. Root cause is psutil not seeing child process trees, same blind spot that zeroes
`hover_plots`' `cpu_time`. **Use the on-disk `*_bytes` manifest columns or `ls`-measured sizes.**

### 3a. RETRACTION — the "uncompressed stream" explanation of `io_out` is wrong

An earlier revision of this file recorded unit 22's conclusion that `io_out` measures the
**uncompressed** text volume written into the compression stream. **That is not what it measures,
and this retraction supersedes it.** Verified against the Snakemake 9.16.3 source in the
`run_snakemake9` env, `snakemake/benchmark.py`:

```
:326        io_in, io_out = 0, 0                       # zeroed EVERY poll
:337        for proc in chain((self.main,), self.main.children(recursive=True)):
:361-363        ioinfo = proc.io_counters(); io_in += ioinfo.read_bytes; io_out += ioinfo.write_bytes
:399        self.bench_record.max_rss = max(self.bench_record.max_rss or 0, rss)   # <- max
:404-405    self.bench_record.io_in = io_in                                        # <- plain ASSIGN
```

Two consequences, and they are not symmetric:

**(a) The reported number is a LAST-POLL SNAPSHOT, not a total and not a max.** Unlike `max_rss`,
`io_in`/`io_out` are overwritten each poll. Only processes **still alive at the final poll** are
summed; every child that already exited contributes nothing. For a `shell:` rule that is a pipeline
of short-lived binaries, the last poll may catch only the tail process. This is the true cause of
`io_out = 0.00` on jobs that demonstrably wrote megabytes (units 14, 16) — the same psutil
child-process blind spot that zeroes `hover_plots`' `cpu_time` and crashes the monitor on
`generate_chrom_sizes_bed_file`.

**(b) `read_bytes`/`write_bytes` are BLOCK-DEVICE counters** (`/proc/<pid>/io`), so page-cache hits
do not count at all. Independent second source of undercount, and it explains why the shortfall is
worst for files an upstream rule has just written — unit 8 measured `io_in` at 9.9–52.6% of the
fragment file, degrading with size; unit 21 saw it on a warm cache from upstream `arc_e2g`.

**What this means for the "write amplification" findings.** Because `write_bytes` counts real device
writes, pipe traffic and in-process gzip **cannot** inflate it. So a large `io_out` indicates
genuine device I/O:

- **Unit 9's ~27 GB `io_out` for `bedSplitSort.sh` is most likely REAL**, not an artifact — the
  script genuinely writes 25 uncompressed per-chromosome BED files to Lustre. An earlier revision of
  this file wrongly talked that finding down. Unit 9's round-trip I/O concern stands.
- Units 20 and 22's 5.5–15× excess over the *compressed* output is now **unexplained again** and
  needs re-checking: `data.table::fwrite` compressing in-process should produce device writes close
  to the compressed size. Possible real causes are uncompressed temp files or R spill. Flag for the
  verify pass; do not carry unit 22's stated mechanism.

**MECHANISM FOUND — unit 6 closes this out, and it confirms the retraction above.** Unit 6
instrumented `/proc/self/io` per stage and found `arc_e2g`'s 7–8× `io_out` amplification is
**86–88% real physical temp-file spill**: both readers decompress the entire `.gz` to a temp file
and then mmap it back. Measured `wchar` of 2,805.2 MB and 1,516.8 MB match the uncompressed input
sizes *exactly*, while the `fwrite` step's own `wchar` is 630.1 MB — i.e. correctly compressed. The
identity `io_out ≈ ABC_uncomp + Kendall_uncomp + arc_bytes` holds to ±1% across all five valid
clusters, and the spill goes to `/tmp` (verified local XFS RAID, not tmpfs).

So the sequence of positions resolves cleanly: unit 22's "uncompressed bytes into the compression
stream" was wrong (pipes and in-process gzip never reach the device), the retraction's reasoning was
right, and the true mechanism is decompress-to-temp-then-mmap. Large `io_out` on the R rules is
**real I/O to a real filesystem**, and Phase 2 should treat it as such — including for units 20 and
22, whose 5.5–15× excess most likely has the same cause.

**Revised guidance (supersedes the blanket "unusable"):** `io_in`/`io_out` are unreliable for any
multi-subprocess `shell:` rule and systematically undercount cached reads. For single long-lived
process rules they do measure device I/O — which is a *different quantity* from the logical
read/write volume most reports wanted, and must not be compared against on-disk file sizes as
though it were the same thing. Prefer the on-disk `*_bytes` manifest columns or `ls`-measured sizes
for volume; use `io_*` only as weak evidence about device traffic, never as a headline number.

## 3b. CROSS-CUTTING: serialisation, not algorithm, dominates the big-table rules

Two independent Tier A units now converge on the same conclusion, by different methods. This is
shaping up to be the largest single cross-cutting theme in Phase 1, and Phase 2 should treat it as
a property of the pipeline rather than of any one rule.

| Unit | Rule | Algorithm share | Serialisation share | Method |
|---|---|---|---|---|
| 2 | `run_e2g_qnorm` | qnorm sort < 1% | **gzip ≥ 52% (multiome) / ≥ 47% (scATAC) of CPU** | measured compress vs decompress rates on real slices |
| 3 | `create_predictions` | **3%** (15.1–18.0 s of a 549–634 s job) | **96%** = `pd.concat` + 3× `to_csv(compression="gzip")` | the script's own per-chromosome self-timer, plus a Slurm microbenchmark on 500k real rows |

Unit 3 decomposed the write further: **52.5% pandas `%.6f` float formatting, 47.5% zlib level 9.**
Both units independently found the same root cause — **pandas passes an empty `compression_args`, so
`GzipFile` defaults to level 9.**

**CORRECTION — the compression asymmetry is NOT a pipeline constant, and "~40×" was my
over-generalisation of unit 2's number.** Unit 2 measured ~15–21 MB/s compressing versus ~194 MB/s
decompressing on pandas/level-9 output. I propagated that ratio into unit 6's brief as a general
prior; unit 6 measured its own and found something quite different:

| | compress | decompress | asymmetry |
|---|---|---|---|
| unit 2 — pandas `to_csv`, gzip level 9 | 15–21 MB/s | ~194 MB/s | ~10–13× |
| unit 6 — R `data.table::fwrite` to `.gz` | **44.1 MB/s** | **244 MB/s** | **5.5×** |

Unit 6 cross-checked its deflate rate against `gzip -6` on identical bytes (40.2 MB/s), confirming
`fwrite` uses in-process zlib with no subprocess. **Applying unit 2's rate to an R rule overstates
its compression cost by roughly 2×.** Compression share must be measured per rule, not inherited:
unit 6's own deflate cost is 28% of the job (81.41 s of 287.9–357.5 s), reproduced to 0.07% by
differencing `fwrite` to `.gz` against `fwrite` to plain `.tsv` — real and the largest single line
item by 2.1×, but well short of unit 2's ≥52%. The theme (serialisation over algorithm) holds
across all three units; the *magnitude* does not transfer.

Unit 3's microbenchmark (35.14 s/M rows) agrees with its production-fitted constant (37.50 s/M) to
6.3%, which is the strongest independent validation of a cost constant anywhere in Phase 1.

**Implication for Phase 2:** for every rule that writes a `n_cand_pairs`-scale gzipped table, the
cost function is dominated by output bytes and float formatting, not by the rule's nominal
algorithm. That is why `feat_bytes` and `arc_bytes` track timing better than row or column counts
across units 2, 15, 21 and 22. Unit 6 (`arc_e2g`, ~500–600 MB output, still running) is the natural
next confirmation.

## 3c. `create_predictions` writes ~20% redundant bytes, and both modalities are replicates

Unit 3 verified per-column that `activity_base_enh` is bit-identical to `activity_base`, and that
`ABC.Score`/`powerlaw.Score` and their numerators are bit-identical (400,000/400,000 rows checked),
because `compute_score` is **called twice with identical arguments** on the powerlaw path
(`predictor.py:62-67` versus `:69-74`) — one fully redundant 14M-row hash groupby. `enh_idx`/
`gene_idx` also leak into the 282 MB output because their `drop` sits inside HiC-only functions.

Separately, and important for the graph: **`create_predictions` is not two modality arms but true
n=2 replicates.** `Expression` is all-NaN in `GeneList.txt` in *both* modalities — the ABC path
never sees RNA — so outputs are bit-identical for 5 of 6 clusters. Phase 2 should charge one cost
function for this node in both graphs.

`getVariantOverlap.py` **does** execute (imported `predict.py:9`, called unconditionally at `:298`,
~5–8 s), but its output is referenced nowhere in the repo, is not a declared Snakemake output, and
its `.tmp` is never unlinked. Factual dependency note only.

## 3d. The tagAlign is traversed SEVEN times across two rules — the pipeline's largest single cost

Units 4 and 5 together account for the fragment-stage cost, and they reach the same machinery.

**Unit 4 (`create_neighborhoods`, the 112-minute leader) — resolved.** 97.4–99.6% of wall in all 12
jobs is **six full single-threaded streaming passes over the tagAlign**. `processCellType` calls
`count_features_for_bed` three times (`neighborhoods.py:126`, `:136`, `:280`) — same tagAlign, three
different target BEDs — and each call does two independent traversals:

- **Pass A** `count_tagalign:448-460` — `bedtools intersect -u | bedtools coverage -counts -sorted | awk` (~61%)
- **Pass B** `count_total:567` → `count_tagalign_total:674-684` — `zcat | grep -E | wc -l` (~39%)

**Pass B computes the identical scalar three times with no memoisation — 26% of total wall is
provably redundant recomputation.**

The 290 MB / 2 hour apparent paradox is resolved: the fragment stream never enters Python. Python
only ever holds the `j_cand_elems` ≈ 157k and `k_genes` = 20,531 frames, so memory stays flat while
CPU is O(`n_frag`). There is **no** `iterrows`/`.apply`/per-element loop on this path —
`count_bam:440-443` and `count_bigwig:467-484` do contain them, but `run_count_reads` dispatches on
filename and never reaches them for a tagAlign input.

**No `O(n_frag · j)` cross term exists.** Unit 4 proved this two ways: the three passes cost the same
to ±5% despite pass 3's `-a` BED having 7.7× the rows, and a direct 1-core experiment gave
`-a`=candidateRegions 74.60 s versus `-a`=GeneList.bed 73.99 s — 0.8% apart.

**Combined with unit 5's one further pass, the pipeline traverses each tagAlign 7 times** (4 bedtools
+ 3 zcat). For `telohaec_crispri` that is ~136 GB of text decoded on one core in
`create_neighborhoods` alone.

**Phase 2 must treat `t_high` as the LIVE case for `telohaec_crispri`, not a pessimum.** Unit 4's
`t_high` predicts 6780 s against measured 6725 s and 6740 s — the observed outcome, twice, because
both telohaec jobs drew the slow CPU generation.

Note also `create_neighborhoods` **hard-codes `mem_mb = 32*1000`** (`neighborhoods.smk:27`) — it does
not use `determine_mem_mb`, unlike the rules in §"over-provisioning". Against a true 1.33 GB cgroup
peak that is 24× over-provisioned.

## 4. A cost-driving table whose size is NOT in the size manifest

`scE2G_predictions.tsv.gz` — the unthresholded output of `run_e2g_qnorm` — is the dominant read for
units 11, 12 and 15, and the manifest has no column for it. Measured directly by unit 15:

| | multiome | scATAC | ratio |
|---|---|---|---:|
| on-disk size | 742.7–889.2 MB | 414.0–476.5 MB | **1.84×** |
| column count | 27 | 20 | 1.35× |
| row count | `n_cand_pairs` | `n_cand_pairs` | 1.00× |

The manifest's `feat_bytes` ratio (~2.69×) is a **poor proxy** for it, and `n_feat_cols` (24/18)
describes a different table (`genomewide_features.tsv.gz`). Phase 2 should carry this table's size
as its own variable for units 11, 12 and 15.

## 5. Where the modality difference actually comes from

For every rule downstream of scoring, row counts are *bit-identical* between modalities; only table
*width* differs. So modality effects are width/byte effects, not row effects — and scATAC is
consistently **faster** despite calling ~2× more links, because its tables are narrower. Confirmed
independently by units 2, 11, 12, 15 and 16.

The one exception is `add_external_features` (unit 21), where the ~15× gap is a
presence/absence effect: the merge happens in multiome and is skipped entirely in scATAC.

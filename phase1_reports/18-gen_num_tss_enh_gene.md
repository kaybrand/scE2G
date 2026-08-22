# Unit 18 — `gen_num_tss_enh_gene.py`

## 1. Unit identity

- Script: `ENCODE_rE2G/workflow/scripts/feature_tables/gen_num_tss_enh_gene.py` (57 lines).
- Rule: `generate_num_tss_enh_gene`, `ENCODE_rE2G/workflow/rules/genomewide_features.smk:91-113`.
- Modalities: **both**, 6 jobs each (one per cluster), 12 measured jobs total. No `threads:`
  declared (runs at 1). `resources: mem_mb = partial(ABC.determine_mem_mb, min_gb=32)` — a
  32 GB floor against an observed peak `max_rss` of 5.0–6.2 GB across all 12 jobs (§5/§6).
- Sibling scripts in the same directory (owned by other units, not re-derived here):
  `gen_num_candidate_enh_gene.py` (#17) and `gen_num_sum_nearby_enhancers.py` (#19).

## 2. What the code actually does

`main()` (`:46-53`) reads the whole predictions file with no column restriction:

```python
pred_df = pd.read_csv(abc_predictions, sep="\t")   # :47 — no usecols
```

`EnhancerPredictionsAllPutative.tsv.gz` has **29 columns** (verified: `zcat ... | head -1 | tr
'\t' '\n' | wc -l` → 29), but `determine_num_tss_enh_gene()` only ever touches 7 of them
(`chr, start, end, distance, TargetGeneTSS, name, TargetGene`). Sibling script #17
(`gen_num_candidate_enh_gene.py:46-47`) restricts to 6 columns via `usecols=[...]` on the same
input file; this script does not. See §6.

`determine_num_tss_enh_gene()` (`:6-37`), over the full `n_cand_pairs`-row frame:

1. `:10-11` — `midpoint = (start+end)/2`, `new_end = midpoint + distance`. Vectorized arithmetic,
   O(`n_cand_pairs`).
2. `:14-20` — boolean mask `downstream_enh` (gene TSS upstream of the enhancer midpoint) with a
   `.loc` reassignment of `start`/`new_end` for that subset. Vectorized, O(`n_cand_pairs`), no
   loop.
3. `:23-25` — `pred_df[["chr","start","new_end","name","TargetGene"]].to_csv(extended_enhancers, ...)`
   writes **one row per row of `pred_df`**, i.e. exactly `n_cand_pairs` rows, not a deduplicated
   set of candidate elements. Verified directly on disk: `extendedEnhancerRegions.txt` for
   `thp1_1` has 10,173,078 lines (1 header + 10,173,077 data rows), matching
   `size_manifest_multiome.tsv`'s `n_cand_pairs=10173077` for that cluster exactly.
4. `:29-35` — a single `os.system()` call running a shell heredoc:
   ```
   printf "name\tgene\tcount\n" > out_file;
   bedtools intersect -a extended_enhancers -b ref_gene_tss -wa -c | cut -f4,5,6 >> out_file
   ```
   `-a` = the file just written in step 3 (`n_cand_pairs` rows). `-b` = `ref_gene_tss` =
   `config['gene_TSS500']`. No `-sorted` flag is passed.

`ref_gene_tss`, resolved at run time from `igvf10_multiome/config/ENCODE_rE2G_config.yml`, points
at `.../scE2G/resources/genome_annotations/CollapsedGeneBound.hg38.GENCODEv43GeneSymbol.TSS500bp.bed`
(the known "`encode_re2g_dir` points elsewhere" gotcha — this run loads the reference from the
main `scE2G` checkout, not `scE2G_optimize`). Verified: **20,532 lines total (1 header + 20,531
data rows) — exactly `k_genes` (20,531, constant per the manifest on every cluster in both
modalities)**. Also verified this file is **not globally sorted** (`sort -k1,1 -k2,2n -c` reports
"disorder" at line 51) — consistent with the script never requesting `-sorted` and confirming
bedtools cannot silently be using the sorted sweep path here even if it wanted to.

## 3. Size variable(s)

- **`n_cand_pairs`** (10,173,077–11,447,564 across the 6 clusters, 1.13× — Tier 2): the row count
  of `pred_df`/`abc_predictions`, and therefore of the bedtools `-a` input `extended_enhancers`.
  Confirmed by direct `wc -l` on the intermediate, not inferred.
- **`k_genes`** (exactly 20,531 on all 6 clusters, both modalities — Tier 3, constant by
  construction): the row count of `ref_gene_tss`, the bedtools `-b` input. Confirmed by `wc -l`
  on the actual reference file resolved at run time.
- **Cross term**: `bedtools intersect -a extended_enhancers(n_cand_pairs rows) -b
  ref_gene_tss(k_genes rows)` is a genuine join between these two variables — see §4 for the
  class this join actually has.

**Correcting the unit brief's "expected drivers" note.** The brief for this unit names
`j_cand_elems` (~156K–159K, from `EnhancerList.txt`) as a co-driver "joined onto `n_cand_pairs`
rows." I read the code and checked the rule's inputs (`genomewide_features.smk:92-93`): this
script never opens `EnhancerList.txt`, and neither the script nor the rule reference `j_cand_elems`
anywhere. The bedtools `-a` file (`extended_enhancers`) has `n_cand_pairs` rows, verified above by
line count — not `j_cand_elems` rows. The genuine Tier-3 co-driver present in this script's own
code is `k_genes` (via `ref_gene_tss`), not `j_cand_elems`. The brief's illustrative arithmetic
("~3.2×10⁹ operations" for a naive `O(j·k)` pass) is `j_cand_elems × k_genes` ≈
157,000 × 20,531 ≈ 3.22×10⁹ — that number describes a different hypothetical script (one whose
`-a` side is the deduplicated candidate-element list), not this one. The analogous number for what
this script actually runs is `n_cand_pairs × k_genes` ≈ 11.1M × 20,531 ≈ **2.28×10¹¹** for a naive
nested-loop equivalent, versus the interval-join's actual `n_cand_pairs × log2(k_genes)` ≈
11.1M × 14.33 ≈ **1.59×10⁸** (§4) — the same qualitative point the brief makes (interval join vs.
cross join is orders of magnitude apart), just with the correct operand and an even larger gap
once corrected. I'm flagging this explicitly rather than silently substituting the brief's numbers,
per the "derive from the code" rule.

## 4. Asymptotic class derived from code

**O(`n_cand_pairs`) [pandas + I/O, §2 items 1–3] + O(`n_cand_pairs`·log2(`k_genes`) +
`k_genes`·log2(`k_genes`)) [bedtools interval join, §2 item 4]**, and because `n_cand_pairs` ≫
`k_genes` at every measured point, the join term dominates: overall
**O(`n_cand_pairs`·log2(`k_genes`))**.

Justification for the join term, `gen_num_tss_enh_gene.py:31` (`bedtools intersect -a
extended_enhancers -b ref_gene_tss -wa -c`, no `-sorted`):

- bedtools' documented default behavior without `-sorted` is to load file B into memory and build
  a searchable (interval-tree-like) structure per chromosome from it, then stream file A and query
  the structure once per interval — as opposed to the `-sorted` linear merge-sweep, which is not
  used here (confirmed neither input is pre-sorted; see §2).
  [Bedtools overview docs](https://bedtools.readthedocs.io/en/latest/content/overview.html);
  interval-tree query complexity of `O(log N + m)` per query against a tree of `N` intervals with
  `m` matches is standard for this class of structure (see also the
  [Bedtk paper's complexity discussion](https://academic.oup.com/bioinformatics/article-pdf/37/9/1315/50359653/btaa827.pdf)).
- Here `-b = ref_gene_tss` = `k_genes` rows (small, constant), `-a = extended_enhancers` =
  `n_cand_pairs` rows (large, mildly variable). So: tree build over `-b` is
  O(`k_genes`·log2(`k_genes`)) ≈ 20,531 × 14.33 ≈ 2.94×10⁵ — negligible. Querying `-a` against the
  tree is O(`n_cand_pairs`·log2(`k_genes`) + total matches).
- Because the command uses `-wa -c` (count mode, `:31`), the reported output is **exactly one row
  per `-a` interval regardless of how many TSS windows it overlaps** — verified directly:
  `NumTSSEnhGene.tsv` for `thp1_1` has exactly the same line count as `extendedEnhancerRegions.txt`
  (10,173,078, header included). So the "total matches" term in the query cost does not inflate
  output size or downstream I/O even though individual enhancer windows can span megabases
  (observed directly: `extendedEnhancerRegions.txt`'s first row for `thp1_1` is
  `chr1  10012  4655044  intergenic|chr1:10012-10512  AJAP1` — a ~4.6 Mb window).
- This is the interval-join class, not a naive cross join: a hypothetical "for each of
  `n_cand_pairs` rows, scan all `k_genes` TSS windows" implementation would be
  O(`n_cand_pairs`·`k_genes`) ≈ 2.28×10¹¹ operations at these sizes. bedtools' interval tree turns
  the `k_genes` side of the join from a linear scan into a `log2(k_genes)` ≈ 14.33 lookup — about
  a **1,433×** reduction in the join term.

No nested loop, no per-row `.apply`/`.iterrows()`, and no pandas merge on a non-unique key
anywhere in this script. The pandas portion (`:10-25`) is fully vectorized; the actual join is
delegated entirely to bedtools' compiled interval-tree implementation via a single `os.system()`
call.

## 5. Benchmark cross-check

All 12 measured jobs (from `benchmarks_{multiome,scatac}.tsv`, rule `generate_num_tss_enh_gene`):

| Modality | Cluster | `n_cand_pairs` | wall s | cpu_time s | cpu/wall | io_in MB | io_out MB | max_rss MB |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| multiome | jurkat | 11,036,283 | 141.69 | 73.72 | 0.520 | 343.26 | 1021.37 | 5878.6 |
| multiome | jurkat_pma_cd3_4hr | 10,854,698 | 87.08 | 73.85 | 0.848 | 8.07 | 947.18 | 5037.4 |
| multiome | k562_crispri | 11,447,564 | 107.62 | 85.24 | 0.792 | 17.97 | 1130.51 | 5192.4 |
| multiome | telohaec_crispri | 11,106,687 | 129.70 | 56.69 | 0.437 | 416.93 | 793.63 | 6015.5 |
| multiome | thp1_1 | 10,173,077 | 74.89 | 45.17 | 0.603 | 45.70 | 641.52 | 5427.8 |
| multiome | thp1_2 | 10,278,325 | 117.24 | 63.36 | 0.541 | 99.89 | 879.23 | 5516.5 |
| scATAC | jurkat | 11,036,283 | 82.14 | 74.97 | 0.913 | 276.00 | 1037.11 | 5871.8 |
| scATAC | jurkat_pma_cd3_4hr | 10,854,698 | 85.81 | 72.25 | 0.842 | 5.65 | 965.09 | 5069.8 |
| scATAC | k562_crispri | 11,447,564 | 87.21 | 74.57 | 0.855 | 292.16 | 1010.91 | 6190.8 |
| scATAC | telohaec_crispri | 11,106,687 | 210.61 | 138.20 | 0.656 | 430.59 | 1009.98 | 6021.3 |
| scATAC | thp1_1 | 10,173,077 | 126.76 | 53.05 | 0.419 | 13.95 | 752.94 | 5448.2 |
| scATAC | thp1_2 | 10,278,325 | 128.95 | 53.03 | 0.411 | 18.86 | 736.91 | 5060.5 |

`n_cand_pairs` is **identical between multiome and scATAC for the same cluster** (verified:
same six values in both `size_manifest_*.tsv` files) — see §9. `k_genes` = 20,531 on every one of
these 12 jobs; it never takes a second value anywhere in Phase 0.

**`n_cand_pairs` spans only 1.13× (10.17M–11.45M) — essentially the same narrow range the
briefing warns is unfittable for Tier 2 variables, and `k_genes` spans 1.00× (never varies at
all).** With 12 points and two variables that barely move, this benchmark set is a weak check on
the O(`n_cand_pairs`·log2(`k_genes`)) class derived in §4, not a confirmation of it: e.g. `k562_
crispri` has the largest `n_cand_pairs` of the six clusters (11,447,564) but is not the slowest job
in either modality (107.62 s multiome, 87.21 s scATAC — both mid-pack), while `telohaec_crispri`
(11,106,687, third-largest) produced the single slowest job of all 12 (210.61 s, scATAC). Per-row
wall time (`wall_s / n_cand_pairs`) ranges from 7.36×10⁻⁶ s/row (`thp1_1`, multiome) to
1.90×10⁻⁵ s/row (`telohaec_crispri`, scATAC) — a **2.6× spread in the "constant"** despite only a
1.13× spread in `n_cand_pairs`, i.e. more of the variance in wall time is unexplained noise than is
explained by the one variable that actually moves.

`cpu_time`/wall ratios (0.41–0.91, mostly 0.4–0.66) show a substantial and inconsistent fraction of
wall time (9–59%) is not active single-core compute — I/O wait, subprocess fork/exec (`os.system`
→ shell → `printf`/`bedtools`/`cut`), or Slurm/Lustre variance, not the O(`n_cand_pairs`·log2(`k_
genes`)) compute itself.

**One point stands out as a likely outlier rather than a real scaling signal**: scATAC
`telohaec_crispri` has wall=210.61 s and cpu_time=138.20 s, while its multiome twin — same cluster,
same `n_cand_pairs` (11,106,687), same `ref_gene_tss` — has wall=129.70 s and cpu_time=56.69 s.
Both wall *and* cpu_time are elevated together on an identical workload, which points to something
external (node contention, a noisy neighbor, cold Lustre metadata) rather than the algorithm
running differently on identical input. This is treated as noise in §7/§10, not folded into the
class-bound calibration.

Net: **the benchmarks cannot confirm the log2(`k_genes`) term, and only weakly constrain the
linear-in-`n_cand_pairs` term** — consistent with the brief's central warning, even though the
specific driver identified in §3 (`k_genes`, not `j_cand_elems`) differs from what the brief named.

## 6. Concrete overhead sources

- **Unfiltered column read.** `pd.read_csv(abc_predictions, sep="\t")` at `:47` has no `usecols`,
  so pandas parses and materializes all 29 columns of `EnhancerPredictionsAllPutative.tsv.gz` even
  though only 7 are used downstream. Sibling script `gen_num_candidate_enh_gene.py:46-47` restricts
  to 6 columns on the *same input file*; this script does not. As with unit #17's finding on
  `usecols`, this mainly costs memory and dtype-conversion work rather than I/O time (the gzip
  stream must still be fully decompressed and tokenized either way), but the memory cost is real:
  `max_rss` is 5.0–6.2 GB against a 32 GB declared floor, and the input file itself is only
  281.5–335.3 MB compressed (`n_cand_pairs_bytes`).
- **No streaming between the pandas stage and the bedtools stage.** `pred_df` must be fully
  materialized and flushed to `extendedEnhancerRegions.txt` (`:23-25`, 638–722 MB observed on disk,
  e.g. `thp1_1`=638,454,897 bytes, `telohaec_crispri`=694,392,580 bytes, `jurkat`=690,931,461
  bytes) before the `bedtools` subprocess launched at `:29-35` can read any of it back. This is a
  full serialize-to-disk / deserialize-from-disk round trip through the filesystem for data that
  never leaves the same host, rather than a pipe from the Python process directly into bedtools.
- **Declared vs. realized memory.** `resources: mem_mb=partial(ABC.determine_mem_mb, min_gb=32)`
  (`genomewide_features.smk:102`) reserves a 32 GB floor; observed peak `max_rss` across all 12
  jobs is 5.0–6.2 GB — roughly **5–6× over-provisioned**, in the same pattern already flagged
  project-wide for `compute_kendall`'s stale 63 GB floor.
- **No threading.** No `threads:` declared → Snakemake schedules this at 1 core; `bedtools
  intersect` itself is also single-threaded, so there is no parallelism to lose here (unlike
  `frag_to_norm_bigWig`/`frag_to_tagAlign`, whose declared thread counts are aspirational).
- **`os.system()` subprocess overhead.** The shell heredoc at `:29-34` spawns a shell plus three
  more processes in a pipe (`printf`, `bedtools`, `cut`). At a 75–210 s job length this overhead is
  small in relative terms (unlike the sub-second checkpoint rules elsewhere in this project where
  interpreter/conda startup dominates), but it is a real, unnecessary extra fork/exec chain — the
  header line could be written directly from Python and the column selection done with pandas
  instead of `cut`.
- **I/O volume roughly tracks the real write, unlike unit #17's rule.** `io_out` (641–1,130 MB per
  job) is in the same order of magnitude as the sum of the two files this rule actually writes per
  cluster (`extendedEnhancerRegions.txt` 638–722 MB + `NumTSSEnhGene.tsv` 425–462 MB ≈ 1.06–1.18
  GB combined) — the benchmark's `io_out` counter is a reasonably faithful (if not exact) proxy
  for write volume here.
- **`io_in` is not a reliable proxy for read volume.** E.g. `jurkat_pma_cd3_4hr` reports
  `io_in`=8.07 MB (multiome) / 5.65 MB (scATAC) despite reading a ~281–335 MB compressed predictions
  file plus a ~1.7 MB `gene_TSS500` reference — almost certainly a page-cache effect (the file was
  recently touched by an upstream rule in the same job graph and much of the "read" was served from
  RAM, not Lustre), not evidence that the read was actually free. Don't use `io_in` at face value
  for this rule; use `n_cand_pairs_bytes` (the compressed input size, confirmed to match the actual
  `EnhancerPredictionsAllPutative.tsv.gz` file size on disk byte-for-byte for `thp1_1`) instead.

## 7. Three functions

Class per §4 is O(`n_cand_pairs`·log2(`k_genes`)) from the code, dominated at these sizes by the
linear term in `n_cand_pairs` because `log2(k_genes)` ≈ 14.33 is a small, and — critically —
**entirely unvarying** multiplier in this dataset (§3, §10). All three functions below are
calibrated on wall time (`s` column) so they sit on a single consistent basis and form a genuine
envelope; they intentionally do **not** try to also bound the Slurm/I/O noise discussed in §5/§6,
which is a separate, unmodeled source of variance layered on top.

```
t_est(n_cand_pairs)              = 1.0637e-5 * n_cand_pairs                              [seconds]
t_low(n_cand_pairs)              = 7.363e-6  * n_cand_pairs                              [seconds]
t_high(n_cand_pairs, k_genes)    = 8.960e-7  * n_cand_pairs * log2(k_genes)              [seconds]
```

Calibration (both modalities pooled, 12 jobs; §9 finds no class or size difference between them):

- `t_low`'s coefficient is the smallest observed `wall_s / n_cand_pairs` across the 12 jobs
  (multiome `thp1_1`: 74.89 s / 10,173,077 = 7.363e-6 s/row) — a pure-linear, best-case-constant
  lower bound on the class.
- `t_high`'s coefficient uses the largest *clean* (non-outlier; see §5) observed
  `wall_s / n_cand_pairs` (multiome `jurkat`: 141.69 s / 11,036,283 = 1.2836e-5 s/row), divided by
  `log2(20,531) = 14.3255` to isolate a per-(row·log₂-gene) rate: `1.2836e-5 / 14.3255 = 8.960e-7`.
  This is the worst-case class bound that takes the interval-tree's `log2(k_genes)` query cost at
  face value rather than folding it into the linear constant.
- `t_est`'s coefficient is the mean of all 12 `wall_s / n_cand_pairs` ratios (1.0637e-5 s/row),
  the best single point-estimate rate given that the log-factor cannot be separated from the
  linear constant at a fixed `k_genes` (§10) — i.e. `t_est` treats the class as effectively
  O(`n_cand_pairs`) for point-estimate purposes, with any true `log2(k_genes)` dependence folded
  into the fitted constant.

These are **not an exponent fit to the 6/12 points** — §5 shows the data cannot support that. They
pin a linear floor and an n·log2(k) ceiling to the cleanest observed extremes, with the class
itself (linear vs. n·log₂k) coming from reading the bedtools join in §4, not from curve-fitting.

## 8. Evaluated at mean / min / max

| Point | Cluster | `n_cand_pairs` | `k_genes` | `t_low` | `t_est` | `t_high` |
|---|---|---:|---:|---:|---:|---:|
| mean (across 6 clusters) | — | 10,816,106 | 20,531 | 79.6 s | 115.1 s | 138.8 s |
| min | `thp1_1` | 10,173,077 | 20,531 | 74.9 s | 108.2 s | 130.6 s |
| max (designated) | `telohaec_crispri` | 11,106,687 | 20,531 | 81.8 s | 118.2 s | 142.6 s |

(Note: `telohaec_crispri` is the project's designated "max" cluster overall, but is not the
cluster with the numerically largest `n_cand_pairs` among the six — that is `k562_crispri` at
11,447,564. Evaluated at `k562_crispri` instead: `t_low`=84.3 s, `t_est`=121.8 s, `t_high`=146.9 s.)

For reference, the actual measured wall times at the two named clusters: `thp1_1` 74.89 s
(multiome) / 126.76 s (scATAC) — the multiome value sits almost exactly on `t_low`, the scATAC
value sits between `t_est` and `t_high` on the *same* `n_cand_pairs`; `telohaec_crispri` 129.70 s
(multiome, between `t_est` and `t_high`) / 210.61 s (scATAC, **above `t_high`** by ~48%, the
flagged outlier from §5). The spread across modalities at fixed `n_cand_pairs` and fixed `k_genes`
is itself evidence of how much non-algorithmic noise this envelope cannot capture.

## 9. Does the class or constant differ between modalities?

**Neither.** `n_cand_pairs` is identical between multiome and scATAC for the same cluster
(verified: `size_manifest_multiome.tsv` and `size_manifest_scatac.tsv` give the same six
`n_cand_pairs` values, and `k_genes`=20,531 in both), because `EnhancerPredictionsAllPutative.tsv.gz`
and `gene_TSS500` are both produced/configured independently of the ARC/Kendall branch that
`checkpoint features_required` decides per cluster. This rule is not among the briefing's four
rules whose cost differs materially by modality (`add_external_features`, `run_e2g_qnorm`,
`gen_final_features`, `get_stats_per_model_per_cluster`) — and the data agrees: the observed
wall-time differences between the two runs of the same cluster (e.g. `telohaec_crispri` 129.70 s
multiome vs. 210.61 s scATAC; `thp1_1` 74.89 s vs. 126.76 s) run in **both directions** across the
six clusters, with no consistent sign — the signature of run-to-run scheduling/I/O noise on an
identical workload, not a modality-dependent algorithm.

## 10. What the measurement cannot tell us

- **`k_genes` is exactly 20,531 on every one of the 12 measured jobs, in both modalities — it never
  takes a second value anywhere in Phase 0.** No amount of re-analysis of this dataset can
  distinguish a genuine O(log2(`k_genes`)) or O(`k_genes`) term in this script's cost from a fixed
  constant folded entirely into the linear-in-`n_cand_pairs` coefficient. This is not the
  usual "six points, narrow range" caveat — it is that the second variable in the join has
  literally one observed value in the entire corpus. The `log2(k_genes)` term in `t_high` (§7) is
  a code-derived structural assumption about how bedtools' interval tree scales, not something
  any Phase 0 measurement could confirm or refute even in principle. Phase 2 should treat `t_est`
  (which has no explicit `k_genes` dependence) as the honest point estimate, and `t_high` as an
  assumption to stress-test only if a future dataset changes the reference gene annotation.
- **`n_cand_pairs` spans only 1.13×** (10.17M–11.45M) — narrow enough, combined with the
  `cpu_time`/wall ratios of 0.41–0.91 and the 2.6× spread in per-row wall time despite the 1.13×
  spread in `n_cand_pairs` (§5), that the linear-vs-n·log₂k distinction in the class is not
  verifiable from these 12 points either. A flat line fits about as well as any curved one would
  at this range.
- **One job (scATAC, `telohaec_crispri`) is a likely outlier** — both wall (210.61 s) and cpu_time
  (138.20 s) are elevated relative to its multiome twin on identical input (129.70 s / 56.69 s).
  With n=1 per (cluster, modality) combination, there is no way to confirm this is node
  contention/noise rather than a real (if unexplained) cost driver; it is excluded from the class
  calibration in §7 on the working assumption that it is noise, but that assumption is unverifiable
  with this data.
- **Cannot apportion the 9–59% of wall time beyond `cpu_time`** between Slurm scheduling latency,
  Lustre I/O wait on the ~638–722 MB `extendedEnhancerRegions.txt` write and ~425–462 MB
  `NumTSSEnhGene.tsv` write, subprocess fork/exec overhead, and conda/interpreter startup — the
  benchmark instrumentation gives wall time, `cpu_time`, RSS, and coarse I/O counters, not a
  breakdown of wait states.
- **`io_in` cannot be trusted at face value for this rule** (§6) — it is inconsistent with the
  known compressed input size by more than an order of magnitude in some jobs, almost certainly a
  page-cache artifact from upstream rules in the same job graph having recently touched the file.
- **Cannot extrapolate beyond the observed ~10.17M–11.45M `n_cand_pairs` / fixed-20,531 `k_genes`
  range** — there is no cluster in this corpus large enough, and no alternate gene annotation
  small or large enough, to reveal whether the interval-join's `log2(k_genes)` term would ever
  become visible against the O(`n_cand_pairs`) linear terms even if `k_genes` did vary.

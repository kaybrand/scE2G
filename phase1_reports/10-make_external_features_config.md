# Unit 10 — `make_external_features_config`

## 1. Unit identity

- Script: `workflow/scripts/feature_computation/format_external_features_config_sc.R` (47 lines)
- Rule: `make_external_features_config`, `workflow/rules/add_external_features.smk:37-52`
- Modalities: **both**. 6 jobs each (one per cluster: jurkat, jurkat_pma_cd3_4hr, k562_crispri,
  telohaec_crispri, thp1_1, thp1_2). 12 measured jobs total.
- Upstream context (owned by unit #26, not re-derived here): `checkpoint features_required`
  (`add_external_features.smk:1-22`) reads each cluster's `feature_table.tsv` and writes
  `to_generate.txt` containing exactly `"ARC"` or `"Neither"` — verified on disk: `"ARC"` for
  all six multiome clusters, `"Neither"` for all six scATAC clusters (per `PHASE1_BRIEFING.md` §2).

**Naming trap avoided:** the briefing's cross-modality cost table (§2) lists an
`add_external_features` row at 12.6–21.3 min (multiome) vs 0.9–1.3 min (scATAC), ratio ~15×.
That row belongs to a *different, similarly-named* rule — `rule add_external_features` at
`ENCODE_rE2G/workflow/rules/genomewide_features.smk:179` — which I confirmed against
`benchmarks_multiome.tsv`/`benchmarks_scatac.tsv`: multiome 756.96–1275.25 s (12.6–21.3 min),
scATAC 52.53–75.87 s (0.9–1.3 min), exact match. **My rule, `make_external_features_config`, is
not that row.** Its own measured times are two to three orders of magnitude smaller (below) and
are not the ~15× signal described in the briefing.

## 2. What the code actually does

`format_external_features_config_sc.R`:

1. `dataset_config = fread(snakemake@params$dataset_config)` — reads
   `configs/tables/igvf10_cell_clusters_{multiome,scatac}.tsv`, a **6-row, 9-column** static
   table (`format_external_features_config_sc.R:7`). Trivial and identical cost every job.
2. Filters that table to the current cluster (`:12`) and checks whether it has an
   `external_features_config` column. For this project's config tables it does not — I verified
   `configs/tables/igvf10_cell_clusters_multiome.tsv` has 9 columns and no
   `external_features_config` column — so the entire branch at `:13-24` (which would `fread` a
   second file and do `file.exists()` stats in a loop) **never executes**. `efc` is always
   initialized empty at `:26-29`.
3. `new_feature = snakemake@input[1]` (`:6`) is a **file path string**, not file content.
   `feature_inputs` in the rule is populated by `features_to_generate()`
   (`add_external_features.smk:26-34`): for multiome it resolves to
   `.../ARC/EnhancerPredictionsAllPutative_ARC.tsv.gz` (a ~400–630 MB file per the size manifest);
   for scATAC it resolves to `RESULTS_DIR` itself (a directory, since `to_generate=="Neither"`).
4. The script does **`grepl()` string matching on that path** (`:40,43`) — it never opens or
   reads the ARC/Kendall file's content. If the path matches
   `EnhancerPredictionsAllPutative_ARC.tsv.gz`, it appends a fixed 5-row data frame (`ARC_rows`,
   `:32-36`) whose `source_file` column is set to that same unread path string. If it matches
   `Pairs.Kendall.tsv.gz` it appends a 1-row frame (`Kendall_row`, `:38`). Otherwise `efc` stays
   empty (0 rows).
5. `fwrite(efc, ...)` (`:47`) writes the tiny result.

**Dead-branch note (code-reading only, not a cost factor):** the checkpoint at
`add_external_features.smk:18-20` sets `final_val = "ARC"` whenever *either* the ARC or Kendall
flag is true, and never emits a standalone `"Kendall"` value. Consequently the
`elif val == "Kendall"` branch in `features_to_generate()` (`add_external_features.smk:29-30`) and
the corresponding `grepl(..., "Pairs.Kendall.tsv.gz", ...)` branch in this script (`:40-42`) are
currently unreachable given how the checkpoint computes `to_generate`. This is a control-flow
observation, not a fix proposal (out of scope), and it does not change the cost analysis: the
live code path is strictly binary — 5 rows (ARC, multiome) or 0 rows (Neither, scATAC).

**Confirmed on disk** — output sizes and content match this analysis exactly:
`igvf10_multiome/*/external_features_config.tsv` are 868–928 bytes, 5 data rows, with
`source_file` literally containing the unread ARC file's path;
`igvf10_scatac/*/external_features_config.tsv` are all exactly 60 bytes (header row only, 0 data
rows).

## 3. Size variables

**None of the manifest's size variables (`n_frag`, `n_umi`, `n_cells`, `n_peaks`,
`n_cand_pairs`, `j_cand_elems`, `k_genes`, `arc_bytes`, etc.) govern this script's cost.** The
script reads a fixed 6×9 table and a path string; it never opens the large file whose size those
variables describe. The only "size" that varies at all is the number of *rows this script itself
constructs* — 5 (ARC) or 0 (Neither) — which is a hard-coded constant set by the branch, not a
function of any manifest variable, and is far too small to matter (a handful of short strings).

## 4. Asymptotic class derived from code

**O(1)** in every manifest variable. There is no loop, join, sort, or merge over any
cluster-scale data structure. The only iteration in the whole script is the
`for (i in 1:nrow(efc))` at `:17`, which is dead code for this project (never entered, since
`external_features_config` is not a column in the cell-clusters table) and even if entered would
iterate over ≤ a handful of externally-supplied rows, not over cluster data. Cost is dominated by
R/`data.table`/`dplyr` interpreter and library load, plus Slurm dispatch — not by algorithmic
work over the declared inputs.

## 5. Benchmark cross-check

All 12 measured jobs (`s` = wall seconds, `cpu_time` = seconds):

| Modality | Cluster | s | cpu_time | io_in (MB) | io_out (MB) |
|---|---|---:|---:|---:|---:|
| multiome | jurkat | 8.69 | 0.68 | 55.99 | 0.06 |
| multiome | jurkat_pma_cd3_4hr | 0.88 | 0.18 | 0.00 | 0.00 |
| multiome | k562_crispri | 0.91 | 0.18 | 0.00 | 0.00 |
| multiome | telohaec_crispri | 1.79 | 0.29 | 1.78 | 0.00 |
| multiome | thp1_1 | 0.89 | 0.18 | 0.00 | 0.00 |
| multiome | thp1_2 | 0.96 | 0.16 | 0.00 | 0.00 |
| scATAC | jurkat | 10.46 | 0.54 | 37.45 | 0.00 |
| scATAC | jurkat_pma_cd3_4hr | 0.93 | 0.17 | 0.00 | 0.00 |
| scATAC | k562_crispri | 0.91 | 0.16 | 0.00 | 0.00 |
| scATAC | telohaec_crispri | 10.76 | 0.64 | 7.39 | 0.00 |
| scATAC | thp1_1 | 6.86 | 0.43 | 50.91 | 0.00 |
| scATAC | thp1_2 | 10.78 | 0.63 | 5.86 | 0.00 |

Pooled (12 jobs): mean wall 4.65 s, min 0.88 s, max 10.78 s. `cpu_time` never exceeds 0.68 s in
any job — actual compute is sub-second everywhere; everything above ~1 s of wall time is
overhead external to the script's own work.

**This confirms the code-derived O(1) class,** and the cross-check is unambiguous because the
observed spread is not correlated with any size variable at all — it isn't just "too narrow to
tell," it is actively inconsistent with a size dependence:
- The largest ARC input file this rule is declared to depend on varies ~1.5× across multiome
  clusters (26.2M–394.8M `n_frag`, `arc_bytes` up to 630 MB at `thp1_1` alone), yet multiome wall
  times (0.88–8.69 s) don't track cluster size at all — `telohaec_crispri` (the largest cluster,
  394.8M `n_frag`) is the *fastest* non-trivial multiome job (1.79 s), while `jurkat` (a
  mid-sized cluster, 75.3M `n_frag`) is the *slowest* (8.69 s).
- scATAC jobs, which write **zero data rows** (60-byte header-only output, confirmed on disk),
  are on average slower (mean 6.95 s) than multiome jobs, which write 5 rows (868–928 bytes,
  mean 2.35 s). If cost scaled with output size or with the branch taken, this would run
  backwards.
- `io_in` (0–56 MB) does not correspond to any declared input or output of this rule (outputs are
  <1 KB; the only genuinely-read input is the 6-row config table). It plausibly reflects
  cold-cache reads of R/`data.table`/`dplyr` shared-library files from the network filesystem
  when the conda env isn't already warm on that node — filesystem/scheduling noise, not
  algorithmic I/O.

**cpu_time/s ratio warning, per the briefing:** these are sub-second-to-low-single-digit-second
jobs; `cpu_time/s` here runs from ~0.078 (jurkat, multiome) up to ~5–19× the other direction
(`s/cpu_time`, e.g. scATAC jurkat 10.46/0.54 ≈ 19.4). This is exactly the sampling-noise pattern
the briefing warns about for sub-second rules — it is **not evidence of parallelism** (nothing in
the code launches threads or subprocesses), and it is not evidence of any real 15–20× slowdown
either; it reflects the benchmark sampler's poor resolution for jobs this short, compounded by
Slurm scheduling latency and page-cache state. **Conclusion: cost is dominated by R interpreter
and library startup (plus scheduling/cache variance), not by input size — exactly the expected
finding for a config-formatting rule.**

## 6. Concrete overhead sources

- **R/library startup**: `library(data.table)` and `library(dplyr)` loads, every job, every time
  — likely the largest single fixed cost, though not separable from Slurm/conda overhead with
  this instrumentation.
- **Conda environment activation** (`envs/sc_e2g.yml`): dispatched via Snakemake's `conda:`
  directive; environment is reused (`--conda-prefix`), but shared-object loading from a networked
  filesystem is still subject to cold-cache reads, visible as the 0–56 MB `io_in` spread that
  doesn't correspond to any real input/output of the script.
- **Slurm scheduling latency**: for jobs whose own compute (`cpu_time`) is ≤0.68 s, wall time
  is dominated by everything Snakemake/Slurm does around the job, not inside it.
- **No redundant passes, no serialization of large data**: the only file this rule is declared to
  depend on (the multi-hundred-MB ARC file) is never opened by the script — it is referenced only
  by its path string. This is worth flagging precisely because it means the DAG-level dependency
  (correctness) is completely decoupled from this rule's runtime cost.
- **Threading**: no `threads:` declared in the rule (`add_external_features.smk:37-52`); none
  needed — there is nothing to parallelize.

## 7. Three functions

Since the code is O(1) in every manifest variable and the benchmarks actively contradict a size
dependence (§5), these are **constants**, not functions of `n`, with a spread that reflects
interpreter/library startup plus Slurm/filesystem noise rather than any exponent uncertainty:

```
t_est(n)  = 3.0   seconds   (for all n; ≈ pooled mean of 4.65 s, rounded down toward the
                             cpu_time-dominated cases which are more representative of a
                             warm-cache node)
t_low(n)  = 0.8   seconds   (for all n; near the observed floor of 0.88–0.91 s, i.e. R +
                             conda-env startup with a warm filesystem cache)
t_high(n) = 15.0  seconds   (for all n; above the observed ceiling of 10.78 s, to leave margin
                             for a colder cache or heavier Slurm contention than was sampled —
                             this is single-measurement data per cluster, not replicated)
```

These bound the **observed operating envelope**, not a complexity class — there is no class
uncertainty here because the class is O(1) with high confidence from the code (§4); the
uncertainty is purely in how much interpreter-startup/scheduling noise a given job happens to hit.

## 8. Evaluated at mean / min / max cluster

Because `t_est`, `t_low`, `t_high` do not depend on `n`, they evaluate identically at every
point:

| Point | cluster | t_low | t_est | t_high |
|---|---|---:|---:|---:|
| mean (across 6 clusters) | — | 0.8 s | 3.0 s | 15.0 s |
| min | `thp1_1` | 0.8 s | 3.0 s | 15.0 s |
| max | `telohaec_crispri` | 0.8 s | 3.0 s | 15.0 s |

(For reference, the actually-measured values at these two named clusters were: `thp1_1` 0.89 s
multiome / 6.86 s scATAC; `telohaec_crispri` 1.79 s multiome / 10.76 s scATAC — both inside the
`[t_low, t_high]` envelope above, and neither at an extreme consistent with "min cluster is
cheapest, max cluster is priciest.")

## 9. Does the class differ between modalities?

**The class is O(1) in both; only a negligible constant differs, and it is not the constant a
size-driven model would predict.** The code does take a genuinely different branch by modality
(§2): multiome always hits the `to_generate=="ARC"` path and writes 5 rows referencing the
(unread) ARC file; scATAC always hits `"Neither"` and writes 0 rows. But the actual work
difference between "append a 5-row literal data.frame" and "do nothing" is a handful of R
vector operations — unmeasurable against interpreter-startup noise. Empirically the two
modalities' wall-time distributions overlap almost entirely (multiome 0.88–8.69 s; scATAC
0.91–10.78 s) and, if anything, scATAC (which does *less* work per job) runs slightly slower on
average — the opposite of what a work-proportional model predicts. This is the same rule the
briefing's four-rules-differ-materially table intends to flag by "`to_generate` differs by
modality," but for *this specific script* the flagged mechanism (ARC vs Neither) does not
translate into a measurable cost difference; the ~15× multiome/scATAC gap actually observed in
the benchmarks belongs to the differently-named `rule add_external_features`
(`ENCODE_rE2G/workflow/rules/genomewide_features.smk:179`), not to `make_external_features_config`.

## 10. What the measurement cannot tell us

- **Cannot separate "R+library startup" from "conda activation" from "Slurm dispatch latency"**
  with this instrumentation — the benchmark only has wall time, `cpu_time`, and RSS/IO, all
  measured *inside* the Snakemake job wrapper. All three plausibly contribute to the 0.88–10.78 s
  spread; the data cannot apportion it.
- **Cannot explain the specific 0–56 MB `io_in` variation** beyond "filesystem cache state" — this
  is a hypothesis, not a measured mechanism, since no file this rule reads or writes is anywhere
  near that size.
- **Six points per modality, one measurement each (n=1 per cluster)** — there is no replication to
  distinguish "this cluster happened to hit a cold cache" from "this cluster is systematically
  slower." The apparent modality difference in §9 could partly be an artifact of when in the
  overall pipeline run scATAC vs multiome jobs happened to execute (warm vs cold shared conda env
  on the node), not a property of the rule itself.
- **Cannot rule out cost at a much larger scale that isn't O(1)** — the code reads a fixed 6-row
  config table regardless of how many clusters exist project-wide, so this conclusion is scoped to
  "per-job, at the current cell_clusters table size (6 rows)." If the number of clusters or the
  columns checked ever grew by orders of magnitude, the `dplyr::filter` and column-name check at
  `:12-13` would still be O(rows in that table), which was not exercised here.
- **The dead Kendall branch (§2) is a code-reading observation, not something the benchmarks can
  confirm or refute** — no job in either modality took that branch, so there is no data on its
  cost even if it were reachable.

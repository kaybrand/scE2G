# Phase 0 measurements — IGVF10 Multiome

Committed backup of the Phase 0 measurement data, so it survives the `$SCRATCH` purge
(90 days from last content write). **Phase 1 should read the live data under `$SCRATCH`,
not these copies** — these exist as the durable record.

| File | Rows | What it is |
|---|---|---|
| `igvf10_multiome_benchmarks.tsv` | 165 | One row per JOB: wall time, peak RSS, I/O, CPU time |
| `igvf10_multiome_size_manifest.tsv` | 6 | One row per cluster: the candidate size variables |
| `igvf10_multiome_checksums.tsv` | 66 | Row counts + md5 of *decompressed* bytes, per output |

Live source (until purged):
`/scratch/users/kaybrand/scE2G_optimize_results/igvf10_multiome/` — per-job files at
`benchmarks/<rule>/<wildcards joined by ~>.tsv`, plus `qc_plots/predictions_qc_report.html`
and the post-run job-level DAG in `../analysis/dag_multiome.{dot,svg}`.

**Units:** `s` = seconds; `max_rss`/`max_vms`/`max_uss`/`max_pss` = MB; `io_in`/`io_out` = MB;
`cpu_time` = seconds. Wall time includes Slurm scheduling latency.

## Run provenance

- Config `configs/igvf10_multiome_config.yaml`, model `multiome_powerlaw_v3`, 6 clusters
- Snakemake 9.16.3, `--profile profiles/measure`, non-preemptible partitions, `retries: 0`
- **30 distinct rules, 165 jobs, 0 rule errors.** `to_generate == ARC` on all six clusters,
  confirming the per-cluster modality branch resolved as intended
- Correctness gate: the toy chr22 run reproduced the committed expected outputs exactly, so
  the `benchmark:` instrumentation does not perturb results
- scATAC half is **not yet run**; this is the Multiome record only

## Headline results

Total wall time per cluster (sum over that cluster's jobs, not critical path):

| Cluster | wall min | n_frag | n_cells |
|---|---:|---:|---:|
| telohaec_crispri | 397.8 | 394,847,593 | 15,555 |
| jurkat | 147.0 | 75,342,313 | 8,239 |
| k562_crispri | 143.3 | 86,424,398 | 5,109 |
| jurkat_pma_cd3_4hr | 106.8 | 55,803,917 | 5,124 |
| thp1_1 | 83.8 | 26,231,026 | 2,676 |
| thp1_2 | 82.0 | 26,292,230 | 2,593 |

Slowest single jobs:

| Rule | Cluster | wall min | peak GB |
|---|---|---:|---:|
| `create_neighborhoods` | telohaec_crispri | 112.1 | 0.28 |
| `frag_to_norm_bigWig` | telohaec_crispri | 81.4 | 24.53 |
| `call_macs_peaks` | telohaec_crispri | 48.4 | 13.29 |
| `frag_to_norm_bigWig` | k562_crispri | 31.6 | 7.05 |
| `frag_to_tagAlign` | telohaec_crispri | 25.9 | 4.65 |

Summed across all clusters: `frag_to_norm_bigWig` 174.9 min, `create_neighborhoods` 171.6,
`call_macs_peaks` 91.1, `add_external_features` 91.0, `create_predictions` 58.2.

## Three findings that should shape Phase 1 and Phase 2

### 1. `compute_kendall` is no longer the bottleneck

It does not appear in the top twelve jobs: 14.5 min on the largest cluster, 3.7 min on the
smallest. This reflects the Kendall optimization work already merged. The prior expectation
that it dominates — encoded in its 63 GB memory floor and 12-hour runtime allowance, both now
oversized — is out of date. `create_neighborhoods` (ABC) now holds the title.

### 2. Only TWO size variables actually vary across these six clusters

| Variable | Range | Ratio |
|---|---|---:|
| `n_frag` | 26.2 M – 394.8 M | **15.1×** |
| `n_cells` | 2,593 – 15,555 | **6.0×** |
| `n_umi` | 14.5 M – 174.2 M | 12.0× |
| `j_elem` (candidate elements) | 156,421 – 158,549 | 1.014× |
| `k_genes` | 20,532 – 20,532 | **1.000×** |
| `n_links` (E2G candidates) | 10.2 M – 11.4 M | 1.13× |
| `n_pairs` (Kendall pairs) | 9.7 M – 11.3 M | 1.17× |

Everything downstream of candidate-region calling has **essentially constant problem size**.
This is a consequence of the pipeline design — the candidate element set is capped and the gene
annotation is shared — not of cluster choice.

**Consequence: the six benchmark points give no scaling leverage for any rule whose cost is
driven by `n_links`, `j_elem` or `k_genes`.** Those rules span a 1.0–1.2× range in problem size,
which cannot distinguish O(n) from O(n²). For them the complexity class must be derived from
the code alone, and Phase 1 agents must say so explicitly rather than fitting an exponent to
noise. Only fragment-driven rules (`process_fragment_file`, `frag_to_tagAlign`,
`frag_to_norm_bigWig`, `call_macs_peaks`, `create_neighborhoods`) and cell-driven rules
(`generate_atac_matrix`, `compute_kendall`) have a usable dynamic range here.

### 3. `create_neighborhoods` is CPU-bound and single-threaded

From `benchmarks/create_neighborhoods/telohaec_crispri.tsv`: `s=6725.41`, `cpu_time=6572.08`,
`mean_load=96.31`, `io_in=1069.96` MB, `max_rss=290` MB. That is **97.7% CPU utilisation on one
core for nearly two hours**, with modest I/O and small memory — not the I/O-bound profile the
low memory first suggests. Every ABC rule declares no `threads:`, so it runs at 1 core by
default.

Contrast `frag_to_norm_bigWig`, which is the other heavy rule but a different shape entirely:
24.5 GB peak RSS, driven by `bedtools genomecov` plus an in-memory `sort`. Note it exists only
because `make_IGV_tracks: True`, and **nothing downstream consumes its output** — so despite
being the largest total consumer it may carry large slack and sit off the critical path. Phase 2
decides that; total time is not the same as critical-path time.

## Caveats

- Wall time includes Slurm queue/scheduling latency, so short jobs are overstated relative to
  their compute cost. Use `cpu_time` to separate the two.
- One rule is uninstrumented by design: `abc_generate_chrom_sizes_bed_file` (~0.03 s) crashes
  under Snakemake's psutil benchmark monitor. See `profiles/measure/config.yaml`.
- Training-only rules are not instrumented; they never load on the prediction path.
- These are single measurements, not replicates. The scATAC run will re-measure the shared
  ABC/fragment stack and provide a second observation of those rules.

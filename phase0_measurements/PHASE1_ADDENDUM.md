# Phase 1 launcher addendum

Read this *after* `PHASE1_BRIEFING.md`. It records facts verified at launch time (2026-08-21)
that the briefing does not carry, or that correct it.

## A. Rule names on disk vs. names in the plan

The plan document prefixes some ABC rules with `abc_`. The **benchmark tables and the `.smk`
files do not.** Use the unprefixed names:

| plan says | actual rule name |
|---|---|
| `abc_call_macs_peaks` | `call_macs_peaks` |
| `abc_sort_narrowpeaks` | `sort_narrowpeaks` |
| `abc_generate_chrom_sizes_bed_file` | `generate_chrom_sizes_bed_file` (no benchmark — by design) |

`process_fragment_file`, `features_required` and `basic_features_required` are **not** plain
`rule` statements: `process_fragment_file` is an indented conditional rule
(`workflow/rules/frag_to_tagAlign.smk:79`), and the other two are `checkpoint` statements
(`workflow/rules/add_external_features.smk:1`, `ENCODE_rE2G/workflow/rules/genomewide_features.smk:3`).

## B. Three columns in `size_manifest_*.tsv` are NOT sizes — do not use them

`n_frag_bytes` (9–10), `n_umi_bytes` (9–10) and `n_cells_bytes` (5–6) are the byte sizes of the
*count text files* (`fragment_count.txt`, `umi_count.txt`, `cell_count.txt`) — i.e. the number of
digits in the number. They are measurement artifacts, not data volumes. If you need fragment-file
volume, take it from the cluster table / the fragment file itself, not from `n_frag_bytes`.

The genuinely useful `*_bytes` columns are the ones the briefing §3 names:
`arc_bytes`, `actonly_bytes`, `feat_bytes`, `kendall_scored_bytes`, `enh_list_bytes`,
plus `n_cand_pairs_bytes`, `n_kendall_pairs_bytes`, `n_peaks_bytes`, `j_cand_elems_bytes`,
`k_genes_bytes`, `n_called_bytes`.

## C. `n_called_all` / `n_called_distal` differ by modality — the briefing quotes multiome only

Briefing §3 Tier 3 gives `n_called_all` 59,438–61,284 and `n_called_distal` 46,155–47,560.
**Those are the multiome numbers.** scATAC is roughly 2× larger:

| | multiome | scATAC |
|---|---|---|
| `n_called_all` | 59,438 – 61,284 | 105,880 – 120,113 |
| `n_called_distal` | 46,155 – 47,560 | 89,681 – 104,046 |
| `n_feat_cols` | 24 | **18** |

So the scATAC model calls ~2× as many links from the same candidate set, and its feature table
is 18 columns wide rather than 24. If your rule reads the thresholded predictions or the final
feature table, **the modality changes your `n` by ~2×, not just the constant.** Read the size
manifest for the modality you are characterising.

The columns that are genuinely **empty in `size_manifest_scatac.tsv`** are exactly the ones for
quantities that do not exist on the scATAC path:

```
n_umi   n_umi_bytes   n_cells   n_cells_bytes
n_kendall_pairs   n_kendall_pairs_bytes
arc_bytes   kendall_scored_bytes
```

**Correction (2026-08-21):** an earlier revision of this section wrongly listed `actonly_bytes`
among the empty columns. It is **populated** for all six scATAC clusters (189,828,751 –
223,454,599 bytes) and verified against the on-disk `ActivityOnly_features.tsv.gz` sizes. The
empty neighbours are `arc_bytes` and `kendall_scored_bytes`. `activity_only_features` runs
upstream of the ARC/Kendall branch point, so it exists and is measured in both modalities.

## D. How to pull your rule's benchmark rows

```bash
cd /scratch/users/kaybrand/scE2G_optimize_results/analysis
head -1 benchmarks_multiome.tsv
awk -F'\t' -v r=YOUR_RULE 'NR==1 || $2==r' benchmarks_multiome.tsv
awk -F'\t' -v r=YOUR_RULE 'NR==1 || $2==r' benchmarks_scatac.tsv
```

Columns: `modality  rule  cluster  wildcards  s  h:m:s  max_rss  max_vms  max_uss  max_pss
io_in  io_out  mean_load  cpu_time`. Memory and I/O are in MB; `s` and `cpu_time` in seconds.

Per-job source of truth, if you want the raw file:
`/scratch/users/kaybrand/scE2G_optimize_results/igvf10_{multiome,scatac}/benchmarks/<rule>/<wildcards ~ joined>.tsv`

Retained intermediates for shape inspection (`--notemp` was set):
`/scratch/users/kaybrand/scE2G_optimize_results/igvf10_{multiome,scatac}/<cluster>/`

## E. Job counts confirmed at launch

30 rules × mostly 6 clusters = 165 multiome jobs; 26 rules = 141 scATAC jobs. `plot_stats`,
`hover_plots` and `save_reference_configs` have **1 job each** (aggregating rules, unwildcarded),
so you get n=1 per modality for those — not 6.

## F. Both DAG exports were regenerated at launch — the old multiome one was unusable

Phase 2 input, corrected 2026-08-21. This matters only to Phase 2, not to Phase 1 agents.

`dag_multiome.dot` as inherited had **125 nodes** — exactly the *pre-checkpoint-resolution* DAG the
plan warned about, not the post-run one. It was missing all 42 jobs that sit behind the two
checkpoints: the four multiome-only rules (`make_kendall_pairs`, `generate_atac_matrix`,
`compute_kendall`, `arc_e2g` = 24 jobs) and the three `generate_num_*` rules (18 jobs). Snakemake
cannot see past an unevaluated checkpoint, so a DAG exported before resolution is structurally
incomplete. The old file is kept as `dag_multiome.PRE_CHECKPOINT.dot.bak`.

Both were re-exported after the runs completed, with checkpoints resolved:

| | nodes | edges | reconciliation against measured jobs |
|---|---:|---:|---|
| `dag_multiome.dot` | 155 | 303 | 155 − `all` − `generate_chrom_sizes_bed_file` = 153 benchmarkable; + 12 resolved checkpoint jobs = **165** ✓ |
| `dag_scatac.dot` | 131 | 243 | 131 − `all` − `generate_chrom_sizes_bed_file` = 129; + 12 = **141** ✓ |

Both now reconcile exactly. Rendered to `.svg` alongside. Two parsing notes for Phase 2:

- **The DAG uses `abc_`-prefixed rule names** (`abc_call_macs_peaks`, `abc_create_neighborhoods`,
  `abc_make_candidate_regions`, `abc_sort_narrowpeaks`, `abc_create_predictions`,
  `abc_generate_chrom_sizes_bed_file`) because the ABC workflow is included as a Snakemake module.
  The **benchmark tables and the `.smk` files use the unprefixed names** (§A). Phase 2 must map
  between them.
- **Node labels carry wildcards only where the wildcard is first introduced** — Snakemake's dot
  export omits wildcards inherited from a job's dependencies. So `run_e2g_qnorm` shows
  `cluster:`/`model_name:` but `gen_final_features` and `abc_create_predictions` show none. Cluster
  identity must be **propagated along edges** from the labelled leaves; label text alone does not
  identify a job.
- The two `checkpoint` node types (`features_required`, `basic_features_required`) do **not** appear
  in either resolved DAG, though both are benchmarked rules with 6 jobs each per modality. Phase 2
  must insert them manually as the barrier nodes they are — see unit #26.

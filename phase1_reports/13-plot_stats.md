# Unit 13 — `plot_stats` (Tier B)

## 1. Unit identity

- Script: `workflow/scripts/prediction_qc/plot_all_qc_stats.R` (218 lines)
- Rule: `plot_stats`, `workflow/rules/sc_predictions.smk:185-205`
- Modality: **both** (multiome and scATAC; same script, same rule, no modality branch in the
  Snakemake rule — branching happens *inside* the R script based on data content)
- Job count: **1 job per modality** (unwildcarded aggregating rule, confirmed by
  `PHASE1_ADDENDUM.md §E` and directly from the benchmark files below — this is a fan-in, not a
  per-cluster rule)

```
$ find .../igvf10_multiome/benchmarks/plot_stats -type f
.../igvf10_multiome/benchmarks/plot_stats/all.tsv
$ find .../igvf10_scatac/benchmarks/plot_stats -type f
.../igvf10_scatac/benchmarks/plot_stats/all.tsv
```

Single file named `all.tsv` in each modality — one job, `wildcards` column empty in
`benchmarks_{multiome,scatac}.tsv` (`cluster` and `wildcards` are both blank for this rule; only
`modality`, `rule="plot_stats"`, and the timing/resource columns are populated).

## 2. Upstream dependency — this rule is a fan-in barrier

`input.stats_files` (`sc_predictions.smk:186-187`):

```python
biosample_model_threshold = list(zip(BIOSAMPLE_DF["biosample"], BIOSAMPLE_DF["model_dir_base"], BIOSAMPLE_DF["model_threshold"]))
rule plot_stats:
    input:
        stats_files = [os.path.join(RESULTS_DIR, biosample, model_name, f"scE2G_predictions_threshold{threshold}_stats.tsv")
                        for biosample, model_name, threshold in biosample_model_threshold]
```

Each entry is produced by `rule get_stats_per_model_per_cluster` (`sc_predictions.smk:164-182`,
one job per `(cluster, model_name, threshold)`). In this config there is one model per cluster,
so `len(biosample_model_threshold) == 6` and `plot_stats` literally cannot start until **all six**
`get_stats_per_model_per_cluster` jobs have written their output — confirmed on disk (six
`*_stats.tsv` files per modality, listed under §5 below). `get_stats_per_model_per_cluster` in
turn reads `scE2G_predictions.tsv.gz` (full) and the thresholded predictions table per cluster —
i.e., the tail end of the whole per-cluster ABC/e2g pipeline for that cluster.

**Phase 2 implication:** `plot_stats` is a synchronization barrier over all 6 clusters. Its own
cost is small (see §7), but because it waits on the slowest of the six upstream
`get_stats_per_model_per_cluster` jobs (each of which waits on that cluster's entire prediction
pipeline), a cheap fan-in node like this can still sit on the critical path — its start time is
`max` over 6 upstream completion times, not a sum, and Phase 2's graph needs that "wait for all
6" edge represented explicitly, not just "cost ≈ 12s."

## 3. Input volume — confirmed small

Per the addendum note for this unit: input is the per-cluster **stats files**, not the
prediction tables. Confirmed on disk — every input file is a few hundred bytes:

```
igvf10_multiome/jurkat/multiome_powerlaw_v3/…_stats.tsv                    367 B
igvf10_multiome/jurkat_pma_cd3_4hr/multiome_powerlaw_v3/…_stats.tsv        379 B
igvf10_multiome/k562_crispri/multiome_powerlaw_v3/…_stats.tsv              373 B
igvf10_multiome/telohaec_crispri/multiome_powerlaw_v3/…_stats.tsv          379 B
igvf10_multiome/thp1_1/multiome_powerlaw_v3/…_stats.tsv                    367 B
igvf10_multiome/thp1_2/multiome_powerlaw_v3/…_stats.tsv                    367 B
                                                                    total ≈ 2.2 KB

igvf10_scatac/jurkat/scATAC_powerlaw_v3/…_stats.tsv                        352 B
igvf10_scatac/jurkat_pma_cd3_4hr/scATAC_powerlaw_v3/…_stats.tsv            364 B
igvf10_scatac/k562_crispri/scATAC_powerlaw_v3/…_stats.tsv                  360 B
igvf10_scatac/telohaec_crispri/scATAC_powerlaw_v3/…_stats.tsv              364 B
igvf10_scatac/thp1_1/scATAC_powerlaw_v3/…_stats.tsv                        352 B
igvf10_scatac/thp1_2/scATAC_powerlaw_v3/…_stats.tsv                        352 B
                                                                    total ≈ 2.1 KB
```

This matches the measured `io_in` of 6.58 MB (multiome) / 13.51 MB (scatac) only loosely — the
bulk of `io_in` is R/conda environment loading (interpreter, shared libraries, package files for
`plyr`, `dplyr`, `tidyr`, `data.table`, `stringr`, `ggplot2`, `ggpubr`), not the ~2 KB of actual
stats data. `io_out` is 0.15 MB / 0.37 MB, matching the small PDF+TSV outputs (see §7).

## 4. What the code actually does

`plot_all_qc_stats.R:161-217` (main block):

1. **Read + merge** (`:161-163`): `lapply(stats_files, fread) %>% rbindlist()` — reads all 6
   tiny per-cluster stats files and stacks them into one data frame of **exactly `nrow = C`
   rows** (`C` = number of clusters/jobs feeding in, `C=6` measured), then `fwrite`s the merged
   table to `all_stats`.
2. **Modality branch on data content, not a Snakemake wildcard** (`:193-204`): if
   `sum(stats$umi_count) > 0` (true for multiome, false for scATAC since scATAC has no RNA and
   `umi_count` is always 0 — consistent with `PHASE1_ADDENDUM.md §C`, "`n_umi` … empty in
   `size_manifest_scatac.tsv`"), the script takes a richer branch: `x_vars` grows from
   `c("fragments_total")` to `c("fragments_total","umi_count","cell_count")`, `gene_vars` grows
   from 2 to 3 entries, and `plot_dataset_stats()` (2 scatter panels) runs instead of
   `plot_fragment_distribution()` (1 violin panel).
3. **Four calls to `plot_scatter_set()`** (`:210-213`, function body `:97-126`): nested loop
   `for y in y_vars: for x in x_vars: ggplot(...)`, `length(y_vars) * length(x_vars)` subplots
   per call, combined with `ggarrange` and written to one PDF per call via `ggsave`.
4. **One call to `plot_violin_metrics()`** (`:217`, function body `:128-156`): one `ggplot` per
   metric in `violin_label_key` (`label_key[6:length(label_key)]`, always **9** metrics — this
   slice is on the static `label_key` vector defined at `:170-175`, unaffected by the modality
   branch), combined into one grid PDF.

None of these operate on anything bigger than `C` rows. There is no join against the prediction
tables, no per-row nested loop, no sort of a large table — `plot_stats` only ever sees the tiny
merged stats data frame.

### Subplot count is a fixed code constant per modality, not a function of C

| Call | y_vars | x_vars (multiome) | x_vars (scATAC) | subplots (multiome) | subplots (scATAC) |
|---|---|---:|---:|---:|---:|
| `plot_dataset_stats` / `plot_fragment_distribution` | — | — | — | 2 | 1 |
| `plot_scatter_set(enh_vars)` | 2 | 3 | 1 | 6 | 2 |
| `plot_scatter_set(eg_vars)` | 2 | 3 | 1 | 6 | 2 |
| `plot_scatter_set(gene_vars)` | 3 (multiome adds `n_genes_not_expressed`) / 2 | 3 | 1 | 9 | 2 |
| `plot_scatter_set(dist_vars)` | 2 | 3 | 1 | 6 | 2 |
| `plot_violin_metrics` | 9 (fixed) | — | — | 9 | 9 |
| **Total ggplot objects** | | | | **38** | **18** |

Total output artifacts written: 6 `ggsave` PDF calls + 1 `fwrite` TSV = **7 files**, fixed by
code, independent of `C`.

## 5. Size variable

Named per `PHASE1_BRIEFING.md`'s menu: **none of the Tier 1/2/3 per-cluster variables
(`n_frag`, `n_umi`, `n_cells`, `n_peaks`, …) apply here.** This rule's cost is not a function of
any single cluster's size — it is a function of **how many clusters feed into it**, which the
size manifest does not carry as a column (it is the *row count* of `size_manifest_*.tsv`
itself). Per the unit brief, I call this `C` — the number of upstream
`get_stats_per_model_per_cluster` jobs (= number of clusters, since one model per cluster in this
config). **Measured: `C = 6` for both modalities.** There is only one value of `C` anywhere in
the Phase 0 data — every run of this pipeline used the same 6-cluster config, so `C` never
varies across the dataset, not even across the two modality measurements.

## 6. Asymptotic class, derived from code

No operation in `plot_all_qc_stats.R` is worse than linear in `C`:

- `rbindlist(lapply(stats_files, fread))` (`:161-162`) — O(`C`) read + O(`C`) concatenation.
- `dplyr::select(...) %>% distinct()` (`:191`) — hash-based dedup, O(`C`).
- Every `ggplot(stats, aes(...))` call draws `geom_point`/`geom_jitter`/`geom_violin` over a
  data frame of `C` rows — O(`C`) point layout per subplot. `geom_violin`'s density estimate
  (`stat = "area"`/`"count"`) is O(`C log C`) at worst (KDE bandwidth selection sorts the input),
  negligible at `C=6`.
- The number of subplots and output files (§4 table) is a **fixed code constant**, not a
  function of `C` — it depends on the static 14-entry `label_key` and the umi/no-umi branch, both
  independent of how many clusters are present.

**Class: O(`C`)** (linear in the number of clusters), dominated by a large additive constant
(R interpreter start, `conda` env activation, loading `plyr/dplyr/tidyr/data.table/stringr/
ggplot2/ggpubr`, and opening 6 PDF/cairo graphics devices) that does **not** scale with `C` at
all. I found no O(`C²`) or worse structure anywhere in this script — there is no all-pairs join,
no per-cluster nested loop over another per-cluster set, and no cross product between clusters.

## 7. Benchmark cross-check

```
modality  rule       cluster  wildcards  s      h:m:s     max_rss  max_vms  max_uss  max_pss  io_in  io_out  mean_load  cpu_time
multiome  plot_stats           all       11.69  0:00:11   174.00   695.20   171.36   171.47   6.58   0.15    39.71      4.81
scatac    plot_stats           all       11.82  0:00:11   170.54   690.92   167.87   168.06   13.51  0.37    34.10      4.17
```

- **Wall time**: 11.69 s (multiome), 11.82 s (scatac) — essentially identical.
- **`cpu_time`**: 4.81 s (multiome), 4.17 s (scatac).
- **`wall − cpu_time` ≈ 6.9 s / 7.65 s** in both — per the briefing's rule 5 ("wall time includes
  Slurm scheduling latency… use `cpu_time` to separate compute from waiting"), this gap is
  scheduling/conda-activation/filesystem latency, not algorithmic cost.
- **Peak memory**: 174 MB / 171 MB `max_rss` — trivial, no memory pressure, consistent with a
  data frame of 6 rows plus R/ggplot2 machinery in memory.

**This cannot confirm or refute an O(`C`) class, or any class, in `C`.** There is exactly one
value of `C` (6) in the entire Phase 0 dataset — both modality measurements used the same
6-cluster config. This is a stronger limitation than "only 6 points to fit an exponent" (the
situation for Tier 1/2/3 per-cluster variables elsewhere in this project): here the independent
variable itself **never varies** across any measurement we have. Any claim about the
coefficient of `C` is 100% code-derived reasoning, not data-supported, even directionally.

The multiome vs. scATAC comparison at fixed `C=6` is informative about a *different* axis: the
0.64 s `cpu_time` gap (4.81 vs 4.17 s) plausibly reflects the ~2× larger subplot count in
multiome (38 vs 18 ggplot objects, §4) — i.e., evidence that per-subplot rendering cost is real
and non-zero, but it says nothing about how cost changes with `C`, since `C` is 6 in both.
`io_in` differs 6.58 → 13.51 MB (≈2×), consistent with `n_called_all`/`n_called_distal` being
~2× larger in scATAC (`PHASE1_ADDENDUM.md §C`) if any residual input-table caching/paging
differs by modality — but note the actual `stats_files` inputs are ~2 KB either way (§3), so
this `io_in` difference is almost certainly dominated by conda/library-loading I/O variance
between the two job's filesystem cache states, not this rule's own data.

## 8. Concrete overhead sources

- **Fixed R/conda startup dominates.** `cpu_time` is only 4.2–4.8 s total for loading 7 R
  packages (`plyr`, `dplyr`, `tidyr`, `data.table`, `stringr`, `ggplot2`, and — only in this
  script among the QC rules — `ggpubr` for `ggarrange`, `plot_all_qc_stats.R:1-9`) plus all
  rendering. Package loading alone (ggplot2 + ggpubr pull in `grid`, `gtable`, `scales`, etc.)
  plausibly accounts for a large fraction of that.
- **7 output files, each with its own device-open/flush cost** (`ggsave` ×6 for PDFs via a
  cairo/pdf graphics device, `fwrite` ×1 for the TSV) — a fixed per-file constant unrelated to
  `C`.
- **Declared `threads:` — none declared** (no `threads:` key in the rule, `sc_predictions.smk:
  185-205`), so it runs at Slurm's default of 1 core. This is consistent with the project-wide
  finding that ABC/QC rules run single-threaded; there is nothing here that could use more than
  one core anyway (R's single-threaded `ggplot2`/`grid` rendering, no parallel-capable step).
- **`mem_mb=encode_e2g.ABC.determine_mem_mb`** is requested but peak `max_rss` is only ~171-174
  MB — heavily over-provisioned relative to actual use, though at this rule's short duration
  that is not itself a cost driver.
- **Wall-vs-cpu gap (~7 s)** is scheduling/filesystem latency per addendum rule 5, not part of
  the algorithmic cost function below.

## 9. Three functions: `t_est(C)`, `t_low(C)`, `t_high(C)`

Because `C` never varies in the data (§7), I cannot fit an intercept and a slope from
measurements — I can only bound the slope from code reasoning and force the fixed term to
reconcile with the single measured point at `C=6`. These are **wall-clock seconds**, per
modality (the code path taken differs by modality per §4, so the subplot-count constant differs;
I hold that difference in the constant term, not in the class, since the class itself — O(`C`)
— does not differ):

**Multiome** (measured: 11.69 s at C=6):

- `t_est(C) = 9.7 + 0.33 · C`  — assumes ~83% of the measured wall time at C=6 is the
  C-independent fixed cost (R interpreter + package loads + Slurm/conda scheduling latency +
  7 file-open/flush operations), and the remaining ~17% is 38 subplots each doing O(C) point
  rendering, apportioned linearly across clusters.
- `t_low(C) = 11.64 + 0.01 · C`  — lower bound on the class: nearly all cost is the fixed
  constant; each added cluster costs only marginally more (one more `fread` of a ~370 B file
  and one more point added to already-open ggplot layers).
- `t_high(C) = 5.7 + 1.0 · C`  — upper bound on the class: a more pessimistic per-cluster
  slope (e.g., if larger `C` also forced larger jitter/violin density computations, wider PDF
  pages, or bigger legends), still linear, not superlinear, because §6 found no O(`C²`)
  structure in the code to justify a steeper class.

**scATAC** (measured: 11.82 s at C=6):

- `t_est(C) = 10.0 + 0.30 · C`
- `t_low(C) = 11.76 + 0.01 · C`
- `t_high(C) = 6.8 + 0.84 · C`

(scATAC's fixed term and slope are marginally smaller than multiome's, reflecting the ~18 vs 38
subplot-object count from §4 — fewer x_vars/y_vars combinations to render — but this is a
constant-term difference, not a class difference; see §10.)

All three functions are linear in `C`; **only the split between the intercept and the slope is
unconstrained by data** — that split is a modeling choice made here to satisfy the deliverable
format, not a measured result. Phase 2 should treat the intercept (≈6–12 s) as the reliable part
(it reproduces the single measured point) and treat the slope as a wide, code-bounded guess.

## 10. Evaluated values

Per-unit instruction: express in terms of `C` rather than a per-cluster size variable (there is
no per-cluster size variable driving this rule — see §5). The general instruction to evaluate at
min/max *cluster* size (`thp1_1`, `telohaec_crispri`) does not apply to this rule, since its cost
does not depend on which cluster's data flows through it, only on how many clusters there are.
Evaluating instead at the one measured point and at illustrative (unmeasured, pure-extrapolation)
values of `C`:

| `C` | multiome `t_est` | multiome `t_low` | multiome `t_high` | scATAC `t_est` | scATAC `t_low` | scATAC `t_high` |
|---:|---:|---:|---:|---:|---:|---:|
| 3 (illustrative, unmeasured) | 10.7 s | 11.67 s | 8.7 s | 10.9 s | 11.79 s | 9.3 s |
| **6 (measured)** | **11.7 s** | **11.7 s** | **11.7 s** | **11.8 s** | **11.8 s** | **11.8 s** |
| 12 (illustrative, unmeasured) | 13.6 s | 11.76 s | 17.7 s | 13.6 s | 11.88 s | 16.9 s |

At `C=6` all three collapse to the measured value by construction (that is the only anchor
point we have). The spread at `C=3`/`C=12` is purely illustrative of how wide the uncertainty
in the slope is — it is not a prediction Phase 2 should trust quantitatively, only qualitatively
("this rule's cost grows slowly and linearly if cluster count changes, and is dominated by a
~6-12 s fixed floor regardless").

## 11. Both-modality comparison — class vs. constant

**Only the constant differs, not the class.** Both modalities execute the identical `for`-loop
structure over `y_vars × x_vars`; the only modality-dependent effect (§4) is *which* branch of an
`if (sum(stats$umi_count) > 0)` check is taken, which changes the fixed *count* of subplots (38
vs 18) and hence the constant term, but does not introduce or remove any dependency on `C` itself
— the branch is selected by whether `umi_count` is present at all (a data-driven fact of the
modality: scATAC has no RNA channel), not by anything that scales with `C`. Both modalities
remain O(`C`).

## 12. What the measurement cannot tell us

- **Whether the class is really linear in `C`.** Every measurement in Phase 0 used exactly the
  same `C=6` config; there is no second data point at any other `C` to check against. The O(`C`)
  claim in §6 is entirely code-derived and has zero empirical confirmation, not "weak"
  confirmation — this is a stronger limitation than the Tier 3 per-cluster variables discussed in
  the briefing, which at least span a small ratio.
- **The true split between the fixed constant and the per-cluster slope.** §9's numbers are
  reverse-engineered to match the single measured point under an assumed split; a different
  split would fit the same data point equally well. Phase 2 should not treat the specific slope
  values as measured quantities.
- **Whether the wall-vs-cpu gap (~7 s) is stable.** It could reflect Slurm queue contention at
  the time this job happened to run rather than a structural property of this rule; with n=1 per
  modality there is no way to separate "typical scheduling latency for this rule" from "this
  particular run's scheduling latency."
- **Behavior at large `C`.** If a future config had, say, 50 or 100 clusters, `ggarrange`'s
  grid layout, PDF page sizing (`wd`/`ht` in `plot_scatter_set`/`plot_violin_metrics` do scale
  with fixed `nrow`/`ncol` counts, not `C`, so the *page* size is stable, but very dense
  `geom_jitter` overplotting could change rendering time in ways not observable in this data),
  or R's `fread`/`rbindlist` overhead across many tiny files (subprocess/file-open overhead per
  file, not tested here) could all shift the real slope away from anything in §9.
- **Whether one-model-per-cluster holds in general.** This report assumes
  `C == number of clusters` because the observed config has exactly one model per cluster
  (verified from directory listings, §3). If a future config runs multiple models per cluster,
  `C` (= `len(biosample_model_threshold)`, the true driver per `sc_predictions.smk:184-187`)
  would exceed the cluster count, and the formulas above should be read as functions of that job
  count, not literally "number of clusters."

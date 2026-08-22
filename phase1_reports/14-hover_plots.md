# Unit 14 (Tier B): `hover_plots`

## 1. Unit identity

- **Scripts**: `workflow/scripts/prediction_qc/qc_report.Rmd` (907 lines), `workflow/scripts/prediction_qc/qc_plot_tab_template.Rmd` (170 lines; not run directly — expanded and knit as a child document, once per model, from inside `qc_report.Rmd`).
- **Rule**: `hover_plots`, `workflow/rules/sc_predictions.smk:208-225`.
- **Modality**: BOTH (multiome and scATAC), same code path, same conda env (`workflow/envs/sc_e2g.yml`).
- **Job count**: **exactly 1 job per modality** (aggregating, unwildcarded — confirmed in `benchmarks_{multiome,scatac}.tsv`, `cluster` column empty, `wildcards` = `all`). **n=1. No fit is possible; every number below is a two-point (one per modality) comparison, not a scaling curve.**

## 2. What the code actually does

`hover_plots` takes exactly one input, `all_qc_stats.tsv` (`sc_predictions.smk:210`, produced by rule `plot_stats`, `sc_predictions.smk:185-205`), and renders `qc_report.Rmd` via `rmarkdown`/`knitr`→`pandoc` into a single self-contained HTML file.

Inside `qc_report.Rmd`:
- `import_data` (`qc_report.Rmd:20-24`) reads two **small** tables with `readr::read_tsv`: the pipeline's own per-cluster stats (`snakemake@input$qc_stats`) and a static reference file (`snakemake@params$reference_clusters` = `config["qc_reference"]` = `resources/reference_qc_metrics_sheth_qiu_2024.tsv`). **It never opens a prediction table, a candidate-region file, or a fragment/UMI file.** All of Tier 1/2/3's `size_manifest` variables (`n_frag`, `n_umi`, `n_cells`, `n_peaks`, `n_kendall_pairs`, `n_cand_pairs`, `j_cand_elems`, `n_called_all`, `n_called_distal`, `k_genes`) are structurally invisible to this rule — they were already collapsed into per-cluster scalar stats by `get_stats_per_model_per_cluster` (`sc_predictions.smk:164-182`) upstream.
- `create_tabs` (`qc_report.Rmd:882-903`) does `models <- unique(stats$model_name)` and calls `knitr::knit_expand(file = tab_template, model = m)` **once per unique model name**, then knits the concatenated result as a child Rmd (`qc_report.Rmd:902`).
  - **Correction to the measurement brief's presumption**: the brief guessed the tab template is "instantiated once per cluster (6 here)." Reading the code shows this is wrong — it is instantiated once **per model** (`unique(stats$model_name)`), not per cluster. In both benchmarked runs `unique(stats$model_name)` has length **1** (`multiome_powerlaw_v3` / `scATAC_powerlaw_v3` — verified directly: `awk` over `all_qc_stats.tsv` shows a single distinct `model_name` value in each modality). Clusters (`C=6` both modalities) enter only as **rows of data bound into existing plot objects**, not as additional plot objects.
- Each tab template instantiation (`qc_plot_tab_template.Rmd`) calls `create_cluster_plots_row()` (defined `qc_report.Rmd:452-515`) **9 times**, once per property (`n_enh_gene_links`, `n_enh_elements`, `n_genes_with_enh`, `n_genes_active_promoter`, `n_genes_not_expressed`, `mean_genes_per_enh`, `mean_enh_per_gene`, `mean_dist_to_tss`, `mean_enh_width` — lines 24-154 of the template), each producing up to 3 `plotly` subplots (one per available x-threshold column: `fragments_total`, `cell_count`, `umi_count`), plus **1 call** to `create_violin_plot_grid()` (template line 164) over 9 variables. All of these are `plotly` (`library(plotly)`, `qc_report.Rmd:14`) scatter/violin/box traces built with `plot_ly()`/`add_trace()`/`add_boxplot()` (`qc_report.Rmd:245, 272, 305, 327, 558, 586`) and combined with `subplot()`.
- `interactive_table` (`qc_report.Rmd:868-880`) renders the entire (6-row) stats table with `DT::datatable(..., extensions=c("FixedColumns","Scroller"))`.
- Output is `html_document` (front-matter, `qc_report.Rmd:5-11`), which defaults to `self_contained: true` — every JS/CSS dependency (plotly.js, DT/datatables, jQuery, Bootstrap) is base64/text-embedded directly into the one output HTML, and every `plotly`/`DT` widget's data is serialized to an inline `<script type="application/json">` block rather than referencing an external file.

## 3. Size variables

Named exactly from `size_manifest_*.tsv` where they apply, plus two structural variables the rule actually uses that are **not** in the size manifest at all (they live in `all_qc_stats.tsv` row/column structure, not in `size_manifest`):

- **None** of the manifest's Tier 1/2/3 variables (`n_frag`, `n_umi`, `n_cells`, `n_peaks`, `n_kendall_pairs`, `n_cand_pairs`, `j_cand_elems`, `k_genes`, `n_called_all`, `n_called_distal`, `n_feat_cols`) are read by this rule's code, directly or transitively through its one input file. They govern upstream rules (through `get_stats_per_model_per_cluster`), not this one.
- `C` — number of cluster-rows in `all_qc_stats.tsv` (= number of clusters run through the pipeline for this modality). Measured `C = 6` for both modalities (both `all_qc_stats.tsv` files are 6 data rows + 1 header, confirmed by `wc -l`).
- `M` — number of distinct `model_name` values in `all_qc_stats.tsv` (drives the number of knit-expanded template tabs, `qc_report.Rmd:884-895`). Measured `M = 1` for both modalities.
- `R` — number of reference rows contributed by `resources/reference_qc_metrics_sheth_qiu_2024.tsv` for the active model, a **static file checked into the repo**, independent of the run. Verified: 57 rows tagged `multiome_powerlaw_v3`, 57 rows tagged `scATAC_powerlaw_v3` (114 total, filtered per model at `qc_plot_tab_template.Rmd:8-10`). `R = 57` is a repo constant, not something Phase 2 can change by scaling the dataset.
- Output HTML byte size (not a size-manifest column, but explicitly requested): multiome `predictions_qc_report.html` = **5,352,969 bytes** (~5.1 MiB); scATAC = **5,131,788 bytes** (~4.9 MiB). Nearly identical (~4% apart) despite `n_called_all`/`n_called_distal` differing ~2× between modalities upstream — consistent with this rule never reading those tables.

## 4. Asymptotic class derived from code

- **In the size-manifest variables** (`n_frag`, `n_umi`, `n_cand_pairs`, etc.): **O(1)** — the code contains no loop, join, sort, or read over any of them. They are fully absorbed by the upstream `get_stats_per_model_per_cluster` rule before this one ever runs.
- **In `C` (clusters)**: every place `C` appears is as the row count of a `dplyr::filter`/`bind_rows`/`add_trace(data=...)` call feeding a *fixed* number of already-existing plot objects (`create_cluster_plot`, `qc_report.Rmd:156-347`; `create_violin_plot`, `qc_report.Rmd:518-646`). No new widgets are created as `C` grows — more clusters just means more points plotted inside existing traces. This is **O(C)** but on a per-widget basis it is trivial vector arithmetic (`mutate`, `filter`, `min`/`max` over ≤(C+R)≈63 numbers) — negligible next to fixed toolchain cost at any C in the observed or plausible range (tens to low hundreds of clusters).
- **In `M` (models)**: `create_tabs` (`qc_report.Rmd:882-895`) is a `sapply` over `models`, each iteration knit-expanding and then knitting the full tab template (9 `create_cluster_plots_row` calls × up to 3 plotly subplots each, + 1 `create_violin_plot_grid` of 9 panels = up to ~36 plotly objects per model). This is genuinely **O(M)** in both knitr chunk-execution count and in JS-widget serialization/embedding cost, since each model gets its own full set of tabs and plots. `M` is fixed at 1 by the current config (one production model per modality: `multiome_powerlaw_v3`, `scATAC_powerlaw_v3`) — if the config were changed to also run e.g. a `megamap_v3` model in the same batch, `M` would become 2 and this rule's cost is the one place in the pipeline where that shows up as a genuine multiplicative *scaling* term, not just a constant.
- **Fixed toolchain term `T_fixed`**: R interpreter + package load (`plotly`, `dplyr`, `DT`), `rmarkdown::render` → `knitr` → `pandoc`, and critically the `self_contained: true` embedding of the JS/CSS dependency bundle. Verified directly on the multiome output: of 5,352,969 total bytes, 4,899,906 (91.5%) sit inside `<script>` tags; a single script block is **3,665,850 bytes** (68.5% of the whole file) — almost certainly the base64/text-embedded `plotly.js` library, which pandoc/htmlwidgets includes **once per document regardless of how many plotly widgets or models exist**. This is a pure constant, not a function of `C`, `M`, or any manifest variable.

Net structural form: **`t ≈ T_fixed + M · T_tab(C, R) + O(C + R)` (negligible last term), i.e. O(1) in every size-manifest variable, O(M) in the number of models, and only trivially O(C) in the number of clusters.**

## 5. Benchmark cross-check

From `benchmarks_multiome.tsv` / `benchmarks_scatac.tsv` (`awk -F'\t' -v r=hover_plots 'NR==1||$2==r' ...`):

| modality | wall `s` | `h:m:s` | max_rss (MB) | io_in (MB) | io_out (MB) | mean_load | cpu_time |
|---|---:|---:|---:|---:|---:|---:|---:|
| multiome | 165.00 | 0:02:44 | 116.52 | 81.16 | 0.00 | 0.21 | 0.00 |
| scatac | 161.14 | 0:02:41 | 118.33 | 82.82 | 0.06 | 0.23 | 0.00 |

- **n=1 per modality, and both points have identical `C=6`, `M=1`.** There is no within-modality variation to fit at all — this is a stricter version of the "six points, cannot fit an exponent" warning: here it's **two points, and they don't even vary in the variables that structurally matter (`C`, `M`)**. The only thing crossing modalities lets us check is whether the *constant* differs, and whether upstream-only variables (`n_called_all` etc., ~2× different between modalities) leak through. They do not: wall time differs by only 2.4% (165.00 vs 161.14 s) despite `n_called_all` differing ~2× — **this is exactly what the code predicts** (the rule never reads that table) and is a genuine (weak) confirmation of the "O(1) in manifest variables" claim, not proof by itself.
- **`cpu_time` reads 0.00 for both** — per the stated caveat, psutil does not capture the knitr→pandoc/Rscript child-process tree for an Rmd render. **Do not read this as "the rule uses no CPU."** `mean_load` (0.21–0.23) is consistent with the *parent* process observed by psutil spending most of its 165 s idle/waiting while a child process (the actual R session doing the knit, and pandoc) does the work — the parent's own CPU usage is genuinely near zero, but that is not the same claim as "the job used no CPU."
- **`io_out` reads ≈0** (0.00 MB multiome, 0.06 MB scatac) despite a >5 MB file being written to disk. This is the same child-process blind spot: the actual HTML write happens inside the child R/pandoc process, not the parent psutil is attached to. **Use the on-disk file size (5.35 MB / 4.9 MB), not `io_out`, as the I/O-volume ground truth for this rule.**
- `io_in` (~81-83 MB) is plausible for R package/library loading (R, plotly, dplyr, DT, rmarkdown and their transitive shared libraries) plus the tiny input TSVs — dwarfing the actual data payload (`all_qc_stats.tsv` is ~1 KB).
- No disagreement to flag beyond the expected "code and data agree that manifest variables don't matter here" — there is no code-predicted scaling that the (nonexistent) within-modality spread could contradict.

## 6. Concrete overhead sources

- **Toolchain startup + package load**: conda env activation, R interpreter start, loading `plotly`/`dplyr`/`DT` and their dependencies (`qc_report.Rmd:13-17`) — fixed, paid once per job regardless of data.
- **Self-contained HTML embedding**: `pandoc`'s `self_contained: true` assembly base64/text-embeds the entire JS/CSS dependency graph into the output. Directly measured: one ~3.67 MB script block (almost certainly `plotly.js`) is 68.5% of the 5.35 MB multiome output; total `<script>` content is 91.5% of the file. **This embedding is the single largest overhead source in the rule, and it is a fixed toolchain cost, not a function of any input size** — it would be identical even if `all_qc_stats.tsv` had 1 row or 100.
- **Widget/data-payload duplication**: 19 `application/json` widget-data blocks were found in the multiome HTML (plot + table payloads); these scale with `M` (one full set of ~36 plotly objects per model) but are individually small since each carries ≤(C+R)≈63 data points.
- **Two sequential in-process renders**: `qc_report.Rmd` itself, then a dynamically-written child Rmd (`tmp_file`, `qc_report.Rmd:898-903`) containing the concatenation of all per-model tabs, knit via `child=tmp_file`. This is a real (if currently invisible at M=1) per-model repeated-pass cost: each model's tab content is generated by string templating (`knit_expand`) then **re-parsed and re-knit** as Rmd source, i.e., there are effectively two knitr passes layered (outer document + inner child) rather than one, adding fixed per-render overhead on top of the per-chunk cost.
- **Declared vs. realized threading**: rule declares no `threads:` (`sc_predictions.smk:208-225` has no `threads:` line) → runs at the Snakemake default of 1. Consistent with the observed `mean_load` ≈ 0.2 for the (psutil-visible) parent process; the actual R/pandoc child process is very unlikely to be multi-threaded either — `rmarkdown::render`/`knitr`/`pandoc` are single-threaded by default and nothing in the code spawns parallel workers.
- **Not I/O-bound in any way this measurement can see**: `io_out≈0` is an artifact (see §5), and `io_in` (~82 MB) is dominated by package loading, not data.

## 7. `t_est`, `t_low`, `t_high`

Variables: `C` = cluster count in `all_qc_stats.tsv` (measured 6, both modalities), `M` = distinct `model_name` count (measured 1, both modalities), `R` = reference rows for the active model (repo constant, 57). All times in seconds.

Because `n=1` per modality and both observations share identical `C` and `M`, **these functions cannot be fit from data — only the additive constant `T_fixed` is anchored by measurement (163 s average, 2 points). The coefficients on `M` and `C` below are code-structure-derived bounds, not calibrated fits**, since no observation exists at `M≠1` or `C≠6` to calibrate against.

```
t_est(C, M)  = 140  + M * 23             +  0.05 * (C + R)
t_low(C, M)  = 120  + M * 12             +  0.01 * (C + R)
t_high(C, M) = 150  + (M ** 1.3) * 45    +  0.20 * (C + R)
```

Rationale for the constants:
- `T_fixed ≈ 120–150 s`: R/toolchain startup + pandoc self-contained assembly, bulk of both measured totals (165.00, 161.14 s), estimated from the fact that ≥68% of the output HTML's bytes (the embedded plotly.js library) are a pure per-render constant unrelated to `M` or `C`.
- `M` term: `t_est` uses linear (23 s/tab — a guess allocating the residual after `T_fixed` across the ~36 widgets/model observed); `t_low` assumes the per-model cost is smaller and purely additive (12 s/tab, no super-linear pandoc growth); `t_high` allows for super-linear growth in `M` (exponent 1.3) because a much longer single self-contained document (many tabs × ~36 widgets each, all base64-embedded) plausibly makes pandoc's whole-document HTML assembly (TOC generation, cross-reference resolution) cost more than strictly linear per added tab — **this is a hypothesis about pandoc internals, not something the code or the C=1-point data can confirm.**
- `C + R` term: bounded by generic per-row plotly JSON-serialization cost (single-digit to tens of ms/row); at `C=6, R=57` this is ≤4 s under any of the three formulas — genuinely negligible at current scale, included only so the formula is honest about not being a pure constant.

**Evaluated at the only measured point (C=6, M=1, R=57)**:

| | t_est | t_low | t_high |
|---|---:|---:|---:|
| multiome / scatac (C=6, M=1, R=57) | 166.2 s | 132.6 s | 208.6 s |

Measured values (165.00 s multiome, 161.14 s scatac) fall inside `[t_low, t_high]` and close to `t_est`, as expected since `t_est`'s constant was set to reproduce them.

**On "mean/min/max cluster" evaluation**: the brief asks for evaluation at the mean cluster size and at `thp1_1` (min) / `telohaec_crispri` (max). **This does not apply to `hover_plots`.** The rule is not parameterized by any individual cluster's size at all — it is a single job that always ingests all 6 clusters' already-reduced stats simultaneously via `plot_stats`'s fan-in output. There is no "run at cluster X's size" version of this job; `C` is the *count* of clusters, not any cluster's `n_frag`/`n_umi`/etc. Evaluating "at thp1_1" is a category error for this rule. The one number that matters, `C=6`, is identical across every plausible per-cluster substitution.

## 8. Modality (BOTH): class or only the constant?

**Only the constant differs, and even that barely.** The code path is identical for multiome and scATAC (same `.Rmd` files, same rule). The two size-manifest variables this rule's code structurally touches (`C`, `M`) are equal across modalities in this run (`C=6`, `M=1` both). The only modality-specific inputs are the reference file's per-model row count (`R=57` for both `multiome_powerlaw_v3` and `scATAC_powerlaw_v3`, coincidentally equal) and the label text ("multiome"/"scATAC" strings). Measured wall time differs by 2.4% (165.00 vs 161.14 s), well within what you'd expect from Slurm scheduling noise and R package cache warmth — **not** a scaling signal. There is no evidence of, and no code path for, a class difference between modalities for this rule.

## 9. What the measurement cannot tell us

- **Cannot fit any exponent for anything.** n=1 per modality, and the two points share identical `C` and `M` — there is zero within-run variation in the variables the code actually uses. Even the weaker "six points, 1.01–1.17× spread" situation other Tier-3 rules have does not apply here; this rule has *no* spread at all in its own structural variables.
- **Cannot separate `T_fixed` from `T_tab` (the per-model term).** Since `M=1` in every observation, any split between "fixed render cost" and "cost per model tab" is an assumption from script-byte-fraction reasoning (§4, §6), not a measurement. If a future config run ever sets `M=2` (e.g., adds a megamap model alongside powerlaw in the same batch), that would be the first real data point on the `M` term.
- **`cpu_time=0.00` must not be read as "no CPU used."** psutil does not attach to the knitr/pandoc/Rscript child-process tree for this Rmd render (stated caveat, confirmed by `io_out≈0` despite a >5 MB file being written — the same blind spot affects both CPU and I/O accounting for child processes). Wall time (165.00 / 161.14 s) is the only trustworthy total-cost number here; it cannot be decomposed into "compute" vs. "waiting" from this benchmark data.
- **Cannot rule out super-linear pandoc-side scaling in `M`** (the `t_high` exponent of 1.3) or confirm it — it is a plausible mechanism (whole-document reprocessing for TOC/cross-refs in a much longer self-contained document) with zero supporting or contradicting data at `M=1`.
- **Cannot attribute the 68.5%-of-file plotly.js block definitively without inspecting its literal contents** — I inferred it from its size and position (a single monolithic script, no other JS library of comparable size is plausible here) rather than by parsing/identifying it byte-for-byte.
- **Fan-in / critical-path note** (per the instructions' explicit ask): `hover_plots` waits only on `plot_stats`'s `all_qc_stats.tsv` (`sc_predictions.smk:210`), but `plot_stats` (`sc_predictions.smk:185-205`) is itself a fan-in over all 6 `get_stats_per_model_per_cluster` jobs (`sc_predictions.smk:164-182`, one per cluster×model×threshold), each of which requires that cluster's full `scE2G_predictions.tsv.gz` — i.e., the entire per-cluster ABC/e2g pipeline for **every** cluster. So although `hover_plots` itself is cheap (~165 s wall, and reads only a ~1 KB file), it is a **second-level barrier**: it cannot start until the slowest of all 6 clusters' full pipelines (and then `plot_stats`, 11.7–11.8 s) has finished. A cheap node here can still sit on the wall-clock critical path purely by being the last thing that has to run after the slowest upstream branch — this rule's own cost (~163 s) is a rounding error next to what it waits on (per the briefing, e.g. `create_neighborhoods` alone is 112 min on the slowest cluster), but it is the very last node in the DAG, so its wait time, not its own compute, is what Phase 2 should care about.

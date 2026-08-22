# Unit 12 — `get_stats_per_model_per_cluster` / `get_stats_per_cluster.R`

## 1. Unit identity

- Script: `workflow/scripts/prediction_qc/get_stats_per_cluster.R` (78 lines).
- Rule: `get_stats_per_model_per_cluster`, `workflow/rules/sc_predictions.smk:164-182`.
- Modality: **both**. Wildcards: `{cluster}`, `{model_name}`, `{threshold}`. In this dataset
  `model_name` and `threshold` are constant per modality (`multiome_powerlaw_v3`/`0.177`,
  `scATAC_powerlaw_v3`/`0.174`), so the wildcard that actually varies is `{cluster}` — 6 jobs
  per modality, 12 total.
- No `threads:` declared → runs at 1 core (`sc_predictions.smk:164-182` has only
  `resources: mem_mb=32000`, no `threads:` line).
- Benchmark: `bench("get_stats_per_model_per_cluster", "cluster", "model_name", "threshold")`.

## 2. What the code actually does

Two independent reads, then independent summary passes over each — **no join between them**:

- `workflow/scripts/prediction_qc/get_stats_per_cluster.R:52` — `pred_full <- fread(snakemake@input$pred_full)`.
  `pred_full` is `{cluster}/{model_name}/scE2G_predictions.tsv.gz` — the **full, unthresholded**
  predictions table (all candidate enhancer–gene pairs, all columns). No `select=` is passed to
  `fread`, so every column is parsed even though only 2 are ever used.
- `:53` — `pred_thresholded <- fread(snakemake@input$pred_thresholded) %>% dplyr::filter(class != "promoter")`.
  This is **the line the QC-report "pairs vs links" trap comes from**: `fread` reads the *entire*
  thresholded file (all `n_called_all` rows) before the `dplyr::filter` drops
  `class == "promoter"` rows down to `n_called_distal`. So for this specific read+filter step,
  cost scales with `n_called_all`, not the smaller `n_called_distal` that survives the filter.
- `get_stats_from_pred()` (`:10-46`), called at `:57`:
  - `:13` `n_genes_active_promoter` — `pred_full %>% select(TargetGene) %>% distinct() %>% nrow()` — full hash-distinct over `pred_full` (`n_cand_pairs` rows).
  - `:14-15` `n_genes_low_TPM` — `pred_full %>% filter(score_column == 0) %>% select(TargetGene) %>% distinct()` — another full scan+distinct over `pred_full`.
  - `:18` `enh_gene <- pred_thresholded %>% select(chr,start,end,TargetGene)` — operates on the **already-filtered**, `n_called_distal`-sized frame.
  - `:21-37` `n_enh`, `n_gene`, `n_links`, `n_gene_per_enh`, `n_enh_per_gene`, `mean_distance`, `mean_size` — all single linear/hash passes (`distinct()`, `group_by()+tally()`, `pull()+mean()`, `mutate()`) over `enh_gene`, i.e. over `n_called_distal` rows.
- `:59-77` — trivial metadata: wildcard values, three 1-line count files (`fragment_count`,
  `cell_count`, `umi_count`) read with `readLines`, plus a `file.info()$isdir` check for
  checkpoint placeholder dirs (`get_count_file`, `sc_predictions.smk:156-162`) — O(1).
- `:77` `fwrite(res, ...)` — writes a 1-row, 9-column TSV. O(1).

No joins, no nested loops, no sort-merge anywhere in this script. Every operation is a single
pass (scan, hash-distinct, or hash-group) over one of two tables.

## 3. Size variables (named as in `size_manifest_*.tsv`)

| Variable | Role in this script | Range (multiome) | Range (scatac) |
|---|---|---|---|
| `n_cand_pairs` | rows of `pred_full` (verified empirically, see §5) — read at `:52`, scanned twice at `:13-15` | 10.17M–11.45M (1.13×) | **identical row counts** to multiome (see §5) |
| `n_called_all` | rows of `pred_thresholded` **before** the `:53` filter — read at `:53` | 59,438–61,284 | 105,880–120,113 |
| `n_called_distal` | rows of `pred_thresholded` **after** the `:53` filter — everything at `:18-37` | 46,155–47,560 | 89,681–104,046 |

None of these is a manifest column for it, but the *width* of `pred_full`/`pred_thresholded`
also matters and differs by modality: **27 columns (multiome) vs 20 columns (scatac)**,
measured directly from the file headers (not in `size_manifest_*.tsv`; see §5). This is the
key fact that resolves the modality "contradiction" flagged for this unit.

## 4. Asymptotic class derived from code

For a single job:

```
T = O(n_cand_pairs · c_pf)        # :52 fread(pred_full), all c_pf columns
  + O(n_cand_pairs)               # :13, :14-15 two full-table distinct()/filter() passes over pred_full
  + O(n_called_all · c_pf)        # :53 fread(pred_thresholded), all c_pf columns, pre-filter
  + O(n_called_distal)            # :18-37 select/distinct/group_by+tally/mutate on the post-filter frame
  + O(1)                          # R/conda startup, package loads, count-file reads, fwrite of a 1-row output
```

where `c_pf` is the modality-dependent column count of the predictions tables (27 multiome /
20 scatac). `dplyr::distinct()` and `group_by()+tally()` are hash-based single passes, not
sorts or joins, so the class is **linear**, not `n log n` or `n²` — there is no complexity-class
ambiguity here (unlike Tier-3-variable rules elsewhere in this project). The two dominant terms
are the two `fread()` calls; `n_cand_pairs · c_pf` is 3–4 orders of magnitude larger than
`n_called_all · c_pf` (see §5), so in practice `T ≈ O(n_cand_pairs · c_pf)`.

## 5. Benchmark cross-check — and the resolved "contradiction"

### 5a. The row-count identity holds *across modalities*, not just within one

`n_cand_pairs` in `size_manifest_multiome.tsv` and `size_manifest_scatac.tsv` is **byte-for-byte
identical per cluster** (e.g. `jurkat`: 11,036,283 both; `thp1_1`: 10,173,077 both). I verified
this is also the row count of `pred_full` directly: `zcat` on
`igvf10_multiome/thp1_1/multiome_powerlaw_v3/scE2G_predictions.tsv.gz` gives 10,173,078 lines
(header + `n_cand_pairs`), and the scATAC counterpart gives the **same** 10,173,078. So the file
this script spends most of its time on (`pred_full`) has the *same row count* in both
modalities — the ARC/Activity-only split does not change the row count, only the column count.

### 5b. Column counts and byte volumes, measured from the files directly (not in the manifest)

```
pred_full header columns:  multiome = 27   scatac = 20
```
Multiome's 7 extra columns are `normalizedATAC_enh`, `Kendall`, `ARC.E2G.Score`,
`RNA_meanLogNorm`, `RNA_pseudobulkTPM`, `RNA_percentCellsDetected`,
`E2G.Score.qnorm.ignoreTPM` — exactly the ARC/RNA features that only exist on the multiome
path. `pred_thresholded` has the same per-modality column count as `pred_full` (verified:
27 for multiome `thp1_1`, 20 for scatac `thp1_1`).

Uncompressed byte size of `pred_full` (via `gzip -l`, all 6 clusters):

| | multiome (bytes) | scatac (bytes) |
|---|---:|---:|
| mean | 3.48×10⁹ | 2.21×10⁹ |
| min (thp1_1) | 3.23×10⁹ | 2.03×10⁹ |
| max (telohaec_crispri) | 3.64×10⁹ | 2.33×10⁹ |

Ratio ≈ **1.57×**, driven by column count (27/20 = 1.35×) plus wider numeric columns among the
extra 7. `pred_thresholded` is tiny by comparison — 4.7–4.9 MB compressed (multiome) vs
4.9–5.5 MB compressed (scatac) — three orders of magnitude smaller than `pred_full`.

### 5c. This directly resolves the apparent contradiction flagged for this unit

The task brief for this unit says: *"`n` is ~2× larger in scATAC but time is ~1.2× smaller —
investigate."* The "~2× larger `n`" refers to `n_called_all`/`n_called_distal` (the thresholded
file). But that file is not the cost driver — it is ~200× smaller by row count and ~150-180×
smaller by bytes than `pred_full`, which the script also reads in full (`:52`) and scans twice
(`:13-15`). `pred_full` has the **same row count** in both modalities but is **1.57× more
bytes in multiome** (27 vs 20 columns). Multiome being slower despite scATAC's larger
`n_called_all`/`n_called_distal` is exactly what you'd expect if `pred_full`, not the thresholded
file, dominates cost — **this is not a contradiction once the dominant read is identified
correctly; it is consistent with `pred_full` bytes being the real driver.**

### 5d. Numeric cross-check

Benchmark rows (`benchmarks_multiome.tsv` / `benchmarks_scatac.tsv`, `s` = wall, `cpu_time` in
seconds):

| cluster | multiome wall | multiome cpu | scatac wall | scatac cpu |
|---|---:|---:|---:|---:|
| jurkat_pma_cd3_4hr | 27.53 | 14.53 | 19.64 | 14.95 |
| jurkat | 28.82 | 13.40 | 20.93 | 12.94 |
| k562_crispri | 28.71 | 14.69 | 20.11 | 14.92 |
| telohaec_crispri | 29.19 | 14.04 | 25.68 | 6.20 |
| thp1_1 | 27.77 | 12.39 | 20.30 | 11.93 |
| thp1_2 | 25.23 | 15.01 | 27.17 | 7.00 |
| **mean** | **27.88** | **14.01** | **22.31** | **11.32** |

Mean wall ratio multiome/scatac = 27.88/22.31 = **1.25×** — matches the ~1.2× noted in the
briefing, and is in the *same direction and rough order of magnitude* as the 1.57× byte ratio
of `pred_full` (not exact, because a sizeable fixed cost is present too — see §7 fit).
`n_called_all`/`n_called_distal` (the ~2× variable) demonstrably does **not** set the direction
of this difference.

Within a modality, `n_cand_pairs` spans only 1.13×, `n_called_all`/`n_called_distal` span
1.01–1.17×, and cpu_time shows **no consistent monotonic trend** against any of them
(e.g. multiome `thp1_2`, `n_cand_pairs`=10.28M, the *second-smallest* value, has the
*largest* cpu_time, 15.01 s). **Within a modality the six points cannot confirm the linear
class** — the code justifies it, but the range is too narrow and the noise (see §6) too large.
The only clean signal in this dataset is the **cross-modality, same-row-count** comparison
in §5c, which the briefing explicitly licenses ("where a rule appears in both modalities with
the same inputs, you have n=2 per cluster; use it").

## 6. Concrete overhead sources

- **Two full reads of large-ish files, one of them fully unnecessary in width.** `:52`
  `fread(pred_full)` parses all 27 (multiome) / 20 (scatac) columns of a 742 MB–889 MB (multiome)
  / 414 MB–476 MB (scatac) gzip file, decompressing to 3.2–3.7 GB / 2.0–2.35 GB, purely to
  compute two `distinct(TargetGene)` counts (`:13-15`) that need only 2 columns
  (`TargetGene`, the score column). No `select=` is used, so `fread` reads every column.
- **Fixed per-job overhead dominates at this scale.** Fitting `wall ≈ s0 + k·c_pf·n_cand_pairs`
  to the two modality means (§7) gives `s0 ≈ 6.4 s` — i.e. roughly a quarter of every job's wall
  time is a size-independent floor (conda env activation, R interpreter start, package loads,
  the Snakemake `script:` wrapper, Slurm dispatch), not proportional work.
- **Unused package loads.** `:1-8` loads `plyr`, `tidyr`, and `ggplot2` in addition to
  `dplyr`/`data.table`/`stringr`; grepping the script shows no `plyr::`, no `tidyr` verb
  (`pivot_*`, `separate`, `unite`), no `ggplot2`/`str_*` call anywhere in this file. Those three
  package loads cost real wall time at every one of the 12 jobs for zero functional benefit in
  this script.
- **`io_in`/`io_out` in the benchmark table do not reflect the real I/O volume, and this is
  worth flagging rather than smoothing over.** `io_in` for these jobs is 0.95–30.6 MB
  (multiome) / 0.2–304 MB (scatac) — two to three orders of magnitude smaller than the true
  ~750–890 MB (multiome) / ~414–476 MB (scatac) compressed bytes actually read from
  `pred_full`+`pred_thresholded`. Most likely explanation: the files are already in the node's
  page cache by the time this job runs (produced or read moments earlier in the same DAG
  branch), so `psutil`'s I/O counters see little fresh block-device traffic — the *parsing* cost
  is still real and paid in CPU time, but it is not "I/O-bound" in the sense the `io_in` column
  would suggest. `io_out` is even more puzzling: 2851–3414 MB (multiome) / 1497–2246 MB
  (scatac) of *writes*, for a job whose only output is a 1-row, 9-column TSV. This cannot be the
  real output volume and is a measurement artifact of this rule's benchmark row — I flag it but
  cannot explain its source from the benchmark table alone (see §8).
- **Declared vs. realized parallelism.** No `threads:` declared, so Snakemake reserves 1 core.
  `mean_load` is 44.6–58.8% (multiome) and 23.9–75.3% (scatac) — consistent with a
  mostly-single-threaded job, with data.table's `fread` plausibly using more than one core
  briefly during decompression/parsing (its default `getDTthreads()` can exceed 1 even when the
  calling R process is otherwise serial), never approaching full multi-core use.

## 7. `t_est(n)`, `t_low(n)`, `t_high(n)`

Let `N = n_cand_pairs`, `A = n_called_all`, `D = n_called_distal` (as in `size_manifest_*.tsv`),
and `c_pf` = pred_full/pred_thresholded column count — **27 for multiome, 20 for scatac**
(measured directly from file headers; not itself a manifest column).

**Fit.** Using the modality-mean wall times (27.88 s at `c_pf=27`, 22.31 s at `c_pf=20`, both at
mean `N=10,816,106`) as two calibration points for `wall ≈ s0 + k·c_pf·N`:

```
k  = 7.35e-8  s per (row · column)
s0 = 6.39     s
```

```
t_est(N, A, D) = 6.39 + 7.35e-8 · c_pf · (N + A + D)      [seconds]
```
(`A` and `D` are included at the same per-cell rate as `N` for formula simplicity; because
`A, D ≤ ~120,000` vs `N ≈ 10.8M`, they change `t_est` by at most ~2% regardless of their true
rate — the dominant term is `7.35e-8 · c_pf · N`. `D`'s operations are in-memory
`dplyr`/hash ops on an already-parsed frame and are almost certainly cheaper per row than a
cold `fread`, so this slightly *overstates* their contribution — flagged as a simplification,
not a calibrated fact.)

```
t_low(N, A, D) = 19.6      [seconds, constant]
```
Lower bound on the *class*: the hypothesis that none of the size variables contribute
detectably at this row-count range (consistent with the lack of a monotonic trend within a
modality, §5d) — a pure fixed-cost floor, set at the minimum wall time actually observed
(scatac `jurkat_pma_cd3_4hr`, 19.64 s).

```
t_high(N, A, D) = 8 + 1.5e-7 · c_pf · (N + A + D)      [seconds]
```
Upper bound on the *class*: the full code-justified linear form (§4) with a coefficient ~2×
the fitted `k`, to bound the uncertainty in the per-cell rate itself and to remain valid if a
future cluster has a much larger `N` than the 1.13× range seen here.

## 8. Evaluated at mean / min (`thp1_1`) / max (`telohaec_crispri`)

| Modality | Point | N | A | D | c_pf | t_low | t_est | t_high | observed wall |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| multiome | mean (6 clusters) | 10,816,106 | 60,189 | 46,771 | 27 | 19.6 | 28.07 | 52.0 | 27.88 |
| multiome | min (`thp1_1`) | 10,173,077 | 60,435 | 46,734 | 27 | 19.6 | 26.79 | 48.9 | 27.77 |
| multiome | max (`telohaec_crispri`) | 11,106,687 | 60,630 | 47,437 | 27 | 19.6 | 28.65 | 53.4 | 29.19 |
| scatac | mean (6 clusters) | 10,816,106 | 113,430 | 97,199 | 20 | 19.6 | 22.60 | 41.1 | 22.31 |
| scatac | min (`thp1_1`) | 10,173,077 | 105,880 | 89,681 | 20 | 19.6 | 21.63 | 39.3 | 20.30 |
| scatac | max (`telohaec_crispri`) | 11,106,687 | 117,185 | 100,864 | 20 | 19.6 | 23.04 | 42.1 | 25.68 |

`t_est` tracks observed wall time within ~1–3 s (≤ ~12%) at every one of the 6 real points
despite being calibrated from only the 2 modality means — the largest miss is scATAC
`telohaec_crispri` (observed 25.68 s vs. `t_est` 23.04 s), which is also the point with the
anomalously low `cpu_time` (6.20 s), suggesting most of that gap is scheduling/queueing noise
rather than a model error.

## 9. Both modalities — class or constant?

**Only the constant differs, not the class.** Both modalities execute the identical control
flow (§2) over tables with the identical row count `N = n_cand_pairs` (§5a) — the only thing
that changes is `c_pf` (27 vs 20 columns), which scales the per-row cost of the two `fread()`
calls, and `A`/`D` (`n_called_all`/`n_called_distal`), which differ ~2× by modality but
contribute negligibly to the total (§7). The single two-point fit in §7 uses one `s0` and one
`k` for both modalities and reproduces all 6×2 = 12 observed wall times to within ~1–3 s, which
is itself evidence that a single linear class with a modality-dependent width constant is
sufficient — no separate exponent is needed per modality.

## 10. What the measurement cannot tell us

- **Cannot separate `N`, `A`, `D` from the benchmark data.** `N` is literally identical across
  modality (same rows); `A`/`D` differ ~2× by modality but span only 1.01–1.17× *within* a
  modality. Six points per modality cannot fit three coefficients; the class-level statement
  in §4 comes entirely from reading the code, not from regression.
- **Cannot rule out mild super-linear behavior in `fread`/`dplyr` internals** (e.g. hash-table
  resizing, gzip-stream buffering effects) — only one column-width contrast exists in this
  dataset (27 vs 20), which is one data point on the "does cost scale with column count"
  question, not enough to distinguish linear from a slightly worse-than-linear dependence on
  `c_pf`.
- **Cannot attribute the ~2.7 s gap between the wall-based fixed cost (`s0≈6.39s`) and the
  cpu_time-based fixed cost (`≈3.6s`, from the same two-point method applied to `cpu_time`
  means of 14.01/11.32 s)** to a specific source (Slurm dispatch latency vs. conda activation
  vs. cold page cache) without finer-grained timestamps than the benchmark table provides.
- **Cannot explain the `io_out` values** (1.5–3.4 GB per job, for an output file of a few
  hundred bytes) from the benchmark table alone (§6) — flagged as an artifact, not resolved.
- **Cannot explain the two scatac cpu_time outliers** (`telohaec_crispri` 6.20 s,
  `thp1_2` 7.00 s, both well below their own modality's other 4 points and below what `t_est`
  implies) — these are single measurements (n=1), not replicates, so there is no way to tell
  noise from a real effect (e.g. transient node contention) with this data.
- **Cannot confirm the linear class empirically within a modality** — only the cross-modality
  (same-row-count, different-column-count) comparison in §5c/§5d provides a clean, direction-
  consistent signal; everything else in this unit's six-point range is consistent with "the
  class exists per the code, but is swamped by fixed cost and job-to-job noise at this scale."

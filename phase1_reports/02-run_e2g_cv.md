# Phase 1 unit #2 — `run_e2g_qnorm` / `run_e2g_cv.py`

Tier A. Model application and quantile-normalised scoring. Characterisation only — no edits,
no fix proposals, no prioritisation.

---

## 1. Unit identity

| item | value |
|---|---|
| Script | `workflow/scripts/model_application/run_e2g_cv.py` — **169 lines** |
| Rule | `run_e2g_qnorm`, `workflow/rules/sc_predictions.smk:47` (shell rule, `sc_e2g.yml` conda env) |
| Shell invocation | `workflow/rules/sc_predictions.smk:68-77` |
| Wildcards | `cluster`, `model_name` → benchmark stem `<cluster>~<model_name>` |
| Input | `<cluster>/genomewide_features.tsv.gz` (one file, gzipped TSV) |
| Output | `<cluster>/<model_name>/scE2G_predictions.tsv.gz` |
| Modalities | **both** |
| Jobs | 6 multiome (`multiome_powerlaw_v3`) + 6 scATAC (`scATAC_powerlaw_v3`) = 12 measured jobs |
| Declared `threads:` | **none** → 1. `profiles/measure/config.yaml` `cpus_per_task=1`. |
| Declared memory | `mem_mb=encode_e2g.ABC.determine_mem_mb` (`ENCODE_rE2G/ABC/workflow/rules/utils.smk:17-24`) |
| Measured wall | multiome 464.8–517.1 s (7.7–8.6 min); scATAC 233.9–294.9 s (3.9–4.9 min) |
| Measured CPU/wall | multiome 0.86–0.96; scATAC 0.70–0.95 → **single-core, partly I/O-stalled** |

Sibling scripts, checked as instructed:

- **`run_e2g_cv_qnormref.py` (162 lines) is dead.** `grep -rn "run_e2g_cv_qnormref"` over the whole
  working tree returns only two hits, both inside the file itself
  (`run_e2g_cv_qnormref.py:143` and `:154`, commented-out example invocations). No `.smk`, no `.py`,
  no config references it. **Not on the path.** Not characterised further.
- **`workflow/scripts/model_application/get_fill_values.R` is NOT on this rule's path.** It is
  `source()`d only by `merge_features_with_crispr_data_apply.R:8,140`, which is the script of rule
  `overlap_features_crispr_apply` (`workflow/rules/sc_predictions.smk:1-19`) — a *consumer* of this
  rule's output, not a producer of its input. `run_e2g_cv.py` does its own NA handling inline
  (`run_e2g_cv.py:139-140`) and hard-codes fill value `0` for every column. Note in passing: it
  reads the model's `feature_table.tsv` but uses **only** the `feature` column
  (`run_e2g_cv.py:134-135`), so the table's declared per-feature `fill_value` is read and discarded.
  For `multiome_powerlaw_v3` / `scATAC_powerlaw_v3` every declared `fill_value` is `0`, so there is
  no behavioural divergence at present. (Unit #22 owns the ENCODE_rE2G copy.)

---

## 2. What the code actually does

Single-process, single-thread, straight-line pandas. Twelve stages, all but one of them a full
pass over the whole table.

| # | stage | `file:line` | cost shape |
|---|---|---|---|
| 1 | `pd.read_csv(predictions, sep="\t")` — gzip inflate + C tokenise + box | `run_e2g_cv.py:138` | `n · C` text bytes |
| 2 | `.replace([inf,-inf], nan)` — full-frame, returns a **new frame** | `run_e2g_cv.py:139` | `n · C` read + `n · C` write |
| 3 | `.fillna(0)` — full-frame, returns another **new frame** | `run_e2g_cv.py:140` | `n · C` read + `n · C` write |
| 4 | read qnorm reference (10,001 rows) | `run_e2g_cv.py:143-146` | O(1) |
| 5 | `X = df.loc[:, feature_list]` then `np.log(np.abs(X) + eps)` — 4 temporaries of `n × k_feat` | `run_e2g_cv.py:13-14` | `n · k_feat` × 4 passes |
| 6 | `pickle.load` + `model.predict_proba(X)` | `run_e2g_cv.py:16-18` | `n · k_feat` matvec + `n` sigmoid |
| 7 | `df["E2G.Score"] = probs[:,1]` | `run_e2g_cv.py:20` | `n` |
| 8 | **CV loop over chromosomes** — *skipped in every measured job*, see §2.2 | `run_e2g_cv.py:25-43`, gated at `:152` | `n · n_chrom` |
| 9 | `calculate_quantiles` → `Series.rank(method="average")` — **the sort** | `run_e2g_cv.py:46-55`, sort at **`:48`** | `n log n`, once per score column |
| 10 | `interp1d` build (10,001 pts) + evaluate + `.clip(0)` | `run_e2g_cv.py:64-66`, `:71` | `n log 10001` |
| 11 | two `print(... .describe())` calls | `run_e2g_cv.py:72-73` | ~4–6 extra `n` passes, **stdout only** |
| 12 | `filter_by_tpm` — **pure-Python list comprehension over every row** | `run_e2g_cv.py:92-98` | `n` Python-level iterations, **multiome only** |
| 13 | `to_csv(compression="gzip")` — Python float→text + **zlib level 9** | `run_e2g_cv.py:165` | `n · (C+2 or C+3)` text bytes |

with `n` = `n_cand_pairs`, `C` = `n_feat_cols`, `k_feat` = number of model features.

Nothing is ever subset: all 24 (multiome) / 18 (scATAC) input columns are carried through every
copy and re-serialised into the output, even though only 6 features + `chr` + `RNA_pseudobulkTPM`
participate in any computation. `read_csv` is called with no `usecols=` and no `dtype=`
(`run_e2g_cv.py:138`).

Stages 2 and 3 apply to **object** columns as well (`chr`, `name`, `class`, `TargetGene`,
`TargetGeneEnsembl_ID`, `CellType` — 6 columns in *both* modalities). Comparing a string column to
`±inf` is a no-op by construction; the pass is paid anyway.

### 2.1 `k_feat` is 6 in BOTH modalities — the ARC feature does not widen the model

`models/multiome_powerlaw_v3/feature_table.tsv` and `models/scATAC_powerlaw_v3/feature_table.tsv`
each declare exactly **6** features. They differ in one row only:

```
multiome: numTSSEnhGene normalizedATAC_prom numNearbyEnhancers ubiqExpressed numCandidateEnhGene ARC.E2G.Score
scATAC:   numTSSEnhGene normalizedATAC_prom numNearbyEnhancers ubiqExpressed numCandidateEnhGene ABC.Score
```

`ARC.E2G.Score` **substitutes for** `ABC.Score`; it is not added to it. `X` is `n × 6` in both
modalities, so stages 5–7 are byte-for-byte identical in cost.

The model itself is a **`sklearn.linear_model.LogisticRegression`** — verified from
`strings models/multiome_powerlaw_v3/model.pkl` (`sklearn.linear_model._logistic`,
`LogisticRegression`, `solver=lbfgs`, `feature_names_in_` listing the 6 names above, sklearn
1.2.1). The pickle is **956 bytes**: a 6-element `coef_` plus an intercept. `predict_proba` is one
`n × 6` dot product plus a sigmoid — a few seconds at `n ≈ 11 M`. **Model application is not where
the time goes.**

### 2.2 The "CV" in the name is dead in this configuration — folds = 23, sequential, apply-only

`run_e2g_cv.py:152` gates the CV branch on `str(crispr_benchmarking) == "True"`. The rule passes
`config["benchmark_performance"]` (`workflow/rules/sc_predictions.smk:57`), and that is **`False`**
in all four relevant configs:

```
configs/igvf10_config.yaml:72          benchmark_performance: False
configs/igvf10_multiome_config.yaml:61 benchmark_performance: False
configs/igvf10_scatac_config.yaml:64   benchmark_performance: False
config/config.yaml:10                  benchmark_performance: False
```

Confirmed on the output side: `scE2G_predictions.tsv.gz` carries **no** `E2G.Score.cv` or
`E2G.Score.cv.qnorm` column in either modality, and the slurm log for
`thp1_1~multiome_powerlaw_v3` prints exactly **two** `describe()` blocks (`E2G.Score`,
`E2G.Score.qnorm`) — the CV `describe()` never fires. **Every measured second of this rule is the
non-CV path.**

For the record, since the prompt asks:

- **Folds = 23**, one per chromosome. `ls models/*/cv_models` = 23 pickles,
  `model_test_chr{1..22}.pkl` + `model_test_chrX.pkl`, in every one of the four model directories.
  There is no `chrY`, no `chrM`.
- **Sequential.** `for chr in chr_list:` at `run_e2g_cv.py:33`, no parallelism of any kind.
- **The model is applied, never refit.** `pickle.load` at `run_e2g_cv.py:38-39` then
  `predict_proba` at `:40`. Each fold pickle is 956 bytes — a pre-trained coefficient vector.
- **The folds partition the rows** (`chr` is a partition of the table), so `predict_proba` totals
  `n`, not `23 n`. But `run_e2g_cv.py:34` re-scans the *whole* frame per fold:
  `df_enhancers[df_enhancers["chr"] == chr]` builds a **full-width DataFrame slice** purely to read
  `.index.values` off it. That is 23 O(`n`) equality scans on an **object** column plus 23 full-width
  materialisations. **This is where a `n_chrom = 23` multiplier genuinely lives** — in the masking,
  not the scoring.
- `run_e2g_cv.py:29-30` recomputes `X` and its log-transform, duplicating `:13-14` exactly.

So the multiplier that "belongs in the constant" is **23, but it currently multiplies zero.** A
code-derived (unmeasured) estimate of turning `benchmark_performance: True` back on is in §7.4.

### 2.3 The quantile normalisation: one global sort per score column, not per feature

Despite the rule name, **the 24/18 feature columns are never quantile-normalised.** What is
normalised is the model's single scalar output. `calculate_quantiles` (`run_e2g_cv.py:46-55`) is
called on **one** column at a time:

- `run_e2g_cv.py:70` → on `E2G.Score`
- `run_e2g_cv.py:77` → on `E2G.Score.cv` (dead, see §2.2)

Each call does exactly one **global** rank over all `n_cand_pairs` values of that column:

```python
run_e2g_cv.py:48   ranks = score_column_zero_replaced.rank(method="average", na_option="top")
```

`method="average"` forces pandas into `algos.rank_1d`, which argsorts the full column →
**Θ(`n_cand_pairs` · log `n_cand_pairs`), k sorts of n with k = 1** (2 if CV were enabled).
It is a global sort, not per-column, and not `k_feat`-wide.

At `n = 11 M` that is `1.1e7 × 23.4 ≈ 2.6e8` comparison-equivalents on a contiguous float64 array —
**order 1–3 s, under 1 % of the job.** The sort implied by "quantile normalisation" is real, is
correctly `n log n`, and is *not* the bottleneck.

The reference side is O(1): `qnorm_reference.tsv.gz` is **10,001 rows × 2 columns**
(`quantile`, `reference_score`) in every model dir, ~132 KB gzipped. `interp1d(kind="linear")`
build is O(10001 log 10001); evaluation (`run_e2g_cv.py:71`) is a vectorised searchsorted,
Θ(`n` log 10001).

### 2.4 `filter_by_tpm` runs in multiome and short-circuits in scATAC — two independent reasons

`run_e2g_cv.py:86`:

```python
if ("RNA_pseudobulkTPM" not in df_enhancers.columns) or (tpm_threshold == 0):
    return df_enhancers  # don't filter
```

Both disjuncts fire for scATAC:

1. **Column absent.** Verified from the on-disk headers of `thp1_1/genomewide_features.tsv.gz`:
   multiome carries 24 columns including `RNA_pseudobulkTPM` (col 21); scATAC carries 18 and the
   six missing ones are exactly `normalizedATAC_enh`, `RNA_meanLogNorm`, `RNA_pseudobulkTPM`,
   `RNA_percentCellsDetected`, `Kendall`, `ARC.E2G.Score`.
2. **Threshold is zero.** `get_tpm_threshold` (`ENCODE_rE2G/workflow/rules/utils.smk:194-202`) globs
   `tpm_threshold_*` in the model dir. `models/multiome_powerlaw_v3/tpm_threshold_1` → `1.0`;
   `models/scATAC_powerlaw_v3/tpm_threshold_0` → `0.0`.

So in **multiome only**, `run_e2g_cv.py:92-98` executes a **pure-Python list comprehension with
`zip` over two pandas Series**, `n_cand_pairs` ≈ 10.2–11.4 M iterations, each boxing two numpy
scalars into Python objects, building a 10-million-element Python list, which pandas then converts
back to a float64 column. This is the highest *unit-cost* operation in the script by a wide margin,
and it is modality-gated. It also adds a third output column, `E2G.Score.qnorm.ignoreTPM`.

Output widths, verified on disk: **multiome 27 columns** (24 + `E2G.Score` + `.qnorm` +
`.qnorm.ignoreTPM`); **scATAC 20 columns** (18 + `E2G.Score` + `.qnorm`).

---

## 3. Size variables

Named exactly as in `size_manifest_{multiome,scatac}.tsv`.

| variable | role in this rule | multiome range | scATAC range |
|---|---|---|---|
| `n_cand_pairs` | rows of the input and of the output. **Identical across modalities.** | 10,173,077 – 11,447,564 (1.13×) | *same six values* |
| `n_feat_cols` | width of the input frame → width of every copy, and of the output +3/+2 | **24** (const) | **18** (const) |
| `feat_bytes` | gzipped bytes of `genomewide_features.tsv.gz` — the read volume | 494,990,740 – 605,033,216 (1.22×) | 189,044,508 – 222,648,282 (1.18×) |

Confirmed `n = n_cand_pairs` exactly: the slurm log for `thp1_1~multiome_powerlaw_v3` prints
`count 1.017308e+07` for `E2G.Score`, against `n_cand_pairs = 10,173,077`.

Derived quantities used below (not in the manifest, measured here):

| quantity | multiome | scATAC | ratio |
|---|---|---|---|
| uncompressed **input** bytes/row | 248.5 | 148.4 | 1.675× |
| uncompressed **output** bytes/row | 309.8 | 193.0 | 1.605× |
| uncompressed in+out volume, `thp1_1` | 5.68 GB | 3.47 GB | 1.635× |
| `k_feat` (model features) | 6 | 6 | **1.00×** |
| `n_chrom` (CV folds, currently unused) | 23 | 23 | 1.00× |

(bytes/row from `zcat … | head -n 100001 | wc -c` on `thp1_1`.)

`n_called_all` / `n_called_distal` are **irrelevant to this rule** — thresholding happens
downstream in `filter_sc_e2g_predictions` (`workflow/rules/sc_predictions.smk:79`). This rule
writes all `n_cand_pairs` rows. The addendum's warning about scATAC's ~2× `n_called_all` does not
touch this unit.

---

## 4. Asymptotic class derived from code

$$
T \;=\; \Theta\!\big(\texttt{n\_cand\_pairs}\cdot\texttt{n\_feat\_cols}\big)
\;+\;\Theta\!\big(\texttt{n\_cand\_pairs}\cdot\log \texttt{n\_cand\_pairs}\big)
\;+\;\Theta\!\big(\texttt{n\_cand\_pairs}\cdot k_{\text{feat}}\big)
\;+\;O(1)
$$

with, if `benchmark_performance` were `True`, an additional
$\Theta(\texttt{n\_cand\_pairs}\cdot n_{\text{chrom}})$ term, $n_\text{chrom}=23$.

The **dominant** term is the first — the linear-in-volume cross term. Justification, term by term:

**`n_cand_pairs · n_feat_cols` (dominant).**
- `run_e2g_cv.py:138` `read_csv` must inflate and tokenise every cell.
- `run_e2g_cv.py:139` `.replace([inf,-inf], nan)` — one full-frame elementwise pass **plus a full
  frame copy**.
- `run_e2g_cv.py:140` `.fillna(0)` — a second full-frame pass **plus a second full frame copy**.
- `run_e2g_cv.py:165` `to_csv(compression="gzip")` must format and compress every cell of a frame
  `n_feat_cols + 3` (multiome) / `+ 2` (scATAC) wide.

**`n_cand_pairs · log n_cand_pairs`.** The single global `Series.rank` at `run_e2g_cv.py:48`
(§2.3). One sort per score column, `k = 1`. This is the only super-linear operation in the entire
script. Measured share at `n ≈ 11 M`: under ~1 %.

**`n_cand_pairs · k_feat`.** `run_e2g_cv.py:13-14` (4 passes over `n × 6`), `:18` `predict_proba`
(one `n × 6` matvec + `n` sigmoid), plus sklearn's internal `check_array` C-contiguous copy.
`k_feat = 6` in both modalities.

**O(1).** Snakemake wrapper (the slurm log shows the DAG built **twice**, ~12 s before the script
starts), conda env activation off Lustre, `import pandas/numpy/scipy/sklearn`, `qnorm_reference`
(10,001 rows) and `model.pkl` (956 B) reads.

**There is nothing quadratic anywhere in this script.** No nested pass over `n`, no self-join, no
`j_cand_elems²` inner loop, no `k_genes` cross term. That is a positive finding and it is why the
`t_low`/`t_high` band in §7 is narrow: the class band is `[Θ(V), Θ(V log n)]`, not
`[Θ(n), Θ(n²)]`.

### Where the constant actually sits: gzip level 9 on the output

`pandas.to_csv(compression="gzip")` supplies no `compresslevel`, so it falls through to
`gzip.GzipFile`'s default of **9**. Measured rates on real slices of the actual output files
(`/usr/bin/time`, user CPU, so login-node contention is excluded):

| operation | data | user CPU | rate |
|---|---|---|---|
| `gzip -9` multiome-shaped output | 62.17 MB (200 k rows) | 4.10 s | 15.2 MB/s |
| `gzip -9` scATAC-shaped output | 39.06 MB (200 k rows) | 1.82 s | 21.5 MB/s |
| `gzip -6` multiome-shaped output | 62.17 MB | 2.08 s | 29.9 MB/s |
| `zcat` multiome input | 100.85 MB | 0.52 s | **194 MB/s** |

Extrapolating to `thp1_1`:

- multiome output 3.151 GB uncompressed ÷ 15.2 MB/s ≈ **208 s** of the job's 397.9 s `cpu_time`
  → **~52 %**
- scATAC output 1.963 GB ÷ 21.5 MB/s ≈ **92 s** of 195.9 s → **~47 %**
- input decompression: multiome 2.53 GB ÷ 194 MB/s ≈ **13 s** → ~3 %

Decompression is ~40× cheaper per byte than compression. The residual (~190 s multiome,
~104 s scATAC) covers `read_csv` tokenising/boxing, the two full-frame copies, `to_csv`'s
Python-level float formatting, the rank, the two `describe()` calls, and (multiome only) the
`filter_by_tpm` list comprehension. That residual budget is consistent with those items; I cannot
split it further without a profiler run, which is out of scope.

*Caveat:* these rates are from system `gzip` (zlib level 9). Python's `gzip` module also uses zlib
level 9 but writes through a `TextIOWrapper`, so the true in-process cost is **≥** the figures
above, not ≤. Treat 52 % / 47 % as lower bounds on gzip's share.

---

## 5. Benchmark cross-check

### 5.1 The raw numbers (12 jobs)

| modality | cluster | `n_cand_pairs` | `feat_bytes` (MB) | wall s | cpu s | CPU/wall | max_rss MB | io_in MB | io_out MB |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| multiome | thp1_1 | 10,173,077 | 494.99 | 464.76 | 397.93 | 0.86 | 3350 | 328.4 | 689.6 |
| multiome | thp1_2 | 10,278,325 | 505.00 | 473.35 | 426.85 | 0.90 | 2772 | 306.4 | 741.3 |
| multiome | jurkat_pma_cd3_4hr | 10,854,698 | 557.72 | 471.95 | 453.11 | 0.96 | 3619 | 322.4 | 796.3 |
| multiome | jurkat | 11,036,283 | 568.23 | 517.12 | 447.60 | 0.87 | 3796 | 367.1 | 761.3 |
| multiome | telohaec_crispri | 11,106,687 | 579.11 | 503.70 | 462.26 | 0.92 | 3700 | 351.1 | 805.0 |
| multiome | k562_crispri | 11,447,564 | 605.03 | 484.36 | 464.58 | 0.96 | 3718 | 317.6 | 817.1 |
| scatac | thp1_1 | 10,173,077 | 189.04 | 277.57 | 195.94 | **0.71** | 2300 | 224.4 | 378.8 |
| scatac | thp1_2 | 10,278,325 | 191.59 | 279.01 | 195.81 | **0.70** | 2414 | 241.2 | 383.0 |
| scatac | jurkat_pma_cd3_4hr | 10,854,698 | 208.35 | 233.92 | 222.74 | 0.95 | 1736 | 203.1 | 442.6 |
| scatac | jurkat | 11,036,283 | 207.84 | 244.64 | 212.13 | 0.87 | 2590 | **27.5** | 390.7 |
| scatac | telohaec_crispri | 11,106,687 | 211.66 | 294.92 | 229.86 | **0.78** | 2866 | 307.8 | 450.6 |
| scatac | k562_crispri | 11,447,564 | 222.65 | 247.02 | 222.12 | 0.90 | 2328 | 220.5 | 414.5 |

### 5.2 Within a modality: the range is too narrow to fit an exponent

`n_cand_pairs` spans **1.13×** and `feat_bytes` **1.22× / 1.18×**. Per-unit costs:

| | `cpu_time / n_cand_pairs` (µs/pair) | `cpu_time / feat_bytes` (s/MB) |
|---|---|---|
| multiome | 39.1 – 41.7 (CV 2.4 %) | 0.768 – 0.845 (CV 3.4 %) |
| scATAC | 19.1 – 20.7 (CV 3.5 %) | 0.998 – 1.086 (CV 3.2 %) |

Both are flat to within noise, and **both are equally flat.** This is the expected outcome and the
honest statement is: *within a modality, `n_cand_pairs` and `feat_bytes` are near-perfectly
collinear (`feat_bytes ≈ 50 B × n_cand_pairs` multiome, `≈ 19 B ×` scATAC), so these six points
cannot distinguish "linear in rows" from "linear in bytes", and cannot resolve the `log n` factor at
all over a 1.13× span.* The class in §4 is derived from the code; the benchmarks are consistent with
it and confirm nothing beyond linearity.

### 5.3 Across modalities: which variable predicts the ~2× gap?

The prompt asks this directly. Per-cluster `cpu_time` ratios (multiome ÷ scATAC):
2.110, 2.034, 2.092, 2.011, 2.031, 2.180 → **mean 2.076**. Wall ratio: mean **1.848** (noisier,
because scATAC has more I/O stall — see §5.4).

Candidate predictors of that 2.076:

| predictor | multiome/scATAC ratio | verdict |
|---|---:|---|
| `n_cand_pairs` (row count) | **1.000** | **cannot explain any of it.** Row counts are identical. |
| `n_feat_cols` (24/18) | 1.333 | **undershoots by 36 %.** |
| `feat_bytes` (compressed input) | 2.686 | **overshoots by 29 %.** |
| uncompressed input+output volume | 1.635 | undershoots by 21 % |
| `k_feat` (6/6) | 1.000 | irrelevant — the model is the same width |
| **uncompressed volume × measured gzip-9 rate** | **2.27** | **closest; overshoots by 9 %** |

So the answer to *"does `feat_bytes` track the time better than either row count or column count?"*
is: **`feat_bytes` beats row count trivially (row count has zero explanatory power) and beats
column count, but it is not the right variable either — it overshoots by 29 %.** The gap is not a
pure width effect and it is not a pure byte-volume effect. Decomposition:

1. **Row count contributes nothing.** `n_cand_pairs` is bit-identical between the two modalities for
   all six clusters (the branch is decided upstream at `checkpoint features_required`; the candidate
   pair set is shared).
2. **`n_feat_cols` 24→18 is the *mechanism*, but the six absent columns are not average-width
   columns.** They are all float64 score columns (`normalizedATAC_enh`, `RNA_meanLogNorm`,
   `RNA_pseudobulkTPM`, `RNA_percentCellsDetected`, `Kendall`, `ARC.E2G.Score`) printed at full
   repr precision, whereas the 18 shared columns include six low-cardinality *string* columns and
   several small integers. Hence 6/24 of the columns carry 40 % of the uncompressed bytes
   (248.5 → 148.4 B/row). **This is why column count undershoots.**
3. **`feat_bytes` overshoots because it is *compressed* bytes, and the multiome extra columns are
   high-entropy floats that compress far worse than the shared string/integer columns.**
   Compression ratios: multiome 2.53 GB → 495 MB (5.11×); scATAC 1.51 GB → 189 MB (7.99×). So
   `feat_bytes` exaggerates the volume difference by 7.99/5.11 = 1.56×, which is precisely the
   2.686/1.675 discrepancy.
4. **The same entropy difference makes gzip-9 *slower per byte* on multiome** (15.2 vs 21.5 MB/s
   measured, 1.42×). Combined with the 1.605× output-byte difference this gives 2.27× on the
   dominant term — the closest single-mechanism prediction of the observed 2.076.
5. **Plus a modality-gated additive term.** `filter_by_tpm` (§2.4) runs only in multiome:
   ~10.2–11.4 M Python-level loop iterations that scATAC never pays, and one extra output column.

**Bottom line:** neither `n_cand_pairs`, nor `n_feat_cols`, nor `feat_bytes` alone predicts the 2×.
It is *not* a pure width effect. The best two-variable description is a fixed overhead + a per-row
term + a per-`feat_bytes` term, which is what §7 uses.

### 5.4 A two-term fit reproduces all 12 points to within 4.6 % on `cpu_time`

Let `N = n_cand_pairs / 1e6` and `B = feat_bytes / 1e6`. Calibrating on the two modality means
(the only contrast that identifies the split — see §6/§10):

```
cpu_est(N, B) = 5.0 + 6.70·N + 0.661·B        [seconds of CPU]
```

| modality | cluster | predicted | observed `cpu_time` | error |
|---|---|---:|---:|---:|
| multiome | jurkat | 454.4 | 447.60 | +1.5 % |
| multiome | jurkat_pma_cd3_4hr | 446.3 | 453.11 | −1.5 % |
| multiome | k562_crispri | 481.6 | 464.58 | +3.7 % |
| multiome | telohaec_crispri | 462.1 | 462.26 | −0.0 % |
| multiome | thp1_1 | 400.4 | 397.93 | +0.6 % |
| multiome | thp1_2 | 407.7 | 426.85 | −4.5 % |
| scatac | jurkat | 216.3 | 212.13 | +2.0 % |
| scatac | jurkat_pma_cd3_4hr | 215.4 | 222.74 | −3.3 % |
| scatac | k562_crispri | 228.9 | 222.12 | +3.1 % |
| scatac | telohaec_crispri | 219.3 | 229.86 | −4.6 % |
| scatac | thp1_1 | 198.1 | 195.94 | +1.1 % |
| scatac | thp1_2 | 200.5 | 195.81 | +2.4 % |

12/12 within 4.6 %. **This is not evidence that the split between the `N` and `B` terms is
correct** — see §10. It is evidence that the job is linear in volume at this scale with a stable
per-unit constant, and that the modality difference is fully captured by `feat_bytes` plus a
constant.

The **wall**-time residuals are worse (±13 %) and they are *anti-correlated between modalities*:
`telohaec_crispri`, `thp1_1` and `thp1_2` are slow in scATAC but on-model in multiome. Those are
exactly the three scATAC jobs with CPU/wall of 0.70, 0.71 and 0.78. **The wall-time residual tracks
CPU/wall, not size** — it is Lustre/node contention, not algorithm. Phase 2 should be aware that
this rule's wall time carries ±13 % of scheduling noise on top of a very tight compute model.

### 5.5 Memory

`max_rss / (N · n_feat_cols)` = 8.9 – 14.3, mean **12.8 MB per (million rows × column)**
≈ 12.8 bytes per cell — consistent with float64 storage plus partially-sampled copies from
`run_e2g_cv.py:139-140`. Predictive model: `max_rss ≈ 12.8 · N · n_feat_cols` MB. The spread is
wide because peak RSS depends on when psutil samples relative to the two transient full-frame
copies.

`determine_mem_mb` (`ENCODE_rE2G/ABC/workflow/rules/utils.smk:17-24`) requests
`max(4 × 8 × input_MiB, 8000)` capped at `max_memory_allocation_mb` (250,000). Verified in the
slurm log: `thp1_1` multiome got `mem_mb=15106` (= 32 × 472 MiB) against a measured 3350 MB peak —
a **4.5× over-request**. Across the 12 jobs the over-request is **3.2× – 5.2×**. Because both the
request and the peak scale with `n_cand_pairs · n_feat_cols`, that headroom ratio is
scale-invariant: **there is no memory cliff as `n` grows**, which is why §7's `t_high` contains no
cliff term. (The over-request is a scheduling/queueing cost, not a runtime cost.)

### 5.6 I/O counters — `io_in` is not usable as read volume

`io_out` tracks the output file well (`thp1_1` multiome: `io_out` 689.6 MB vs a 742.7 MB file;
scATAC 378.8 vs 414.0 MB). `io_in` does **not**: multiome `io_in` is 306–367 MB against 495–605 MB
inputs, and `scatac/jurkat` reports **27.5 MB against a 208 MB input**. Snakemake's benchmark
`io_in` comes from psutil's block-device counters, so page-cache hits are invisible. Use
`feat_bytes` for read volume, never `io_in`. (This also means `io_in` is contaminated by whatever
else was running on the node.)

---

## 6. Concrete overhead sources

Ordered by measured or code-derived size.

1. **gzip level 9 on the output — ≥52 % (multiome) / ≥47 % (scATAC) of CPU time.**
   `run_e2g_cv.py:165` passes `compression="gzip"` with no `compresslevel`; `gzip.GzipFile` defaults
   to 9. Measured 15.2 MB/s (multiome-shaped) and 21.5 MB/s (scATAC-shaped) of uncompressed input,
   against 194 MB/s for decompression. Single-threaded by construction.
2. **Eight to ten full passes over the `n × n_feat_cols` frame.** Enumerated in §2. Two of them
   (`run_e2g_cv.py:139`, `:140`) each materialise a **complete extra copy** of the frame and each
   also walks the six object (string) columns comparing them to `±inf` / filling NA — provably
   useless work on 25 % (multiome) / 33 % (scATAC) of the columns.
3. **A pure-Python per-row loop, multiome only.** `run_e2g_cv.py:92-98`: a list comprehension over
   `zip` of two pandas Series, `n_cand_pairs` ≈ 10.2–11.4 M iterations, materialising a
   10-million-element Python list. This is the single most expensive *per element* operation in the
   script and is one of the two mechanisms behind the modality gap (§5.3).
4. **No column pruning.** `run_e2g_cv.py:138` reads all 24/18 columns with no `usecols=`; all are
   carried through both copies and re-serialised. Only 6 features + `chr` + `RNA_pseudobulkTPM` are
   used computationally; the other 16/11 are pure pass-through. Because the output is a superset of
   the input, the ~500 MB input is effectively re-encoded from scratch into a ~740 MB output.
5. **`to_csv` float formatting.** Python-level `repr`-precision formatting of `n × 27` (multiome) /
   `n × 20` (scATAC) cells. No `float_format=` is set, so every score is written at full
   repr precision (≈ 17 significant digits for the float columns), which both costs formatting time
   *and* inflates the byte volume that gzip-9 then has to compress. Not separately measurable here,
   but it is the largest single item in the ~190 s / ~104 s residual of §4.
6. **Two `describe()` calls purely for stdout.** `run_e2g_cv.py:72-73`. Each computes count, mean,
   std, min, three quantiles and max → ~2–3 extra full passes over `n` each (the quantiles go
   through `np.partition`). Their entire output is the six lines visible in the slurm log.
7. **`np.log(np.abs(X) + epsilon)` allocates three temporaries** of `n × 6` float64
   (`run_e2g_cv.py:14`, and again at `:30` in the dead CV path) = ~530 MB of transient allocation per
   temporary at `n = 11 M`, with no in-place ops.
8. **Fixed per-job overhead, ~30–45 s.** The slurm log shows Snakemake **building the DAG twice**
   (remote wrapper + inner localrule; 19:02:09 → 19:02:21 = 12 s before the script begins), conda
   env activation off Lustre (`../scE2G/.snakemake/conda/c1ae1ceebb…`), and
   `import pandas/numpy/scipy/sklearn`. Recall this is the *reported* difference between wall
   (mean 485.9 / 262.9 s) and CPU (mean 442.1 / 213.1 s): 43.8 s and 49.8 s respectively.
9. **Declared vs realised threading: no gap to report — there is nothing declared.** The rule
   carries no `threads:`, `profiles/measure/config.yaml` sets `cpus_per_task=1`, and the log confirms
   `Provided cores: 1`. Measured CPU/wall 0.70–0.96 ≤ 1 throughout. The shortfall below 1.0 is
   **Lustre read/write stall**, not lost parallelism — and it is systematically worse in scATAC
   (0.70–0.95, mean 0.82) than multiome (0.86–0.96, mean 0.91), because scATAC does less CPU work per
   byte moved. Everything in the script is single-threaded pandas/zlib; there is no thread flag being
   ignored anywhere.
10. **Duplicated read when both models run on one cluster.** The rule fans out on `{model_name}`
    (`workflow/rules/sc_predictions.smk:65`) while the input depends only on `{cluster}`
    (`:49`). The measured configs assign one model per cluster (6 jobs), but the parent
    `configs/igvf10_config.yaml` assigns both, which would read the same 495–605 MB
    `genomewide_features.tsv.gz` twice and write two ~740 MB outputs. Job count, not per-job cost.
11. **Redundant work inside the dead CV branch**, listed for completeness: `run_e2g_cv.py:29-30`
    recomputes `X` and its log-transform identical to `:13-14`; `:34` builds a full-width DataFrame
    slice per chromosome only to discard everything but `.index.values`.

---

## 7. The three functions

Variables, named exactly as in the manifest:

```
N = n_cand_pairs / 1e6        # millions of candidate enhancer-gene pairs (10.17 - 11.45)
B = feat_bytes  / 1e6         # MB of gzipped genomewide_features.tsv.gz
                              #   multiome 495 - 605, scATAC 189 - 223
```

`n_feat_cols` does **not** appear explicitly: it enters through `B`, which is the modality-carrying
variable (§5.3 shows `n_feat_cols` alone mispredicts the modality ratio by 36 %). **One formula
covers both modalities** — see §9.

All three return **wall seconds** (the critical-path currency). They differ in **class**, not in
point estimate, and they are not blended.

### 7.1 `t_est` — the derived class, Θ(volume)

```
t_est(N, B) = 33.0 + 9.00*N + 0.644*B
```

Θ(`n_cand_pairs` · `n_feat_cols`), i.e. linear in volume. The `n log n` rank
(`run_e2g_cv.py:48`) is folded into the `9.00*N` term because its measured share at `N ≈ 11` is
under 1 %. Reproduces both modality wall-time means to within 0.2 % and all 12 individual jobs to
within 13 % (residual = Lustre contention, §5.4).

### 7.2 `t_low` — same class, fastest constants

```
t_low(N, B) = 20.0 + 7.60*N + 0.548*B
```

This is *the same asymptotic class as* `t_est`, with constants set ~15 % below the calibration.
**That is deliberate and it is the correct lower bound**: the script provably must read every cell
of the input once (`run_e2g_cv.py:138`) and write every cell of the output once (`:165`), so
Θ(`n_cand_pairs` · `n_feat_cols`) is a hard **lower** bound on the class, not merely an estimate.
There is no sub-linear regime available. The only slack is the constant: a warm page cache, an
uncontended node, and the fastest gzip-9 rate observed. `t_low` sits below all 12 observations
(8–13 % below the fastest job in each modality).

### 7.3 `t_high` — the log factor made explicit, Θ(volume · log n)

```
t_high(N, B) = 33.0 + (9.00*N + 0.644*B) * log2(N * 1e6) / 23.366
```

(`log2` = base-2 logarithm; `23.366 = log2(1.08162e7)`, the mean `n_cand_pairs`, so the factor is
exactly 1 at the calibration point and the three functions are anchored consistently.)

This promotes the single `Series.rank` sort at `run_e2g_cv.py:48` and the `np.partition` calls
inside the two `describe()` calls at `:72-73` from "negligible" to "governing", i.e. it asks what
happens if the whole job behaves like the sort rather than like the stream. It is the honest
**upper** bound because **nothing in the script is worse than `n log n`** (§4): there is no nested
pass over `n`, no join, no self-comparison, no `j_cand_elems²`, no `k_genes` cross term. I checked
for these explicitly and they are absent.

`t_high` deliberately contains **no memory-cliff term**: §5.5 shows the `mem_mb` request and the
peak RSS both scale as `n_cand_pairs · n_feat_cols`, so the 3.2–5.2× headroom is scale-invariant.

The resulting band is narrow — 1.36× at 10× the measured scale (§8). **That narrowness is a finding,
not a hedge**: this script's cost is structurally determined and there is no exponent ambiguity to
resolve. Contrast this with a rule containing an O(`j²`) loop, where the band would span orders of
magnitude.

### 7.4 Separate multiplier: `benchmark_performance: True`

**Not covered by the band above, and not measured.** All 12 measured jobs ran with
`benchmark_performance: False` (§2.2). If Phase 2 needs to cost a CRISPR-benchmarking run, the code
at `run_e2g_cv.py:25-43`, `:76-80` and `:100-110` adds:

- Θ(`n_cand_pairs` · `n_chrom`) with `n_chrom = 23` — 23 O(`n`) equality scans on an **object**
  column plus 23 full-width DataFrame materialisations (`run_e2g_cv.py:34`), of which only
  `.index.values` is used
- one `np.unique` over an object column of `n` (`:31`) — a sort with Python-string comparisons
- a duplicate `n × 6` log-transform (`:29-30`)
- a **second** Θ(`n log n`) rank (`:77`)
- 23 scatter `.loc[]` assignments totalling `n` (`:41`)
- multiome only: a **second** `n`-iteration Python list comprehension (`:104-110`)
- output width 27 → 30 (multiome) / 20 → 22 (scATAC), inflating the dominant gzip-9 term ~1.15×

Code-derived estimate, **wide bounds, unmeasured**: `× 1.5 – 1.9` on multiome, `× 1.6 – 2.0` on
scATAC. Flagged as a guess; I have no measurement of this path at any scale, and the object-column
mask cost in particular could be off by 2× in either direction.

---

## 8. Evaluated

Instructions ask for mean, min (`thp1_1`) and max (`telohaec_crispri`). **Note:
`telohaec_crispri` is NOT the maximum for this rule** — it is the max in `n_frag`, but this rule's
variables peak at `k562_crispri` (`n_cand_pairs` 11,447,564; `feat_bytes` 605.0 / 222.6 MB).
`k562_crispri` is included as the true max.

### Multiome

| point | `N` | `B` | `t_low` | `t_est` | `t_high` | observed wall |
|---|---:|---:|---:|---:|---:|---:|
| mean of 6 | 10.8162 | 551.68 | 404.5 s (6.74 min) | **485.6 s (8.09 min)** | 485.6 s (8.09 min) | 485.9 s (8.10 min) |
| `thp1_1` (min) | 10.1731 | 494.99 | 368.6 s (6.14 min) | **443.3 s (7.39 min)** | 441.8 s (7.36 min) | 464.8 s (7.75 min) |
| `telohaec_crispri` | 11.1067 | 579.11 | 421.8 s (7.03 min) | **505.9 s (8.43 min)** | 506.7 s (8.44 min) | 503.7 s (8.40 min) |
| `k562_crispri` (true max) | 11.4476 | 605.03 | 437.7 s (7.29 min) | **525.7 s (8.76 min)** | 527.4 s (8.79 min) | 484.4 s (8.07 min) |

### scATAC

| point | `N` | `B` | `t_low` | `t_est` | `t_high` | observed wall |
|---|---:|---:|---:|---:|---:|---:|
| mean of 6 | 10.8162 | 205.19 | 214.6 s (3.58 min) | **262.5 s (4.37 min)** | 262.5 s (4.37 min) | 262.8 s (4.38 min) |
| `thp1_1` (min) | 10.1731 | 189.04 | 200.9 s (3.35 min) | **246.3 s (4.11 min)** | 245.5 s (4.09 min) | 277.6 s (4.63 min) |
| `telohaec_crispri` | 11.1067 | 211.66 | 220.4 s (3.67 min) | **269.3 s (4.49 min)** | 269.7 s (4.50 min) | 294.9 s (4.92 min) |
| `k562_crispri` (true max) | 11.4476 | 222.65 | 228.1 s (3.80 min) | **279.4 s (4.66 min)** | 280.3 s (4.67 min) | 247.0 s (4.12 min) |

At the measured scale `t_est` and `t_high` are indistinguishable (`log2(n)/23.366 ≈ 1.00 ± 0.004`
over `n_cand_pairs` 10.17–11.45 M) — which is exactly the point of §5.2. They separate only at
scales the data does not reach. **At 10× the measured `n`** (multiome, `N = 108.16`, `B = 5517`):
`t_low` 3865 s (64 min), `t_est` 4559 s (76 min), `t_high` 5203 s (87 min) — a 1.36× band.

---

## 9. Both modalities: class or constant?

**The class is identical. Only the constants differ — but not by a single multiplier, so a
"constant" is the wrong mental model.**

Same in both:
- `n_cand_pairs` — **bit-identical** for all six clusters (the modality branch is decided upstream at
  `checkpoint features_required`; both DAGs share the same candidate pair set)
- `k_feat = 6` — `ARC.E2G.Score` **substitutes for** `ABC.Score` in the feature table, it is not
  added (§2.1). `X` is `n × 6` in both. The model is a 956-byte `LogisticRegression` in both.
- `n_chrom = 23` fold pickles present in both model dirs (both unused)
- the single global rank at `run_e2g_cv.py:48`, `k = 1` sort of `n`
- the 6 object columns walked pointlessly by `:139-140`
- gzip level 9 on the output

Different in scATAC:
- `n_feat_cols` 18 vs 24 → **a genuine cross term**, not a scalar. It changes the frame width, the
  copy volume at `:139-140`, the parse volume at `:138`, and the output width at `:165`.
- Uncompressed bytes/row 148.4 vs 248.5 in, 193.0 vs 309.8 out — the 6 missing columns are all
  high-entropy float64, so they carry ~40 % of the bytes despite being 25 % of the columns.
- Compressibility: 7.99× vs 5.11×, hence `feat_bytes` ratio 2.686 while uncompressed volume ratio is
  only 1.635×. And gzip-9 is 1.42× **slower per byte** on the multiome data.
- **`filter_by_tpm` (`run_e2g_cv.py:85-112`) does not run at all in scATAC** — it returns at `:86`
  because `RNA_pseudobulkTPM` is absent *and* `tpm_threshold` is 0 (`models/scATAC_powerlaw_v3/
  tpm_threshold_0`). Multiome pays a ~10–11 M-iteration pure-Python loop that scATAC skips.
  Strictly this is O(`n`) vs O(1) — a term that exists in one modality and not the other — but since
  it is dominated by the Θ(volume) term, the **overall class is unchanged**.
- Output width 20 vs 27 columns.

**Is the 2× a pure width effect?** No. Width (`n_feat_cols`, 1.333×) explains under half of it.
The measured `cpu_time` ratio is 2.076; the closest single mechanism is uncompressed output volume ×
measured gzip-9 rate (1.605 × 1.42 = 2.27), with the multiome-only `filter_by_tpm` loop and the
identical fixed overhead pulling in opposite directions. `feat_bytes` (2.686×) tracks the time
better than either row count (1.000×, zero explanatory power) or column count (1.333×), but it
overshoots by 29 % because it is compressed bytes and the modalities compress differently.

The single formula in §7 handles both modalities correctly because `feat_bytes` is measured per
modality; there is no separate multiome/scATAC calibration.

---

## 10. What the measurement cannot tell us

1. **It cannot separate `n_cand_pairs` from `feat_bytes` within a modality.** They are collinear at
   `feat_bytes ≈ 50 B × n_cand_pairs` (multiome) and `≈ 19 B ×` (scATAC), with a coefficient of
   variation under 3.5 % in each. The `6.70*N` / `0.661*B` split in §5.4 rests **entirely on the
   two-modality contrast** — effectively two data points for two parameters, plus an assumed 5 s
   interpreter-startup constant. A different assumed startup constant moves the split without
   changing any of the 12 predictions. **Do not read the split as physically identified.** The
   *sum* is well-determined; the *apportionment* is not.
2. **It cannot confirm the `log n` factor.** `n_cand_pairs` spans 1.13×, so `log2(n)` spans
   23.278 → 23.449, i.e. **1.007×**. The rank at `run_e2g_cv.py:48` is derived from the code
   (`method="average"` forces an argsort in `algos.rank_1d`); the benchmarks cannot see it and never
   will at this range.
3. **It cannot confirm linearity in `n_feat_cols`.** `n_feat_cols` takes exactly **two** values,
   24 and 18. Two points determine a line by construction. The claim "linear in `n_feat_cols`" comes
   from the code (`read_csv`/`replace`/`fillna`/`to_csv` are all elementwise over the frame), not
   from the data.
4. **It has measured nothing at all about the CV path.** `benchmark_performance: False` in all four
   configs, so `run_e2g_cv.py:25-43`, `:76-80` and `:100-110` — roughly 40 % of the script's lines,
   including the only place `n_chrom = 23` appears — contributed **zero** seconds to every one of the
   12 jobs. §7.4 is code-reading with wide bounds, not measurement.
5. **The gzip-9 attribution (~52 % / ~47 %) is an extrapolation from a 200,000-row slice, using
   system `gzip`, run on a contended login node.** I used user CPU time to remove contention, and
   Python's `gzip` module uses the same zlib level 9, but Python writes through a `TextIOWrapper`,
   so the real in-process figure is ≥ mine. I could not profile the actual job (no code edits
   permitted, and re-running would consume cluster time). Treat 52 %/47 % as **lower bounds** on
   gzip's share, with maybe ±10 percentage points of uncertainty.
6. **I cannot split the ~190 s (multiome) / ~104 s (scATAC) non-gzip residual.** `read_csv`
   tokenising, the two full-frame copies, `to_csv` float formatting, the rank, the two
   `describe()` calls and the multiome `filter_by_tpm` loop all live inside it. My budget estimates
   for each are consistent with the total but are **not measurements**. In particular I cannot say
   whether `filter_by_tpm` costs 10 s or 80 s.
7. **The wall-time model carries ±13 % of irreducible scheduling noise.** These are single
   measurements, not replicates. The three worst wall residuals (`scatac/thp1_1`, `scatac/thp1_2`,
   `scatac/telohaec_crispri`) are exactly the three jobs with CPU/wall ≤ 0.78 — Lustre contention,
   not size. The `cpu_time` model is far tighter (≤ 4.6 %) and is what Phase 2 should use if it can
   model queueing separately.
8. **`io_in` is unusable as read volume** (§5.6). `scatac/jurkat` reports 27.5 MB against a 208 MB
   input — a page-cache hit. Any conclusion drawn from `io_in` for this rule is wrong.
9. **`max_rss` is sampled, and this script has two transient full-frame copies.** The 8.9–14.3
   MB/(Mrow·col) spread is psutil sampling luck relative to `run_e2g_cv.py:139-140`, not real
   variation. `multiome/thp1_2` reports 2772 MB where `multiome/jurkat` (only 1.07× more rows)
   reports 3796 MB.
10. **Nothing here says whether this rule is on the critical path or whether it matters.** Its
    output is consumed by `filter_sc_e2g_predictions` (`workflow/rules/sc_predictions.smk:79-81`),
    `element_and_gene_summaries` (`:140`), `overlap_features_crispr_apply` (`:5`) and the rule at
    `:166`. That is Phase 2's arithmetic, not mine.

---

### Appendix — sources

- Benchmarks: `/scratch/users/kaybrand/scE2G_optimize_results/analysis/benchmarks_{multiome,scatac}.tsv`,
  rows with `rule == "run_e2g_qnorm"` (6 + 6)
- Sizes: `.../analysis/size_manifest_{multiome,scatac}.tsv`
- Slurm log quoted for `mem_mb`, DAG double-build and the two `describe()` blocks:
  `/scratch/users/kaybrand/scE2G_optimize_results/igvf10_multiome/slurm_logs/rule_run_e2g_qnorm/thp1_1_multiome_powerlaw_v3/39922797.log`
- On-disk shapes (headers, byte widths, output sizes):
  `/scratch/users/kaybrand/scE2G_optimize_results/igvf10_{multiome,scatac}/thp1_1/`
- gzip/zcat rate measurements: `/usr/bin/time` on 200,000-row and 400,000-row slices of the real
  `scE2G_predictions.tsv.gz` and `genomewide_features.tsv.gz`. **Leftover: `/tmp/kb_qnorm_probe`
  (~101 MB of sampled TSV rows) on the login node — the cleanup `rm` was blocked by the permission
  system; please remove it manually.**
- Toy `chr22` run excluded throughout, per the briefing.

# Unit 16 — `process_model_output.py`

## 1. Unit identity

- Script: `ENCODE_rE2G/workflow/scripts/model_application/process_model_output.py` (36 lines,
  the local submodule copy at `/oak/.../scE2G_optimize/ENCODE_rE2G/`). Verified this is the copy
  actually executed: the rule's `scripts_dir` param resolves via `os.path.join(config["encode_re2g_dir"],
  encode_e2g.SCRIPTS_DIR)`, and `encode_e2g.SCRIPTS_DIR` (`ENCODE_rE2G/workflow/Snakefile:37`) is
  itself already absolute (`os.path.abspath(config["E2G_DIR_PATH"])`-derived), so Python's
  `os.path.join` discards the first argument and the effective path is exactly
  `encode_e2g.SCRIPTS_DIR`. The igvf10 configs (`configs/igvf10_multiome_config.yaml:81`,
  `configs/igvf10_scatac_config.yaml:84`) point `encode_re2g_dir` at this checkout's
  `ENCODE_rE2G/`, not the main `scE2G` checkout — consistent with the known
  `sce2g-encode-re2g-dir-gotcha`, but for these two config files the "gotcha" happens to resolve
  to the local copy anyway.
- Rule: `write_sc_e2g_predictions_bedpe`, `workflow/rules/sc_predictions.smk:105-125` — defined
  directly in the outer `scE2G_optimize` Snakefile stack, not inside the `encode_e2g` module.
- Modalities: **both**, 6 jobs each, 12 measured jobs total. No `threads:` declared (runs at 1).
  `resources: mem_mb = encode_e2g.ABC.determine_mem_mb` (`sc_predictions.smk:117-118`).

## 2. What the code actually does

The entire script (36 lines) is reproduced in control-flow order:

1. `:31` — `pred_thresh = pd.read_csv(predictions_file, sep="\t")`. `predictions_file` is
   `scE2G_predictions_threshold{threshold}.tsv.gz` (`sc_predictions.smk:107`), the thresholded
   predictions file whose row count the manifest calls `n_called_all`. Pandas auto-detects the
   `.gz` extension and decompresses inline; there is no separate decompression pass. This is a
   single-pass parse: O(`n_called_all` · `n_cols`), where `n_cols` is the file's column count.
   Verified on disk (not in the manifest, which only tracks `n_feat_cols` for a *different* file,
   `genomewide_features.tsv.gz`): the multiome threshold file has **27 columns**, the scATAC file
   has **20 columns** (`zcat … | head -1 | awk -F'\t' '{print NF}'` on `k562_crispri`'s file in
   both modalities).
2. `:8` — `pred = pred.drop_duplicates()`, called with no `subset=`, so it hashes every column of
   every row. Pandas implements `duplicated()`/`drop_duplicates()` via hash-based grouping
   (`factorize`/hashtable), average-case **O(n)**, not a sort — this is standard, documented
   pandas behavior, not something visible in this repo's code, so I flag it as inferred from
   pandas' implementation rather than read directly. Empirically, on every file I checked
   (`k562_crispri` multiome, `thp1_2` multiome [max `n_called_all`], `k562_crispri` scATAC [max
   `n_called_all`]), **output row count exactly equals input row count minus header** — zero
   duplicates are ever actually dropped, so this call is a full-cost no-op in every observed job.
3. `:10-21` — eight column assignments building `towrite`: five are direct vectorized numeric/str
   copies (`chr1`, `x1`, `x2`, `chr2`, `y1`, `y2`, `score`), two are constant fills (`strand1`,
   `strand2`), and one (`:18`, `towrite["name"] = pred["TargetGene"] + "_" + pred["name"]`) is an
   elementwise Python-object string concatenation — still O(n), single pass, but with a larger
   per-row constant than the numeric copies since it touches Python `str` objects rather than a
   contiguous numeric buffer.
4. `:23` — `towrite.to_csv(outfile, header=False, index=False, sep="\t")`: a single-pass,
   uncompressed plain-text write, O(`n_called_all`).

**No `.sort_values` call anywhere in the file**, and no `subprocess`/`gzip`/`bgzip`/`tabix`
imports or calls. The rule's shell block (`sc_predictions.smk:119-125`) invokes only this one
Python script — no post-processing shell command follows it. I confirmed on disk that the output
`.bedpe` files *are* in genomic block order (grouped by `chr1`, ascending `start` within each
chromosome, verified on `k562_crispri` multiome) — but that ordering is **inherited from the
upstream thresholded-predictions file**, not produced by any sort in this script. This directly
answers the unit brief's question: **no sort, and no bgzip/tabix indexing** — a single O(n) pass
end to end.

## 3. Size variable(s)

- **Primary: `n_called_all`** — Tier 3, near-constant within a modality. Multiome: 59,438–61,284
  (1.031×). scATAC: 105,880–120,113 (1.134×). Verified directly: `zcat` row count minus header on
  `k562_crispri`'s multiome threshold file = 59,438, matching the manifest exactly.
- **`n_called_bytes`** (gzip-compressed size of `predictions_file`, governs the read/decompress
  cost): 4,739,099–4,886,331 bytes multiome (4.74–4.89 MB), 4,864,744–5,471,944 bytes scATAC
  (4.86–5.47 MB) — matches the brief's stated ranges.
- **Column count** (not a named manifest column, observed directly): 27 (multiome) vs. 20
  (scATAC) for `predictions_file`. This is a per-modality *constant* multiplier on the per-row
  parse/hash cost, not a second independent variable and not the same thing as the manifest's
  `n_feat_cols` (24/18), which describes `genomewide_features.tsv.gz`, a file this script never
  opens.
- No cross term: the script touches no other manifest variable (`j_cand_elems`, `k_genes`,
  `n_cand_pairs`, etc. do not appear). This is a genuinely single-variable, single-pass script.

## 4. Asymptotic class derived from code

**O(`n_called_all`)** — every operation in the script (`:31` read+decompress+parse, `:8`
hash-based dedup, `:12-21` eight vectorized/elementwise column builds, `:23` write) is a single
linear pass over the rows. There is no sort, no join, no nested loop, no groupby. This is, by
design, the simplest class encountered in this project's rule set — confirmed by reading the
entire 36-line file rather than assuming a sort/dedup step implies more than O(n). The only
per-modality constant is the column count (27 vs. 20), which scales the per-row cost of step 1
but not the exponent.

## 5. Benchmark cross-check

All 12 measured jobs (`s` = wall seconds, `cpu_time` = seconds):

| Modality | Cluster | `n_called_all` | `s` | `cpu_time` | cpu/s | `io_in` MB | `io_out` MB |
|---|---|---:|---:|---:|---:|---:|---:|
| multiome | jurkat_pma_cd3_4hr | 59,650 | 1.19 | 0.66 | 0.55 | 4.53 | 4.00 |
| multiome | jurkat | 59,698 | 1.34 | 0.55 | 0.41 | 4.63 | 0.00 |
| multiome | k562_crispri | 59,438 | 2.57 | 0.75 | 0.29 | 7.87 | 6.21 |
| multiome | telohaec_crispri | 60,630 | 2.75 | 0.70 | 0.25 | 11.85 | 0.00 |
| multiome | thp1_1 | 60,435 | 1.28 | 0.60 | 0.47 | 4.70 | 0.00 |
| multiome | thp1_2 | 61,284 | 5.44 | 1.86 | 0.34 | 26.66 | 15.06 |
| scATAC | jurkat_pma_cd3_4hr | 114,015 | 1.58 | 1.00 | 0.63 | 5.07 | 11.91 |
| scATAC | jurkat | 116,138 | 1.40 | 0.66 | 0.47 | 5.05 | 0.00 |
| scATAC | k562_crispri | 120,113 | 1.55 | 1.04 | 0.67 | 5.23 | 12.56 |
| scATAC | telohaec_crispri | 117,185 | 3.24 | 2.15 | 0.66 | 5.21 | 23.06 |
| scATAC | thp1_1 | 105,880 | 1.52 | 0.53 | 0.35 | 4.80 | 0.00 |
| scATAC | thp1_2 | 107,248 | 1.59 | 0.95 | 0.60 | 4.70 | 11.22 |

**Within a modality, `n_called_all` cannot be correlated with anything — range is 1.03×
(multiome) / 1.13× (scATAC), squarely Tier 3.** Confirming this rather than assuming it: the
multiome cluster with the *smallest* `n_called_all` (`k562_crispri`, 59,438) has the *second
highest* wall time (2.57 s), while `jurkat_pma_cd3_4hr` (59,650, near the low end) has the
lowest wall time (1.19 s) — no monotonic pattern.

**Cross-modality (the only real ~2× lever available) gives a weak but sign-consistent signal in
`cpu_time`, and none at all in wall time `s`.** Mean `n_called_all` is 60,189 (multiome) vs.
113,430 (scATAC) — a 1.88× increase. Mean `cpu_time` rises from 0.853 s to 1.055 s (+24%) —
positive, but far sub-proportional to the 88% rise in `n`. Mean wall time `s` actually *falls*,
2.428 s (multiome) → 1.813 s (scATAC), the wrong sign for an n-driven model — this is an artifact
of two multiome outliers (`k562_crispri` 2.57 s despite the *smallest* multiome `n`;
`thp1_2` 5.44 s, discussed below), not evidence against O(n).

A pooled OLS of `cpu_time` on `n_called_all` across all 12 jobs gives slope ≈ +4.73×10⁻⁶ s/row,
intercept ≈ 0.543 s, but **R² ≈ 0.065** — i.e. `n_called_all` explains only ~6.5% of the variance
in `cpu_time`. The sign agrees with the code-derived O(n) class; the fit itself is not
statistically meaningful and should not be read as confirmation.

**The cleanest actual O(n) confirmation is in output *bytes*, not time.** The `.bedpe` output
file size scales almost perfectly linearly with `n_called_all`, at essentially the same
per-row rate in *both* modalities:

| Modality | bytes/row (min–max across 6 clusters) |
|---|---|
| multiome | 109.32–109.49 |
| scATAC | 109.30–109.62 |

This is expected — the code performs one `to_csv` write per row with a fixed 10-column schema —
but it is a genuinely clean linear signal, unlike the noisy timing columns, and it is
modality-independent because the output schema doesn't depend on the input's column count.

**`cpu_time`/`s` ratios are 0.25–0.67**, meaning 33–75% of wall time is not measured as active
CPU, at total wall times of only 1.2–5.4 s. This is squarely the "sub-second/small-job" regime the
briefing warns about: startup and wait dominate, and the timing signal for the O(n) algorithmic
content is close to invisible at this scale.

**`thp1_2` (multiome) is a flagged outlier, not a fit point**: `s`=5.44, `cpu_time`=1.86 (both the
highest of the 12 jobs), `n_called_all`=61,284 (the largest in multiome, so the *direction* is
consistent), but `io_in`=26.66 MB against a 4.89 MB actual input file — a ~5.5× read
amplification that has no algorithmic explanation in this script. This is far more likely
filesystem/cold-cache or conda-environment shared-library-loading noise than a `n`-driven cost,
and I have not folded it into the calibration below (see §6, §10).

## 6. Concrete overhead sources

- **Heavy shared conda env for a 3-import script.** `workflow/envs/sc_e2g.yml` (the env this rule
  uses, `sc_predictions.smk:115-116`) bundles R (Seurat, Signac, tidyverse, rmarkdown, doParallel,
  …), scikit-learn, matplotlib, seaborn, anndata, `fast_kendall_sc`, etc. — dozens of packages —
  even though `process_model_output.py` imports only `click`, `pandas`, `os`. Environment
  activation and Python's site-packages path resolution pay a cost proportional to the whole env,
  not to what this script actually uses, which plausibly explains most of the intercept in §5.
- **`io_out` under-reports the real write volume.** Several jobs report `io_out`=0.00 MB
  (`jurkat` multiome, `telohaec_crispri` multiome, `thp1_1` both modalities, `jurkat` scATAC)
  despite writing a 6.5–13.2 MB `.bedpe` file (verified with `ls -la` on-disk: 6,508,288–6,709,835
  bytes multiome, 11,604,530–13,163,564 bytes scATAC). This is the same psutil-buffered-write
  undercount pattern noted for other rules in this project — use the on-disk output size, not
  `io_out`, for this rule's write volume.
- **`io_in`=26.66 MB for `thp1_2`** against a 4.89 MB actual input file (§5) — a concrete,
  unexplained-by-this-script I/O cost, most plausibly conda/shared-library loading or Lustre
  cold-cache reads rather than anything scaling with `n_called_all`.
- **No compression on write.** The 4.7–5.5 MB gzip input expands to a 6.5–13.2 MB plain-text
  `.bedpe` output — purely from dropping compression, independent of the 27/20 → 10 column
  reduction (fewer columns, but no gzip, nets a larger file).
- **Full double-precision float serialization.** Observed score values in the output are written
  at full precision (e.g. `0.1812275828053579`, 17 significant digits) — `to_csv` at `:23` passes
  no `float_format`, so pandas serializes the score column via its default (near-full-precision)
  float-to-string path. This doesn't change the O(n) class but is a real, verifiable contributor
  to per-row output bytes and write time (each score field averages far more characters than a
  rounded value would need).
- **No threading**: no `threads:` declared on the rule (`sc_predictions.smk:105-118`); consistent
  with the project-wide finding that ABC/sc_e2g rules run at 1 regardless of declared resources.
- **No redundant passes**: input is read exactly once (`:31`), output written exactly once
  (`:23`); no intermediate temp files or re-reads.
- **Zero duplicates ever actually removed** (§2 item 2) — the `drop_duplicates()` call at `:8`
  pays its full O(n) hashing cost on every job but never changes the row count in any file I
  checked; it is real, verified-necessary-in-principle defensive code, not dead weight in terms
  of correctness, but it is 100% overhead in terms of observed effect on this data.

## 7. Three functions

The class itself is **not ambiguous** here — unlike Tier-3-driven rules elsewhere in this
project where `t_low`/`t_high` bracket competing *exponents* (e.g. O(n) vs. O(n log n) vs.
O(n²)), every operation in this 36-line script is demonstrably a single linear pass (§2, §4),
with no sort, join, or nested loop to create genuine class uncertainty. `t_low` and `t_high`
therefore bracket the **per-row constant / fixed-overhead uncertainty** — a legitimate unknown,
given the R²≈0.065 pooled fit in §5 — not the exponent.

Basis: wall-clock seconds (`s`), because that is the actual pipeline cost. Calibration route:
fit `cpu_time` (the better-behaved signal) against `n_called_all` using the two modality means as
anchors — mean multiome (n=60,189, cpu_time=0.853 s) and mean scATAC (n=113,430, cpu_time=1.055 s)
— giving slope 3.79×10⁻⁶ s/row and intercept 0.625 s; then add the median non-CPU wall overhead
across all 12 jobs (median of `s − cpu_time` = 0.765 s) as a roughly n-independent addend, since
§5 shows no scaling signal in that residual either.

```
t_est(n_called_all)  = 1.4  + 3.8e-6 * n_called_all   [seconds, wall-clock]
t_low(n_called_all)  = 0.9  + 1.0e-6 * n_called_all   [seconds]
t_high(n_called_all) = 1.3  + 2.5e-5 * n_called_all   [seconds]
```

- `t_low`'s slope (1.0e-6, ~26% of `t_est`'s) represents the possibility that the per-row cost is
  much smaller than the anchored estimate — the within-modality data (1.03×–1.13× range) cannot
  rule this out. It sits at or below every one of the 12 observed wall times (e.g. at
  `n`=59,438: 0.96 s, below the smallest observed multiome wall of 1.19 s).
- `t_high`'s slope (2.5e-5, ~6.6× `t_est`'s) represents a pessimistic per-row constant — e.g. if
  the elementwise string concatenation at `:18` and the full-row hash in `drop_duplicates` at `:8`
  dominate rather than being amortized the way pandas' vectorized numeric ops are. It covers 11 of
  the 12 observed jobs (e.g. at `n`=117,185: 4.23 s vs. observed 3.24 s). **It deliberately does
  not cover the `thp1_2` outlier** (observed 5.44 s vs. `t_high`=2.83 s at n=61,284) — per §5, that
  job's `io_in` anomaly (5.5× its own input file size) is flagged as exogenous noise, not
  something the class bound should be stretched to absorb.

## 8. Evaluated at mean / min / max

| Point | Modality | `n_called_all` | `t_low` | `t_est` | `t_high` |
|---|---|---:|---:|---:|---:|
| mean (6 clusters) | multiome | 60,189 | 0.96 s | 1.63 s | 2.80 s |
| `thp1_1` | multiome | 60,435 | 0.96 s | 1.63 s | 2.81 s |
| `telohaec_crispri` | multiome | 60,630 | 0.96 s | 1.63 s | 2.82 s |
| mean (6 clusters) | scATAC | 113,430 | 1.01 s | 1.83 s | 4.14 s |
| `thp1_1` | scATAC | 105,880 | 1.01 s | 1.80 s | 3.95 s |
| `telohaec_crispri` | scATAC | 117,185 | 1.02 s | 1.85 s | 4.23 s |

(As §5 notes, `thp1_1` and `telohaec_crispri` are *not* the min/max clusters for `n_called_all`
itself within either modality — those are `k562_crispri`/`thp1_2` for multiome and
`thp1_1`/`k562_crispri` for scATAC — but `thp1_1`/`telohaec_crispri` are reported here for
consistency with the other 25 units, per the shared instructions. All three functions barely move
within a modality because `n_called_all` is genuinely near-constant there.)

## 9. Does the class differ between modalities?

**No — same class, and the dominant per-row constant (output writing) is empirically the same
too.** The output `.bedpe` bytes/row is 109.3–109.6 across *all 12* jobs regardless of modality
(§5) — expected, since the 10-column output schema doesn't depend on the input's column count.
The input-side constant differs modestly (27 vs. 20 input columns, `n_called_bytes` 4.7–4.9 MB
multiome vs. 4.9–5.5 MB scATAC despite scATAC's `n_called_all` being ~2× larger — i.e. a smaller
per-row compressed footprint on the scATAC side, consistent with fewer columns), but this
difference is swamped by fixed process/environment startup at these 1–5 s wall-clock scales, so
the timing data cannot cleanly separate a modality-dependent constant from noise. What genuinely
differs ~2× between modalities, per the addendum, is `n_called_all` itself (the scATAC model calls
roughly twice as many rows from the same upstream candidate set) — not the constant and not the
class.

## 10. What the measurement cannot tell us

- **Cannot confirm O(`n_called_all`) from within-modality data alone.** The variable spans only
  1.031× (multiome) / 1.134× (scATAC) — Tier 3, exactly the regime the briefing says six (here,
  twelve) points cannot fit. The class in §4 comes from reading all 36 lines of the script, not
  from these numbers.
- **The weak 2-point cross-modality comparison (cpu_time +24% against n +88%) is compatible with
  several stories** — a genuinely small per-row constant, fixed startup dominating both regimes,
  or some mix — and 2 aggregate points cannot distinguish them. The pooled 12-point regression's
  R²≈0.065 says the same thing from a different angle: essentially no explanatory power.
- **Cannot apportion wall time between Python/pandas/click import, conda-environment shared-
  library loading, actual I/O wait, and Slurm/filesystem scheduling noise.** The benchmark table
  gives wall time, `cpu_time`, RSS, and coarse I/O counters — not an import-time or syscall trace.
  Isolating the ~0.5–2.2 s of non-algorithmic overhead precisely would need direct profiling
  (e.g. `python -X importtime`, `py-spy`), which is out of scope for characterization from
  existing artifacts.
- **Cannot fully trust `io_out` for this rule** — it reports 0.00 MB for several jobs that
  demonstrably wrote 6.5–13.2 MB of output (verified via `ls -la`); use on-disk output size
  instead (§6).
- **Cannot explain the `thp1_2` `io_in` anomaly (26.66 MB against a 4.89 MB input file) from the
  benchmark table alone** — flagged as likely filesystem/environment noise, not confirmed.
- **Single (n=1) measurement per cluster per modality** — no replicate exists to separate "this
  job hit a slow node/cold cache" from "this cluster is systematically different," which matters
  here more than usual given how small and noise-dominated these jobs are.
- **Cannot say anything about behavior at `n_called_all` values well outside the observed
  ~59k–120k range** (e.g. a much less stringent threshold producing millions of called rows). The
  O(n) class in §4 is a code-derived asymptotic statement, not validated at larger scale — though
  for this particular script (no sort, no join, no nested loop of any kind), there is little
  reason to expect the class itself to change with scale, only the relative weight of the fixed
  startup overhead vs. the linear term (which would shrink as `n` grows).

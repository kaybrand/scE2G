# Unit 23 — `make_biosample_feature_table`

## 1. Unit identity

- Script: `ENCODE_rE2G/workflow/scripts/feature_tables/combine_feature_tables_apply.R` (33 lines).
- Rule: `make_biosample_feature_table`, `ENCODE_rE2G/workflow/rules/predictions.smk:2-18`.
- Modality: **both** (multiome and scATAC). Benchmark wildcard: `biosample` (= cluster). 6 jobs
  per modality, 12 total.
- Declared `resources: mem_mb=4*1000`; no `threads:` (runs at 1, consistent with every other ABC
  rule in this project).

## 2. What the code actually does

`combine_feature_tables_apply.R:6-33`:

1. `model_dirs` (line 6) comes from `snakemake@params$model_dirs`, which is
   `BIOSAMPLE_DF.loc[BIOSAMPLE_DF['biosample']==wildcards.biosample]['model_dir'].to_list()`
   (`predictions.smk:8`) — i.e. the list of trained-model directories assigned to this biosample,
   split on `" "`.
2. Lines 8-11: for each model dir, `fread(model_dir/feature_table.tsv)`; the first iteration
   assigns `df`, every subsequent one does `df = rbind(df, ft)`.
3. Lines 14-30: **sc-E2G-specific patch.** If any row already has `feature == "ARC.E2G.Score"`
   or `feature == "Kendall"`, append a fixed 7-row block (`ARC_rows`) describing the ARC-derived
   features (`RNA_meanLogNorm`, `RNA_pseudobulkTPM`, `RNA_percentCellsDetected`, `Kendall`,
   `ARC.E2G.Score`, `ABC.Score`, `normalizedATAC_enh`), dropping the `ABC.Score` row from that
   block if `df` already has one (avoids a duplicate ABC.Score/powerlaw.Score row).
4. Line 32: `dplyr::distinct(df)` — a single dedup pass over the whole (tiny) table.
5. Line 33: `fwrite(df, snakemake@output$biosample_features, sep="\t")`.

**The rule's declared `input:` is `config["ABC_BIOSAMPLES"]` (`predictions.smk:4`), but the
script never reads it** — `grep snakemake@` in the script shows only `snakemake@params$model_dirs`
and `snakemake@output$biosample_features` (verified: no `snakemake@input` reference anywhere).
`config["ABC_BIOSAMPLES"]` is `tmp/config_abc_biosamples.tsv`, the file rewritten at every
Snakemake parse by `make_biosample_config()` (top-level `workflow/Snakefile:49`, function at
`workflow/rules/utils.smk:20`) — the landmine flagged in my briefing. Declaring it as `input:`
without reading it means this rule is marked out-of-date on every invocation regardless of
`rerun-triggers`. **This does bear on the measured cost**: every one of the 12 slurm logs for
this rule shows `reason: Forced execution` and `jobid: 0`, confirming each measurement is an
isolated, single-rule Snakemake invocation rather than a job embedded in a longer batched DAG run
— relevant context for interpreting the "constant" below as a fresh cold start each time, not
something amortized across a run. (No fix proposed; this is the landmine's owner's territory, not
mine.)

## 3. Size variable(s)

**None of the `size_manifest_*.tsv` variables govern this rule.** Its only real input is the
*model's* `feature_table.tsv` (a feature-definition table shipped with the trained model, e.g.
`models/multiome_powerlaw_v3/feature_table.tsv`), not any per-cluster data product. Confirmed on
disk:

| | multiome (`models/multiome_powerlaw_v3/feature_table.tsv`) | scATAC (`models/scATAC_powerlaw_v3/feature_table.tsv`) |
|---|---:|---:|
| rows | 6 | 6 |
| bytes | 461 | 453 |

identical for **all six clusters** in a modality, because `config_biosamples.tsv` assigns every
cluster the same `model_dir` (verified in
`/scratch/.../igvf10_{multiome,scatac}/config/expanded_biosample_config.tsv`: `model_dir` column
is `.../models/multiome_powerlaw_v3` for all 6 multiome rows, `.../models/scATAC_powerlaw_v3` for
all 6 scATAC rows — no commas, i.e. one model dir per biosample). Output `feature_table.tsv` is
likewise **byte-identical across clusters within a modality**: 818 bytes (multiome, 13 lines) vs.
442 bytes (scATAC, 7 lines), for all six clusters each.

The one variable the code genuinely has is **`M`, the number of model directories assigned to a
biosample** (`length(model_dirs)`, the loop bound at line 8) — **not a column in
`size_manifest_*.tsv`**. In this dataset `M = 1` for all 12 jobs (both modalities), so it is
Tier-3-like (worse than Tier 3: it doesn't even vary by 1.01×, it is exactly constant). I name it
`M` explicitly per the instructions' guidance not to force a manifest variable where the code has
none.

## 4. Asymptotic class derived from code

Let `M` = number of model dirs for a biosample, `r` = rows in one model's `feature_table.tsv`
(constant 6 here), `k` = size of the fixed ARC block (constant 7 rows, only added conditionally).

- Loop `combine_feature_tables_apply.R:8-11`: `M` calls to `fread` on an `r`-row table → `O(M·r)`
  I/O + parse, trivially small since `r` is a small constant.
- **The accumulation pattern is the classic R `rbind`-in-a-loop antipattern**: `df = rbind(df, ft)`
  (line 10) reallocates and copies the entire accumulated frame on every iteration after the
  first. If R's `rbind.data.frame` does not pre-allocate, the total cost of the loop is
  `Σ_{i=1}^{M} O(i·r) = O(M²·r)` — quadratic in `M`. If it were implemented as a single
  `rbindlist(list_of_tables)` after the loop (it is not — it is done inside), the same operation
  would be `O(M·r)`, linear. **The code as written is the quadratic pattern**, so `O(M²·r)` is
  the honest upper bound on that step; `O(M·r)` is the honest lower bound (in case the current R
  data.table/dplyr environment's `rbind` amortizes better than naive copy — I did not verify this
  empirically since it would require an artificial multi-model config to exercise `M > 1`).
- ARC-row logic (lines 14-30): membership test over `df$feature` (`O(M·r)` elements, tiny),
  possible one more `rbind` of a constant `k=7`-row block → `O(k)`, independent of `M`.
- `dplyr::distinct(df)` (line 32): one hash-based dedup pass over the final table, `O(M·r + k)`.
- `fwrite` (line 33): `O(M·r + k)`.

**Overall class: `O(M·r)` best case / `O(M²·r)` worst case, both dominated by an additive
constant `k=O(1)` for the ARC block.** With `r` and `k` fixed and `M=1` throughout this dataset,
every one of these terms evaluates to the *same tiny number of rows (≤13)* regardless of cluster,
modality's ARC branch aside. This is a case exactly like the `j_cand_elems` warning in the
briefing: the class-in-`M` matters in principle (an ensemble config with many models per biosample
would expose the `O(M²)` term), but it is **completely unconfirmable from these six clusters**
because `M` never varies.

## 5. Benchmark cross-check

```
multiome  s: 1.09 1.07 1.10 1.10 1.06 1.03   (jurkat, jurkat_pma_cd3_4hr, k562_crispri, telohaec_crispri, thp1_1, thp1_2)
scatac    s: 13.62 13.63 13.63 13.63 13.63 13.63
```
(`awk -F'\t' -v r=make_biosample_feature_table 'NR==1||$2==r' benchmarks_{multiome,scatac}.tsv`)

Two findings, and they disagree with each other in an informative way:

1. **Within each modality, wall time is flat to measurement noise** (multiome: 1.03–1.10 s,
   1.07× spread; scATAC: 13.62–13.63 s, ~1.001× spread) even though `n_frag` spans 15.1×, `n_umi`
   12.0×, and `n_cells` 6.0× across these same six clusters. This is not "too narrow a range to
   tell" (as with Tier 3 variables) — it is a **zero-effect finding**, exactly as expected from
   the code: none of `n_frag`/`n_umi`/`n_cells`/`n_peaks`/etc. is ever touched by this script.
   `cpu_time` (0.60–0.67 s multiome, 0.60–0.67 s scatac) is likewise flat and small.
   `mean_load` values of ~46–58 (multiome) vs. ~4–5 (scatac) are just `cpu_time/wall×100` and
   reflect the wall-time gap below, not parallelism (no `threads:` declared; single R process).
2. **Between modalities, there is an unexplained ~13× constant gap** (1.07 s mean vs. 13.63 s
   mean) that the code cannot explain: the script is byte-identical in both runs, the input
   `feature_table.tsv` is nearly identical in size (461 vs. 453 bytes), `M=1` in both, and the
   conda env hash is literally the same directory
   (`../scE2G/.snakemake/conda/49efddd2545177ca66d0194982a8c760_`, confirmed in both slurm logs).
   `io_in` is reported as 0.00 MB for every multiome job but 6.8–57 MB for every scATAC job —
   far larger than the <1 KB actual data file, so this is almost certainly page-ins for R/conda
   shared libraries (`library(data.table)`, `library(dplyr)`), not the rule's real input. **I
   cannot determine why the scATAC run's cold-start cost was ~13× higher** — candidates are
   filesystem/page-cache state at the time each modality's batch ran (multiome: 2026-08-19
   17:12, scATAC: 2026-08-21 00:02 — different sessions, ~36 h apart), or contention among the six
   concurrently-launched sibling jobs reading the same conda-prefix libraries over Lustre/NFS at
   job start (all 6 scATAC jobs started within the same second per their slurm log timestamps).
   Both are guesses; I did not find code or config that would make the modality itself the cause.
   This is a **single measurement per modality**, not a replicate (per Briefing §4.6), so the gap
   cannot be resolved further from this data.

**Conclusion: this rule's cost is dominated by R interpreter + conda environment startup, not by
its algorithm or its inputs — exactly the "constant `t(n)` dominated by R startup" case the
launch instructions anticipated. Do not read a real scaling law into the 13× modality gap; it is
environmental, not a property of `M`, `r`, or `k`.**

## 6. Concrete overhead sources

- **I/O volume**: real data read is tiny — one `feature_table.tsv` per job, 461 or 453 bytes.
  Real data written is 818 or 442 bytes. Both are far below any threshold where I/O volume could
  plausibly cost seconds.
- **Declared-vs-realised threading**: no `threads:` declared → 1 core, matches every other ABC/
  ENCODE_rE2G rule surveyed in this project; not a special finding for this rule.
- **Serialization/deserialization**: `fread`/`fwrite` on ≤13-row tables — negligible compute,
  but each call still pays R's per-call function/package dispatch overhead, part of the constant.
- **Redundant/unused input declaration**: `config["ABC_BIOSAMPLES"]` is declared as this rule's
  sole `input:` (`predictions.smk:4`) but is never opened by the script (§2). It exists purely to
  create a dependency edge, and because it is rewritten at every parse by `make_biosample_config()`
  (landmine, not mine to fix), this rule always shows `reason: Forced execution` — meaning the
  "constant" measured here is a cold-start cost paid on every single invocation of this rule,
  with no batching/warm-cache benefit to expect in production runs either, since the same
  landmine will force it every time Snakemake reparses.
- **Repeated passes over the same data**: none beyond the `rbind`-in-loop pattern noted in §4,
  which is unexercised at `M=1`.
- **Subprocess/interpreter startup**: this is the actual bottleneck — conda env activation +
  `library(data.table)` + `library(dplyr)` — visible directly in both slurm logs before any of
  the script's own logic runs, and consistent with `io_in` values on the scATAC side being
  50-100× larger than the actual input file.

## 7. `t_est(M)`, `t_low(M)`, `t_high(M)`

Because the measured cost is startup-dominated and `M` cannot be varied in this dataset, these
functions separate a **measured, modality-specific startup constant `C`** from a **code-derived,
unconfirmed `M`-dependent term**. `r=6` (rows in one model's feature table), `k=7` (ARC block
rows), both constants in this dataset.

- `C_multiome = 1.07` s (mean of the six measured values 1.09, 1.07, 1.10, 1.10, 1.06, 1.03)
- `C_scatac = 13.63` s (mean of 13.62, 13.63, 13.63, 13.63, 13.63, 13.63)

```
t_est(M)  = C + 0.001·M·r        seconds     (guessed marginal fread+rbind cost of 1 ms per row;
                                               NOT measured — M never varies in this dataset)

t_low(M)  = C + 0.0005·M·r       seconds     (class lower bound: O(M·r), i.e. rbind amortizes
                                               and does not recopy the growing frame each time)

t_high(M) = C + 0.0005·M²·r      seconds     (class upper bound: O(M²·r), i.e. the naive
                                               rbind-in-loop antipattern at line 10 recopies the
                                               whole accumulated frame every iteration)
```

where `C = C_multiome` or `C_scatac` depending on modality, per §5. The `0.001`/`0.0005`
coefficients are order-of-magnitude guesses for R's per-row `fread`/`rbind` cost on tables this
small (tens of rows), not fits to data — flagged explicitly as guesses because no measurement in
this dataset exercises `M>1`. **The important claim is the exponent bound (`M` vs. `M²`), not
these coefficients**: at `M=1` (the only value ever observed), all three collapse to
`t ≈ C + O(10⁻³)` s — i.e., the `M`-term is unconditionally negligible next to `C` for any `M`
this pipeline would plausibly use (a handful of ensembled models, not thousands), so in practice
`t(M) ≈ C` is the operative model regardless of which class is correct.

## 8. Evaluated at mean / min (`thp1_1`) / max (`telohaec_crispri`) cluster size

Because this rule has **no dependence on cluster identity or size** — only on modality and the
constant `M=1` — all three cluster points evaluate to the *same* number within a modality; that
flatness is itself the finding, not an artifact of rounding.

| Modality | Cluster | `t_est` | `t_low` | `t_high` | measured `s` |
|---|---|---:|---:|---:|---:|
| multiome | mean-cluster (n/a, no size dependence) | 1.076 | 1.0725 | 1.0725 | 1.075 (mean) |
| multiome | thp1_1 (min `n_frag`/`n_umi`/`n_cells`) | 1.076 | 1.0725 | 1.0725 | 1.06 |
| multiome | telohaec_crispri (max) | 1.076 | 1.0725 | 1.0725 | 1.10 |
| scatac | mean-cluster | 13.636 | 13.6325 | 13.6325 | 13.628 (mean) |
| scatac | thp1_1 (min) | 13.636 | 13.6325 | 13.6325 | 13.63 |
| scatac | telohaec_crispri (max) | 13.636 | 13.6325 | 13.6325 | 13.63 |

(`t_est`/`t_low`/`t_high` computed at `M=1, r=6` → `M·r=6`, `M²·r=6`; `0.001×6=0.006`,
`0.0005×6=0.003`.)

## 9. Both modalities: class or constant?

**Only the constant differs, and the code path is almost identical.** The script executed is
byte-identical in both modalities. The only branch that differs is lines 14-30: multiome's model
feature table contains `ARC.E2G.Score`, so the 7-row ARC block is appended (one extra `O(k)`
`rbind`, `k=7` constant); scATAC's model feature table contains `ABC.Score` directly and neither
`ARC.E2G.Score` nor `Kendall`, so that block is skipped entirely. This is a fixed `O(1)` amount of
extra work for multiome, not a class difference. The dominant, order-of-magnitude difference
between modalities (13.63 s vs. 1.07 s) is the unexplained startup-constant gap from §5, which has
nothing to do with this branch (it shows up as a flat per-modality wall time, not a per-row
difference of ~7 rows).

## 10. What the measurement cannot tell us

- **Whether the `rbind`-in-loop at line 10 is `O(M·r)` or `O(M²·r)`.** `M=1` in all 12 measured
  jobs (both modalities have exactly one model dir per biosample — confirmed from
  `expanded_biosample_config.tsv`'s `model_dir` column, no commas). Testing this would require an
  ensemble-model config with `M>1`, which does not exist in this dataset.
- **Why the scATAC run's cold-start cost (13.63 s) is ~13× the multiome run's (1.07 s)**, given
  the identical script, nearly identical input size, identical conda env path, and `M=1` in both.
  Candidate explanations (filesystem cache state at run time, contention among six
  simultaneously-launched sibling jobs sharing a conda prefix over Lustre) are plausible but
  unverified; this is one measurement per modality, and the difference could equally be an
  artifact of when/how each batch was launched rather than anything about the rule itself.
- **Whether any real production config ever sets `M>1`** (ensembles of models per biosample). If
  it does, the `O(M²)` risk in §4 becomes relevant; this project's data cannot say whether it
  ever occurs.
- **Whether `config["ABC_BIOSAMPLES"]` being declared-but-unread as `input:` (§2, §6) has any
  cost consequence beyond forcing re-execution.** I observed the `Forced execution` /
  `jobid: 0` pattern in all 12 logs and connected it to the known landmine, but I did not attempt
  to measure what a non-forced, DAG-batched invocation of this rule would cost, since producing
  one would require changing the pipeline's rerun-trigger behavior — out of scope here.
- Whether the checkpoint's downstream ARC/Neither decision (§ mechanism below) would ever flip
  *within* a modality for a different model assignment; in this dataset it is modality-determined
  (by which features the assigned model uses), not cluster-determined, for all 12 jobs observed.

## Mechanism note for Phase 2 (not my rule, but I own the table that feeds it)

This rule's only externally consequential output is the **union feature table**
(`{biosample}/feature_table.tsv`, §3: 818 bytes/13 rows for multiome, 442 bytes/7 rows for
scATAC, constant across all six clusters in a modality). `checkpoint features_required`
(`workflow/rules/add_external_features.smk:1-22`) reads exactly this file, checks for
`"ARC.E2G.Score"` or `"Kendall"` in column 2 or 3 of any row (`add_external_features.smk:14-17`),
and writes `to_generate.txt` = `"ARC"` if either is present, else `"Neither"`
(`add_external_features.smk:18-22`). Verified on disk (§ above): `to_generate.txt` = `ARC` for
all six multiome clusters, `Neither` for all six scATAC clusters — driven entirely by which model
(`multiome_powerlaw_v3` vs. `scATAC_powerlaw_v3`) was assigned, not by anything cluster-specific.
`features_to_generate()` (`add_external_features.smk:26-34`) then routes to the
`make_kendall_pairs → generate_atac_matrix → compute_kendall → arc_e2g` chain (24 jobs: 4 rules ×
6 multiome clusters) when `ARC`, or short-circuits to nothing extra when `Neither`. This rule
(`make_biosample_feature_table`) is therefore the node whose ~1-14 s of near-fixed cost gates a
~15×-cost-differential subgraph (`add_external_features`: 12.6–21.3 min multiome vs. 0.9–1.3 min
scATAC, per the briefing) plus the entire 24-job ARC chain — a large topology decision made by an
essentially free node. Unit #26 (`features_required` itself) owns the checkpoint's own cost;
this is only the mechanism by which this table's *content* (not its cost) determines the DAG
shape Phase 2 needs to reason about.

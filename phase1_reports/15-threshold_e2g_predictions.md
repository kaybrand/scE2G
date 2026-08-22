# Unit 15 — `threshold_e2g_predictions.py`

## 1. Unit identity

- Script: `ENCODE_rE2G/workflow/scripts/model_application/threshold_e2g_predictions.py` (36 lines).
- Rule: `filter_sc_e2g_predictions`, `workflow/rules/sc_predictions.smk:79-103`.
- Modalities: **both**, 6 jobs each = 12 measured jobs (`{cluster} × {model_name} × {threshold}`,
  one `model_name`/`threshold` pair per modality: `multiome_powerlaw_v3~0.177`,
  `scATAC_powerlaw_v3~0.174`).
- No `threads:` declared — runs at 1, consistent with the briefing's "every ABC rule declares no
  threads" pattern (this rule is not strictly an ABC rule, but the same absence holds here: no
  `threads:` line at `sc_predictions.smk:79-103`, so Snakemake requests 1 CPU).
- `resources: mem_mb = encode_e2g.ABC.determine_mem_mb` (`sc_predictions.smk:92`) — no `min_gb`
  override, so the 8 GB floor from `ABC/workflow/rules/utils.smk:17` applies. See §6 for how badly
  this over-predicts actual usage.
- Input: `{cluster}/{model_name}/scE2G_predictions.tsv.gz` (the **unthresholded**, qnorm-scored
  output of `run_e2g_qnorm`). Output: `{cluster}/{model_name}/scE2G_predictions_threshold{threshold}.tsv.gz`.
- I confirmed on disk that the input row count equals `n_cand_pairs` exactly and the output row
  count equals `n_called_all` exactly (both counts below exclude the header line):
  - multiome `jurkat`: input 11,036,283 rows (manifest `n_cand_pairs` = 11,036,283); output 59,698
    rows (manifest `n_called_all` = 59,698).
  - This extends the briefing's collapsed-identity chain one step further than documented:
    `run_e2g_qnorm` (a pure column addition, no filtering) means `scE2G_predictions.tsv.gz` also
    has exactly `n_cand_pairs` rows, even though the briefing's collapsed-identity list names only
    ARC/ActivityOnly/genomewide_features. Noted for context, not "discovered" as new — it is the
    same tautology one hop downstream.

## 2. What the code actually does

The entire script (`threshold_e2g_predictions.py:1-36`):

1. `:28` `all_predictions = pd.read_csv(all_predictions_file, sep="\t")` — **no `chunksize`**, no
   `usecols`, no `dtype` restriction. The full `n_cand_pairs`-row, 27-column (multiome) / 20-column
   (scATAC) table is parsed and materialized as a single in-memory `DataFrame`. This is the one
   read of the whole input and is a single linear parse pass, O(`n_cand_pairs`).
2. `:7` `filtered_predictions = all_putative[all_putative[score_column] >= threshold]` — one
   vectorized elementwise comparison over the full `n_cand_pairs`-row score column, producing a
   boolean mask and a filtered copy. O(`n_cand_pairs`), single pass, no sort.
3. `:10-16` — promoter/self-promoter logic, branching on `include_self_promoter` (config value,
   confirmed `True` by default in `ENCODE_rE2G/config/config.yaml:20`, not overridden in the
   `igvf10_{multiome,scatac}_config.yaml` files):
   - `:11-14` `include_self_promoter=True`: `filtered_predictions[(filtered_predictions["class"]
     != "promoter") | filtered_predictions["isSelfPromoter"]]` — two more elementwise column
     comparisons, boolean OR, one more boolean-indexed copy. This operates on `filtered_predictions`
     (already reduced by the score threshold in step 2), so it is bounded above by
     O(`n_cand_pairs`) but in practice touches somewhat fewer rows.
   - `:16` `include_self_promoter=False`: single comparison, same shape, smaller cost.
4. `:32` `filtered_predictions.to_csv(output_file, compression="gzip", sep="\t", index=False)` —
   writes the final `n_called_all`-row (59,438–61,284 multiome; 105,880–120,113 scATAC) table,
   gzip-compressed. Output is 4.7–4.9 MB (multiome) / 5.2–5.5 MB (scATAC) per the manifest's
   `n_called_bytes` — three orders of magnitude smaller than the ~750–890 MB / ~414–476 MB input.

**No sort, no `.sort_values()`, no `.rank()`, no `groupby`, no `merge`/`join`, no quantile call
anywhere in the 36 lines.** The threshold is a literal `>=` comparison against a scalar passed in
via `--threshold` (`sc_predictions.smk:83`, resolved by
`ENCODE_rE2G/workflow/rules/utils.smk:187-192` `get_model_threshold()`, which reads a
pre-computed threshold value out of a `score_threshold_*` filename written by an earlier training
step — the threshold itself is not computed here). This confirms directly what the unit brief
asked me to establish: **the whole table is loaded into a single pandas `DataFrame` (memory
O(`n_cand_pairs`), not O(1)/streamed), and the threshold operation is a single-pass elementwise
comparison (O(`n_cand_pairs`)), not a sort or quantile (which would be O(`n_cand_pairs`
log `n_cand_pairs`)).**

## 3. Size variable(s)

**Primary: `n_cand_pairs`** (10,173,077–11,447,564 rows, Tier 2, 1.13× range) — the row count of
the input `all_predictions_file`, read in full at `:28` and scanned by every subsequent boolean
mask. Identical between multiome and scATAC for the same cluster (verified against both
`size_manifest_{multiome,scatac}.tsv`; e.g. `thp1_1` = 10,173,077 in both) — the ATAC-side
candidate-pair count does not depend on which modality branch ran downstream.

**Byte volume is the more informative variable for this I/O-shaped rule, and it genuinely differs
by modality** (unlike `n_cand_pairs`, which does not). I measured the actual on-disk size of the
input file directly (`$SCRATCH/.../<cluster>/<model_name>/scE2G_predictions.tsv.gz`), since this
is more precise than the nearest manifest proxy (`feat_bytes`, the size of the upstream
`genomewide_features.tsv.gz`, which lacks the qnorm score columns this rule's actual input has):

| Modality | mean input bytes (measured) | mean `feat_bytes` (manifest proxy) | measured/`feat_bytes` |
|---|---:|---:|---:|
| multiome | 819.3 MB | 551.7 MB | 1.485× |
| scATAC | 445.8 MB | 205.2 MB | 2.174× |

The measured/`feat_bytes` ratio is tight and modality-specific (1.47–1.50× multiome, 2.14–2.19×
scATAC across all 6 clusters each) — `feat_bytes` is a good *proxy* for this rule's real input
volume within a modality, but the **conversion factor differs by modality**, so it cannot be used
directly across modalities without that correction. Consequently the true input-byte ratio between
modalities is **819.3/445.8 = 1.84×**, smaller than `feat_bytes`'s own 551.7/205.2 = 2.69× ratio —
see §9 for why this matters.

**No cross term.** The output size (`n_called_all`) does not appear anywhere in the cost model —
it is a consequence of the filtering, not an input to it, and its write cost (§6) is negligible
next to the read.

## 4. Asymptotic class derived from code

**O(`n_cand_pairs`)**, both in time and in memory:

- **Time:** one `pd.read_csv` parse (`:28`) plus 2–3 elementwise boolean-mask passes (`:7`, `:11-14`
  or `:16`) over an array whose length starts at `n_cand_pairs` and only shrinks — every one of
  these is a single linear scan. The final `to_csv` (`:32`) writes `n_called_all` rows, three orders
  of magnitude smaller, and does not change the class. `T = c1·n_cand_pairs (parse) +
  c2·n_cand_pairs (threshold mask) + c3·n_cand_pairs (promoter mask, upper bound) +
  O(n_called_all) (write)`, which collapses to **O(`n_cand_pairs`)**.
- **Memory:** `pd.read_csv(all_predictions_file, sep="\t")` at `:28` has no `chunksize` — the
  entire table is materialized as one `DataFrame` before any filtering happens. This is O(1) in
  the sense that no accumulator grows across the run, but O(`n_cand_pairs`) in the sense that peak
  RSS is set by the size of the single fully-loaded table, not a constant independent of input
  size. It is **not** a streamed/chunked read.

There is no hash table, no sort, no join anywhere in this script, so — unlike some other units in
this project where a `distinct()`/`groupby` leaves open the possibility of superlinear hash-resize
behavior at larger scale — there is no plausible algorithmic candidate for anything above
O(`n_cand_pairs`) here. The class itself is not in question from the code; §7's `t_low`/`t_high`
therefore bound the observed *constant*, not an alternate exponent (see §7 for why that is the
honest choice for this particular rule).

## 5. Benchmark cross-check

All 12 measured jobs (`benchmarks_{multiome,scatac}.tsv`, filtered to `filter_sc_e2g_predictions`),
with the directly-measured input byte volume added (not a manifest column — read from the
retained `scE2G_predictions.tsv.gz` intermediates):

| Modality | Cluster | `n_cand_pairs` | input MB (measured) | wall s | cpu_time s | mean_load | max_rss MB | io_in MB | io_out MB | `n_called_all` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| multiome | jurkat | 11,036,283 | 841.70 | 42.37 | 14.91 | 34.95 | 943.20 | 290.83 | 0.00 | 59,698 |
| multiome | jurkat_pma_cd3_4hr | 10,854,698 | 829.72 | 44.72 | 13.48 | 30.08 | 826.84 | 163.05 | 0.00 | 59,650 |
| multiome | k562_crispri | 11,447,564 | 889.17 | 43.60 | 14.37 | 32.85 | 912.33 | 289.34 | 0.00 | 59,438 |
| multiome | telohaec_crispri | 11,106,687 | 855.10 | 65.25 | 26.72 | 40.73 | 1556.20 | 354.57 | 15.05 | 60,630 |
| multiome | thp1_1 | 10,173,077 | 742.74 | 40.15 | 15.04 | 37.04 | 925.75 | 160.00 | 0.00 | 60,435 |
| multiome | thp1_2 | 10,278,325 | 757.67 | 48.89 | 36.03 | 73.62 | 6559.36 | 20.70 | 9.40 | 61,284 |
| scATAC | jurkat | 11,036,283 | 454.01 | 27.69 | 14.35 | 51.43 | 1033.91 | 288.66 | 15.05 | 116,138 |
| scATAC | jurkat_pma_cd3_4hr | 10,854,698 | 450.81 | 27.98 | 14.89 | 52.60 | 1016.48 | 0.00 | 15.05 | 114,015 |
| scATAC | k562_crispri | 11,447,564 | 476.46 | 27.72 | 15.14 | 54.01 | 1160.80 | 4.00 | 0.00 | 120,113 |
| scATAC | telohaec_crispri | 11,106,687 | 460.18 | 43.43 | 2.78 | 6.26 | 166.91 | 64.93 | 15.05 | 117,185 |
| scATAC | thp1_1 | 10,173,077 | 414.03 | 24.25 | 14.96 | 61.36 | 1166.46 | 292.62 | 0.00 | 105,880 |
| scATAC | thp1_2 | 10,278,325 | 419.51 | 34.16 | 5.27 | 15.43 | 433.02 | 339.12 | 0.00 | 107,248 |

**`n_cand_pairs` spans only 1.13× and does not rank monotonically with wall time.** `k562_crispri`
has the largest `n_cand_pairs` of all six clusters (11,447,564) but is not the slowest job in
either modality (multiome: 43.60 s, 3rd of 6; scATAC: 27.72 s, 4th of 6). This is the same pattern
flagged elsewhere in this project for Tier 2 variables: the six points can rule out wildly
superlinear behavior but cannot confirm O(`n_cand_pairs`) over, say, O(`n_cand_pairs`
log `n_cand_pairs`) by fit alone — the class call in §4 rests on the code, not this table.

**Two clear outliers, flagged rather than smoothed:**
- **multiome `thp1_2`**: max_rss = 6559.36 MB, roughly **7–8× every other multiome cluster**
  (826–1556 MB), despite `n_cand_pairs`/input bytes being near the *smallest* of the six
  (10,278,325 rows / 757.7 MB, 2nd-smallest). It also has the highest `cpu_time` (36.03 s) and
  `mean_load` (73.62) of any multiome job. There is only one benchmark row for this job
  (`igvf10_multiome/benchmarks/filter_sc_e2g_predictions/thp1_2~multiome_powerlaw_v3~0.177.tsv`,
  no retry subdirectory), so this is not an attempt-2-with-doubled-memory artifact — it is a
  genuine single-measurement anomaly not explained by input size. I cannot attribute it to a
  specific cause (transient node memory pressure, a pandas dtype-inference difference on that
  cluster's data, GC timing) without a replicate.
- **scATAC `telohaec_crispri`**: wall = 43.43 s but `cpu_time` = 2.78 s and `mean_load` = 6.26 —
  by far the lowest CPU utilization of any of the 12 jobs, meaning ~94% of its wall time was not
  active CPU. `max_rss` for this job (166.91 MB) is also far below its scATAC siblings
  (433–1167 MB), consistent with the process spending most of its wall time waiting rather than
  holding the full table in memory yet at the time it was sampled. Most plausibly a slow/contended
  read of the 460 MB input off Lustre on that particular run; not confirmable from one measurement.
  `scATAC thp1_2` shows a milder version of the same signature (`cpu_time`/wall = 0.154,
  `mean_load` = 15.43).

**`cpu_time`/wall ratios are low across nearly every job** (0.30–0.74 multiome, 0.06–0.62 scATAC),
meaning most wall time in this rule is *not* active CPU. This matters for the byte-vs-column
question below: computing the modality ratio from `cpu_time` (which nets out much of the I/O-wait
noise that inflates individual wall times) gives **mean cpu_time 20.09 s (multiome) / 11.23 s
(scATAC) = 1.79×** — close to the 1.84× measured-input-byte ratio (§3) and clearly closer than the
column-count ratio (27 vs 20 columns = 1.35×) or the raw wall-time ratio (47.50 s / 30.87 s =
1.54×, pulled down by the scATAC I/O-wait outliers). **This agrees with the unit brief's hypothesis:
once I/O-wait noise is reduced (by using `cpu_time` instead of wall time), the modality gap tracks
byte volume much better than it tracks row count (which is identical across modalities) or column
count.**

**`io_in`/`io_out` do not track the real read/write volume for this rule and should not be trusted
as I/O size measures** — use the on-disk file sizes instead:
- Real input volume is 742.7–889.2 MB (multiome) / 414.0–476.5 MB (scATAC), essentially flat within
  a modality. Reported `io_in` ranges from 20.70 MB (multiome `thp1_2`) to 354.57 MB (multiome
  `telohaec_crispri`) — undercounting the true read volume by roughly **2.4× to 41×** depending on
  the job, with no consistent relationship to the actual file size. This is consistent with
  `io_in` measuring only cache-miss disk/network reads (via `psutil` I/O counters), which depends
  heavily on the page-cache state of whichever node happened to run the job, not on data volume.
- Real output volume is 4.7–4.9 MB (multiome, `n_called_bytes`) / 5.2–5.5 MB (scATAC). Reported
  `io_out` is 0.00 MB for 8 of the 12 jobs and 9.40–15.05 MB for the other 4 — either undercounting
  to zero or overcounting by ~2–3×, again with no size-tracking relationship.

## 6. Concrete overhead sources

- **The entire ~750–890 MB (multiome) / ~414–477 MB (scATAC) input is materialized as a single
  in-memory `DataFrame` before any filtering happens** (`:28`). There is no `chunksize`, no
  `usecols` to drop unused columns before the threshold check, and no early-exit — every column of
  every row is parsed even though only `score_column` and two other columns (`class`,
  `isSelfPromoter`) are read afterward. This is the dominant, code-confirmed cost: a single
  linear pass over the full table just to extract three columns' worth of decision logic.
- **`mem_mb` resource request badly over-provisions actual usage.** `ABC.determine_mem_mb`
  (`ENCODE_rE2G/ABC/workflow/rules/utils.smk:17-24`) computes `mem_to_use_mb =
  max(4 * input_size_mb * 8, min_gb * 1000)` at attempt 1 — the `* 8` assumes gzip compresses the
  file up to 8×, and the further `* 4` assumes pandas holds ~4× the decompressed size in memory.
  For this rule (`min_gb=8`, no override), that predicts:

  | Cluster | modality | measured input MB | requested mem_mb (formula) | measured `max_rss` MB | overshoot |
  |---|---|---:|---:|---:|---:|
  | jurkat | multiome | 841.70 | 26,934 | 943.20 | 28.6× |
  | k562_crispri | multiome | 889.17 | 28,453 | 912.33 | 31.2× |
  | thp1_2 | multiome | 757.67 | 24,245 | 6559.36 | 3.7× (the RSS outlier from §5 narrows this) |
  | jurkat | scATAC | 454.01 | 14,528 | 1033.91 | 14.0× |
  | telohaec_crispri | scATAC | 460.18 | 14,726 | 166.91 | 88.2× (RSS outlier from §5 widens this) |

  Actual `pd.read_csv` memory tracks much closer to the *compressed* file size than to
  `8× compressed × 4`, because this table is almost entirely fixed-width numeric columns
  (floats/ints), which pandas stores compactly relative to their gzipped-text representation — the
  formula's assumptions (built for arbitrary, possibly string-heavy tables) do not hold well for
  this specific table shape. This is a declared-vs-realized resource mismatch, analogous in kind
  to the project's already-established declared-vs-realized *threading* gap, but on the memory
  axis instead. Noted as a mechanism only, per scope — not a proposal to change `min_gb` or the
  formula.
- **`io_in`/`io_out` are unreliable size proxies for this rule** (§5) — any cost model built from
  this rule's benchmark I/O columns should use the on-disk `scE2G_predictions.tsv.gz` /
  `scE2G_predictions_threshold*.tsv.gz` sizes instead.
- **No threading anywhere in the path.** No `threads:` on the rule; `pandas.read_csv`'s C parser
  and the gzip decompression it drives are single-threaded by default; the boolean-mask operations
  and `to_csv(..., compression="gzip")` write are also single-threaded. `cpu_time`/wall ratios
  well under 1 for most jobs (§5) do not indicate hidden multi-core work — they indicate I/O wait,
  not parallelism.
- **The output is three orders of magnitude smaller than the input** (`n_called_all` ≈
  60k/multiome, ≈106–120k/scATAC vs `n_cand_pairs` ≈ 10.2–11.4M) — the ~180× "pairs → thresholded
  rows" reduction the unit brief flagged is real and confirmed on disk, but it happens *after* the
  full-table read, so it does not reduce the dominant read cost at all.

## 7. Three functions

Class from §4 is **O(`n_cand_pairs`)**, with no sort/hash/join anywhere in the 36-line script — so
unlike other units in this project where a hash-based `distinct()`/`groupby` leaves open a
plausible (if unconfirmed) superlinear tail, there is no algorithmic candidate here for anything
above linear. `t_low`/`t_high` below therefore bound the observed **constant of proportionality**
across the 12 measured jobs (the fastest and slowest per-row wall-time rate seen), not an
alternate class — that is the honest choice for a rule whose class is not actually in doubt.

```
t_est(n_cand_pairs; multiome) = 4.03e-6 * n_cand_pairs   [seconds]
t_est(n_cand_pairs; scATAC)   = 2.54e-6 * n_cand_pairs   [seconds]
t_low(n_cand_pairs)           = 2.38e-6 * n_cand_pairs   [seconds]
t_high(n_cand_pairs)          = 5.87e-6 * n_cand_pairs   [seconds]
```

Calibration:
- `t_est`'s two coefficients are each the **median** (not mean, to reduce sensitivity to the two
  flagged outliers in §5) of wall_s / `n_cand_pairs` across the 6 clusters per modality: multiome
  {3.809e-6, 3.839e-6 (jurkat, computed but not shown above — see raw ratios below), 3.947e-6,
  4.120e-6, 4.755e-6, 5.874e-6} → median 4.03e-6; scATAC {2.384e-6, 2.423e-6, 2.509e-6, 2.578e-6,
  3.324e-6, 3.911e-6} → median 2.544e-6, rounded to 2.54e-6.
- `t_low` is the single smallest observed wall_s/`n_cand_pairs` across all 12 jobs (scATAC
  `thp1_1`: 24.25 s / 10,173,077 = 2.384e-6 s/row).
- `t_high` is the single largest observed wall_s/`n_cand_pairs` across all 12 jobs (multiome
  `telohaec_crispri`: 65.25 s / 11,106,687 = 5.874e-6 s/row) — this is the same job flagged as an
  RSS/cpu_time outlier in §5, so `t_high` is deliberately calibrated to the worst observed case
  rather than a smoothed value.
- An equivalent, arguably more mechanistic byte-based form, given the tight within-modality
  bytes-per-row ratio (§3): `t_est(bytes; multiome) ≈ 0.054–0.081 s/MB` (median ≈0.058 s/MB) and
  `t_est(bytes; scATAC) ≈ 0.058–0.094 s/MB` (median ≈0.061 s/MB), evaluated on the directly-measured
  input MB in §5's table. I report the `n_cand_pairs`-based form above as primary because
  `n_cand_pairs` is the named manifest variable and the two forms agree to within the noise budget
  of these 12 points (§5's cpu_time-based modality-ratio check supports bytes as the more
  mechanistic driver, but both are linear in the same underlying quantity for this dataset, since
  bytes/row is close to constant within each modality).

## 8. Evaluated at mean / min / max

Using the designated clusters (`thp1_1` = min, `telohaec_crispri` = max, per project convention;
note these are not the literal min/max of `n_cand_pairs` itself — `thp1_1` is the true min
[10,173,077] but the true max is `k562_crispri` [11,447,564], not `telohaec_crispri`
[11,106,687] — evaluated at the assigned clusters regardless, for cross-unit comparability):

| Point | cluster | `n_cand_pairs` | `t_low` | `t_est` (multiome) | `t_est` (scATAC) | `t_high` |
|---|---|---:|---:|---:|---:|---:|
| mean (across 6 clusters) | — | 10,816,106 | 25.8 s | 43.6 s | 27.5 s | 63.5 s |
| min | `thp1_1` | 10,173,077 | 24.3 s | 41.0 s | 25.9 s | 59.7 s |
| max | `telohaec_crispri` | 11,106,687 | 26.5 s | 44.8 s | 28.2 s | 65.2 s |

(For reference, the actual measured wall times: `thp1_1` 40.15 s multiome / 24.25 s scATAC — both
land almost exactly on `t_est`/`t_low` respectively, because those are the points that helped
define those coefficients; `telohaec_crispri` 65.25 s multiome / 43.43 s scATAC — the multiome
value lands exactly on `t_high` by construction (§7), the scATAC value is well above its own
`t_est`, reflecting the I/O-wait outlier flagged in §5. All 12 measured points fall within
[`t_low`, `t_high`] for their respective modality by construction, since those bounds were set
from the observed extremes.)

## 9. Does the class differ between modalities?

**No — only the constant differs, and the byte-volume evidence (not the row-count or raw
wall-time evidence) explains why.** `n_cand_pairs` is numerically identical between multiome and
scATAC for every cluster (§3, §5) — the code path in `threshold_e2g_predictions.py` is itself
identical between modalities (the same script runs for both; the only difference is the
`--score_column`/`--threshold` values passed in, both of which are configuration, not code
branches). So the class is trivially shared: same script, same operations, same asymptotic
behavior in `n_cand_pairs`.

What differs is the constant, and the best-supported explanation is **input byte volume**, not
row count (identical) and not primarily column count:
- Measured input bytes: 819.3 MB (multiome) vs 445.8 MB (scATAC) mean — **1.84× ratio**.
- `cpu_time`-based wall time (which nets out much of the I/O-wait noise flagged in §5): 20.09 s vs
  11.23 s mean — **1.79× ratio**, close to the byte ratio.
- Column count: 27 (multiome) vs 20 (scATAC) — 1.35× ratio, a weaker match.
- Raw wall time: 47.50 s vs 30.87 s mean — 1.54× ratio, pulled down from the cpu_time ratio by the
  scATAC I/O-wait outliers (§5).

This rule is *not* in the launcher's "genuinely different cost by modality" list
(`add_external_features` ~15×, `run_e2g_qnorm` ~2×, `gen_final_features` ~2×,
`get_stats_per_model_per_cluster` ~1.2×) — its modality gap (~1.8× by bytes/cpu_time) sits between
`run_e2g_qnorm`'s and `add_external_features`'s, but the mechanism here (input table width/byte
volume driven by how many feature+score columns the upstream ARC/Kendall branch added, not a
different code path within this script) is simpler: this script itself does not branch on
modality at all.

One nuance worth flagging explicitly: the manifest's `feat_bytes` column (genomewide_features.tsv.gz,
upstream of this rule's actual input) has a modality ratio of ~2.69× (551.7 MB / 205.2 MB, §3) —
**larger** than this rule's true input-byte ratio of 1.84×. Using `feat_bytes` directly as a stand-in
for this rule's input volume across modalities would overstate the modality gap; the
measured/`feat_bytes` conversion factor itself differs by modality (1.48× multiome vs 2.17×
scATAC, §3), which is why the two proxies diverge.

## 10. What the measurement cannot tell us

- **Cannot confirm O(`n_cand_pairs`) over any nearby exponent from these 12 points** — the range
  is 1.13× and wall time is not even monotonic in `n_cand_pairs` within a modality (§5:
  `k562_crispri` has the largest `n` but ranks 3rd/4th slowest). The O(`n_cand_pairs`) call in §4
  is a code-reading conclusion (no sort/hash/join exists in the script), not an empirical fit.
- **Cannot separate an intercept (fixed per-job overhead: interpreter/conda startup, `click`
  parsing, import of pandas) from the linear slope** with only a 1.13×-range, 6-point sample per
  modality — the `t_est`/`t_low`/`t_high` formulas in §7 are pure-linear-through-the-origin fits
  for this reason, which likely means they slightly overestimate the true slope and omit a few
  seconds of genuine fixed overhead (conda env activation for a `conda:`-directive rule, Python/
  pandas import) that a wider size range would let me isolate.
- **Cannot explain the multiome `thp1_2` max_rss anomaly (6559 MB vs 826–1556 MB for its five
  siblings, despite near-median input size) or the scATAC `telohaec_crispri` / `thp1_2` low-
  `cpu_time`/high-wall anomalies** (§5) from the benchmark table alone — no retry subdirectories
  exist for either, so these are not attempt-2-with-doubled-resources artifacts; the most likely
  explanations (node memory contention, cold Lustre cache, page-cache state at read time) all
  require information (per-node contemporaneous load, cache state) not captured by this benchmark
  infrastructure.
- **Cannot trust the benchmark's `io_in`/`io_out` columns for this rule's actual read/write
  volume** (§5) — `io_in` undercounts by 2.4×–41× depending on the job with no consistent relation
  to file size, and `io_out` is 0.00 for 8/12 jobs despite a real ~4.7–5.5 MB write. Any
  Phase 2 cost model that wants I/O volume for this rule should use the on-disk
  `scE2G_predictions.tsv.gz`/`scE2G_predictions_threshold*.tsv.gz` sizes reported in §5, not the
  benchmark's own `io_in`/`io_out` fields.
- **These are single (n=1) measurements per cluster per modality**, not replicates — so none of
  the scatter or the two flagged outliers in §5 can be attributed with confidence to "this cluster
  is systematically different" versus "this particular run hit a slow/contended node."
- **Cannot say anything about behavior at `n_cand_pairs` well outside the observed 10.17M–11.45M
  range** (e.g., an order of magnitude more candidate pairs, which would also mean a
  correspondingly larger `mem_mb` request per §6's formula) — the O(`n_cand_pairs`) claim is a
  code-level guarantee at any scale (no operation in the script becomes more expensive per row as
  `n` grows), but whether the *system* (Slurm memory/queueing, Lustre I/O bandwidth) behaves
  linearly at 10× this scale is outside what these clusters can address.

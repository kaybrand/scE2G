# Unit 17 — `gen_num_candidate_enh_gene.py`

## 1. Unit identity

- Script: `ENCODE_rE2G/workflow/scripts/feature_tables/gen_num_candidate_enh_gene.py` (55 lines).
- Rule: `generate_num_candidate_enh_gene`, `ENCODE_rE2G/workflow/rules/genomewide_features.smk:70-88`.
- Modalities: **both**, 6 jobs each (one per cluster), 12 measured jobs total. No `threads:`
  declared (runs at 1, consistent with the briefing's "every ABC rule declares no threads"
  cross-cutting finding). `resources: mem_mb = partial(ABC.determine_mem_mb, min_gb=16)`.
- Sibling scripts in the same family (owned by units #18/#19, not re-derived here, but see the
  factual note below): `gen_num_tss_enh_gene.py` and `gen_num_sum_nearby_enhancers.py`, both in
  the same directory.

## 2. What the code actually does

Single input, no cross-join with any element or gene list:

```python
pred_df = pd.read_csv(abc_predictions, sep="\t", compression="gzip",
                      usecols=["name", "chr", "start", "end", "TargetGene", "TargetGeneTSS"])
```
(`gen_num_candidate_enh_gene.py:46-47`) reads `Predictions/EnhancerPredictionsAllPutative.tsv.gz`
— the file whose row count *is* `n_cand_pairs` by the manifest's own definition (Tier 2 table:
"`n_cand_pairs` | candidate enhancer–gene pairs; `EnhancerPredictionsAllPutative.tsv.gz`"). I
confirmed the byte sizes match exactly: `n_cand_pairs_bytes` in `size_manifest_multiome.tsv`
(e.g. `281522213` for `thp1_1`) equals the on-disk size of that cluster's
`EnhancerPredictionsAllPutative.tsv.gz` (`281522213` bytes, verified with `ls -la`). This is
already a table of (enhancer, gene) **pairs** — the pairing against genes happened upstream, in
ABC's own candidate-region/prediction code, not in this script. **This script never opens
`EnhancerList.txt` (`j_cand_elems`) or `GeneList.txt` (`k_genes`)** — neither appears anywhere in
the file. See §3 for why this matters against the unit's briefed "expected drivers."

`determine_num_candidate_enh_gene()` (`:6-39`) then, over the full `n_cand_pairs`-row frame:

1. `:8` — computes `midpoint = (start+end)/2`, O(`n_cand_pairs`).
2. `:9` — stashes `_orig_idx = range(len(df))` to allow restoring input order later.
3. `:10` — **`df.sort_values(by=["TargetGene","midpoint"])`** — a full two-key sort of all
   `n_cand_pairs` rows. This is the first of three sorts and the one that does the real work:
   pandas' default sort algorithm (`quicksort`) on a two-key `sort_values` is
   O(`n_cand_pairs`·log(`n_cand_pairs`)).
4. `:12-13` — `is_downstream`/`is_upstream` boolean masks, O(`n_cand_pairs`).
5. `:16-21` — for the upstream subset (~half the rows): `df[is_upstream]` (boolean-index copy),
   **`.sort_values("midpoint", ascending=False)`** (second sort, O(m·log m), m≈`n_cand_pairs`/2),
   then **`.groupby("TargetGene").cumcount()`** — pandas groupby is hash-based, so this is O(m),
   not a nested loop.
6. `:24-28` — for the downstream subset: `.groupby("TargetGene").cumcount()` only (no extra sort
   needed because the global sort at `:10` already ordered by `midpoint` ascending within each
   gene) — O(d), d≈`n_cand_pairs`/2.
7. `:31-32` — `fillna(0)` + `astype("int")`, O(`n_cand_pairs`).
8. `:34` — **`df.sort_values("_orig_idx")`** — third sort, restoring input row order,
   O(`n_cand_pairs`·log(`n_cand_pairs`)).
9. `:35-39` — `to_csv` of 3 columns (`name`, `TargetGene`, `NumCandidateEnhGene`) for all
   `n_cand_pairs` rows, uncompressed. Confirmed on disk: output files are 433.6–487.6 MB
   (`igvf10_{multiome,scatac}/<cluster>/new_features/NumCandidateEnhGene.tsv`), identical byte-for-
   byte between modalities for the same cluster (e.g. `thp1_1`: 433,609,152 bytes in both).

**No nested loop, no cross join, and no interval/windowed overlap anywhere in this script** —
the per-gene ranking is a groupby+cumcount over an already-materialized pairs table, not a
gene×element scan. This distinguishes it from its sibling `gen_num_tss_enh_gene.py`, which *does*
run a `bedtools intersect` between the full predictions rows and a reference-TSS file (a genuine
interval join with a `k_genes`-scale second operand), and from `gen_num_sum_nearby_enhancers.py`,
which runs `bedtools slop`+`intersect` directly on `EnhancerList.txt` (`j_cand_elems` rows).
**Factual note for units #18/#19:** those two scripts are the ones whose code plausibly produces
an actual `j`/`k`-scale operation; this script (#17) does not, despite sharing the `pred_df =
pd.read_csv(...); pred_df["midpoint"] = ...` idiom with `gen_num_tss_enh_gene.py`.

## 3. Size variable(s)

**Primary: `n_cand_pairs`** (10,173,077–11,447,564 rows across the 6 clusters; identical between
multiome and scATAC for the same cluster, since the ABC `EnhancerPredictionsAllPutative.tsv.gz`
this script reads does not depend on the ARC/Kendall branch). This is a Tier 2 variable in the
briefing (range 1.13×), not the Tier 3 `j_cand_elems`/`k_genes` the unit brief names as "expected
drivers."

**Correcting the brief's expected-drivers note:** the brief for this unit states the expected
drivers are `j_cand_elems` and `k_genes`, "with the output joined onto `n_cand_pairs` rows." I
read the code and found no operation in this script over `j_cand_elems` or `k_genes` at all —
they are not inputs. The output row count equals `n_cand_pairs` not because of a join the script
performs, but because the script does zero filtering on an already-`n_cand_pairs`-row input (the
pairing of elements to genes happened upstream in ABC). `n_cand_pairs` is the one and only size
variable this script's control flow touches. I am flagging this discrepancy explicitly rather
than forcing the briefed variables onto code that doesn't use them — per the top-line instruction
to derive the class from the code, not from priors.

There is no cross term: this script does not multiply or join `n_cand_pairs` against any other
manifest variable.

## 4. Asymptotic class derived from code

**O(`n_cand_pairs` · log(`n_cand_pairs`))**, driven by three explicit `sort_values` calls over
all or half of the pairs table (`:10` full sort, `:18` ~half-table sort, `:34` full sort to
restore order — cited in §2 items 3, 5, 8). The two `groupby(...).cumcount()` calls (`:19`, `:26`)
and the `read_csv`/`to_csv` I/O are each O(`n_cand_pairs`) and do not change the class; the
repeated sorts are the asymptotic ceiling. There is no O(`n_cand_pairs`²) or O(`n_cand_pairs` ·
`k_genes`) term anywhere — I looked specifically (per the unit brief's instruction) for a nested
loop, cross join, or non-unique-key merge and found none; the groupby is hash-based (O(n)), not a
per-group nested scan.

## 5. Benchmark cross-check

All 12 measured jobs:

| Modality | Cluster | `n_cand_pairs` | wall s | cpu_time s | io_in MB | io_out MB |
|---|---|---:|---:|---:|---:|---:|
| multiome | jurkat | 11,036,283 | 97.41 | 27.38 | 383.40 | 16.19 |
| multiome | jurkat_pma_cd3_4hr | 10,854,698 | 45.03 | 14.75 | 284.02 | 15.79 |
| multiome | k562_crispri | 11,447,564 | 67.08 | 25.07 | 41.36 | 16.18 |
| multiome | telohaec_crispri | 11,106,687 | 95.22 | 30.29 | 55.94 | 15.87 |
| multiome | thp1_1 | 10,173,077 | 97.08 | 19.48 | 362.60 | 16.14 |
| multiome | thp1_2 | 10,278,325 | 83.12 | 31.53 | 316.38 | 60.17 |
| scATAC | jurkat | 11,036,283 | 84.88 | 32.19 | 368.99 | 0.00 |
| scATAC | jurkat_pma_cd3_4hr | 10,854,698 | 46.18 | 42.19 | 289.41 | 435.67 |
| scATAC | k562_crispri | 11,447,564 | 44.44 | 14.01 | 289.16 | 0.00 |
| scATAC | telohaec_crispri | 11,106,687 | 127.57 | 47.53 | 384.43 | 16.19 |
| scATAC | thp1_1 | 10,173,077 | 93.75 | 22.98 | 279.04 | 16.05 |
| scATAC | thp1_2 | 10,278,325 | 93.73 | 22.62 | 306.27 | 16.07 |

**The `n_cand_pairs` range across clusters is only 1.13× (10.17M–11.45M), and the 12 points do
not show any monotonic relationship with wall time or `cpu_time`.** Sorting all 12 rows by
`n_cand_pairs`: the smallest-`n` cluster (`thp1_1`, 10.17M) has wall times of 97.08 s (multiome)
and 93.75 s (scATAC) — near the top of the observed range — while the largest-`n` cluster
(`k562_crispri`, 11.45M) has wall times of 67.08 s and 44.44 s — near the bottom. This is the
opposite of what an increasing function of `n_cand_pairs` predicts, and it is not a subtle effect:
`k562_crispri`'s scATAC job (44.44 s) is the fastest of all 12 despite `k562_crispri` having the
*largest* `n_cand_pairs` of any cluster. **The data actively cannot confirm the O(n log n) class**
— consistent with the briefing's general warning that Tier 2 variables, spanning barely more than
Tier 3, are also unconfirmable from six (here, twelve) points, especially once other noise
sources are larger than the size effect.

`cpu_time`/wall ratios are low and inconsistent (0.20–0.42) across all 12 jobs, meaning 58–80% of
wall time is not active CPU work. Combined with the lack of size-correlation, this says the
measured wall time is dominated by I/O wait, filesystem/Slurm variance, and conda/interpreter
startup, not by the O(n log n) sort cost this size regime would predict if it dominated. At
`n_cand_pairs` ≈ 10–11M, pandas' vectorized C-level sorts are fast in absolute terms (likely low
single-digit seconds each); the dominant, noisy cost is almost certainly the ~280–490 MB of actual
I/O (see §6), not the algorithmic sorting/grouping work — i.e. **the code-derived class and the
measured behavior disagree at this scale**, and the likely explanation is I/O dominance, not that
the class is wrong.

## 6. Concrete overhead sources

- **I/O volume asymmetry between the benchmark's `io_out` and the real output size.** `io_out`
  in the benchmark table is only 15–60 MB per job, but the actual output file
  (`new_features/NumCandidateEnhGene.tsv`) is 433.6–487.6 MB on disk (verified with `ls -la`,
  identical between modalities for a given cluster, e.g. `thp1_1` = 433,609,152 bytes both). The
  psutil-based `io_out` counter evidently undercounts buffered writes that hadn't reached the
  block device by the time the job's I/O counters were sampled — the on-disk file size is the
  more trustworthy volume figure for this rule, and it means the write-side I/O is roughly
  **10–30× larger than the benchmark table's `io_out` column suggests.**
- **Reading only 6 of 29 columns does not save decompression/tokenization cost.** `usecols` at
  `:46-47` selects `name, chr, start, end, TargetGene, TargetGeneTSS` out of 29 total columns in
  `EnhancerPredictionsAllPutative.tsv.gz`, but the file must still be fully gzip-decompressed and
  every row fully tokenized by the C parser before unused columns are dropped — `usecols` saves
  memory and dtype-conversion work, not I/O or parse-scan time. `io_in` (41–385 MB, noisy) is
  consistent with reading the full compressed file (`n_cand_pairs_bytes`, 281.5–335.3 MB).
- **Output is uncompressed while input is gzip.** The script writes plain-text `.tsv` (`:35-39`,
  no `compression=` argument), so the 281–335 MB compressed input expands to a 433–488 MB
  uncompressed output — a ~1.5× byte-volume increase purely from dropping compression on write,
  independent of the 23-column reduction (29→3 columns in the output).
- **Three full/half-table sorts (§2, §4)** are the only non-linear cost in the script; at this
  `n_cand_pairs` scale they are likely fast in absolute wall-clock terms relative to I/O, but they
  are also the term that would eventually dominate at much larger `n` (see §4 — no other operation
  in the script is worse than O(n)).
- **No threading**: no `threads:` declared on the rule (`genomewide_features.smk:70-88`); pandas
  sorting/groupby here is single-threaded by default and nothing in the script parallelizes
  across chromosomes or genes.
- **No redundant re-reads**: the input file is read exactly once; there is no repeated pass over
  `EnhancerPredictionsAllPutative.tsv.gz` beyond the single `read_csv` call.

## 7. Three functions

Class per §4 is O(`n_cand_pairs`·log(`n_cand_pairs`)) from the code; the benchmarks (§5) show no
detectable scaling signal and are dominated by non-algorithmic noise at this `n` range. `t_low`
and `t_high` therefore bound the **class itself** — a linear floor (if the sorts are cheap enough
relative to I/O at this scale that the process behaves close to O(n)) and an n·log(n) ceiling
(the actual work the three `sort_values` calls impose) — rather than being point-estimate
uncertainty bands.

```
t_est(n_cand_pairs)  = 3.214e-7 * n_cand_pairs * log2(n_cand_pairs)     [seconds]
t_low(n_cand_pairs)  = 3.88e-6  * n_cand_pairs                          [seconds]
t_high(n_cand_pairs) = 4.909e-7 * n_cand_pairs * log2(n_cand_pairs)     [seconds]
```

Calibration (both modalities pooled, since §9 finds no class or constant difference between
them):
- `t_low`'s coefficient is the minimum observed wall-time/`n_cand_pairs` rate across the 12 jobs
  (scATAC `k562_crispri`: 44.44 s / 11,447,564 = 3.88e-6 s/row), used as a linear lower bound.
- `t_high`'s coefficient is fit so `t_high` passes through the maximum observed point (scATAC
  `telohaec_crispri`: 127.57 s at `n_cand_pairs`=11,106,687, log2(n)=23.402): 127.57 /
  (11,106,687 × 23.402) = 4.909e-7.
- `t_est`'s coefficient is fit through the pooled mean (mean wall time across all 12 jobs = 81.29
  s at mean `n_cand_pairs` = 10,816,106, log2(n)=23.375): 81.29 / (10,816,106 × 23.375) = 3.214e-7.

These are **not derived from an exponent fit to the 6/12 points** (§5 shows that would be fitting
noise) — they are calibrated only to pin the linear and n·log(n) curves to the observed
central/extreme wall-clock values, with the class itself (n vs n·log n) coming from the code
reading in §4, not from the data.

## 8. Evaluated at mean / min / max

| Point | cluster | `n_cand_pairs` | t_low | t_est | t_high |
|---|---|---:|---:|---:|---:|
| mean (across 6 clusters) | — | 10,816,106 | 42.0 s | 81.3 s | 124.1 s |
| min | `thp1_1` | 10,173,077 | 39.5 s | 76.1 s | 116.3 s |
| max | `telohaec_crispri` | 11,106,687 | 43.1 s | 82.2 s | 127.6 s |

(For reference, the actual measured wall times at these two named clusters: `thp1_1` 97.08 s
multiome / 93.75 s scATAC — both *above* `t_high` at the min point, underscoring that `thp1_1`
is not simply "the fast one" despite having the smallest `n_cand_pairs`; `telohaec_crispri`
95.22 s multiome / 127.57 s scATAC — the scATAC value sits almost exactly at `t_high`, the
multiome value well below it. This scatter, not a clean fit, is the expected and correctly
reported outcome given §5.)

## 9. Does the class differ between modalities?

**No — neither the class nor, as far as can be told, the constant differs.** `n_cand_pairs` is
identical between multiome and scATAC for the same cluster (verified: `size_manifest_multiome.tsv`
and `size_manifest_scatac.tsv` give the same six `n_cand_pairs` values), because
`EnhancerPredictionsAllPutative.tsv.gz` — the sole input to this script — is produced by the
shared ABC candidate-calling stack and does not depend on the ARC/Kendall branch that
`checkpoint features_required` decides. This rule is not among the briefing's four rules whose
cost differs materially by modality (`add_external_features`, `run_e2g_qnorm`,
`gen_final_features`, `get_stats_per_model_per_cluster`), and the benchmarks are consistent with
that: multiome and scATAC wall times for the same cluster differ by up to ~2× in either direction
(e.g. `k562_crispri` 67.08 s multiome vs 44.44 s scATAC; `telohaec_crispri` 95.22 s multiome vs
127.57 s scATAC) with no consistent sign — the pattern one would expect from run-to-run
scheduling/I/O noise on an identical workload, not from a modality-dependent algorithm.

## 10. What the measurement cannot tell us

- **Cannot confirm the O(n log n) class, or even distinguish it from O(n), from these 12 points.**
  `n_cand_pairs` spans only 1.13×, and the observed wall/cpu times are not monotonic in
  `n_cand_pairs` at all (§5) — the six-cluster (here twelve-job) design cannot separate an
  n·log(n) sort cost from a flat I/O-dominated cost at this scale. `t_low`/`t_high` in §7 bound
  the class as read from the code, not as confirmed by data.
- **Cannot apportion the ~58–80% of wall time that is not `cpu_time`** between Slurm scheduling
  latency, cold-cache Lustre reads of the ~280–380 MB compressed input, conda-environment
  shared-library loading, and Python/pandas interpreter startup — the benchmark instrumentation
  gives wall time, `cpu_time`, RSS, and coarse I/O counters, not a breakdown of wait states.
- **Cannot fully trust the `io_out` column for this rule** — it reports 15–60 MB while the actual
  written file is 433.6–487.6 MB (§6); some unknown fraction of the write is happening outside the
  window or mechanism psutil's counters capture (buffered pages not yet flushed, or a different
  I/O accounting path on the underlying Lustre filesystem). Any cost model that relies on
  `io_out` for this rule should use the on-disk output size instead.
- **These are single (n=1) measurements per cluster per modality** — 12 points total, not
  replicated — so no per-cluster variance estimate exists to separate "this job hit a slow
  node/cold cache" from "this cluster is systematically more expensive." The scatter described in
  §5 and §9 is consistent with pure measurement noise but cannot be proven to be only that.
- **Cannot say anything about behavior far outside the observed 10.17M–11.45M `n_cand_pairs`
  range** — e.g. a hypothetical cluster with 10× more candidate pairs. The O(n log n) ceiling in
  §4 is a code-derived asymptotic statement, not an extrapolation validated by any observed point
  beyond this narrow band.

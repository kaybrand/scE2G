# Phase 0 findings — briefing for Phase 1 agents

Read this before touching anything. It tells you where the measurements are, which size
variables exist and what they mean, and — most importantly — **what these six clusters can
and cannot establish**. Getting that last part wrong is the main way a Phase 1 report can be
confidently misleading.

Phase 0 is **complete for both modalities**. Multiome: 165 jobs, 30 rules, 0 errors. scATAC:
141 jobs, 26 rules, 0 errors. Correctness gate passed — the toy chr22 run reproduced the
committed expected outputs exactly, so the `benchmark:` instrumentation does not perturb
results.

---

## 1. Where the data is

Live (authoritative; `$SCRATCH` purges after 90 days of no content writes):

```
/scratch/users/kaybrand/scE2G_optimize_results/
├── igvf10_multiome/              # multiome_powerlaw_v3, 6 clusters
│   ├── benchmarks/<rule>/<wildcards joined by ~>.tsv     ← one file per JOB
│   ├── slurm_logs/rule_<rule>/<wildcards>/<jobid>.log
│   ├── qc_plots/predictions_qc_report.html
│   └── <cluster>/…                                        ← all intermediates (--notemp)
├── igvf10_scatac/                # scATAC_powerlaw_v3, same layout
└── analysis/
    ├── benchmarks_{multiome,scatac}.tsv      ← START HERE
    ├── size_manifest_{multiome,scatac}.tsv   ← the size variables
    ├── checksums_{multiome,scatac}.tsv
    └── dag_multiome.{dot,svg}, rulegraph_multiome.{dot,svg}
```

A committed backup of the multiome tables lives in the repo at `phase0_measurements/`.
**Use the `$SCRATCH` copies** — they are current.

Benchmark columns: `s` (wall seconds), `h:m:s`, `max_rss`/`max_vms`/`max_uss`/`max_pss` (MB),
`io_in`/`io_out` (MB), `mean_load`, `cpu_time` (seconds).

Intermediates are retained, so you can inspect the actual shape of your rule's inputs and
outputs rather than inferring them.

---

## 2. The two DAG topologies

The Multiome/scATAC branch is decided **per cluster**, at `checkpoint features_required`,
which reads the cluster's union feature table. Verified on disk: `to_generate == "ARC"` for
all six multiome clusters, `"Neither"` for all six scATAC clusters.

Exactly **four rules exist only in Multiome**:

```
make_kendall_pairs → generate_atac_matrix → compute_kendall → arc_e2g
```

Nothing exists in scATAC that is absent from Multiome. 165 − 141 = 24 jobs = 4 rules × 6
clusters.

**Four further rules run in both but cost materially different amounts**, because the ARC
feature is present or absent. This is a real signal, not noise — if your rule is in this list,
model the two modalities separately:

| Rule | Multiome (min) | scATAC (min) | ratio |
|---|---:|---:|---:|
| `add_external_features` | 12.6 – 21.3 | 0.9 – 1.3 | **~15×** |
| `run_e2g_qnorm` | 7.7 – 8.6 | 3.9 – 4.9 | ~2× |
| `gen_final_features` | 1.7 – 2.2 | 0.8 – 0.9 | ~2× |
| `get_stats_per_model_per_cluster` | 0.4 – 0.5 | 0.3 – 0.5 | ~1.2× |

---

## 3. The size-variable menu

One row per (modality, cluster) in `size_manifest_*.tsv`. Pick the variable(s) that actually
govern *your* script. Cross terms are expected and encouraged — do not force a single `n`.

### Tier 1 — genuinely variable (exponents are fittable here)

| Variable | Definition / source | Range | Ratio |
|---|---|---|---:|
| `n_frag` | ATAC fragments after chr filtering; `fragment_count.txt` | 26.2M – 394.8M | **15.1×** |
| `n_umi` | total RNA UMIs; `umi_count.txt` (multiome only) | 14.5M – 174.2M | **12.0×** |
| `n_cells` | cells in cluster; `cell_count.txt` (multiome only) | 2,593 – 15,555 | **6.0×** |

**Caveat:** these three are strongly correlated — deeper datasets have more of everything.
With six points you generally *cannot* separate them. Claiming O(`n_frag`) over O(`n_umi`)
from benchmarks alone is overreaching; say which one the *code* uses.

### Tier 2 — mildly variable

| Variable | Definition / source | Range | Ratio |
|---|---|---|---:|
| `n_peaks` | raw MACS2 peaks, pre-cap; `macs2_peaks.narrowPeak.sorted` | 432,039 – 595,372 | 1.38× |
| `n_kendall_pairs` | candidate (peak, gene) pairs to correlate; `Kendall/Pairs.tsv.gz` | 9.69M – 11.32M | 1.17× |
| `n_cand_pairs` | candidate enhancer–gene pairs; `EnhancerPredictionsAllPutative.tsv.gz` | 10.17M – 11.45M | 1.13× |

### Tier 3 — near-constant, but DO NOT dismiss as constant overhead

| Variable | Definition / source | Range |
|---|---|---|
| `j_cand_elems` | candidate elements post-cap; header-free `*.candidateRegions.bed` | 156,420 – 158,548 |
| `k_genes` | genes considered; `Neighborhoods/GeneList.txt` | 20,532 (constant) |
| `n_called_all` | all rows of the thresholded predictions | 59,438 – 61,284 |
| `n_called_distal` | same, excluding `class == "promoter"` | 46,155 – 47,560 |
| `n_feat_cols` | columns of `genomewide_features.tsv.gz` | 24 (constant) |

`k_genes` is constant *by construction* — it is the human gene annotation. `j_cand_elems` is
near-constant because the candidate set is capped by `nStrongestPeaks`.

**This does not make them unimportant.** The complexity class *in* `j_cand_elems` still
matters enormously: an O(j²) inner loop over ~157,000 elements is ~2.5×10¹⁰ operations, while
O(j log j) is ~2.8×10⁶ — four orders of magnitude, and exactly the kind of finding that
justifies rewriting a script. Derive the class in `j` from the code and report it. Separately
note that the six clusters cannot *empirically confirm* it, because `j` spans only 1.01×.
Those are two different claims and both belong in your report.

### Collapsed variables — do not go looking for these

These row counts are **identities**, because the stages between them are pure column-additions
with no filtering. They were collapsed to avoid an agent "discovering" a relationship that is
tautological:

```
n_cand_pairs    == rows of ARC == rows of ActivityOnly == rows of genomewide_features
n_kendall_pairs == rows of Pairs.Kendall.tsv.gz
j_cand_elems    == rows of EnhancerList.txt − 1   (header offset only)
```

`compute_kendall` adds 4 columns and drops zero rows. The `size_manifest` carries
`row_identities_ok`, re-verified each run, so if a future dataset breaks one we find out.

**What is NOT redundant is bytes.** Identical row counts, very different widths (thp1_1):

| Table | cols | bytes |
|---|---:|---:|
| `ActivityOnly_features.tsv.gz` | 19 | 203 MB |
| `EnhancerPredictionsAllPutative.tsv.gz` | 29 | 282 MB |
| `genomewide_features.tsv.gz` | 24 | 495 MB |
| `EnhancerPredictionsAllPutative_ARC.tsv.gz` | 36 | **630 MB** |

A 3.1× I/O spread at constant `n`. If your rule is I/O-shaped, use the `*_bytes` columns
(`arc_bytes`, `actonly_bytes`, `feat_bytes`, `kendall_scored_bytes`, `enh_list_bytes`), not the
row count.

### One naming trap

"Pairs" and "links" are different things and were conflated once already:

```
n_cand_pairs     ~10.8M   candidate enhancer–gene pairs
      ↓  model scoring + threshold          ~180× reduction
n_called_all      ~60.4k  every row of the thresholded file  (COST-relevant)
      ↓  drop class == "promoter"           −23%
n_called_distal   ~46.7k  called links      (what the QC report shows)
```

The QC report's link count comes from `get_stats_per_cluster.R:53`, which filters
`class != "promoter"`. Use `n_called_all` for cost (rules read the whole file),
`n_called_distal` when you mean biological links.

---

## 4. How to use the benchmarks — methodological rules

1. **Derive the complexity class from the CODE first.** The six benchmark points are a
   *secondary* check on the exponent, never the primary source.
2. **Six points, spanning 1.01×–1.17× for most variables, cannot fit an exponent.** If your
   rule is driven by a Tier 3 variable, say plainly: "class derived from code; benchmarks
   cannot confirm because the variable spans only N×." Do not fit to noise.
3. **Flag disagreement, don't smooth it.** If the code says O(n²) and the timings look linear,
   report the conflict and hypothesise why (caching, I/O dominance, a constant that dominates
   at this scale).
4. **Exclude the toy chr22 run entirely.** At 49 KB of fragments — ~85,000× smaller than
   `telohaec_crispri` — it measures interpreter and conda startup, not algorithm. It exists
   only as a correctness gate.
5. **Wall time includes Slurm scheduling latency.** For sub-minute rules this dominates. Use
   `cpu_time` to separate compute from waiting.
6. **These are single measurements, not replicates** — except for the shared ABC/fragment
   stack, which was measured independently in both runs. Where a rule appears in both
   modalities with the same inputs, you have n=2 per cluster; use it.

---

## 5. Cross-cutting findings you should know before starting

### Declared threads are largely fictional

| Rule | declared `threads:` | measured CPU/wall |
|---|---:|---:|
| `frag_to_norm_bigWig` | 16 | **0.99** |
| `frag_to_tagAlign` | 8 | **0.92** |
| `process_fragment_file` | 8 | **0.94** |

All three reserve 8–16 CPUs from Slurm and use approximately **one**. Every ABC rule declares
no `threads:` at all and therefore runs at 1. `compute_kendall` requests `config["threads"]`,
currently **1**, passed through to `compute_kendall.py --threads`. If your rule is in this
list, the gap between declared and realised parallelism is a first-class finding — identify
*why* (pipe serialisation, a tool that ignores its thread flag, a `sort --parallel` that is
not the bottleneck).

### The current cost leaders

Slowest single jobs (multiome), and note the two profiles are completely different in shape:

| Rule | Cluster | wall min | peak GB | CPU/wall |
|---|---|---:|---:|---:|
| `create_neighborhoods` | telohaec_crispri | 112.1 | 0.28 | 0.98 |
| `frag_to_norm_bigWig` | telohaec_crispri | 81.4 | **24.53** | 0.99 |
| `call_macs_peaks` | telohaec_crispri | 48.4 | 13.29 | 0.96 |
| `add_external_features` | jurkat | 21.3 | 20.48 | 0.97 |
| `frag_to_tagAlign` | telohaec_crispri | 25.9 | 4.65 | 0.92 |

`create_neighborhoods` is **CPU-bound on a single core** (290 MB RSS, ~1 GB read, 97.7% of one
core for two hours) — not I/O-bound as the low memory might suggest.

### `compute_kendall` is no longer the bottleneck

14.5 min on the largest cluster; absent from the top twelve. This reflects optimization work
already merged. Its declared 63 GB memory floor and 12-hour runtime allowance are now heavily
oversized. Do not assume it dominates.

### `frag_to_norm_bigWig` may be off the critical path entirely

It is the largest total consumer (174.9 min summed) but exists only because
`make_IGV_tracks: True`, and **nothing downstream consumes its output**. Characterise it
normally; Phase 2 decides whether it matters. Total time ≠ critical-path time.

### Measurement artifacts to ignore

- Sub-second rules (`features_required`, `basic_features_required`, `save_reference_configs`,
  `make_external_features_config`) report `cpu_time/s` ratios of 3–14. That is sampling noise
  on jobs too short to measure, not real parallelism. Treat them as ~0 cost.
- `hover_plots` reports `cpu_time` 0.00 — not captured for the Rmd render. Use wall time.
- `abc_generate_chrom_sizes_bed_file` has **no benchmark by design**: Snakemake's psutil
  monitor crashes on it (the awk exits before psutil can attach). ~0.03 s; ignore.

---

## 6. What each agent must deliver

For your script unit:

1. **Size variable(s)**, named from the manifest above. Cross terms where the code has them
   (e.g. `O(n_cand_pairs · k_genes)`). Do not force a single `n`.
2. **Asymptotic class, derived primarily from reading the code**, with the specific loops,
   joins, sorts or merges that justify it. Cite file:line.
3. **A benchmark cross-check** — does the observed scaling agree? State explicitly if the
   variable's range is too narrow to tell, and flag any disagreement.
4. **Concrete overhead sources**: I/O volume, serialisation/deserialisation, redundant writes,
   declared-vs-realised threading, repeated passes over the same data.
5. **Three functions of your variable(s)**: `t_est(n)`, `t_low(n)`, `t_high(n)` — best-guess
   class calibrated to the measurements, plus bounds on the *class itself*.
   **Do not blend them.** PERT's `(O + 4M + P)/6` assumes point-estimate human uncertainty
   about a fixed task; ours is uncertainty about an exponent, and averaging an O(n) with an
   O(n²) guess produces a number wrong at every scale. Phase 2 needs the functions.
6. All three evaluated at the mean cluster size, **for reporting convenience only** — the
   functions are the deliverable.
7. If the rule appears in both modalities, say whether the class differs or only the constant.

## 7. Out of scope for Phase 1

Characterisation only. **No code edits, no fix proposals, no prioritisation, no speculation
about what should be optimized first.** Phase 2 computes the critical path and decides what
matters; a plausible-sounding recommendation made without the path arithmetic is exactly the
error this whole exercise exists to avoid.

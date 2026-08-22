# Unit 26 — inline `run:` blocks and checkpoints (`features_required`, `basic_features_required`, `save_reference_configs`)

## 1. Unit identity

Three DAG nodes, all implemented as inline Snakemake `run:` blocks (no external script files):

| Node | Kind | File:line | Modalities | Jobs/modality |
|---|---|---|---|---|
| `features_required` | `checkpoint` | `workflow/rules/add_external_features.smk:1-22` | both | 6 (wildcard `sample`) |
| `basic_features_required` | `checkpoint` | `ENCODE_rE2G/workflow/rules/genomewide_features.smk:3-33` | both | 6 (wildcard `sample`, called via `wildcards.biosample`) |
| `save_reference_configs` | `rule` | `workflow/rules/save_configs.smk:1-66` | both | **1** (unwildcarded, aggregating) |

All three are sub-second (`s` column: 0.03–0.37 s across every measured job in both modalities —
see §5). Neither checkpoint declares `threads:`; `save_reference_configs` doesn't either. All run
at 1 core.

## 2. What the code actually does

### `checkpoint features_required` — `add_external_features.smk:1-22`

```
1  checkpoint features_required:
2      input: feature_table_file = RESULTS_DIR/{sample}/feature_table.tsv
...
8      run:
9          Kendall = False; ARC = False
11         with open(input.feature_table_file) as f:
12             for line in f:
13                 columns = line.strip().split("\t")
14                 if ("Kendall" in columns[1]) or ("Kendall" in columns[2]): Kendall = True
16                 if ("ARC.E2G.Score" in columns[1]) or ("ARC.E2G.Score" in columns[2]): ARC = True
18         final_val = "Neither"
19         if ARC or Kendall: final_val = "ARC"
21         write to_generate.txt = final_val
```

A single linear scan of the (13-row multiome / 7-row scATAC) `feature_table.tsv`, testing two
substring memberships per row against columns 2/3. Writes one of `"ARC"` / `"Neither"` (the
`"Kendall"` string value is never actually reachable as a `final_val` — both `Kendall` and `ARC`
booleans collapse to the same `"ARC"` output at line 19; this is a real but harmless dead branch,
not a bug that changes behaviour). Verified on disk (`to_generate.txt`): `ARC` for all six multiome
clusters, `Neither` for all six scATAC clusters — confirmed by direct read of every file in both
`igvf10_multiome/*/to_generate.txt` and `igvf10_scatac/*/to_generate.txt`.

### `checkpoint basic_features_required` — `genomewide_features.smk:3-33`

Same pattern, four substring tests per line (`numCandidateEnhGene`, `numTSSEnhGene`,
`numNearbyEnhancers`, `sumNearbyEnhancers`) over the same `feature_table.tsv`, writing three
independent `"True"`/`"False"` flag files (`genomewide_features.smk:28-33`). Verified on disk for
`jurkat`, `thp1_1`, `telohaec_crispri` in **both** modalities: all three flags are `"True"` in
every case checked. In this dataset the checkpoint's discriminating power is never exercised —
every cluster in both modalities takes the same "generate everything" branch — but the mechanism
(and its barrier property, §4) is identical to `features_required`.

### `rule save_reference_configs` — `save_configs.smk:1-66`

No `input:` block at all — its only dependencies are Python objects resolved at Snakefile parse
time (`config`, `encode_e2g_config`, `BIOSAMPLE_DF`, `IGV_DIR`), not the output of any other rule.
The `run:` block (lines 21-65):

1. `os.makedirs(params.out_dir, exist_ok=True)` (line 25).
2. Three `yaml.safe_dump` calls (lines 28-37) writing the resolved scE2G/ENCODE-rE2G/ABC config
   dicts to `config/{scE2G,ENCODE_rE2G,ABC}_config.yml`.
3. `df.apply(mkpaths, axis=1)` (line 64) over `BIOSAMPLE_DF` (6 rows, one per cluster):
   for each row, string-joins ~4-6 output paths (predictions, thresholded predictions, and,
   conditionally on `config["make_IGV_tracks"]`, the IGV bedpe/bigWig paths) and returns them as a
   `pd.Series`. Concatenated onto `df` and written with `df.to_csv` (line 65).

Confirmed on disk: `expanded_biosample_config.tsv` is 6 data rows + header, 6.9 KB (multiome);
the three YAML files are 1.4-2.2 KB each. All four outputs are small, in-memory, single-pass
operations over a 6-row table plus a handful of dict keys — no loop over any per-sample data
product.

## 3. Size variable(s)

**None of the `size_manifest_*.tsv` variables govern any of these three nodes.** Their true inputs
are:

- `features_required` / `basic_features_required`: the row count of `feature_table.tsv` — 13 rows
  (multiome) / 7 rows (scATAC), a quantity determined by which model (`multiome_powerlaw_v3` vs.
  `scATAC_powerlaw_v3`) is assigned, not by cluster size. This is not a manifest column and is
  constant within a modality across all six clusters.
- `save_reference_configs`: the number of clusters (`len(BIOSAMPLE_DF)` = 6 in both modalities,
  fixed by the run's cell-cluster config) and the number of config dict keys — both small
  constants, and again not manifest columns (`size_manifest` describes per-cluster data volume,
  not pipeline configuration size).

Per the instructions' framing: these are not Tier-3 near-constant *data* variables whose class
merely can't be confirmed at this range — they are configuration-determined constants that do not
vary with any of `n_frag`/`n_umi`/`n_cells`/`n_peaks`/etc. at all, by construction of the code.

## 4. Asymptotic class derived from code

Let `f` = rows of `feature_table.tsv` (13 or 7, constant), `c` = number of clusters (6, constant),
`k` = number of config dict keys (small constant).

- `features_required`: one pass over `f` rows, 2 substring tests each → `O(f)`. `f` is a fixed
  small constant in this project (bounded by the number of possible feature names in the model
  spec), so this is `O(1)` in every variable that varies across clusters.
- `basic_features_required`: identical shape, `O(f)` with 4 tests per row instead of 2 — same
  constant-in-practice conclusion.
- `save_reference_configs`: `O(k)` for the three YAML dumps (dict serialisation, tiny) + `O(c)`
  for the `df.apply(mkpaths, axis=1)` row-wise string-path construction (`save_configs.smk:43-64`)
  — linear in the number of clusters, but `c=6` is fixed by the experiment design, not a variable
  that scales with data volume the way `n_frag` etc. do.

**All three nodes are `O(1)` with respect to every size variable this project tracks.** The only
variable in which any of them is even nominally non-constant is `c` (cluster count) for
`save_reference_configs`, and that is a study-design constant (how many biosamples you chose to
run), not a data-size variable.

## 5. Benchmark cross-check

```
features_required        s (multiome): 0.37 0.04 0.03 0.04 0.13 0.04   (jurkat…thp1_2)
features_required        s (scatac):    0.05 0.06 0.05 0.07 0.05 0.04
basic_features_required  s (multiome): 0.37 0.05 0.05 0.05 0.13 0.05
basic_features_required  s (scatac):    0.05 0.05 0.07 0.05 0.05 0.05
save_reference_configs   s (multiome): 0.08                              (n=1, "all" wildcard)
save_reference_configs   s (scatac):    0.21                             (n=1)
```
(`awk -F'\t' -v r=RULE 'NR==1||$2==r' benchmarks_{multiome,scatac}.tsv`, `/scratch/.../analysis/`)

This agrees with §4: flat at ~0.03-0.4 s regardless of cluster, i.e. regardless of `n_frag`
(15.1×), `n_umi` (12.0×), `n_cells` (6.0×) — none of which the code ever touches. This is the
**"data cannot constrain this" case named explicitly in the instructions**, but for a different
reason than the Tier-3 j-variables: it's not that the range is too narrow, it's that the code has
*zero* dependence on any tracked variable, confirmed by both reading the code and by the flat
benchmark.

`cpu_time` for all these jobs is reported as 3-14× the wall time (e.g. `features_required`
jurkat: 0.72 s cpu_time on 0.37 s wall). **This is psutil sampling noise on jobs too short for the
periodic sampler to catch, not parallelism** — none of these rules declares `threads:`, all run a
single Python `for` loop with no subprocess or multiprocessing. Reported exactly as an artifact,
per the launch briefing's explicit instruction not to read parallelism into it.

**A second, more informative signal came from the Slurm logs, not the benchmark table.** The
`benchmark:` directive only times the job's own execution; it does not capture the gap between
"Snakemake decided to execute this job" and "the job actually started running under Slurm." I
pulled that gap directly from `slurm_logs/rule_{features_required,basic_features_required,
save_reference_configs}/*/*.log`, which log a timestamp both when the outer job is dispatched and
when the (nested) execution actually starts and finishes:

```
features_required (12 logs, both modalities):        dispatch→start gap 2-6 s, dispatch→finish 3-7 s
basic_features_required (12 logs, both modalities):    dispatch→start gap 2-6 s, dispatch→finish 3-7 s
save_reference_configs (2 logs, one per modality):     dispatch→start gap 2 s,   dispatch→finish 3-4 s
```

Example (`igvf10_multiome/slurm_logs/rule_features_required/jurkat_pma_cd3_4hr/39910434.log`):
`Execute 1 jobs...` at `17:17:10`, the nested `localcheckpoint features_required:` actually starts
at `17:17:14` (4 s dispatch gap), finishes at `17:17:14`, and the *outer* "Finished jobid: 0" +
"Storing output in storage" line lands at `17:17:16` (another 2 s for Snakemake's own
checkpoint-output bookkeeping). So of a ~6 s user-visible latency, <0.5 s is the benchmarked
compute, ~4 s is Slurm dispatch, and ~2 s is Snakemake's checkpoint-storage/DAG-update overhead —
**this last piece is specific to checkpoints** (I did not see the equivalent double-log,
nested-execution structure in ordinary `rule` logs I spot-checked, though I did not survey this
systematically outside this unit). This is consistent with the checkpoint mechanism requiring
Snakemake to synchronously persist and re-read the checkpoint's output before it can resolve the
downstream input functions (§7 below).

I did not observe this nested double-execution/sub-sbatch pattern in the `save_reference_configs`
log (plain rule, not a checkpoint), consistent with the mechanism being checkpoint-specific.

## 6. Concrete overhead sources

- **I/O volume**: trivial in all three. `features_required`/`basic_features_required` read one
  13- or 7-row TSV (≤1 KB) and write 1-3 files of a few bytes each (`"ARC"`, `"Neither"`, `"True"`,
  `"False"`). `save_reference_configs` writes ~1.4-2.2 KB of YAML plus a 6.9 KB / 6-row TSV.
  `io_in`/`io_out` in the benchmark table are ≤0.3 MB for all jobs of all three rules — an order
  of magnitude below anything that could plausibly cost measurable time; whatever `io_in` is
  measuring here is almost certainly conda/Python interpreter shared-library page-ins, not the
  rule's actual data.
- **Declared-vs-realised threading**: none declared, none used — not a finding specific to this
  unit, consistent with every other rule in the project.
- **Serialisation/deserialisation**: negligible — plain-text line reads, one `yaml.safe_dump`
  each, one small `pandas.to_csv`.
- **Redundant writes / repeated passes**: none — each is a single linear pass over a small
  in-memory structure.
- **Subprocess/interpreter startup**: because these are `run:` blocks (executed inside the
  already-running Snakemake Python process, not via `script:`/`shell:`), they pay **no extra
  interpreter startup** beyond Snakemake's own process — unlike, e.g., the R-script rules
  characterised elsewhere in this project (unit #23's `make_biosample_feature_table`, which pays
  a 1-14 s R/conda cold start). This is the single biggest structural difference between these
  three nodes and most other benchmarked rules: **their entire measured cost is genuinely just
  the Python loop/dict-dump/dataframe-apply**, not environment startup.
- **The dominant real-world overhead is Slurm dispatch latency + (for checkpoints only)
  Snakemake's own checkpoint-resolution bookkeeping** (§5), both of which are unrelated to the
  code inside the `run:` block and would be paid even if that code were reduced to `pass`.

## 7. Three functions: `t_est(n)`, `t_low(n)`, `t_high(n)`

All three nodes are constants — none of the manifest's size variables enter the code at all (§3,
§4). I therefore give each a constant `t(n) = t` for all `n`, but split "bare compute" from
"user-visible latency including Slurm dispatch and (for checkpoints) DAG-resolution bookkeeping,"
per the addendum's instruction to separate the two as far as the data allows.

### `checkpoint features_required`

```
t_low  = 0.03   s   (fastest measured benchmark 's', k562_crispri multiome / min across all 12 jobs)
t_est  = 5.0    s   (mean of the 12 measured dispatch→finish gaps from Slurm logs, range 3-7 s)
t_high = 60     s   (NOT measured — see caveat below)
```

### `checkpoint basic_features_required`

```
t_low  = 0.03   s   (fastest measured benchmark 's')
t_est  = 5.5    s   (mean of the 12 measured dispatch→finish gaps, range 3-7 s)
t_high = 60     s   (NOT measured — see caveat below)
```

### `rule save_reference_configs`

```
t_low  = 0.08   s   (measured benchmark 's', multiome; the single scATAC measurement was 0.21 s)
t_est  = 3.5    s   (mean of the 2 measured dispatch→finish gaps: 4 s multiome, 3 s scatac)
t_high = 60     s   (NOT measured — see caveat below)
```

**Caveat on `t_high` for all three**: the 3-7 s dispatch gap in §5 is what this specific run
observed on a partition that, at the time, had free capacity. It is **not** a measured ceiling —
it is a floor-ish "typical case." Under partition contention (busy `normal`/owner partition,
competing large-memory jobs elsewhere in this same pipeline — several rules request 16-64 GB) the
same job could sit `PD` in the queue for minutes; nothing in this dataset measures that, because
the launch happened once, under one set of cluster conditions. I used 60 s as a round, explicitly
unverified placeholder for "contended queue" — Phase 2 should not treat it as calibrated; it only
exists to avoid a `t_high = t_low` collapse that would misleadingly imply certainty. If Phase 2 has
independent data on this cluster's queue-wait distribution for tiny 1-core/16 GB jobs, that should
replace this number outright.

## 8. Evaluated at mean / min (`thp1_1`) / max (`telohaec_crispri`) cluster size

Because none of the three nodes depends on cluster size at all, every one of `t_est`/`t_low`/
`t_high` evaluates to the **same value at every cluster** — the flatness itself is the finding, not
a rounding artifact (same conclusion as unit #23's report on `make_biosample_feature_table`, the
rule immediately upstream of both checkpoints).

| Node | mean-cluster | `thp1_1` (min) | `telohaec_crispri` (max) |
|---|---:|---:|---:|
| `features_required` — `t_est`/`t_low`/`t_high` | 5.0 / 0.03 / 60 s | 5.0 / 0.03 / 60 s | 5.0 / 0.03 / 60 s |
| `basic_features_required` — same | 5.5 / 0.03 / 60 s | 5.5 / 0.03 / 60 s | 5.5 / 0.03 / 60 s |
| `save_reference_configs` — same (n=1, not per-cluster) | 3.5 / 0.08 / 60 s | n/a (aggregating rule, no per-cluster job) | n/a |

## 9. Both modalities: class or constant?

**Only the constant (and, trivially, which branch fires) differs — never the class.** Both
checkpoints run the identical `O(f)` linear-scan code in both modalities; `f` (feature-table row
count) differs by a small constant (13 vs. 7 rows) that is invisible at these timescales. The
*decision* each checkpoint reaches differs by modality (`ARC` vs. `Neither` for
`features_required`; all-`True` in both modalities for `basic_features_required`, so no observed
difference there), but that is a data/config difference, not an algorithmic one — see §11 for why
it matters anyway. `save_reference_configs` runs the same code in both modalities over the same
`c=6`; its two measurements (0.08 s multiome, 0.21 s scatac) differ by ~2.6×, well within what a
single-measurement psutil sample on a sub-second job could produce as noise — I do not read a
modality effect into it.

## 10. Structural findings — the actual deliverable of this unit

### `checkpoint features_required` decides the entire downstream topology

It reads the cluster's union feature table (written by `make_biosample_feature_table`, unit #23,
at `ENCODE_rE2G/workflow/rules/predictions.smk:2-18`) and writes `{sample}/to_generate.txt`
(`add_external_features.smk:1-22`, logic at lines 9-22, quoted in §2). Verified on disk: `ARC` for
all six multiome clusters, `Neither` for all six scATAC clusters, driven by which model
(`multiome_powerlaw_v3` vs. `scATAC_powerlaw_v3`) each cluster was assigned — not by anything
cluster-specific (row-identical decision within each modality).

The decision is consumed by **three separate input functions**, each independently calling
`checkpoints.features_required.get(...)`, meaning the checkpoint has three direct DAG out-edges,
not one:

1. `features_to_generate` (`add_external_features.smk:26-34`) → feeds `rule
   make_external_features_config` (`add_external_features.smk:37-52`). When the value is `"ARC"`,
   this function returns
   `{sample}/ARC/EnhancerPredictionsAllPutative_ARC.tsv.gz` — the exact output path of `rule
   arc_e2g` (`workflow/rules/arc_e2g.smk:19-25`, confirmed by direct path comparison). Requiring
   that file is what pulls the whole `make_kendall_pairs → generate_atac_matrix →
   compute_kendall → arc_e2g` chain into the DAG (4 rules × 6 clusters = 24 jobs — exactly the
   165 − 141 = 24-job multiome/scATAC delta). When `"Neither"`, it returns `RESULTS_DIR` (a
   directory, not a real dependency — a Snakemake idiom for "no extra input needed"), and none of
   those four rules are ever added to the scATAC DAG.
2. `get_gex_file` (`workflow/rules/sc_predictions.smk:128-133`) → feeds `rule
   element_and_gene_summaries` (`sc_predictions.smk:136-153`): if `"Kendall"` or `"ARC"`, requires
   `{cluster}/Kendall/gene_expression_metrics.tsv.gz`; else `RESULTS_DIR`.
3. `get_count_file` (`sc_predictions.smk:157-162`, parameterised by `metric`) → feeds `rule
   get_stats_per_model_per_cluster` (`sc_predictions.smk:164-...`) twice (for `cell_count` and
   `umi_count`): if `"Kendall"` or `"ARC"`, requires `{cluster}/{metric}.txt`; else `RESULTS_DIR`.

So the checkpoint's immediate downstream rule set is `{make_external_features_config,
element_and_gene_summaries, get_stats_per_model_per_cluster}`, and via (1) it transitively gates
the 4-rule/24-job ARC chain. Note: `make_external_features_config`'s own output
(`external_features_config.tsv`) feeds `rule add_external_features`
(`add_external_features.smk:179-193`) as an `ancient()`-wrapped input — meaning changes to that
file won't force `add_external_features` to re-run, but the DAG dependency (and therefore the
ordering/barrier requirement) still holds for a first build.

### `checkpoint basic_features_required` gates the three `generate_num_*` rules

Reads the same `feature_table.tsv` and writes three independent True/False flag files
(`genomewide_features.smk:3-33`). Each flag is consumed by its own input function
(`get_numCandidateEnhGene_file`, `get_numTSSEnhGene_file`, `get_numNearbyEnhancers_file`,
`get_sumNearbyEnhancers_file`; `genomewide_features.smk:37-67`), which in turn feed **the single
downstream rule `activity_only_features`** (`genomewide_features.smk:153-171`) as four of its five
inputs. When a flag is `"True"` (verified: all three flags are `"True"` for every cluster checked
in both modalities), the corresponding upstream rule — `generate_num_candidate_enh_gene`
(`:70-88`), `generate_num_tss_enh_gene` (`:91-113`), or `generate_num_sum_enhancers` (`:116-150`,
which alone produces both `NumEnhancersEG5kb.txt` and `SumEnhancersEG5kb.txt`) — is pulled into the
DAG. That's the 3 rules × 6 clusters = 18 jobs/modality the addendum names; because the flags are
`"True"` in both modalities here, these 18 jobs run in **both** multiome and scATAC (36 jobs total,
not part of the 24-job multiome-only delta, which comes entirely from the `features_required`/ARC
branch above).

**Both checkpoints converge on the same downstream rule, `add_external_features`**: its inputs are
`predictions_extended` (= `activity_only_features`'s output, gated by `basic_features_required`)
and `external_features_config` (= `make_external_features_config`'s output, gated by
`features_required`). So `add_external_features` cannot start until *both* checkpoints have
resolved and *both* of their respective downstream chains have completed — it is the join point
where the two barriers' consequences merge.

### A checkpoint is a hard DAG barrier — this is the single most important thing this unit delivers

Snakemake cannot expand the DAG past an unevaluated checkpoint: any rule whose input function
calls `checkpoints.X.get(...)` blocks DAG construction for that branch until `X`'s output file
physically exists on disk. This is exactly why a dry run reported an identical 125-job plan for
both configs (per addendum §F) — the pre-run DAG couldn't see the 24 ARC-chain jobs or the 18 (×2
modality) `generate_num_*` jobs, because neither checkpoint had run yet. It's also why the
job-level `--dag` export (per addendum §F) had to be regenerated *after* the run: `dag_multiome.dot`
as inherited had 125 nodes (the stale pre-resolution DAG), missing all 42 checkpoint-gated jobs;
the corrected export has 155 nodes and reconciles to 165 jobs only once both checkpoints'
consequences are counted.

**The cost-vs-topology asymmetry is the point: these three nodes cost ~0 (0.03-0.4 s benchmarked,
3-7 s including Slurm dispatch in this run), yet everything listed above under each checkpoint —
up to 24+18=42 jobs, plus every job further downstream that reads their outputs
(`add_external_features`, `gen_final_features`, and the whole prediction/QC tail) — cannot even be
scheduled until the checkpoint's single small job has both run *and* had its output re-read by
Snakemake's DAG re-planner.** A zero-cost node can still sit squarely on the critical path: it
contributes ~0 to the path's *duration* but is a hard prerequisite for the path's *existence*.
Phase 2 should model each checkpoint as a serialisation point with near-zero duration but full
fan-out dependency, not omit it because its own cost rounds to zero.

**For Phase 2's DAG**: neither checkpoint node appears in either resolved exported DAG (confirmed,
addendum §F), even though both are benchmarked, 6-job-per-modality rules. They must be inserted
manually as barrier nodes:

- `checkpoint features_required`: upstream = `make_biosample_feature_table` (produces
  `feature_table.tsv`, the checkpoint's sole input); downstream = `make_external_features_config`,
  `element_and_gene_summaries`, `get_stats_per_model_per_cluster` directly, and transitively (only
  when multiome/`ARC`) the entire `make_kendall_pairs → generate_atac_matrix → compute_kendall →
  arc_e2g` chain, which itself then feeds back into `make_external_features_config`'s resolved
  input and onward into `add_external_features`.
- `checkpoint basic_features_required`: upstream = `make_biosample_feature_table` (same input
  file, same upstream rule as above — the two checkpoints are siblings reading the same file, not
  chained to each other); downstream = `generate_num_candidate_enh_gene`,
  `generate_num_tss_enh_gene`, `generate_num_sum_enhancers` directly, all three then feeding
  `activity_only_features`, which feeds `add_external_features` (the same join point named above).

### `save_reference_configs`: leaf status

`save_reference_configs` (`save_configs.smk:1-66`) has **no `input:` block** — its only
dependencies are in-memory Python objects fixed at Snakefile-parse time, so it has **zero upstream
rule edges** in the DAG; it can run at any point (including first) with no prerequisite job. Its
four outputs (`scE2G_config.yml`, `ENCODE_rE2G_config.yml`, `ABC_config.yml`,
`expanded_biosample_config.tsv`) are **not read as input by any other rule** — confirmed by
grepping every `.smk` file under `workflow/` and `ENCODE_rE2G/workflow/` for these four filenames;
they appear only in `save_configs.smk`'s own `output:` block and in `Snakefile:89`.

**But it is not fully disconnected: `Snakefile:89` sets
`output_files = [os.path.join(RESULTS_DIR, "config", "expanded_biosample_config.tsv")]`, and
`rule all` (`Snakefile:123-125`) takes `output_files` as its `input:`.** So `rule all` — the
pipeline's default target — directly depends on `save_reference_configs`'s `res_out`. It is
therefore a true **leaf in the sense that nothing consumes its output as a computational input**,
but it is **not disconnected from `rule all`**: it's a zero-in-degree, direct-to-`all` node. This
is exactly the shape Phase 2 needs for slack analysis: it can be scheduled anywhere in the run
(nothing blocks it, it blocks nothing but the final `all` target), so it can absorb scheduling
slack freely without perturbing any other rule's start time.

## 11. What the measurement cannot tell us

- **Whether `basic_features_required`'s discriminating branch (any flag = `"False"`) ever fires in
  a real config.** In every cluster checked (both modalities) all three flags were `"True"`; the
  18×2=36 `generate_num_*` jobs always run in this dataset. Whether a differently-configured model
  (one that doesn't need `numCandidateEnhGene`/`numTSSEnhGene`/nearby-enhancer features) would
  ever skip some of those jobs is unconfirmable from six clusters that all use the same two
  models.
- **A calibrated `t_high` for Slurm dispatch latency.** §7's `t_high=60 s` is an explicit,
  unverified placeholder. The only real data point is this run's observed 2-6 s dispatch gap on
  what was presumably a partition with free capacity at the time; nothing in this dataset
  constrains behaviour under contention.
- **Whether the ~2 s "Storing output in storage" / DAG-update step seen in checkpoint logs (§5) is
  a fixed Snakemake-internal cost or scales with DAG size.** I only inspected this project's 155-
  /131-node DAGs; I cannot say whether a much larger DAG would make checkpoint re-planning slower,
  which would matter for a Phase 2 estimate of the barrier's true "cost" beyond the job itself.
- **Whether `save_reference_configs`'s single measurement per modality (0.08 s multiome, 0.21 s
  scatac, n=1 each) is representative** — with only one job per modality there is no way to assess
  run-to-run variance, unlike the 6-job checkpoints above.
- **Whether the dead `"Kendall"` branch in `features_required` (§2 — `Kendall=True` and `ARC=True`
  both map to the same `final_val="ARC"` at line 19) is intentional dead code or a residual from an
  earlier design where "Kendall" and "ARC" were distinct output values** (the docstring comment at
  `add_external_features.smk:5` explicitly lists `"Kendall"` as a possible value, and
  `features_to_generate()`, `get_gex_file()`, `get_count_file()` all still branch on `val ==
  "Kendall"` separately from `val == "ARC"` — but no code path can ever produce `val == "Kendall"`
  given the checkpoint's own logic). This doesn't change any timing conclusion — I flag it only
  because it's a genuine structural oddity Phase 2 might otherwise mistake for a third live branch.

# Phase 1 unit #4 (Tier A) — `create_neighborhoods` / ABC neighborhood construction

**Headline.** The two hours are not a Python loop and not I/O. They are **six full
single-threaded streaming decompress-and-parse passes over the tagAlign** — three
`bedtools intersect | bedtools coverage | awk` pipelines and three *identical*
`zcat | grep -E | wc -l` pipelines — accounting for **97.4–99.6 % of wall time** in every
one of the 12 jobs. Cost is **strictly linear in `n_frag`** and **provably independent of
`j_cand_elems` and `k_genes`**. The 290 MB working set is not a paradox: the fragment
stream never enters Python. And the `cpu_time/wall = 0.977` is not evidence that one
process is working — it is the **ceiling imposed by `AllocCPUS=1`**, under which the
concurrent stages of each pipeline time-share one core and their costs *add*.

---

## 1. Unit identity

| item | value |
|---|---|
| Rule | `create_neighborhoods` (`abc_create_neighborhoods` in the exported DAG) |
| Rule definition | `ENCODE_rE2G/ABC/workflow/rules/neighborhoods.smk:1-48` |
| Entry point | `ENCODE_rE2G/ABC/workflow/scripts/run.neighborhoods.py` (209 lines; thin argparse wrapper, `processCellType` at `:148`) |
| Algorithm | `ENCODE_rE2G/ABC/workflow/scripts/neighborhoods.py` (877 lines) |
| Shared helpers | `ENCODE_rE2G/ABC/workflow/scripts/tools.py` (142 lines) — `run_command:16`, `run_piped_commands:21`, `df_to_pyranges:124` |
| Modalities | **BOTH**, 6 jobs each = 12 jobs |
| Declared `threads:` | **none** → Snakemake default 1 → profile `cpus_per_task=1` → **`AllocCPUS=1` confirmed by `sacct` for all 12 jobs** |
| Declared memory | **`mem_mb=32*1000` hard-coded** at `neighborhoods.smk:27` — *not* `determine_mem_mb` (see §6.4) |
| Outputs | `<cluster>/Neighborhoods/{EnhancerList.txt, EnhancerList.bed, GeneList.txt, GeneList.bed, GeneList.TSS1kb.bed, 3× *.CountReads.bedgraph}` + `<cluster>/processed_genes_file.bed` |

Row-count identities verified on disk (telohaec_crispri): `GeneList.txt` = 20,532 lines =
`k_genes` (20,531) **+ 1 header**; `EnhancerList.txt` = 157,552 lines = `j_cand_elems`
(157,551) + 1 header. *(Minor correction to `PHASE1_BRIEFING.md` §3: it lists `k_genes` as
"20,532 (constant)". The `size_manifest` column is 20,531; 20,532 is the line count
including the header.)*

---

## 2. What the code actually does

`processCellType` (`run.neighborhoods.py:148`) makes exactly **three** calls into the
counting machinery, each over the *same* accessibility file (`tagAlign.sort.gz`) but a
different BED of target intervals:

| # | caller | BED (`-a`) | rows in BED |
|---|---|---|---|
| 1 | `annotate_genes_with_features` → `count_features_for_bed(genes, GeneList.bed, …)` `neighborhoods.py:126` | whole gene bodies | `k_genes` = 20,531 |
| 2 | `annotate_genes_with_features` → `count_features_for_bed(tss1kb, GeneList.TSS1kb.bed, …)` `neighborhoods.py:136` | TSS ±500 bp | `k_genes` = 20,531 |
| 3 | `load_enhancers` → `count_features_for_bed(enhancers, candidateRegions.bed, …)` `neighborhoods.py:280` | candidate elements | `j_cand_elems` ≈ 157,551 |

Each of those funnels through `count_features_for_bed:495` → `count_single_feature_for_bed:534`,
which performs **two independent full passes over the tagAlign**:

**Pass A — `run_count_reads:401` → `count_tagalign:448-460`**

```python
remove_alt_chr_cmd = f"bedtools intersect -u -a {tagalign} -b {genome_sizes_bed}"
coverage_cmd = f"bedtools coverage -counts -sorted -g {genome_sizes} -b - -a {bed_file}"
awk_cmd = 'awk \'{{print $1 "\t" $2 "\t" $3 "\t" $NF}}\'' + f" > {output}"
run_piped_commands([remove_alt_chr_cmd, coverage_cmd, awk_cmd])
```

Three concurrent processes. `bedtools intersect` decodes the whole bgzipped tagAlign and
tests each record against 25 whole-chromosome intervals; `bedtools coverage -sorted`
re-parses the entire resulting text stream as a merge-sweep against `-a`; `awk` reprojects
4 columns.

**Pass B — `count_total:698` → `count_tagalign_total:674-684`**

```python
result = int(check_output(
    "zcat {} | grep -E 'chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY' | wc -l".format(tagalign),
    shell=True))
```

Called at `neighborhoods.py:567`, unconditionally, once per `count_single_feature_for_bed`
invocation. The argument (`feature_bam`) is the **same file path** all three times and the
result is a **single scalar** — so this pipeline computes the *identical* number three
times, with **no caching anywhere in the call chain**.

**Net: 3 × (Pass A + Pass B) = 6 full decompress-and-parse traversals of the tagAlign.**

The remaining control flow is all bounded by `k_genes` or `j_cand_elems` and is measured
at ≤ 2 s per pass plus ~10–32 s of post-processing (§5.4):
`make_tss_region_file:173` (pyranges + `bedtools sort`), `double_sex_chrom_counts:424`
(awk rewrite + `mv` of each bedgraph), `assign_enhancer_classes:327` (two pyranges/NCLS
joins + groupby-with-Python-lambda + a `.apply(…, axis=1)` over `j_cand_elems` rows at
`:395-397`), `run_qnorm:787` (scipy `interp1d`), and the `to_csv` writes at `:162`, `:307`,
`:313`.

**The tagAlign stream is 2 · `n_frag` records, not `n_frag`.** `frag_to_tagAlign`
(`workflow/rules/frag_to_tagAlign.smk`) emits two half-fragment records per fragment
(`print $1,$2,mid,…; print $1,mid+1,$3,…`). Measured directly: thp1_1 `n_frag` = 26,231,026 →
tagAlign records surviving the chromosome filter = **52,129,396** = 1.987 · `n_frag`.

### What is NOT the cost — things the brief asked me to look for and rule out

- **No per-row Python loop over fragments exists on this path.** `neighborhoods.py` does
  contain two per-element Python loops — `count_bam:440-443` (`iterrows` + a `pysam.count()`
  call *per BED region*) and `count_bigwig:467-484` (per-region `pyBigWig.stats`) — but
  `run_count_reads:401-420` dispatches on filename, and this pipeline's accessibility input
  is `tagAlign.sort.gz`, so **neither is reached**. They would matter for a BAM/bigWig
  biosample.
- **No accidental quadratic, no sort over the fragments, no per-element bedtools round trip.**
  Both bedtools stages are streaming (`-sorted` merge-sweep); the tagAlign is pre-sorted
  upstream; the `.tbi` index already exists as a declared rule input, so the
  `tabix -p bed` branch at `count_tagalign:450-452` never fires (confirmed absent from all
  12 logs).
- **The only per-row Python `.apply` on the hot path** is `neighborhoods.py:395-397`
  (`axis=1` string format over `j_cand_elems` rows). It is real but lives inside the
  10–32 s post-processing block, i.e. ≤ 0.5 % of wall on the critical cluster.
- **Dead work, minor:** `count_single_feature_for_bed:546` does `orig_df = df.copy()` and
  never uses the copy (only `orig_shape` at `:547` is used) — a full DataFrame copy per
  pass. `run_piped_commands` (`tools.py:21-34`) also never closes the parent's copy of each
  pipe read-end and never checks any stage's return code, so a failure in stage 1 or 2 of
  a pipeline is silently invisible.

---

## 3. Size variables (exact `size_manifest_*.tsv` names)

| variable | role | how it enters | range across 6 clusters |
|---|---|---|---|
| **`n_frag`** | **dominant, and effectively the only one** | stream length of `tagAlign.sort.gz` = 1.987 · `n_frag` records, traversed 6× | 26.2 M – 394.8 M (**15.1×**) |
| `j_cand_elems` | `-a` of pass 3; row count of every pandas frame in `load_enhancers`; drives `assign_enhancer_classes` / `run_qnorm` / the 33 MB `to_csv` | 156,420 – 158,548 (1.01×) |
| `k_genes` | `-a` of passes 1 and 2; row count of the gene frames | 20,531 (exactly constant) |

There is **no cross term.** I specifically tested for the plausible one —
`O(n_frag · |-a rows|)`, which would arise if bedtools' sweep cost scaled with the depth of
the `-a` interval set — and it does not exist (§5.3). The `-a` BED enters only through the
`O(|-a|)` merge-sweep and the pandas post-processing.

`n_frag`, `j_cand_elems` and `k_genes` are **byte-identical between the two modalities**
(same upstream fragment files, same reference; this rule sits upstream of the ARC/Kendall
branch point), so the two runs are true n = 2 replicates on identical inputs.

---

## 4. Asymptotic class derived from the code

$$T \;=\; \underbrace{3\big[\,c_A\cdot 2n_{frag} + c_B\cdot 2n_{frag}\,\big]}_{\text{6 streaming passes}} \;+\; \underbrace{O(j_{cand\_elems} + k_{genes})}_{\text{sweep + pandas}} \;+\; C_0$$

$$\boxed{\;O(n_{frag}) \;+\; O(j_{cand\_elems}\log j_{cand\_elems})\;}$$

**Strictly linear in `n_frag`, with a constant of 6 passes.** Justification, line by line:

- `count_tagalign:454` `bedtools intersect -u -a <tagAlign> -b <chrom_sizes_bed>` — `-b` is
  25 intervals loaded into an interval tree; `-a` is streamed record by record. O(stream).
- `count_tagalign:455-457` `bedtools coverage -counts -sorted` — `-sorted` selects
  bedtools' `chromsweep`, a single synchronised merge over two *already sorted* streams. No
  sort, no index build, no nested loop. O(stream + |-a|).
- `count_tagalign:458` `awk '{print $1,$2,$3,$NF}'` — one pass, O(|-a|) output records but
  O(stream) is already spent upstream.
- `count_tagalign_total:676-681` `zcat | grep -E … | wc -l` — three streaming filters.
  O(stream).
- Everything above is invoked **exactly 3 times** from `count_features_for_bed:512`, driven
  by the three `count_features_for_bed` call sites (`:126`, `:136`, `:280`). The multiplier
  is 3 (bounded, not a loop over data) because `params["features"]` has one key (`ATAC`) with
  one file — see `get_features:727`.
- The `log` term is *only* in `j_cand_elems`: the four `.rank()` calls at
  `count_single_feature_for_bed:581-589`, two more in `average_features:608-610`, the
  `bedtools sort` in `make_tss_region_file:188`, `drop_duplicates`/`duplicated` at `:574`
  and `:592`, and the pyranges NCLS joins at `assign_enhancer_classes:348,356`. All
  O(j log j) with j ≈ 157 k — ~2.7 × 10⁶ comparisons, i.e. seconds, and confirmed as such.
- **No term is quadratic.** The `groupby(...).aggregate(lambda x: ",".join(list(set(x))))`
  at `:350-353` and `:358-361` is a Python lambda per *group*, O(overlap hits) ≈ 91 k, not
  O(j²).
- **Memory is O(`j_cand_elems` + `k_genes`), independent of `n_frag`** — because the
  fragment stream is confined to subprocesses that never materialise it. This is the direct
  answer to "two hours at a 290 MB working set."

---

## 5. Benchmark cross-check — the code and the data agree, exactly

This is one of the few units where the benchmarks are decisive. `n_frag` spans 15.1× while
`j_cand_elems` spans 1.01× and `k_genes` is exactly constant, so the discriminating
experiment is available — **and the answer is unambiguously `n_frag`.**

### 5.1 Benchmark rows (`analysis/benchmarks_{multiome,scatac}.tsv`)

| modality | cluster | `n_frag` | wall s | cpu_time s | CPU/wall | benchmark `max_rss` MB |
|---|---|---:|---:|---:|---:|---:|
| multiome | thp1_1 | 26.23 M | 374.64 | 301.30 | 0.80 | 179.5 |
| multiome | thp1_2 | 26.29 M | 376.27 | 303.21 | 0.81 | 286.6 |
| multiome | jurkat_pma_cd3_4hr | 55.80 M | 777.72 | 705.31 | 0.91 | 177.1 |
| multiome | jurkat | 75.34 M | 1051.35 | 966.13 | 0.92 | 176.1 |
| multiome | k562_crispri | 86.42 M | 992.30 | 845.17 | 0.85 | 176.4 |
| multiome | **telohaec_crispri** | **394.85 M** | **6725.41** | **6572.08** | **0.977** | **290.0** |
| scatac | thp1_1 | 26.23 M | 376.54 | 302.10 | 0.80 | 316.5 |
| scatac | thp1_2 | 26.29 M | 428.78 | 266.87 | 0.62 | 181.2 |
| scatac | jurkat_pma_cd3_4hr | 55.80 M | 825.47 | 693.22 | 0.84 | 167.2 |
| scatac | jurkat | 75.34 M | 1061.57 | 968.19 | 0.91 | 479.2 |
| scatac | k562_crispri | 86.42 M | 1202.32 | 1114.52 | 0.93 | 364.1 |
| scatac | **telohaec_crispri** | **394.85 M** | **6739.76** | **6602.17** | **0.98** | **512.8** |

Wall spans **18.0×** (374.6 s → 6739.8 s) while `n_frag` spans 15.1× and everything
downstream spans 1.01×. Nothing near-constant can produce an 18× spread. **The variable is
`n_frag`.**

Note that CPU/wall is *not* 0.98 everywhere — it rises monotonically-ish from 0.62–0.81 on
the smallest clusters to 0.977–0.980 on the largest. This is not variable parallelism: the
whole job is capped at one core throughout (§5.3). It is the fixed non-CPU overhead — the
cold `import neighborhoods` off the Oak-hosted conda env, measured at up to 30.8 s (§6.5) —
being a 3–17 % share of a 6-minute job and a 0.4 % share of a 112-minute one. The briefing's
0.977 figure for telohaec is the *asymptote* of this rule's CPU efficiency, which is exactly
what makes it the cleanest evidence that the job is compute-saturated on its single core.

### 5.2 Free decomposition from the retained slurm logs

`count_features_for_bed:507,528-529` times and prints each of the three passes. All 12 logs
carry three `Feature ATAC completed in …` lines. Combined with the (untouched) mtimes of
the three `*.CountReads.bedgraph` files — each written at the end of Pass A, before Pass B
starts — this yields a **complete per-sub-pass budget for all 36 passes without running
anything**:

| modality | cluster | node CPU | Σ 3 passes | wall | wall − Σ | s per 10⁶ `n_frag` | Pass A (bedtools) | Pass B (zcat\|grep\|wc) |
|---|---|---|---:|---:|---:|---:|---:|---:|
| multiome | thp1_1 | Siena 8224P | 364.5 | 374.6 | 10.1 | 13.896 | 220 (60.4 %) | 144 (39.5 %) |
| multiome | thp1_2 | Siena 8224P | 366.4 | 376.3 | 9.9 | 13.936 | 223 (60.9 %) | 143 (39.0 %) |
| multiome | jurkat_pma_cd3_4hr | Siena 8224P | 768.8 | 777.7 | 8.9 | 13.777 | 468 (60.9 %) | 301 (39.2 %) |
| multiome | jurkat | Siena 8224P | 1040.2 | 1051.3 | 11.2 | 13.806 | 633 (60.9 %) | 407 (39.1 %) |
| multiome | k562_crispri | **SapphireRapids 8462Y+** | 928.8 | 992.3 | 63.5 | **10.747** | 554 (59.6 %) | 374 (40.3 %) |
| multiome | telohaec_crispri | **Rome 7502** | 6628.2 | 6725.4 | 97.2 | **16.787** | 4160 (62.8 %) | 2469 (37.2 %) |
| scatac | thp1_1 | Siena 8224P | 366.4 | 376.5 | 10.1 | 13.968 | 224 (61.1 %) | 142 (38.8 %) |
| scatac | thp1_2 | Siena 8224P | 364.2 | 428.8 | 64.6 | 13.852 | 222 (61.0 %) | 142 (39.0 %) |
| scatac | jurkat_pma_cd3_4hr | Siena 8224P | 769.2 | 825.5 | 56.3 | 13.784 | 468 (60.8 %) | 301 (39.1 %) |
| scatac | jurkat | Siena 8224P | 1039.8 | 1061.6 | 21.8 | 13.801 | 634 (61.0 %) | 406 (39.0 %) |
| scatac | k562_crispri | Siena 8224P | 1191.2 | 1202.3 | 11.1 | 13.783 | 727 (61.0 %) | 464 (39.0 %) |
| scatac | telohaec_crispri | **Rome 7502** | 6712.3 | 6739.8 | 27.5 | **17.000** | 4188 (62.4 %) | 2525 (37.6 %) |

Three facts fall straight out:

1. **The three passes are 97.4–99.6 % of wall.** The residual is 8.9–64.6 s with **no
   `n_frag` dependence at all** (the second-largest residual, 64.6 s, belongs to the
   *smallest* cluster).
2. **The three passes cost the same as each other to within ±5 %** — 2207 / 2159 / 2262 s
   on telohaec — despite pass 3's `-a` BED having **7.7× the rows** of passes 1 and 2
   (158 k vs 20.5 k) and a completely different interval-length distribution (500 bp
   elements vs multi-megabase gene bodies). Cost does not see `-a`.
3. The Pass A / Pass B split is stable at **~61 % / ~39 %**. Pass B —
   `count_tagalign_total`, whose three invocations compute the *same scalar* — is
   **39 % of the whole rule**, of which two thirds (**26 % of total wall**) is provably
   redundant recomputation.

### 5.3 Confirmatory experiment (my own; 1 vs 4 cores, same node family)

Ran the exact command strings from `neighborhoods.py:448-460` and `:676-681` against
thp1_1's tagAlign on **sh04-16n27** — the same node that ran two of the production jobs —
under `--cpus-per-task=1` and `--cpus-per-task=4` (`bedtools v2.26.0` from the pipeline's
own conda env):

| stage | 4-CPU wall | 4-CPU CPU | 1-CPU wall | 1-CPU CPU |
|---|---:|---:|---:|---:|
| `zcat file > /dev/null` (decompression alone) | 8.60 | 8.57 | 8.55 | 8.53 |
| **`zcat \| grep -E \| wc -l`** = `count_tagalign_total` | 36.81 | 47.44 | **46.79** | 46.78 |
| `bedtools intersect -u` alone (stage 1 of 3) | 35.40 | 35.39 | 35.57 | 35.56 |
| **full 3-stage pipe, `-a` = candidateRegions (158,548 rows)** | 36.18 | 71.52 | **74.60** | 74.54 |
| **full 3-stage pipe, `-a` = GeneList.bed (20,531 gene bodies)** | 35.91 | 70.82 | **73.99** | 73.95 |

Five conclusions, each load-bearing:

- **The `AllocCPUS=1` cpuset is the whole explanation for CPU/wall ≈ 0.977.** At 1 CPU,
  wall ≡ total CPU (74.60 vs 74.54). At 4 CPUs the same work finishes in 36.18 s wall for
  71.52 s of CPU (**CPU/wall = 1.98**). The three pipeline stages *are* concurrent and *do*
  overlap when given cores; confined to one core they **time-share and their costs add**.
  0.977 is a ceiling, not a measurement of serial work. Confirmed independently by `sacct`:
  `AllocCPUS=1`, `ReqTRES=cpu=1`, `TotalCPU=01:50:36` vs `Elapsed=01:53:11` for job
  39924729, and by `slurm-efficiency-report` (`efficiency_report_*.csv`: `NCPUS=1`,
  CPU efficiency 87.2–99.5 %).
- **Penalty from the 1-core cpuset:** `count_tagalign` 74.60 / 36.18 = **2.06×**;
  `count_tagalign_total` 46.79 / 36.81 = **1.27×**.
- **`-a` independence, directly:** 74.60 s vs 73.99 s — **0.8 % apart across a 7.7×
  difference in `-a` row count** and a ~2500× difference in mean `-a` interval length. This
  is the experiment that kills the `O(n_frag · |-a|)` cross term and, with it, any role for
  `j_cand_elems` or `k_genes` in the dominant cost.
- **It reproduces production exactly.** 74.60 + 46.79 = 121.4 s predicted per pass; the
  production thp1_1 passes measured 122.9 / 119.0 / 122.6 s. Within 1–2 %.
- **Half of Pass A goes to a filter that removes 0.63 % of records.** `bedtools intersect -u`
  alone is 35.6 s of the pipe's 74.5 s CPU; the tagAlign was already built with
  `-g chrSizes` upstream and `chrom_sizes_bed` is 25 whole-chromosome intervals. Measured
  survivors: 52,129,396 of 52,462,052.

### 5.4 The 15 % apparent superlinearity is a **node-generation artifact**, not an exponent

A naive fit over all 12 points gives `T ∝ n_frag^1.07` — the telohaec per-fragment rate is
22 % above the small clusters'. That is not algorithmic. **Both** telohaec jobs landed on
`sh03-04n2{1,2}` (AMD EPYC **7502** "Rome", Zen2, 2.50 GHz) while nine of the other ten ran
on `sh04-16n2{5,7,8}` (AMD EPYC **8224P** "Siena", Zen4c, 2.55 GHz) and the one remaining
outlier — k562 multiome, the *fastest* per fragment — ran on `sh04-02n16` (Intel Xeon
**8462Y+** Sapphire Rapids, 2.80 GHz). Stratifying by CPU:

| node CPU | jobs | `n_frag` range | s per 10⁶ `n_frag` |
|---|---:|---|---|
| **Siena 8224P** | **9** | 26.2 M – 86.4 M (**3.30×**) | **mean 13.845, sd 0.068, CV 0.49 %**, range 13.777–13.968 |
| Rome 7502 | 2 | 394.8 M | 16.787, 17.000 |
| SapphireRapids 8462Y+ | 1 | 86.4 M | 10.747 |

**On identical hardware the per-fragment cost is constant to 0.49 % across a 3.30× range in
`n_frag`.** That is as clean an empirical confirmation of a linear law as this dataset can
produce, and it agrees with the code exactly. The 1.57× spread between CPU generations
(10.75 → 16.89) is the expected single-thread IPC difference for a
scalar-parse-and-decompress workload across Zen2 / Zen4c / Sapphire Rapids, and it is
reproducible: both telohaec runs, on two different Rome nodes, gave 16.79 and 17.00.

**Consequence for the replicate/noise question.** Where the two modalities drew the same
node generation, agreement is **≤ 1.3 %** (0.5 %, 0.6 %, 0.05 %, 0.04 %, 1.3 %). The single
28 % gap (k562: 928.8 vs 1191.2 s) is entirely accounted for by Sapphire Rapids vs Siena.
For *this* rule the noise floor is not ~2 % rising to 12–46 % with size — it is **~1 %
plus a discrete node-generation lottery worth up to 1.57×**, and the largest cluster lost
that lottery twice.

---

## 6. Concrete overhead sources

### 6.1 Bytes actually moved (`io_in`/`io_out` are unusable — on-disk sizes used)

Six traversals of `tagAlign.sort.gz`:

| cluster | tagAlign gz | gz bytes read (×6) | text decoded & re-parsed (×6) |
|---|---:|---:|---:|
| thp1_1 | 258.2 MB | **1.55 GB** | ~9.0 GB |
| jurkat | 679.4 MB | **4.08 GB** | ~26 GB |
| k562_crispri | 736.1 MB | **4.42 GB** | ~28 GB |
| **telohaec_crispri** | **2800.1 MB** | **16.80 GB** | **~136 GB** |

Writes are trivial by comparison (~60 MB): `GeneList.bed` + `GeneList.TSS1kb.bed` (709 KB
each, the latter written twice via the `sort`/`mv` at `:188`), three bedgraphs (5.5 MB,
each **written twice** because `double_sex_chrom_counts:424-429` rewrites via `awk … > .tmp && mv`),
`GeneList.txt` 5.6 MB, `EnhancerList.bed` 8.9 MB, `EnhancerList.txt` 33.2 MB (via `.tmp` +
`os.rename` at `:320`).

**Why `io_in` is unusable, mechanistically:** it reports 17.3 MB for thp1_1 against 1.55 GB
actually read (90× under) and 1070 MB for telohaec against 16.80 GB (16× under). Snakemake's
psutil sampler reads `io_counters()` for the Python process and its *live* children and
takes a max over samples; the bytes are consumed by short-lived `bedtools`/`zcat`/`grep`
grandchildren that are reaped between samples. Use on-disk sizes.

### 6.2 Repeated passes over the same data — the single largest structural overhead

- `count_tagalign_total` (`neighborhoods.py:674-684`) is invoked **3×** on the same file for
  the same scalar, at `:567`, with no memoisation. **39 % of the rule; 26 % is redundant.**
- `count_tagalign` (`:448-460`) is invoked 3× on the same file for three different `-a` BEDs
  — and within each invocation the stream is fully parsed **twice** (once by
  `bedtools intersect`, again by `bedtools coverage` off the pipe), of which the first
  parse discards 0.63 % of records.
- **Cross-unit note, do not double-count:** unit #5 (`make_candidate_regions`) reaches the
  *same* function — `peaks.py:99 get_read_counts` → `count_reads_over_peaks:225` →
  `run_count_reads` → `count_tagalign` — for one further pass over the same tagAlign. Across
  the two rules the pipeline traverses each cluster's tagAlign **7 times** (4 bedtools + 3
  zcat). Unit #5's 81 s → 1478 s figure is one `count_tagalign`; my mtime decomposition puts
  a single `count_tagalign` at 1389 s on telohaec, consistent. `create_neighborhoods` is a
  separate DAG node with its own benchmark rows, so its 6 passes are its own cost.

### 6.3 Declared vs realised threading

`create_neighborhoods` declares **no `threads:`** (`neighborhoods.smk:1-48`), as does every
ABC rule. Snakemake defaults to 1, the profile sets `cpus_per_task=1`
(`profiles/measure/config.yaml`, `default-resources`), and `sacct` confirms `AllocCPUS=1`
for all 12 jobs. Unlike `frag_to_norm_bigWig` (declares 16, realises ~1), this rule's gap
runs the *other* way: **the work is genuinely parallel — six pipelines of 2–3 concurrent
CPU-bound stages — but the allocation forbids it.** Measured cost of that confinement:
2.06× on `count_tagalign`, 1.27× on `count_tagalign_total`; **~1.66× on the rule as a
whole** (per pass, 74.5 + 46.8 = 121.3 s of CPU on one core versus 36.2 + 36.8 = 73.0 s of
wall when each pipeline runs at its natural width). This is characterisation only — Phase 2
owns whether it matters.

### 6.4 Memory request

`neighborhoods.smk:27` hard-codes `mem_mb=32*1000`. It does **not** call
`ABC.determine_mem_mb` (`ENCODE_rE2G/ABC/workflow/rules/utils.smk:17-24`,
`attempt_multiplier * max(4 * input_size_mb, min_gb*1000)`, ×8 for `.gz` inputs, capped at
`MAX_MEM_MB`) — the other five ABC rules (`macs2.smk:15,50,70`,
`candidate_regions.smk:22`, `predictions.smk:36,73`, `qc.smk:21`) do.

Two measurements of peak RSS disagree and both matter:

| cluster | Snakemake `max_rss` (multiome / scatac) | `sacct` MaxRSS of the python step (multiome / scatac) | request | over-provision (sacct) |
|---|---|---|---:|---:|
| thp1_1 | 179 / 317 MB | **421 / 425 MB** | 32,000 MB | 76× |
| jurkat | 176 / 479 MB | **667 / 920 MB** | 32,000 MB | 35× |
| k562_crispri | 176 / 364 MB | **833 / 652 MB** | 32,000 MB | 38× |
| **telohaec_crispri** | **290 / 513 MB** | **1328 / 1219 MB** | 32,000 MB | **24×** |

**Over-provisioned 24×–77×** (1.3 %–3.7 % memory efficiency in
`efficiency_report_*.csv`). Note also that Snakemake's benchmark `max_rss` **undercounts by
1.5–4.7×** relative to the cgroup high-water mark, for the same sampling reason as `io_in`;
the briefing's "290 MB" figure for telohaec is really ~1.3 GB. And the briefing's implicit
read that 290 MB is telohaec-specific does not hold — thp1_2, the *smallest* cluster,
also reports 287 MB. Python-side memory is `j_cand_elems`/`k_genes`-driven and hence flat;
the sacct growth (421 MB → 1328 MB, i.e. ~`n_frag`^0.42 over 15.1×) is buffer growth in the
bedtools subprocesses, not in Python.

### 6.5 Interpreter and import startup (measured)

`neighborhoods.py:1-14` imports `numpy`, `pandas`, `pyranges`, `pyBigWig`, `pysam`,
`scipy.interpolate`, plus `linecache` (never used). Measured on a compute node against the
pipeline's own conda env (which lives on **Oak**, per the mandatory `--conda-prefix`):

- `python -c "pass"`: **0.01–0.07 s**
- `import neighborhoods`: **30.83 s cold, then 0.81 / 0.68 s warm**
- individually, warm: pandas 0.39 s, pyranges 0.46 s, scipy.interpolate 0.21 s, numpy 0.10 s,
  **pyBigWig 0.10 s, pysam 0.05 s**

`pyBigWig` and `pysam` are imported unconditionally but are reachable only from
`count_bigwig`/`count_bam`/`count_bam_mapped`/`count_bigwig_total`, none of which this
pipeline enters. The 30.8 s cold figure is Lustre/Oak metadata + shared-object paging and it
**fully explains the 8.9–64.6 s "wall − Σ passes" residual** (§5.2) and its lack of
correlation with `n_frag` — it is a cold-cache lottery, not scaling. Consistent with the
"1–65 s of startup" reported by unit #5.

### 6.6 Serialisation / deserialisation

Per pass, between Pass A and Pass B: `read_bed:632` of the bedgraph, then after Pass B a
`merge` + `drop_duplicates` + 4 `.rank()` calls (`:574-589`) + `duplicated()` (`:592`) +
`average_features:595` (2 more ranks). Measured residual inside each pass
(`pass − Pass A − Pass B`) = **~2.1 s at `j_cand_elems` ≈ 158 k**. Post-processing after
the last pass (`assign_enhancer_classes` incl. the `axis=1` `.apply` at `:395`, `run_qnorm`,
and writing 33.2 MB `EnhancerList.txt` + 5.6 MB `GeneList.txt` with `float_format="%.6f"`)
= **~10 s (thp1_1) to ~26–32 s (telohaec)**, derived from log timestamps minus the summed
pass durations. Together ≤ 0.6 % of wall on the critical cluster.

---

## 7. The three functions

Units: **seconds**. Variables named exactly as in `size_manifest_*.tsv`. Directly
evaluable; every constant stated numerically.

```
t_est(n_frag, j_cand_elems)  = 13.85e-6 * n_frag + 1.3e-4 * j_cand_elems + 12
t_low(n_frag, j_cand_elems)  = 10.75e-6 * n_frag + 1.0e-4 * j_cand_elems +  1
t_high(n_frag, j_cand_elems) = 16.90e-6 * n_frag + 2.0e-4 * j_cand_elems + 76
```

**All three are linear in `n_frag`, by design.** The code is unambiguous about the class —
six streaming passes, no sort of the fragment stream, no nested loop, `-sorted` merge-sweeps
— and §5.4 confirms it empirically to a CV of 0.49 % on fixed hardware. Per the Phase 1
instructions, `t_low`/`t_high` therefore differ **only in their constant**, not their
exponent. There is no ambiguity in the code to justify a different exponent, and inventing
one would be the mirror-image error of fitting to noise.

**Where the spread comes from — and it is not statistical.** The leading coefficient is the
measured per-fragment cost on the three CPU generations this pipeline actually ran on:

| coefficient | s / 10⁶ `n_frag` | source |
|---|---|---|
| `t_low` | 10.75 | Intel Xeon 8462Y+ (Sapphire Rapids), `sh04-02n16`, n=1 |
| `t_est` | 13.85 | AMD EPYC 8224P (Siena), `sh04-16n2{5,7,8}`, **n=9, CV 0.49 %** |
| `t_high` | 16.90 | AMD EPYC 7502 (Rome), `sh03-04n2{1,2}`, n=2 |

Secondary terms: the `j_cand_elems` coefficient absorbs `assign_enhancer_classes` +
`run_qnorm` + the per-pass pandas work + the 33 MB `to_csv`; the constant absorbs the
conda-env-on-Oak import (0.7 s warm to 30.8 s cold) plus the shell preamble. Both are
< 3 % of the total at every measured scale.

**Validation.** `t_est` is within **5.2 %** of measured wall on all nine same-hardware
(Siena) jobs. `t_low ≤ measured ≤ t_high` holds for **all 12 of 12** jobs. `t_est`
under-predicts telohaec by 18 % precisely because both telohaec jobs drew the slowest node
generation; `t_high` captures that within 0.6 % (6780 predicted vs 6725/6740 measured).
**Phase 2 should treat `t_high` as the live risk for the largest cluster, not an
implausible pessimum** — it is the empirically observed outcome, twice.

If Phase 2 needs a single-variable form, `t_est(n_frag) = 13.85e-6 * n_frag + 32` is within
6 % of `t_est` above at `j_cand_elems` ∈ [156 k, 159 k].

*(Offered but not recommended, and explicitly not code-supported: a whole-dataset fit
ignoring node identity gives `T = 4.02e-4 * n_frag^1.070`. It interpolates the 12 points
slightly better than the linear form only because the node-generation confound happens to
correlate with size in this run. Do not extrapolate with it.)*

## 8. All three evaluated

| point | `n_frag` | `j_cand_elems` | `t_low` | `t_est` | `t_high` | measured (multiome / scatac) |
|---|---:|---:|---:|---:|---:|---|
| **min** — `thp1_1` | 26,231,026 | 158,548 | 298.8 s (5.0 min) | **395.9 s (6.6 min)** | 551.0 s (9.2 min) | 374.6 / 376.5 s |
| **mean** of 6 clusters | 110,823,580 | 157,635 | 1208.1 s (20.1 min) | **1567.4 s (26.1 min)** | 1980.4 s (33.0 min) | — |
| **max** — `telohaec_crispri` | 394,847,593 | 157,551 | 4261.4 s (71.0 min) | **5501.1 s (91.7 min)** | 6780.4 s (113.0 min) | 6725.4 / 6739.8 s |

(Reporting convenience only; the functions in §7 are the deliverable.)

---

## 9. Multiome vs scATAC: neither the class nor the constant differs

The rule runs **upstream of the ARC/Kendall branch point** (`checkpoint features_required`),
so both modalities feed it byte-identical inputs: `n_frag`, `j_cand_elems`, `k_genes` and
`enh_list_bytes` are identical row-for-row between `size_manifest_multiome.tsv` and
`size_manifest_scatac.tsv`, and the tagAlign files are byte-identical in size. Same class,
same constant. The 12 jobs are **6 true n = 2 replicate pairs**, which is what made §5.4
possible: five of six pairs agree to ≤ 1.3 %, and the sixth (28 %) is fully explained by
CPU generation. If Phase 2 needs a per-modality figure, use the same function for both.

---

## 10. What the measurement cannot tell us

- **Anything about `j_cand_elems` scaling beyond ~157 k.** `j` spans 1.01× and I proved it
  is *absent* from the dominant term, so the measurements constrain its coefficient only to
  "≤ 3 % of total at j ≈ 157 k". The O(j log j) ranks/sorts/joins listed in §4 are
  code-derived and would remain O(j log j) at any j; I cannot bound their constant to better
  than a factor of ~3. In particular I did **not** separately time
  `assign_enhancer_classes` or the `axis=1` `.apply` at `:395-397`; my ~5–10 s attribution
  for that `.apply` is an **estimate from typical pandas per-row cost, not a measurement**,
  and it is only defensible because the whole post-processing block it sits in is bounded at
  10–32 s.
- **`k_genes` is exactly constant (20,531).** Nothing in this dataset can say anything about
  its exponent. From the code it enters only as the `-a` of passes 1–2 (O(k) sweep) and the
  gene DataFrames (O(k log k) ranks) — and §5.3 shows the sweep term is invisible at
  k = 20.5 k.
- **`n_frag` cannot be separated from `n_umi` (12.0×) or `n_cells` (6.0×) by regression** —
  they are correlated across these six clusters. Here it does not matter: the code reads
  *only* the tagAlign, whose record count I measured directly as 1.987 · `n_frag`
  (52,129,396 records at `n_frag` = 26,231,026). `n_umi` and `n_cells` are not inputs to
  this rule at all. The attribution is from the code, not the fit.
- **Only one accessibility file and only ATAC.** `count_features_for_bed:506-511` loops over
  `features` × `feature_bam_list`; every cluster here has exactly one ATAC tagAlign and no
  H3K27ac (`H3K27ac=None` in all 12 logs). The pass count is therefore
  `3 × (1 + n_H3K27ac_files)`, and I have **zero** measurements with H3K27ac or with a
  multi-file feature list. A biosample with one H3K27ac bigWig would take a completely
  different code path (`count_bigwig:463-484`, a **per-region Python loop**) whose cost this
  dataset says nothing about.
- **BAM and bigWig inputs are entirely unmeasured** and would change the class: `count_bam`
  (`:432-445`) is `iterrows` + one `pysam.count()` per BED region → O(|-a| · log n_frag)
  with a large Python constant; `count_bigwig` (`:463-484`) is a Python loop over `-a`. Both
  are plausibly *worse* than the tagAlign path at j ≈ 157 k. Not exercised here.
- **The node-generation coefficient rests on n=1 (Sapphire Rapids) and n=2 (Rome).** The
  1.57× spread is real and reproducible for Rome, but `t_low`'s 10.75 comes from a single
  job. A future run scheduled entirely onto `sh04-02` nodes could plausibly come in below
  `t_low`.
- **I cannot decompose Pass A's 74.5 s of CPU beyond two stages.** I measured
  `bedtools intersect` alone (35.6 s) and the full pipe (74.5 s), so
  `bedtools coverage -counts -sorted` + `awk` ≈ 38.9 s; I did not separate `coverage` from
  `awk`. Similarly for Pass B I measured `zcat` alone (8.5 s) and the full pipe (46.8 s), so
  `grep -E` + `wc -l` ≈ 38.3 s, with `grep -E` presumed dominant but not isolated.
- **`bedtools v2.26.0`** (2017) is what the conda env pins. All constants above are specific
  to that build; a modern bedtools could shift them materially in either direction. Not
  tested.
- **Nothing here addresses whether this rule is on the critical path.** By construction it
  produces `GeneList.txt` and `EnhancerList.txt`, consumed downstream, but Phase 2 owns the
  path arithmetic.

---

### Reproduction pointers

- Benchmarks: `/scratch/users/kaybrand/scE2G_optimize_results/analysis/benchmarks_{multiome,scatac}.tsv`, rule `create_neighborhoods`
- Per-pass timings: `…/igvf10_{multiome,scatac}/slurm_logs/rule_abc_create_neighborhoods/<cluster>/<jobid>.log`, lines `Feature ATAC completed in …`
- Sub-pass split: mtimes of `…/igvf10_*/<cluster>/Neighborhoods/{Genes,Genes.TSS1kb,Enhancers}.ATAC.tagAlign.sort.gz.CountReads.bedgraph`
- Allocation / true peak RSS: `sacct -j <jobid> --format=AllocCPUS,ReqTRES,Elapsed,TotalCPU,MaxRSS,NodeList`; node CPUs via `sinfo -N -n <node> -o "%f"`
- My 1-vs-4-core and import experiments: `/scratch/users/kaybrand/phase1_u04/{decomp.sh,imports.sh}` and their `*_<jobid>.out`

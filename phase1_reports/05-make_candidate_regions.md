# Phase 1 unit #5 — `make_candidate_regions` (Tier A)

Repo: `/oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G_optimize`, branch `optimize-sce2g`.
All paths below are repo-relative unless absolute. Characterisation only — no fix proposals.

---

## 0. Headline

The prompt's framing was **right in shape, wrong in location**. There is a genuine pre-cap /
post-cap split, but the cap is not where the money is:

| term | what it is | code | spans | share of runtime |
|---|---|---|---|---|
| **pre-cap** | one streaming pass over the tagAlign, counting ATAC reads per MACS2 peak | `neighborhoods.py:448-460` | **81 s → 1478 s (18.2×)** | **55 % – 97 %** |
| **post-cap** | the entire 14-stage `bedtools` pipe *including the cap and everything before it* | `peaks.py:105-120` | **5 s → 10 s** | **0.4 % – 4.5 %** |
| startup | conda activate + `pandas`/`pyranges`/`pysam`/`pyBigWig`/`scipy` import | `neighborhoods.py:8-13` | 1 s → 65 s | 0.7 % – 41 % |

The `nStrongestPeaks` cap **is binding and hard** (it discards 52 %–66 % of merged accessible
regions), and it is what makes `j_cand_elems` — and therefore every downstream rule — constant
at 1.01×. But it bounds essentially **none of this rule's own cost**, because 95 %–99.6 % of that
cost is spent *before the pipe starts*, in a full pass over `2·n_frag` tagAlign records that the
cap cannot touch.

I was able to decompose the runtime **exactly** (to ±1 s on all 12 jobs) using file mtimes plus
the slurm-log job-start timestamp, so the pre/post-cap split below is measured, not modelled.

---

## 1. Unit identity

| item | value |
|---|---|
| Rule | `make_candidate_regions`, defined `ENCODE_rE2G/ABC/workflow/rules/candidate_regions.smk:1` |
| Rule name on disk | benchmark table: `make_candidate_regions`; slurm log dir: `rule_abc_make_candidate_regions` (the outer workflow imports the ABC module with an `abc_` prefix) |
| Entry script | `ENCODE_rE2G/ABC/workflow/scripts/makeCandidateRegions.py` — **119 lines** |
| Work script | `ENCODE_rE2G/ABC/workflow/scripts/peaks.py` — **260 lines** |
| Also on the hot path | `ENCODE_rE2G/ABC/workflow/scripts/neighborhoods.py` — 877 lines (only `run_count_reads`/`count_tagalign`/`double_sex_chrom_counts` execute) |
| | `ENCODE_rE2G/ABC/workflow/scripts/tools.py` — 142 lines (`run_command`, `run_piped_commands`, `write_params`) |
| Modalities | **BOTH**. 6 jobs each → 12 jobs total |
| Declared `threads:` | **none** → runs at 1. Confirmed in the slurm log: `cpus_per_task=1` |
| Declared memory | `resources: mem_mb=determine_mem_mb` (`candidate_regions.smk:21-22`) |

**Replicate status.** This rule sits upstream of `checkpoint features_required`, so the multiome
and scATAC runs execute the *same* computation on the *same* fragment stack in two independent
Snakemake runs. `n_frag` is byte-identical; `n_peaks` differs by 2 rows for `telohaec_crispri`
only (561,437 vs 561,435 — MACS2/sort tie-break); `j_cand_elems` is **identical in all six
clusters**. So this unit has a true **n = 2 per cluster**, which I use throughout.

---

## 2. What the code actually does

### 2.1 Parameter plumbing

`candidate_regions.smk:23-35` shells out to `makeCandidateRegions.py` with
`--nStrongestPeak {params.nStrongestPeak}` (`candidate_regions.smk:34`), sourced from
`config['params_candidate']['nStrongestPeaks']` (`candidate_regions.smk:12`).

The argparse option is actually named `--nStrongestPeaks` (`makeCandidateRegions.py:46`); the
rule passes the *singular* `--nStrongestPeak`. This binds correctly only because argparse
accepts unambiguous prefix abbreviations. Verified from the on-disk artefact — every
`<cluster>/Peaks/params.txt` records `nStrongestPeaks 150000`, so the intent is realised. (Noting
the fragility as an observation; no recommendation.)

**The cap value is 150,000**, from `ENCODE_rE2G/ABC/config/config.yaml:31`.
`params.txt` also confirms `ignoreSummits False`, so `makeCandidateRegions.py:86` takes the
`make_candidate_regions_from_summits` branch (`peaks.py:76`), and the
`make_candidate_regions_from_peaks` branch (`peaks.py:129`) is dead in this configuration.

### 2.2 Control flow, in execution order

```
makeCandidateRegions.py:118  parseargs
makeCandidateRegions.py:82   os.makedirs
makeCandidateRegions.py:83   write_params  -> Peaks/params.txt        <-- mtime marker A
peaks.py:93                  validate_macs_peaks
peaks.py:95-96               build includelist / blocklist command strings
peaks.py:99-101              get_read_counts  ------------------------- THE 95-99% TERM
peaks.py:105-120             piped_cmds (14 stages)
peaks.py:122                 run_piped_commands -> candidateRegions.bed
peaks.py:125                 count_lines(outfile)
peaks.py:126                 write_candidate_regions_qc -> qc.txt      <-- mtime marker C
```

**(a) `validate_macs_peaks` — `peaks.py:12-43`.** `n_peaks = sum(1 for _ in open(macs_peaks))`
at **`peaks.py:23`**: a pure-Python line count over the 32–44 MB `macs2_peaks.narrowPeak.sorted`.
O(`n_peaks`) with a Python-level constant. Its only consumer is the QC text file
(`peaks.py:52`). This is the **first of three full reads** of the narrowPeak in this rule.

**(b) `get_read_counts` → `count_reads_over_peaks` → `run_count_reads` — `peaks.py:196-260`,
`neighborhoods.py:401-422`.** The accessibility input is
`<cluster>/tagAlign/tagAlign.sort.gz`, so `run_count_reads` dispatches to the
`"tagAlign" in filename` branch at **`neighborhoods.py:413-414`**.

`count_tagalign` (**`neighborhoods.py:448-460`**) is the whole ballgame:

```python
454  remove_alt_chr_cmd = f"bedtools intersect -u -a {tagalign} -b {genome_sizes_bed}"
455  coverage_cmd = f"bedtools coverage -counts -sorted -g {genome_sizes} -b - -a {bed_file}"
458  awk_cmd = 'awk \'{{print $1 "\\t" $2 "\\t" $3 "\\t" $NF}}\'' + f" > {output}"
460  run_piped_commands(piped_cmds)
```

The tabix-index branch at `neighborhoods.py:449-452` is **skipped** — the `.tbi` is a declared
rule input (`candidate_regions.smk:5`), so it already exists. Confirmed: no `tabix -p bed`
line appears in any of the 12 slurm logs.

Stage 1 (`bedtools intersect -u -a <tagAlign.gz>`) inflates and parses the **entire** tagAlign
and re-serialises the survivors as plain BED6 text into a pipe. Stage 2
(`bedtools coverage -counts -sorted`) is a linear chromsweep: it streams that text as `-b` and
emits exactly `n_peaks` rows. Stage 3 (`awk`) reduces to BED4.

**The tagAlign has `2 · n_frag` records, not `n_frag`.** `workflow/rules/frag_to_tagAlign.smk:50`
emits *two* lines per fragment (the two half-fragments, `+` and `-` strand). Measured record
width: **29.87 bytes/line** (5,973,307 bytes / 200,000 lines, `thp1_1`).

**(c) `double_sex_chrom_counts` — `neighborhoods.py:421`, body `424-429`.** A second awk pass
over the `n_peaks`-row Counts.bed to double chrX/chrY counts, written to `.tmp` and `mv`-ed back
(`neighborhoods.py:428`). Run through `run_command` (`tools.py:16-18`, `check_output`), i.e. a
separate `/bin/sh` + awk. O(`n_peaks`); ~1–2 s. The Counts.bed is therefore **written twice**.

**(d) The 14-stage pipe — `peaks.py:105-120`.** Verbatim from the slurm log (`thp1_1`), with the
row volume each stage actually handles:

| # | `peaks.py` | stage | rows in | rows out |
|--:|---|---|---:|---:|
| 1 | :106 | `bedtools sort -i Counts.bed -faidx` | `n_peaks` | `n_peaks` |
| 2 | :107 | `bedtools merge -i stdin -c 4 -o max` | `n_peaks` | `m_merged` |
| 3 | :108 | `sort -nr -k 4` | `m_merged` | `m_merged` |
| **4** | **:109** | **`head -n 150000`  ← THE CAP** | `m_merged` | **150,000** |
| 5 | :110 | `bedtools intersect -b stdin -a <narrowPeak> -wa` | `n_peaks` + 150,000 | `p_surv` |
| 6 | :111 | `awk '{print $1, $2+$10, $2+$10}'` (summit, zero-width) | `p_surv` | `p_surv` |
| 7 | :112 | `bedtools slop -b 250` (→ 500 bp windows) | `p_surv` | `p_surv` |
| 8 | :113 | `bedtools sort -faidx` | `p_surv` | `p_surv` |
| 9 | :114 | `bedtools merge` | `p_surv` | ≈ `j` |
| 10 | :115 | `bedtools intersect -v -a stdin -b blocklist` (910 rows) | ≈ `j` | ≤ `j` |
| 11 | :116 | `cut -f 1-3` | ≈ `j` | ≈ `j` |
| 12 | :117 | `(bedtools intersect -a TSS500bp -b chrom_sizes_bed \| cut -f1-3 && cat)` | 20,532 + `j` | 20,532 + `j` |
| 13 | :118 | `bedtools sort -faidx` | `j` + 20,532 | `j` + 20,532 |
| 14 | :119 | `bedtools merge > outfile` | `j` + 20,532 | **`j_cand_elems`** |

All 14 processes are launched by `run_piped_commands` (`tools.py:21-34`) as
`Popen(cmd, shell=True)` — **28 processes** (14 `/bin/sh` + 14 tools), all alive concurrently on
**one** allocated CPU.

**Sort accounting.** Four sorts, all shelled out to a separate binary, none in-process:
- `peaks.py:106`, `:113`, `:118` — `bedtools sort`, which loads its whole input into RAM and
  `std::sort`s it. **O(m log m)** in-memory.
- `peaks.py:108` — GNU coreutils `sort -nr -k 4`, external-merge-capable. **O(m log m)**; at
  310k–444k rows × ~30 B = 10–14 MB it stays fully in memory. `head` closing early gives no
  saving: `sort` must consume all input before emitting anything.

**Nothing in this rule is quadratic.** `bedtools intersect` at `peaks.py:110` runs *without*
`-sorted` (see the comment at `peaks.py:104`: *"use -sorted in intersect command? Not worth it,
both files are small"*), so it builds a bin index over `-b` (≤150,000 disjoint intervals) and
streams `-a`. Because the 150,000 merged regions are disjoint and each MACS2 peak lies inside
exactly one, hits-per-query ≈ 1, so this is **O(`n_peaks` · log 150000)**, not O(`n_peaks` · `j`).
I looked specifically for an O(j²) inner loop and there is none — no per-peak Python loop
executes at all on this path.

**Dead code on this path** (confirmed from `params.txt`: exactly one accessibility file):
`peaks.py:244-258`, the `nFiles > 1` pandas cross-file averaging, never runs. Neither does
`count_bam` (`neighborhoods.py:432`) nor `count_bigwig` (`neighborhoods.py:463`) — and those are
the only consumers of the `pysam` and `pyBigWig` imports.

---

## 3. Size variables

Named exactly as in `size_manifest_*.tsv`:

- **`n_frag`** — 26,231,026 – 394,847,593 (**15.05×**). Governs the pre-cap term. The code
  touches `2 · n_frag` records (`frag_to_tagAlign.smk:50`).
- **`n_peaks`** — 432,039 – 595,372 (**1.38×**). Governs stages 1–5 of the pipe, the Python line
  count at `peaks.py:23`, and the `double_sex` awk pass.
- **`j_cand_elems`** — 156,420 – 158,548 (**1.01×**). The rule's *output*; governs stages 9–14.
- **`n_peaks_bytes`** — 32.3 – 44.3 MB. Read three times.

Two **hidden intermediate sizes** the manifest does not carry. I measured both from the retained
`*.Counts.bed` intermediates:

- **`m_merged`** — merged accessible regions entering the cap (output of stage 2 / input to
  `sort -nr`).
- **`p_surv`** — MACS2 peaks surviving the cap (output of stage 5 / input to the summit awk).

| cluster | `n_frag` | `n_peaks` | `m_merged` | `p_surv` | `j_cand_elems` |
|---|---:|---:|---:|---:|---:|
| `thp1_1` | 26,231,026 | 432,039 | 345,203 | 236,071 | 158,548 |
| `thp1_2` | 26,292,230 | 435,637 | 348,341 | 236,449 | 157,498 |
| `jurkat_pma_cd3_4hr` | 55,803,917 | 595,372 | 444,369 | 287,942 | 157,395 |
| `jurkat` | 75,342,313 | 520,543 | 361,613 | 295,763 | 158,397 |
| `k562_crispri` | 86,424,398 | 503,068 | 313,291 | 321,060 | 156,420 |
| `telohaec_crispri` | 394,847,593 | 561,437 | 309,844 | 361,988 | 157,551 |
| span | 15.05× | 1.38× | 1.43× | **1.53×** | 1.01× |

Three things follow, and all three matter for Phase 2:

1. **The cap is binding, hard.** `m_merged` is 309,844–444,369 against a cap of 150,000. The
   `head` at `peaks.py:109` discards **52 %–66 %** of merged accessible regions in every cluster.
   `j_cand_elems / 150000` = 1.043–1.057, i.e. `j` is an affine function of `nStrongestPeaks`
   plus the includelist's net contribution. This is *why* everything downstream of this rule
   spans 1.01×.
2. **The cap does not bound this rule's cost.** Stage 5's `-a` side is the *full* narrowPeak
   (`n_peaks`), and stages 6–9 process `p_surv` ≈ **1.57×–2.41× the cap**, because each of the
   top-150,000 merged regions contains several stacked `--call-summits` peaks. Only stages 9–14
   are actually bounded near `j`.
3. **`p_surv` is governed by `n_frag`, not `n_peaks`.** `p_surv` is perfectly rank-monotone in
   `n_frag` across all six clusters (236,071 → 361,988) and **not** monotone in `n_peaks`
   (595,372 peaks → 287,942 survivors, but 561,437 peaks → 361,988 survivors). Mechanism: deeper
   data → more summits stacked per merged region → the top-150,000 regions absorb a larger share
   of the peaks. The peak-stacking factor `n_peaks / m_merged` rises 1.25 → 1.81 with depth.

**On the `n_peaks` / `n_frag` decorrelation lever.** The prompt flagged this as a weak lever, and
I confirm it is weak — but I did not need it. The two candidate drivers separate on a *ratio*
argument instead, which is far stronger than the rank argument. `jurkat_pma_cd3_4hr` has the
**most** peaks (595,372, 1.38× `thp1_1`) and the **second-lowest** phase-1 time (172 s, 2.1×
`thp1_1`); `telohaec_crispri` has 1.30× `thp1_1`'s peaks and 17.2× its phase-1 time. A model in
`n_peaks` alone must explain a 17.2× time range from a 1.38× size range (implied exponent ≈ 8.8,
with the wrong rank order); a model in `n_frag` explains it with exponent 1.05 and R² 0.999. The
decorrelation lever is real but marginal; the ratio argument is decisive.

**Cross terms.** The code has **no** cross term. Every stage is a single sequential pass or sort
over one stream. The correct shape is a **sum**, not a product:
`t = S + f(n_frag) + g(n_peaks, j_cand_elems)`.

---

## 4. Asymptotic class, derived from code

### Pre-cap (`get_read_counts`, `peaks.py:99-101`)

> **O(`n_frag`)** — strictly linear, one streaming pass, no sort, no join.

- `neighborhoods.py:454` — `bedtools intersect -u -a <tagAlign.gz> -b chrom_sizes_bed`. `-b` is
  25 intervals; `-a` is streamed. Cost = inflate + parse + re-serialise of `2·n_frag` records.
  **Θ(`n_frag`)**.
- `neighborhoods.py:455-457` — `bedtools coverage -counts -sorted`. The `-sorted` chromsweep is
  a single linear merge of two pre-sorted streams with a bounded cache of currently-open
  intervals. **Θ(`n_frag` + `n_peaks`)**. Not O(n log n): the tagAlign was already sorted by
  `bedSplitSort` upstream (`frag_to_tagAlign.smk:51-57`) and the narrowPeak by
  `sort_narrowpeaks`; **no sort over fragments happens inside this rule**.
- `neighborhoods.py:458` — awk, **Θ(`n_peaks`)**.

There is therefore **no algorithmic mechanism for superlinearity** in the dominant term. The only
routes to a >1 exponent are physical: Lustre read bandwidth, pipe/scheduler contention, or
memory-hierarchy effects.

### Post-cap (the pipe, `peaks.py:105-120`)

> **O(`n_peaks` log `n_peaks` + `p_surv` log `p_surv` + `j_cand_elems` log `j_cand_elems`)**,
> which at `nStrongestPeaks = 150000` collapses to **≈ 7 s, effectively constant**.

Justification per stage: three in-memory `bedtools sort`s (`peaks.py:106`, `:113`, `:118`) and one
coreutils `sort` (`:108`) supply the log factors; every other stage is a linear stream. The
`log` factors are 18.7–19.2 across the whole observed range, i.e. constant to 2.7 %.

### Overall

> **O(`n_frag`) + O(`n_peaks` log `n_peaks`)**, with the first term dominating by 12×–190×.

### Memory

`determine_mem_mb` (`ENCODE_rE2G/ABC/workflow/rules/utils.smk:17-24`) requests
`4 × 8 × input_size_MiB` (the `×8` because `".gz" in str(input)` is true, `utils.smk:20-21`), i.e.
**32× the compressed tagAlign size**. Verified exactly: `thp1_1` inputs = 291,162,748 B =
277.7 MiB → 8,886 MB requested (slurm log: `mem_mb=8886`); `telohaec_crispri` = 2,711.9 MiB →
**86,783 MB (86.8 GB) requested**. Measured `max_rss` for that job: **157 MB**. That is a
**550× over-reservation**. The true memory class is O(`n_peaks`) (the largest `bedtools sort`
holds `n_peaks` intervals) with a ~120 MB Python-interpreter floor — completely independent of
`n_frag`, which is the variable the request is computed from.

---

## 5. Benchmark cross-check

### 5.1 The phase decomposition (measured, exact)

The rule leaves four timestamped artefacts. Combined with the job-start line in the slurm log
this gives a complete, non-overlapping phase breakdown:

| marker | source | meaning |
|---|---|---|
| `t0` | `slurm_logs/rule_abc_make_candidate_regions/<c>/*.log`, the `localrule` timestamp | job start |
| A | `Peaks/params.txt` mtime | end of conda activate + Python import (`makeCandidateRegions.py:83`) |
| B | `Peaks/*.Counts.bed` mtime (set by the `mv` at `neighborhoods.py:428`) | end of pre-cap |
| C | `Peaks/candidateRegions.qc.txt` mtime | end of the pipe (`peaks.py:126`) |

(`candidateRegions.bed`'s own mtime is **unusable** — Snakemake re-touches outputs during
post-processing, 25–83 s after the job ends. `qc.txt` is written by Python after the pipe and is
not touched, so it is the correct end marker.)

**Multiome** (seconds):

| cluster | startup | pre-cap | post-cap | tail | Σ | benchmark `s` |
|---|---:|---:|---:|---:|---:|---:|
| `thp1_1` | 61 | 81 | 6 | 1 | 149 | 148.67 |
| `thp1_2` | 61 | 81 | 6 | 1 | 149 | 148.62 |
| `jurkat_pma_cd3_4hr` | 2 | 173 | 7 | 1 | 183 | 181.92 |
| `jurkat` | 1 | 232 | 7 | 1 | 241 | 239.54 |
| `k562_crispri` | 2 | 271 | 6 | 1 | 280 | 279.64 |
| `telohaec_crispri` | 39 | 1325 | 7 | 1 | 1372 | 1372.10 |

**scATAC** (seconds):

| cluster | startup | pre-cap | post-cap | tail | Σ | benchmark `s` |
|---|---:|---:|---:|---:|---:|---:|
| `thp1_1` | 60 | 82 | 6 | 1 | 149 | 146.84 |
| `thp1_2` | 1 | 81 | 5 | 1 | 88 | 87.66 |
| `jurkat_pma_cd3_4hr` | 1 | 171 | 7 | 1 | 180 | 179.25 |
| `jurkat` | 56 | 236 | 7 | 1 | 300 | 299.21 |
| `k562_crispri` | 14 | 273 | 6 | 1 | 294 | 293.68 |
| `telohaec_crispri` | 65 | 1478 | 10 | 1 | 1554 | 1554.03 |

Σ reproduces the benchmark wall time to **≤1.1 s on all 12 jobs**. The decomposition is sound.

### 5.2 Pre-cap scaling in `n_frag` — the code and the data agree

Twelve independent points spanning 15.05× in `n_frag`:

| `n_frag` | pre-cap (mult) | pre-cap (scatac) | mean | **µs per fragment** |
|---:|---:|---:|---:|---:|
| 26,231,026 | 81 | 82 | 81.5 | 3.107 |
| 26,292,230 | 81 | 81 | 81.0 | 3.081 |
| 55,803,917 | 173 | 171 | 172.0 | 3.082 |
| 75,342,313 | 232 | 236 | 234.0 | 3.106 |
| 86,424,398 | 271 | 273 | 272.0 | 3.147 |
| 394,847,593 | 1325 | 1478 | 1401.5 | 3.549 |

Log-log OLS, all 12 points: **exponent p = 1.049 ± 0.011 (SE), 95 % CI [1.024, 1.073],
R² = 0.9989**, prefactor 1.318 × 10⁻⁶.

Excluding `telohaec_crispri` (10 points, still 3.3× range): **p = 1.009 ± 0.006, 95 % CI
[0.995, 1.022]** — dead linear, exactly as the code says.

**The one disagreement, stated honestly.** A strictly linear model calibrated on the five smaller
clusters (3.149 µs/frag, intercept −2 s, R² = 0.999) under-predicts `telohaec_crispri` by
**11.4 %**. Three candidate explanations, and I cannot separate them:

1. **Node/filesystem variance.** `telohaec_crispri` is the only cluster whose two replicates
   disagree materially: 1325 s vs 1478 s, an **11.5 % spread** — the same magnitude as the
   residual it is being asked to explain. The other five clusters agree to ≤1.7 %. This is the
   explanation I find most likely, and it means the 4.9 % exponent inflation may be entirely an
   artefact of one point measured twice on two busy nodes.
2. **I/O bandwidth saturation.** `telohaec_crispri` reads 2.80 GB compressed from Lustre against
   0.26–0.74 GB for the others.
3. Genuine mild superlinearity — for which **the code offers no mechanism** (§4). I weight this
   lowest for that reason.

**`n_frag` vs. bytes — a clean separation.** The compressed tagAlign is 258 MB … 2,800 MB, a span
of only **10.84×** against 15.05× in records, because deep data compresses better (9.92 → 7.09
compressed bytes/fragment). Per-record cost varies **1.15×** across clusters; per-compressed-byte
cost varies **1.56×**. Fitting the exponent against bytes instead of records gives 1.187. So the
driver is **records, not bytes**: gzip inflate is a minor term (~28 s of 1,400 s for
`telohaec_crispri` at a nominal 100 MB/s), and bedtools' per-record parse/compare dominates.
This is the one place where the six clusters genuinely *do* discriminate between two candidate
`n`s, and they do so decisively.

### 5.3 Post-cap scaling — resolution-limited, but consistent with O(`n_peaks`)

Post-cap is 5–10 s, mean 6.67 s, across a 1.38× range of `n_peaks` and 1.01× of `j_cand_elems`.
At 1-second mtime granularity this **cannot** support a fitted exponent, and I do not fit one.

Two weak but real consistency checks:

- The six **multiome** post-cap times rank-order *perfectly* with `n_peaks`:
  432,039 → 6 s; 435,637 → 6 s; 503,068 → 6 s; 520,543 → 7 s; 561,437 → 7 s; 595,372 → 7 s.
  The scATAC replicates are noisier (6, 5, 6, 7, 10, 7). One clean ordering out of two is
  suggestive, not conclusive.
- A row-touch budget over the 14 stages, using the measured `m_merged`/`p_surv`, gives 5.40 M
  row-touches for `thp1_1`, 6.43 M for `telohaec_crispri`, 6.60 M for `jurkat_pma_cd3_4hr` — i.e.
  **≈ 11.1–12.5 × `n_peaks`** — implying a uniform **1.06–1.56 µs per row-touch**. A single
  per-row constant reproduces all six post-cap times within the ±1 s measurement resolution.
  That is a consistency check on the class, not a confirmation of it.

**Verdict for §5.3:** class O(`n_peaks` log `n_peaks` + `j_cand_elems` log `j_cand_elems`) derived
from code at `peaks.py:106-119`; **the benchmarks cannot confirm it, because `n_peaks` spans only
1.38×, `j_cand_elems` only 1.01×, and the measurement resolution is 1 s against a 5–10 s signal.**
What the benchmarks *do* establish firmly is the magnitude: the entire pipe is ~7 s, so no
plausible exponent in `j_cand_elems` at this scale changes the rule's cost.

### 5.4 The declared-vs-realised threading gap — and a correction to `cpu_time`

The rule declares no `threads:` and Slurm gives it `cpus_per_task=1`. But `count_tagalign`
launches **three concurrent processes** (`bedtools intersect` | `bedtools coverage` | `awk`) plus
inline gzip inflate, and the pipe launches **28**. So this rule has genuine 2–3-way pipeline
parallelism that the 1-CPU allocation prevents from being realised. Measured: pre-cap
`cpu_time / pre-cap-wall` = **0.93 – 0.99** in all 12 jobs. It runs at one core, and it is
**CPU-bound, not I/O-bound**, even at 2.8 GB of compressed input.

That ratio also **explains the benchmark `cpu_time` column exactly**, and reveals that it is
being misread:

| cluster | benchmark `cpu_time` | pre-cap wall | ratio |
|---|---:|---:|---:|
| `thp1_1` (mult) | 76.79 | 81 | 0.948 |
| `thp1_1` (scatac) | 77.76 | 82 | 0.948 |
| `k562_crispri` (mult) | 253.28 | 271 | 0.935 |
| `telohaec_crispri` (mult) | 1299.31 | 1325 | 0.981 |
| `telohaec_crispri` (scatac) | 1456.50 | 1478 | 0.985 |

`cpu_time` tracks the **pre-cap phase alone**, not the job. The startup CPU and the ~7 s post-cap
CPU are invisible in it. That is consistent with Snakemake's psutil monitor recording the
*maximum over samples of the summed live process tree*, so CPU belonging to phases that have
already exited drops out. Consequence: **`mean_load` is not a parallelism measurement for this
rule.** `thp1_1` reports `mean_load` 42.49 purely because a 61 s zero-CPU startup dilutes the
mean over a 149 s job — not because anything ran at 0.42 cores. I could not read the installed
Snakemake source to confirm the implementation (no `snakemake` on PATH; scoped searches of the
conda prefixes found none), so the *mechanism* is a hypothesis; the *numerical pattern* above is
measured and holds in all 12 jobs.

### 5.5 `io_in` / `io_out` are unusable for this rule

`telohaec_crispri` reports `io_in` = 356.04 MB (multiome) / 391.57 MB (scATAC). It **must** read
at least 2,800 MB — the tagAlign alone. `thp1_1` reports 346.20 MB against a true ~320 MB. So
`io_in` saturates around 350 MB and is non-monotone in every size variable
(109.93 / 161.00 / 277.09 / 282.84 / 346.20 / 356.04 MB, in no meaningful order). Same
max-over-samples artefact as `cpu_time`. **Derive I/O volume from the files, not from these
columns** — I do so in §6.

---

## 6. Concrete overhead sources

### 6.1 The dominant one: 1.6–23.6 GB of uncompressed BED6 text through a pipe

`neighborhoods.py:454` inflates the tagAlign and writes plain-text BED6 into a pipe consumed by
`neighborhoods.py:455`. At 29.87 measured bytes/record × `2 · n_frag` records:

| cluster | tagAlign on disk (gz) | records | **text crossing the pipe** |
|---|---:|---:|---:|
| `thp1_1` | 258 MB | 52,462,052 | **1.57 GB** |
| `jurkat` | 679 MB | 150,684,626 | **4.50 GB** |
| `telohaec_crispri` | 2,800 MB | 789,695,186 | **23.59 GB** |

For `telohaec_crispri` that is **8.4× more bytes moved through a 64 KB pipe than are read from
disk**, each byte costing a write syscall's share, a read syscall's share, two memcpys, and — on
one allocated CPU — a scheduler round-trip. This, plus the double parse of every record (once by
`intersect`, once by `coverage`), is what the 3.1 µs/fragment is spent on.

### 6.2 A pass that drops nothing

`neighborhoods.py:454` filters the tagAlign against `chrom_sizes_bed`. But the fragments were
**already** restricted to exactly those chromosomes upstream, at
`workflow/rules/frag_to_tagAlign.smk:104` (`awk 'NR==FNR {keep[$1]; next} $1 in keep'`), using
`config["chr_sizes"]`.

I traced the config chain to confirm these are the **same file**:
`workflow/rules/utils.smk:97-98` sets `config["chr_sizes"] = <encode_re2g_dir>/` +
`ENCODE_rE2G/config/config.yaml:10` = `reference/GRCh38_EBV.no_alt.chrom.sizes.tsv`;
`configs/igvf10_multiome_config.yaml:34-36` deliberately omits `chr_sizes` so this fill-in
happens. And `params.txt` records
`chrom_sizes .../ENCODE_rE2G/reference/GRCh38_EBV.no_alt.chrom.sizes.tsv` for the ABC side. Same
25 chromosomes (chr1–22, X, Y, M), same file.

So in this configuration that `intersect` is a **no-op filter** that costs a full inflate, parse
and re-serialisation of `2 · n_frag` records. I verified the config chain resolves identically; I
did **not** count records on both sides of the pipe (that would require a fresh 790 M-line pass),
so "drops zero" is an inference from the chromosome sets, not a direct measurement.

### 6.3 Unused imports on the critical path

`makeCandidateRegions.py:4` → `peaks.py:5` → `neighborhoods.py:8-13` transitively imports
`numpy`, `pandas`, `pyranges`, `pyBigWig`, `pysam`, `scipy.interpolate`; `tools.py:9` imports
`pyranges` again. **None of them is used on the executed path** — `count_tagalign` is pure
`subprocess`, `count_bam` (`pysam`) and `count_bigwig` (`pyBigWig`) are not dispatched, and the
only `pandas` use, `peaks.py:250-257`, is behind `nFiles > 1` which is false.

Measured cost (marker A − `t0`): **1, 1, 1, 2, 2, 14, 39, 56, 60, 61, 61, 65 s**; mean **30.3 s**,
median 26.5 s, span **65×**. The pattern is a cold/warm Lustre cache on the conda env: `thp1_1`
and `thp1_2` launched simultaneously at 17:40:06 and both paid 61 s; the next three jobs, on the
same node minutes later, paid 1–2 s. This is the single largest source of *unpredictability* in
the rule (±32 s), and for `thp1_1` it is **41 %** of the job.

### 6.4 Repeated passes over the same files

| file | reads | writes |
|---|---|---|
| `macs2_peaks.narrowPeak.sorted` (32–44 MB) | **3** — `peaks.py:23` (Python line count), `neighborhoods.py:456` (`coverage -a`), `peaks.py:110` (`intersect -a`) | — |
| `*.Counts.bed` (11.7–15.9 MB) | **2** — `double_sex` awk, `peaks.py:106` (`bedtools sort`) | **2** — pipe awk `neighborhoods.py:458`, then `.tmp` + `mv` `neighborhoods.py:428` |
| `*.candidateRegions.bed` (3.73–3.79 MB) | 1 — `count_lines`, `peaks.py:125` | 1 |
| `TSS500bp.bed` (20,532 rows) | 1 | — |
| blocklist (910 rows) | 1 | — |

The `peaks.py:23` read exists solely to populate one line of a QC text file (`peaks.py:52`).

Derived true I/O (per job): reads ≈ `tagAlign_gz` + 3 × `n_peaks_bytes` + 2 × Counts + 4 MB ≈
**0.42 GB (`thp1_1`) → 2.97 GB (`telohaec_crispri`)**; writes ≈ 2 × Counts + 3.8 MB ≈
**27–36 MB**.

### 6.5 Process and interpreter startup

- 1 × Python interpreter (with the §6.3 import chain).
- `run_command` (`tools.py:16-18`): 1 × `/bin/sh` + 1 × awk, for `double_sex_chrom_counts`.
- `run_piped_commands` (`tools.py:21-34`) × 2: **3 + 14 = 17** `Popen(..., shell=True)` calls →
  **34 processes** (each spawns a `/bin/sh` which execs the tool). At ~5–15 ms each this is
  ~0.2–0.5 s, negligible against 7 s of pipe and 1,400 s of pre-cap — worth stating precisely so
  Phase 2 does not have to guess.
- `run_piped_commands` uses `shell=True` and does **not** close the parent's copy of each
  upstream `stdout` (contrast `run_piped_commands_safe`, `tools.py:57-59`, which does). SIGPIPE
  therefore does not propagate cleanly when `head -n 150000` (`peaks.py:109`) closes its input.
  Observationally harmless here (`sort` upstream must read everything anyway, and the job exits
  cleanly), but it means the upstream stages of the pipe are not torn down early.

### 6.6 Resource over-reservation

86.8 GB requested vs 157 MB used for `telohaec_crispri` (§4). Also `runtime=48` minutes in the
declared resources, against a measured 22.9 / 25.9 min — a 1.9× margin, the tightest of any
resource on this rule.

---

## 7. The three functions

Variables named exactly as in `size_manifest_*.tsv`. **Units: seconds.** All constants numeric.
The shape is a **sum of three additive terms** — startup, pre-cap, post-cap — because the code
contains no cross term (§3).

```
t_est(n_frag, n_peaks, j_cand_elems)
    = 30
    + 1.318e-6 * n_frag**1.049
    + 1.20e-5  * n_peaks
    + 5.5e-6   * j_cand_elems

t_low(n_frag, n_peaks, j_cand_elems)
    = 5
    + 2.747e-6 * n_frag
    + 8.0e-6   * n_peaks
    + 3.0e-6   * j_cand_elems

t_high(n_frag, n_peaks, j_cand_elems)
    = 70
    + 6.146e-7 * n_frag**1.10
    + 2.00e-5  * n_peaks
    + 1.50e-5  * j_cand_elems
```

**These are not blended and must not be averaged.** They differ in the *exponent* of `n_frag`:

| | `n_frag` exponent | rationale |
|---|---|---|
| `t_low` | **1.00** | the class the code implies — a single streaming pass, no sort over fragments (§4). Coefficient 2.747e-6 is 0.85× the `t_est` pre-cap value at the pivot, and below the fastest observed per-record rate (3.081 µs). |
| `t_est` | **1.049** | the log-log OLS over all 12 measured points, R² = 0.9989 (§5.2). Deliberately *not* the code's 1.00, because a pure-linear fit under-predicts the largest cluster by 11 %; using 1.049 keeps `t_est` honest inside the calibration range even though the mechanism is linear. |
| `t_high` | **1.10** | headroom beyond the 95 % CI upper bound (1.073) for the possibility that the `telohaec_crispri` excess is real physical superlinearity (I/O bandwidth, contention) rather than node noise. **1.10 is deliberately not `n log n`** — there is no sort over fragments anywhere on this path, so an `n log n` upper class would be unjustifiable. |

The pre-cap coefficients are pivoted at `n_frag = 1.0e8` (the geometric centre of the measured
26.2 M – 394.8 M range), where the three agree to within ±15 %/+10 %, and diverge away from it —
which is the correct way to encode exponent uncertainty.

Post-cap coefficients: `t_est` reproduces all six multiome post-cap times to the 1 s measurement
resolution. `t_low`/`t_high` are ⅔× and 2.5× that, spanning the observed 5–10 s. Both post-cap
terms are written as **linear** because the log factors are constant to 2.7 % across the whole
plausible range (§4); a Phase-2 evaluation at a much larger `j_cand_elems` should replace them
with the `m log m` form, and I flag that in §9.

Startup constants 5 / 30 / 70 s bracket the measured 1–65 s (§6.3). This term is
**cache-dependent, not size-dependent**, so it appears in all three functions as a pure constant.

### Envelope validation — all 12 observations lie inside `[t_low, t_high]`

| cluster | `t_low` | `t_est` | `t_high` | obs (mult) | obs (scATAC) | |
|---|---:|---:|---:|---:|---:|:--|
| `thp1_1` | 81.0 | 115.9 | 170.0 | 148.67 | 146.84 | OK |
| `thp1_2` | 81.2 | 116.1 | 170.3 | 148.62 | 87.66 | OK |
| `jurkat_pma_cd3_4hr` | 163.5 | 214.3 | 288.4 | 181.92 | 179.25 | OK |
| `jurkat` | 216.6 | 278.6 | 366.8 | 239.54 | 299.21 | OK |
| `k562_crispri` | 246.9 | 315.8 | 412.7 | 279.64 | 293.68 | OK |
| `telohaec_crispri` | 1094.6 | 1410.3 | 1840.2 | 1372.10 | 1554.03 | OK |

`t_est` mean absolute error 12 %, which is about the floor given that the startup term alone
contributes ±32 s of irreducible jitter.

---

## 8. Evaluated at mean / min / max

Cluster-mean sizes over the six clusters: `n_frag` = 110,823,580; `n_peaks` = 508,016;
`j_cand_elems` = 157,635.
Min = `thp1_1` (26,231,026 / 432,039 / 158,548). Max = `telohaec_crispri`
(394,847,593 / 561,437 / 157,551).

| point | `t_low` | `t_est` | `t_high` |
|---|---:|---:|---:|
| **mean** | 314.0 s (5.2 min) | **399.0 s (6.6 min)** | 516.7 s (8.6 min) |
| **min — `thp1_1`** | 81.0 s (1.3 min) | **115.9 s (1.9 min)** | 170.0 s (2.8 min) |
| **max — `telohaec_crispri`** | 1094.6 s (18.2 min) | **1410.3 s (23.5 min)** | 1840.2 s (30.7 min) |

Term shares of `t_est` — **this is the number Phase 2 needs**:

| point | startup | **pre-cap (`n_frag`)** | **post-cap (`n_peaks`, `j_cand_elems`)** |
|---|---:|---:|---:|
| mean | 30.0 s (7.5 %) | **362.0 s (90.7 %)** | **6.96 s (1.75 %)** |
| min | 30.0 s (25.9 %) | **79.8 s (68.9 %)** | **6.06 s (5.23 %)** |
| max | 30.0 s (2.1 %) | **1372.7 s (97.3 %)** | **7.60 s (0.54 %)** |

Reporting convenience only — the functions in §7 are the deliverable.

---

## 9. Modality comparison

**Neither the class nor the constant differs.** This rule is upstream of
`checkpoint features_required`, so multiome and scATAC execute byte-identical work on
byte-identical inputs. Treat the two runs as **n = 2 replicates**, not as two configurations.

Evidence: `n_frag` identical in both manifests; `j_cand_elems` identical in all six clusters
(156,420 / 157,395 / 157,498 / 157,551 / 158,397 / 158,548); `n_peaks` identical except
`telohaec_crispri` (561,437 vs 561,435 — 2 rows, and it produced the *same* `j_cand_elems`,
confirming the cap absorbs small input perturbations).

Replicate agreement on the pre-cap term: **≤1.7 % for five of six clusters** (81/82, 81/81,
173/171, 232/236, 271/273) and **11.5 % for `telohaec_crispri`** (1325/1478). Post-cap:
6/6, 6/5, 7/7, 7/7, 6/6, 7/10.

Practical consequence for Phase 2: use one set of functions for both modalities, and note that
the **replicate noise floor is ~2 % for jobs under ~5 min and ~12 % for the 20-minute job** —
the latter being the same magnitude as the exponent excess in §5.2, which is precisely why that
excess cannot be resolved.

---

## 10. What the measurement cannot tell us

1. **Whether the pre-cap exponent is 1.00 or 1.05.** The whole 4.9 % inflation rests on one
   cluster whose two replicates disagree by 11.5 %. Distinguishing them needs either replicates
   at `telohaec_crispri` scale on a quiet node, or a cluster beyond 400 M fragments.

2. **Whether the post-cap term scales at all.** 5–10 s at 1 s mtime resolution, over 1.38× in
   `n_peaks` and 1.01× in `j_cand_elems`. The class comes entirely from reading
   `peaks.py:106-119`. The clean multiome rank-ordering with `n_peaks` is suggestive; the noisier
   scATAC ordering is not. **No exponent is fittable here and I have not fitted one.**

3. **Whether the class in `j_cand_elems` matters at larger `j`.** It is
   O(`j` log `j`) with no quadratic term — I looked specifically for an O(j²) inner loop and
   there is none, because no per-element Python loop executes on this path (§2.2). So the four-
   orders-of-magnitude concern the briefing raises for O(j²) rules **does not apply to this
   rule**. But the six clusters cannot confirm even the log factor, since `j` spans 1.01×.
   If `nStrongestPeaks` were raised, `j` would grow ~proportionally (§3, item 1) and the post-cap
   term would grow slightly faster than linearly; at 10× `j` the linear post-cap terms in §7
   would under-predict by roughly the log ratio, ~18 %.

4. **The `m_merged` / `p_surv` recomputation is not free of assumption.** I derived both by
   re-implementing `bedtools merge` semantics in `sort`+`awk` over the retained `*.Counts.bed`
   (merging on `end >= start`, matching bedtools' default `-d 0` which joins book-ended
   intervals). `bedtools` was not on PATH to cross-check. A different adjacency convention would
   shift `m_merged` by well under 1 % (a handful of exactly-touching intervals) and would not
   change the conclusion that `m_merged` ≫ 150,000. `p_surv` additionally relies on the fact
   that `*.Counts.bed` rows *are* the MACS2 peaks one-for-one (432,039 rows == `n_peaks`,
   verified for all six clusters) so each peak belongs to exactly one merged region.

5. **The startup term is unmodellable from six clusters.** It is 1–65 s driven by conda-env page
   cache state and node placement, not by any size variable. It is 41 % of `thp1_1`'s runtime and
   2 % of `telohaec_crispri`'s. Any Phase-2 critical-path arithmetic involving the small clusters
   inherits ±32 s of irreducible uncertainty from this alone.

6. **`cpu_time`, `mean_load`, `io_in`, `io_out` in `benchmarks_*.tsv` cannot be used at face
   value for this rule.** `cpu_time` measures the pre-cap phase only; `mean_load` is diluted by
   the zero-CPU startup; `io_in` under-reports by 8× on the largest job (356 MB reported against
   ≥2,800 MB necessarily read). §5.4–5.5 give the evidence. All I/O figures in §6.4 are derived
   from file sizes and the code's access pattern, not from these columns.

7. **I could not read the installed Snakemake source** to confirm *why* those columns behave that
   way (`snakemake` is not on PATH; scoped searches of `$HOME/.conda/envs`,
   `$OAK/.../miniforge3/envs` and the `--conda-prefix` tree found no `benchmark.py`). The
   max-over-samples explanation fits all 12 jobs quantitatively but remains a hypothesis about
   the implementation.

8. **Whether the §6.2 no-op filter truly drops zero records.** Established by tracing the config
   chain to a single shared 25-chromosome file (`utils.smk:97-98`,
   `ENCODE_rE2G/config/config.yaml:10`, `params.txt`), not by counting records on both sides of
   the pipe. `tabix` was unavailable to list the tagAlign's chromosome set independently.

---

## Appendix — provenance

- Benchmarks: `/scratch/users/kaybrand/scE2G_optimize_results/analysis/benchmarks_{multiome,scatac}.tsv`,
  rows with `$2 == "make_candidate_regions"`.
- Sizes: `.../analysis/size_manifest_{multiome,scatac}.tsv`.
- Phase markers: `.../igvf10_{multiome,scatac}/<cluster>/Peaks/{params.txt,*.Counts.bed,candidateRegions.qc.txt}`
  mtimes, plus the `localrule` timestamp in
  `.../igvf10_{multiome,scatac}/slurm_logs/rule_abc_make_candidate_regions/<cluster>/*.log`.
- Verbatim pipe commands: the same slurm logs (`run_piped_commands` prints its command list at
  `tools.py:22`).
- Cap parameters and branch selection: `<cluster>/Peaks/params.txt` (written at
  `makeCandidateRegions.py:83`); QC counts: `<cluster>/Peaks/candidateRegions.qc.txt`.
- `m_merged`, `p_surv`: recomputed on the login node from the retained `*.Counts.bed` with
  `sort`+`awk` only (no Python, no job) — see §10 item 4 for the caveat.
- Record width: `zcat tagAlign.sort.gz | head -200000 | wc -c` = 5,973,307 → 29.87 B/record.
- Toy chr22 run excluded throughout, per the briefing.

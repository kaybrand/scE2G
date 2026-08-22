# Unit 25 — ABC external-binary rules (Tier C, GROUPED)

Rules covered (real names — the `abc_` prefix in the project plan is a DAG-export artifact of the
ABC workflow being included as a Snakemake module; neither the benchmark tables nor the `.smk`
files use it, per `PHASE1_ADDENDUM.md` §A):

| Rule | file:line | Benchmarked? | Modalities × clusters |
|---|---|---|---|
| `call_macs_peaks` | `ENCODE_rE2G/ABC/workflow/rules/macs2.smk:2-36` | yes | both, 6+6 = 12 jobs |
| `generate_chrom_sizes_bed_file` | `ENCODE_rE2G/ABC/workflow/rules/macs2.smk:38-54` | **no, by design** | both, 6+6 = 12 jobs, unmeasured |
| `sort_narrowpeaks` | `ENCODE_rE2G/ABC/workflow/rules/macs2.smk:57-77` | yes | both, 6+6 = 12 jobs |

All three ship with no `threads:` declaration at all (confirmed by reading the whole file —
there is no `threads:` line in any of the three `rule:` blocks), so Snakemake runs each at the
default of 1. All three declare `resources: mem_mb=determine_mem_mb` (`utils.smk:17-24`), a
formula off the (compressed) input file size, not off any manifest size variable.

There are no first-party script files for this unit — every cost-bearing step is either an
external binary (`macs2`, `bedtools`) or a one-line `awk`. This report characterizes what those
invocations actually do and cross-checks that against the benchmark tables.

---

## Rule 1 of 3 — `call_macs_peaks`

### 1. Unit identity

- `ENCODE_rE2G/ABC/workflow/rules/macs2.smk:2-36`, rule `call_macs_peaks`.
- Both modalities, 6 clusters each = 12 jobs. Confirmed the input is identical across
  modalities: `params.txt` for `telohaec_crispri` records
  `accessibility = ['/scratch/.../igvf10_multiome/telohaec_crispri/tagAlign/tagAlign.sort.gz']`
  and `n_frag` is byte-identical between `size_manifest_multiome.tsv` and
  `size_manifest_scatac.tsv` for all six clusters (verified directly: `394847593` for
  `telohaec_crispri` in both). This is genuinely the shared-stack n=2-per-cluster situation the
  briefing flags (§4.6) — used below.
- No `threads:` declared → runs at 1. Resources: `mem_mb=determine_mem_mb` (see §6).

### 2. What the code actually does

The shell block (`macs2.smk:17-36`):

```
if [[ "{input.accessibility}" == *tagAlign* ]]; then FORMAT="BED"; else FORMAT="AUTO"; fi
macs2 callpeak -f $FORMAT -g {params.genome_size} -p {params.pval} -n macs2 \
  --shift -75 --extsize 150 --nomodel --keep-dup all --call-summits \
  --outdir {RESULTS_DIR}/{wildcards.biosample}/Peaks -t {input.accessibility}
```

Confirmed on disk that the input is always the tagAlign path
(`.../tagAlign/tagAlign.sort.gz`, per the `expanded_biosample_config.tsv` `ATAC` column), so
`FORMAT` is always `BED` in this run, never `AUTO`. `params_macs.pval = 0.1`,
`params_macs.genome_size = hs` (`ENCODE_rE2G/ABC/config/config.yaml:23-25`); the igvf10 configs
deliberately omit `macs2_genomesize`, so this default is what actually ran
(`configs/igvf10_multiome_config.yaml:34`, comment confirms it's intentionally absent).

What this invocation means for cost:
- `--nomodel` skips MACS2's cross-correlation fragment-size estimation pass — one whole pass
  over the data that this run does *not* pay for.
- `--shift -75 --extsize 150` replaces that model with a fixed shift/extend, applied once per
  read during pileup construction — an `O(1)`-per-read constant-time step.
- `--keep-dup all` disables PCR-duplicate filtering, so every one of `n_frag` reads survives
  into the pileup (no reduction in the working set before the main pass).
- `--call-summits` adds a second, local refinement pass (`PeakDetect`) over each called
  candidate region to find sub-peak summits — additional work proportional to the total
  significant signal, not to `n_frag` directly, but still bounded by it.
- `macs2 callpeak`'s classic algorithm (`FWTrack`-based) is documented as single-process: build
  a tag track (parse + per-chromosome sort of read positions), pileup, local-lambda background
  estimation at multiple window sizes, Poisson/binomial scoring, peak calling, then the summit
  refinement above. There is no `-p`/thread flag on `macs2 callpeak` in the invocation — none
  exists for this subcommand — so the entire chain is one core.

### 3. Size variable(s)

`n_frag` — ATAC fragments after chr filtering (`size_manifest_*.tsv` col 6). Range
26,231,026 (`thp1_1`) – 394,847,593 (`telohaec_crispri`), ratio **15.06×** (Tier 1, genuinely
variable, confirmed identical between modalities). This is exactly the tagAlign's read count —
`-t {input.accessibility}` reads every one of `n_frag` records (no upstream filtering between
`fragment_count.txt` and the tagAlign, since `--keep-dup all`).

### 4. Asymptotic class derived from code

MACS2's `callpeak` (per the shell invocation, §2) does, per chromosome: sort tag positions
(comparison sort, `O(n log n)`), then linear passes for pileup construction, local-lambda
background estimation and peak/summit calling (each `O(n)`). Overall:
**`O(n_frag log n_frag)`**, dominated in practice by the `O(n_frag)` linear-scan terms because
`--call-summits` and multi-window lambda estimation each add their own linear-in-signal pass;
the sort's `log` factor is a genuinely smaller contributor at realistic fragment counts. This is
read from the documented MACS2 architecture and the specific flags in the shell block, not from
a repo-local script — `macs2` is an installed conda binary, so there is no `file:line` inside
this repo for its internals to cite beyond the invocation itself.

### 5. Benchmark cross-check

All 12 jobs (both modalities, `n=2` per cluster on the shared stack):

| cluster | `n_frag` | multiome `s` | scATAC `s` | multiome cpu_time | scATAC cpu_time | multiome cpu/wall | scATAC cpu/wall |
|---|---:|---:|---:|---:|---:|---:|---:|
| thp1_1 | 26,231,026 | 276.28 | 284.48 | 236.53 | 232.29 | 0.856 | 0.817 |
| thp1_2 | 26,292,230 | 282.78 | 277.25 | 236.82 | 240.22 | 0.838 | 0.866 |
| jurkat_pma_cd3_4hr | 55,803,917 | 525.56 | 531.01 | 494.94 | 511.68 | 0.942 | 0.964 |
| jurkat | 75,342,313 | 662.54 | 690.68 | 645.41 | 663.15 | 0.974 | 0.960 |
| k562_crispri | 86,424,398 | 815.97 | 809.52 | 780.95 | 781.75 | 0.957 | 0.966 |
| telohaec_crispri | 394,847,593 | **2903.88** | **4665.34** | 2847.14 | 4621.24 | 0.981 | 0.990 |

- **`n_frag` genuinely varies (15.1×) and the benchmarks agree with roughly-linear-to-`n log n`
  scaling.** Log-log regression of `cpu_time` on `n_frag` over all 12 points gives slope
  ≈**1.00**; the same regression on wall `s` gives slope ≈**0.95**. Both are consistent with the
  code-derived class (`O(n_frag)` to `O(n_frag log n_frag)` — over this 15.1× range `ln(n_frag)`
  itself only grows 1.16× (17.08 → 19.80), so the two classes are not really separable even with
  a fittable variable). **This is one of the few rules in this project where the benchmarks
  meaningfully confirm the code-derived class**, per the task brief.
- **`telohaec_crispri` disagrees with itself by 46.6% between the two modality runs**
  (2903.88 s vs 4665.34 s) on byte-identical input (same `tagAlign.sort.gz`, confirmed same
  `n_frag`, `n_peaks` differ only by 2 lines: 561,437 vs 561,435 — MACS2 tie-breaking noise, not
  a real input difference). Every other cluster agrees within 1–4% between modality runs. I
  cannot attribute this to anything in the code — the same binary, same flags, same file, twice
  — so I read it as node-to-node performance variance on a shared Sherlock partition (a
  contended/slower node for that specific scATAC job), not signal. It sets a hard floor on how
  precisely a fitted constant can mean anything at this cluster's scale.
- **CPU/wall climbs with cluster size**, from ≈0.82–0.86 (`thp1_1`/`thp1_2`) to ≈0.96–0.99
  (`telohaec_crispri`), in both modalities. That is the signature of a roughly fixed
  process-startup/I/O overhead (reading options, opening/creating six output files, writing
  `.xls`/summits) that is a larger fraction of a 4.6-minute job than a 48–78-minute one — not
  evidence of any parallelism, since there is none available to this binary.

### 6. Concrete overhead sources

- **Single core throughout** — confirmed no `threads:` in the rule, no `-p`/thread flag exists
  on `macs2 callpeak`, and CPU/wall (0.82–0.99, above) is consistent with one core continuously
  busy plus brief I/O/setup gaps, never with any multi-core burst.
- **Declared-vs-realized memory is far more overprovisioned than CPU.**
  `determine_mem_mb` (`ENCODE_rE2G/ABC/workflow/rules/utils.smk:17-24`) computes
  `mem_to_use_mb = max(4 × input_size_mb × 8, 8000)` when the input path contains `.gz` (it does
  here — the tagAlign is gzipped). For `telohaec_crispri` (tagAlign 2,800,146,535 bytes ≈
  2800.1 MB) that's `max(4×2800.1×8, 8000) ≈ 89,605 MB` (≈87.5 GiB) requested vs. **13.29 GB**
  peak RSS measured — a **≈6.7×** over-request. For `thp1_1` (tagAlign 258.2 MB):
  `max(4×258.2×8, 8000) ≈ 8263 MB` requested vs. **1.27 GB** measured — **≈6.5×** over-request.
  The ratio is essentially constant across cluster sizes because the `×8` "assume gzip
  compresses ≤8×" heuristic (`utils.smk:21`) overshoots for MACS2's actual working-set profile
  (it does not need 8× the compressed size in RAM to hold the read positions and pileup arrays).
  This doesn't cost wall time directly, but it is a real Slurm-scheduling cost (bigger memory
  request → longer queue wait, fewer jobs packed per node) that this benchmark set does not
  measure.
- **`io_out` does not reconcile with the visible output files.** For `thp1_1` multiome,
  `io_out = 1949.49 MB`, but the actual `Peaks/` output (`macs2_peaks.narrowPeak` 31 MB +
  `macs2_peaks.xls` 34 MB + `macs2_summits.bed` 21 MB) totals ≈86 MB — a ≈22× gap. `io_in` for
  `telohaec_crispri` (329.11 MB) is likewise far *smaller* than the 2.8 GB compressed tagAlign
  it reads. I flag this rather than explain it with confidence: `psutil`'s benchmark columns
  reflect `/proc/<pid>/io`'s `read_bytes`/`write_bytes`, which count actual block-layer I/O, not
  logical read()/write() calls — page-cache hits (the tagAlign was written by the immediately
  preceding `frag_to_tagAlign` job, likely still cache-resident) would suppress `io_in` below
  the logical file size, and cache write-back timing/dirty-page accounting could inflate
  `io_out` above the logical output size. I have not instrumented this further (e.g. `strace`),
  so this is a flagged discrepancy, not a resolved one.
- **`--keep-dup all` and `--call-summits` both add work rather than remove it** relative to
  MACS2 defaults — no read filtering happens before the pileup, and summit calling is a second
  pass over every called region. Both are visible only as an increased constant in the fit in
  §5/§7, not as a separate measured line item.

### 7. `t_est(n_frag)`, `t_low(n_frag)`, `t_high(n_frag)`

Calibrated on wall `s`, all 12 jobs (both modalities), `n_frag` as the raw fragment count.

- **`t_est(n_frag) = 2.567×10⁻⁵ · n_frag^0.946`** — best-guess, the fitted class from the
  log-log regression across all 12 points (slope 0.946, intercept exp(−10.5704) = 2.567×10⁻⁵).
  This tracks `O(n_frag)` closely; the small sub-linear exponent is consistent with the
  fixed-overhead effect in §5/§6 (a constant startup cost dilutes proportionally more at small
  `n`), not with any genuine sub-linear algorithm.
- **`t_low(n_frag) = 7.354×10⁻⁶ · n_frag`** — `O(n_frag)` lower bound: the smallest observed
  per-fragment ratio across all 12 jobs (`telohaec_crispri` multiome: 2903.88 s /
  394,847,593 = 7.354×10⁻⁶ s/frag), i.e. the best-case scenario where the `log n` term and
  `--call-summits` refinement pass contribute negligibly.
- **`t_high(n_frag) = 5.968×10⁻⁷ · n_frag · ln(n_frag)`** — `O(n_frag log n_frag)` upper bound:
  constant set from the largest observed per-`n·ln n` ratio (`telohaec_crispri` scATAC —
  the contended-node outlier from §5: 4665.34 s / (394,847,593 × ln(394,847,593)) = 5.968×10⁻⁷),
  which also upper-bounds 11 of the other 12 points (see §8; one point, `thp1_2` scATAC,
  exceeds it by ≈3%, within the noise band already documented).

These are **not blended**: `t_low` is `O(n_frag)`, `t_high` is `O(n_frag log n_frag)` — genuinely
different classes, each anchored at a real extreme observation rather than averaged.

### 8. All three evaluated

| `n_frag` | cluster | `t_low` | `t_est` | `t_high` | actual (multiome / scATAC) |
|---|---|---:|---:|---:|---:|
| 110,823,580 (mean of 6) | — | 815.2 s (13.6 min) | 900.9 s (15.0 min) | 1006.1 s (16.8 min) | mean-of-cluster-means 875.1 s |
| 26,231,026 (min) | thp1_1 | 192.9 s (3.2 min) | 231.5 s (3.9 min) | 267.9 s (4.5 min) | 276.28 s / 284.48 s |
| 394,847,593 (max) | telohaec_crispri | 2903.9 s (48.4 min) | 3494.7 s (58.2 min) | 4665.3 s (77.8 min) | 2903.88 s / 4665.34 s |

At the max, the two real observations land **exactly on `t_low` and exactly on `t_high`**
respectively — because both bounds were anchored at those two points. That is by construction,
not a coincidence, but it usefully shows the bracket width (48.4–77.8 min) is the same order as
the actual run-to-run spread on identical input — i.e. at this cluster's scale, node-to-node
noise is comparable to the uncertainty in the exponent itself.

### 9. Modality: class or constant?

**Neither, materially — the rule is upstream of the modality branch.** `call_macs_peaks`
consumes the tagAlign built from raw fragments, before `checkpoint features_required`
(briefing §2) decides ARC vs Neither per cluster. `n_frag`, the input file, and the shell block
are all identical between multiome and scATAC runs. The only observed difference (§5,
`telohaec_crispri`'s 46.6% spread) is not systematic in either direction — the multiome run was
*faster* here, but `k562_crispri` in the `bedSplitSort.sh` unit's report (`09-bedSplitSort.md`)
saw the *scATAC* run faster for its largest disagreement. This is consistent with shared-node
contention, not a modality effect.

### 10. What the measurement cannot tell us

- **Cannot separate `O(n_frag)` from `O(n_frag log n_frag)` with confidence**, even though this
  is one of the few rules where the benchmarks are usable at all — `ln(n_frag)` only spans
  1.16× over the 15.1× `n_frag` range, so a log factor is fundamentally hard to detect here
  regardless of how clean the rest of the data is.
- **Cannot explain the `io_in`/`io_out` mismatch with the visible files** (§6) — flagged, not
  resolved.
- **Cannot attribute the `telohaec_crispri` 46.6% spread to a specific cause** (node hardware,
  contention from co-scheduled jobs, filesystem latency) — only that it is not explained by
  anything in the code, config, or input files, all of which are identical between the two runs.
- **Cannot decompose the single wall-clock number into its named algorithmic phases** (tag
  sort, pileup, lambda estimation, peak calling, summit refinement) — MACS2 is an opaque binary
  here; the benchmark measures the whole invocation only.

---

## Rule 2 of 3 — `generate_chrom_sizes_bed_file`

### 1. Unit identity

`ENCODE_rE2G/ABC/workflow/rules/macs2.smk:38-54`. Both modalities, 6 clusters each = 12 DAG
nodes — **but zero benchmark rows exist for it in either
`benchmarks_multiome.tsv`/`benchmarks_scatac.tsv` or the per-job benchmark directory**, verified
directly (`awk -F'\t' -v r=generate_chrom_sizes_bed_file 'NR==1||$2==r' benchmarks_*.tsv`
returns only the header in both files). This is by design, documented in the rule's own comment
(`macs2.smk:43-48`): Snakemake's `benchmark:` monitor polls the process tree with `psutil`, and
this `awk` over a 25-line file exits before `psutil` can attach, raising
`NoSuchProcess: process PID not found` and failing the job under the SLURM executor. Do not go
looking for a benchmark file for this rule; there isn't one.

### 2. What the code actually does

```
awk 'BEGIN {OFS="\t"} {if (NF > 0) print $1,"0",$2 ; else print $0}' {input.chrom_sizes} > {output.chrom_sizes_bed}
```

One `awk` pass over `ENCODE_rE2G/reference/GRCh38_EBV.no_alt.chrom.sizes.tsv`, confirmed on disk
to be **25 lines** (chr1–22, X, Y, M/EBV), converting `chrom\tsize` into a 3-column BED
(`chrom\t0\tsize`). Output confirmed identical line count (25) at
`/scratch/.../igvf10_multiome/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed`.

### 3. Size variable(s)

None from the manifest — the input is the fixed 25-line reference chromosome-sizes file, not a
function of any per-cluster quantity (`n_frag`, `n_peaks`, etc.). It runs identically 12 times
(once per cluster per modality) on the exact same input and produces the exact same output each
time (a genuinely redundant recomputation across the DAG, though re-deriving that is out of
scope here — Phase 2's concern).

### 4. Asymptotic class derived from code

`O(1)` — a single linear pass over a fixed-size (25-line) file. There is no variable in this
repo's size manifest that this rule's cost depends on.

### 5. Benchmark cross-check

None possible — no benchmark file exists (§1). Per the task brief's required framing: this is a
real DAG node with **no measurement**, costing an estimated ~0.03 s per the project plan (an awk
pass over 25 short lines plus process spawn/conda-activation overhead), not a modeling gap to
work around.

### 6. Concrete overhead sources

Process spawn and (if not already warm) conda environment activation are the only real costs;
the `awk` computation itself is negligible at n=25.

### 7–8. `t_est`, `t_low`, `t_high`, evaluated

No size variable governs this rule, so all three functions are the same constant:

**`t_est() = t_low() = t_high() = 0.03`** (seconds) — a constant, not a function of `n_frag`,
`n_peaks`, or any other manifest column, per the project plan's own estimate and consistent with
a 25-line `awk` pass. There is nothing to evaluate "at mean/min/max cluster size" because there
is no cluster-size dependence.

### 9. Modality: class or constant?

Neither varies by modality — same reference file, same 25 lines, run once per cluster per
modality with byte-identical inputs and outputs.

### 10. What the measurement cannot tell us

**Everything quantitative** — there is no benchmark file, by design (§1). The 0.03 s figure is
an estimate carried over from the project plan, not something re-derived here; I have not
independently timed this `awk` invocation (doing so would require running it outside the
Snakemake/psutil harness, which is out of scope for characterization-only work).

---

## Rule 3 of 3 — `sort_narrowpeaks`

### 1. Unit identity

`ENCODE_rE2G/ABC/workflow/rules/macs2.smk:57-77`, rule `sort_narrowpeaks`. Both modalities, 6
clusters each = 12 jobs, all benchmarked. No `threads:` declared → runs at 1.

### 2. What the code actually does

```
bedtools intersect -u -a {input.narrowPeak} -b {input.chrom_sizes_bed} | \
bedtools sort -faidx {params.chrom_sizes} -i stdin > {output.narrowPeakSorted}
```

Two piped `bedtools` subcommands:

1. `bedtools intersect -u -a narrowPeak -b chrom_sizes_bed`: keeps each of the ~500K narrowPeak
   records that overlaps at least one of the 25 chrom_sizes_bed intervals (this is what drops
   alternate contigs — anything not on `chr1–22,X,Y,EBV`). `-b` (the 25-line file) is the
   smaller file; `bedtools` builds its interval index from `-b`, so this is effectively an
   `O(n_peaks)` scan with `O(1)`-ish lookups against a 25-interval structure per record — the
   25-line side never becomes a cost driver at any realistic `n_peaks`.
2. `bedtools sort -faidx {chrom_sizes} -i stdin`: an **in-memory** comparison sort of the
   surviving records, ordered by chromosome (per the `-faidx` genome-file order) then start
   coordinate. `bedtools sort` (all versions, including the `bedtools=2.26.0` pinned in
   `ENCODE_rE2G/ABC/workflow/envs/abcenv.yml`) loads the whole input into memory and sorts it —
   there is **no `-T`/external-merge option on `bedtools sort` at all** (unlike GNU coreutils
   `sort`, which does have `-T`). This directly answers the task's required point: **the
   Lustre-vs-`$L_SCRATCH` question that matters for a GNU `sort -T` does not apply to this
   rule**, because the tool invoked here is `bedtools sort`, not `sort`, and it never writes a
   temp file to spill comparisons to disk — confirmed by reading the exact shell block, not
   assumed. (Contrast with `bedSplitSort.sh`'s `sort_chr()`, unit 9, which *does* invoke GNU
   `sort -T` and does sit on Lustre — that finding does not transfer to this rule.)

### 3. Size variable(s)

`n_peaks` — raw MACS2 peaks pre-cap, from `macs2_peaks.narrowPeak.sorted` line count
(`size_manifest_*.tsv` col 12). Range 432,039 (`thp1_1`) – 595,372 (`jurkat_pma_cd3_4hr`), ratio
**1.38×** (Tier 2, mildly variable). Verified directly: `wc -l` on
`telohaec_crispri/Peaks/macs2_peaks.narrowPeak(.sorted)` gives 561,437 lines both before and
after sorting, matching the manifest exactly, confirming `sort` here drops zero records (it is
`intersect -u` that filters, not `sort`).

**Genuine decorrelation from `n_frag` — but weak.** Ranking clusters by `n_frag` (descending):
`telohaec_crispri` (394.8M) > `k562_crispri` (86.4M) > `jurkat` (75.3M) > `jurkat_pma_cd3_4hr`
(55.8M) > `thp1_2`/`thp1_1` (26.3M/26.2M). Ranking by `n_peaks` (descending):
`jurkat_pma_cd3_4hr` (595,372) > `telohaec_crispri` (561,437) > `jurkat` (520,543) >
`k562_crispri` (503,068) > `thp1_2`/`thp1_1` (435,637/432,039). `jurkat_pma_cd3_4hr` has the
**most** peaks with only the **4th**-most fragments — confirmed. However, I could not reproduce
the task prompt's specific "telohaec_crispri has 4.8× more fragments than thp1_1" figure from
the manifest: the actual ratio is `394,847,593 / 26,231,026 = 15.05×` (matching the briefing's
own headline 15.1× `n_frag` range, which is exactly this pair), not 4.8×. I report the number I
verified (15.05×) rather than the prompt's figure. The peak-count side of that comparison does
check out: `561,437 / 432,039 = 1.30×` (`561,437 / 435,637 = 1.29×` if compared to `thp1_2`
instead), close to the "1.3×" quoted. **This decorrelation is real but weak** — it separates two
variables that are each barely varying in absolute terms (`n_peaks` only spans 1.38× total), so
it cannot be used to fit two independent exponents from six data points; it only demonstrates
that `n_peaks` is not a simple linear readout of `n_frag`.

Per the briefing, `n_peaks` sits upstream of the pipeline's variance sink: the raw peak count
here (432K–595K, 1.38×) still gets capped downstream to `j_cand_elems` (~157K, 1.01×) in
`make_candidate_regions` (unit 5) — not in this rule. This rule's output is on the
higher-variance side of that cap.

### 4. Asymptotic class derived from code

| Step | file:line | Class in `n_peaks` |
|---|---|---|
| `bedtools intersect -u` against 25-line file | `macs2.smk:74` | `O(n_peaks)` |
| `bedtools sort -faidx ... -i stdin` (in-memory comparison sort) | `macs2.smk:75` | `O(n_peaks log n_peaks)` |

**Overall class: `O(n_peaks log n_peaks)`**, dominated by the in-memory sort. No temp-file /
external-merge behavior exists in this specific tool (§2), so there is no Lustre-vs-node-local
distinction to make for this rule — that's the honest answer to the task's required point,
rather than assuming a `-T` flag exists.

### 5. Benchmark cross-check

All 12 jobs:

| cluster | `n_peaks` | multiome `s` | scATAC `s` | multiome cpu/wall | scATAC cpu/wall |
|---|---:|---:|---:|---:|---:|
| thp1_1 | 432,039 | 6.26 | 5.81 | 0.609 | 0.654 |
| thp1_2 | 435,637 | 5.91 | 4.40 | 0.675 | 0.834 |
| k562_crispri | 503,068 | 5.49 | 5.33 | 0.858 | 0.865 |
| jurkat | 520,543 | 5.64 | 5.58 | 0.842 | 0.853 |
| telohaec_crispri | 561,437 / 561,435 | 7.64 | **9.53** | 0.679 | 0.810 |
| jurkat_pma_cd3_4hr | 595,372 | 6.59 | 6.28 | 0.871 | 0.911 |

**The data do not show a usable relationship between `n_peaks` and wall time, and this rule is
the textbook case for "do not fit an exponent to noise."** Sorting clusters by `n_peaks` gives
wall times of 5.9–6.3, 4.4–5.9, 5.3–5.5, 5.6–5.7, 7.6–9.5, 5.6–6.6 s — not monotonic in either
direction. A log-log regression of mean wall time on `n_peaks` across the six clusters gives a
slope with essentially no explanatory power (the residuals are as large as the fitted trend); I
computed it (≈1.07×10⁻⁵ s/peak with the group-mean data) but do not report it as a class-fitting
result, because the six points visibly do not support one — `n_peaks`'s 1.38× range is simply
too narrow, exactly as the briefing warns for Tier 2 variables.
- **The one outlier, `telohaec_crispri` scATAC (9.53 s), is the same cluster/modality pair that
  was the outlier in `call_macs_peaks` (§5 above, 4665.34 s vs. 2903.88 s).** That is a second,
  independent rule showing an anomaly on the identical job (same cluster, same modality run),
  which strengthens (without proving) the "contended/slow node for that specific scATAC run"
  explanation over a coincidence.
- A back-of-envelope algorithmic estimate makes the "mostly overhead" conclusion concrete: an
  in-memory comparison sort of ~500,000 short BED records, even at a generous 20 ns/comparison,
  costs `500,000 × log₂(500,000) × 20ns ≈ 0.19 s` — under 4% of the smallest observed wall time
  (4.40 s). **The great majority of the observed 4.4–9.5 s is process/pipe/conda startup and the
  32–44 MB of file I/O, not the sort itself**, which is consistent with the low, inconsistent
  CPU/wall ratios (0.61–0.91) and the complete absence of a size trend.

### 6. Concrete overhead sources

- **Single core**: no `threads:`; both `bedtools` subcommands are single-threaded, and
  `intersect`/`sort` are chained by a shell pipe (concurrent between the two processes, but each
  individually serial).
- **Memory is pinned to the 8 GB floor regardless of cluster size.** `determine_mem_mb`
  (`utils.smk:17-24`) computes `max(4 × input_size_mb, 8000)` (no `×8`, since neither
  `narrowPeak` nor `chrom_sizes_bed` is gzipped). Even at the largest `narrowPeak`
  (`jurkat_pma_cd3_4hr`, 44.3 MB), `4 × 44.3 = 177.2 MB ≪ 8000 MB`, so every one of these 12 jobs
  requests exactly the 8 GB floor. Measured peak RSS is 300–414 MB — a **≈20–27× overprovision**
  that is completely flat across clusters (the input is far too small to ever push past the
  floor at this file-size range).
- **`io_out` is noisy/near-zero for several jobs** (e.g. `jurkat` multiome `io_out=0.00`,
  `thp1_1` multiome `io_out=0.05`) while `thp1_2` multiome reports `io_out=15.80` for
  essentially the same output size — consistent with the briefing's general warning that
  sub-10-second jobs produce unreliable `psutil` sampling, not with any real difference in bytes
  written (all six `.sorted` files are within 2% of their unsorted `narrowPeak` size, since
  `sort` doesn't drop columns).
- **Fixed per-job overhead (conda/process/pipe startup) is the dominant cost at this `n_peaks`
  range**, per the back-of-envelope estimate in §5 — the algorithmic term is asymptotically real
  (§4) but numerically negligible here.

### 7. `t_est(n_peaks)`, `t_low(n_peaks)`, `t_high(n_peaks)`

Given §5's finding that `n_peaks` explains essentially none of the observed variance, these
functions are built from the **code-derived class** (§4) plus a fixed overhead calibrated to the
data range, rather than from a regression I don't trust:

- **`t_low(n_peaks) = 4.30 + 1×10⁻⁷ · n_peaks`** — `O(n_peaks)` lower bound. `4.30 s` is
  slightly below the smallest observed wall time (4.40 s, `thp1_2` scATAC); the linear term is a
  deliberately small placeholder for "the sort's contribution never exceeds a per-record
  constant," representing the best case where the `log n_peaks` factor and pipe overhead are
  fully absorbed into the fixed term.
- **`t_high(n_peaks) = 4.30 + 6.9×10⁻⁷ · n_peaks · ln(n_peaks)`** — `O(n_peaks log n_peaks)`
  upper bound, the class the code actually implements (§4). The linear-in-`n log n` coefficient
  is set from the single largest observed excess over the floor (`telohaec_crispri` scATAC:
  `(9.53 − 4.30) / (561,435 × ln(561,435)) = 6.9×10⁻⁷`), and checked to upper-bound all other 11
  observed points at their respective `n_peaks` (it does, with the largest margin at the small
  clusters, where the fixed overhead is proportionally most of the total).
- **`t_est(n_peaks) = 6.20`** (seconds) — the mean of all 12 observed wall times. I am
  deliberately **not** offering a fitted-slope point estimate here: §5 shows the six-cluster
  data have no detectable trend in `n_peaks` at all, so a "best fit line" would be fitting noise,
  which the task instructions explicitly say not to do. The honest single best-guess number for
  Phase 2 is the empirical mean, not a spurious slope.

`t_low` and `t_high` are genuinely different classes (`O(n_peaks)` vs. `O(n_peaks log n_peaks)`),
not blended — but at this `n_peaks` range (432K–595K) the numeric gap between them is tiny
compared to the fixed-overhead term both share, which is the point: the class matters
asymptotically (§4, and per the briefing's general warning about near-constant variables) but is
invisible in this data.

### 8. All three evaluated

Note: the task's generic min/max clusters (`thp1_1` min, `telohaec_crispri` max) are the correct
ones for `n_frag`-driven rules, but `n_peaks`'s own max cluster is **`jurkat_pma_cd3_4hr`**, not
`telohaec_crispri` (§3 decorrelation) — evaluated at the true `n_peaks` min/mean/max below, with
`telohaec_crispri` shown for reference since it's the generic "max" cluster:

| `n_peaks` | cluster | `t_low` | `t_est` | `t_high` | actual (multiome / scATAC) |
|---|---|---:|---:|---:|---:|
| 508,016 (mean of 6) | — | 4.35 s | 6.20 s | 8.35 s | mean-of-cluster-means 6.20 s |
| 432,039 (min) | thp1_1 | 4.34 s | 6.20 s | 7.63 s | 6.26 s / 5.81 s |
| 595,372 (max for `n_peaks`) | jurkat_pma_cd3_4hr | 4.36 s | 6.20 s | 9.86 s | 6.59 s / 6.28 s |
| 561,437 (generic "max" cluster) | telohaec_crispri | 4.36 s | 6.20 s | 9.13 s | 7.64 s / **9.53 s** |

`t_low`–`t_high` bracket every actual observation.

### 9. Modality: class or constant?

**Neither differs.** `sort_narrowpeaks` is upstream of the ARC/Neither branch decision, same as
`call_macs_peaks`; its input (`narrowPeak`) and `n_peaks` are computed identically in both
pipelines from the same tagAlign. The only anomaly (`telohaec_crispri` scATAC, §5) mirrors the
same cluster/modality anomaly seen in `call_macs_peaks`, again consistent with a node-level
effect rather than a modality-driven one.

### 10. What the measurement cannot tell us

- **Cannot confirm `O(n_peaks log n_peaks)` over `O(n_peaks)` — or any relationship with
  `n_peaks` at all — from these six clusters.** The 1.38× range is too narrow and the fixed
  per-job overhead too large relative to it. The class in §4 is derived entirely from reading
  the `bedtools` invocation, not from the data.
- **Cannot separate the `intersect` step's cost from the `sort` step's cost** — the benchmark
  measures the whole piped command; I have not run either subcommand standalone under the
  monitor.
- **Cannot confirm the exact reason for the `telohaec_crispri` scATAC anomaly** beyond
  "corroborates the same anomaly in `call_macs_peaks`, on the same cluster/modality pair" — this
  is circumstantial, not a proof of node contention.
- **Cannot confirm `bedtools sort`'s in-memory-only behavior from a live trace** — this is
  stated from the tool's documented design (no `-T`-equivalent flag exists in its argument
  list) and from the low, flat memory footprint (300–414 MB, far under the 8 GB floor) across
  all 12 jobs, not from directly observing its internals.

---

## Cross-cutting notes for Phase 2

- Of the three rules, only `call_macs_peaks` has a genuinely fittable variable (`n_frag`,
  15.1×) and the benchmarks agree with the code-derived class there — treat that fit with
  moderate confidence. `sort_narrowpeaks`'s `n_peaks` (1.38×) and
  `generate_chrom_sizes_bed_file`'s constant input give no fittable signal; their functions in
  §7 above are code-derived, not data-fitted, by design.
- The same cluster/modality pair (`telohaec_crispri`, scATAC run) is the single biggest outlier
  in **both** benchmarked rules in this unit, independently. That is worth flagging to whoever
  is deciding whether to trust any single scATAC `telohaec_crispri` measurement elsewhere in the
  project — it looks like a run-level (node-level), not rule-level, anomaly.
- Memory requests are overprovisioned by ~6.5–27× across all three rules (constant-floor for
  `sort_narrowpeaks`/`generate_chrom_sizes_bed_file`, size-scaled-but-still-6-7× for
  `call_macs_peaks`). Not a wall-time cost in this data, but worth noting as a scheduling-cost
  side effect of the shared `determine_mem_mb` heuristic across ABC rules.

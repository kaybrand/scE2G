# Unit 24 — Fragment shell pipelines (Tier C, GROUPED)

## 1. Unit identity

Three rules, no first-party script — all logic lives in each rule's own `shell:` block, in
`workflow/rules/frag_to_tagAlign.smk` (168 lines):

| Rule | Definition | Modalities | Jobs | Declared `threads:` |
|---|---|---|---:|---:|
| `process_fragment_file` | `frag_to_tagAlign.smk:79-108` (indented conditional rule, see §2a) | both | 12 (6 clusters × 2) | 8 |
| `frag_to_tagAlign` | `frag_to_tagAlign.smk:9-62` | both | 12 (6 × 2) | 8 |
| `frag_to_norm_bigWig` | `frag_to_tagAlign.smk:138-167` | both | 12 (6 × 2) | 16 |

A fourth rule, `frag_to_bigWig` (`frag_to_tagAlign.smk:111-136`, unnormalized bigWig), shares
`frag_to_norm_bigWig`'s shape but is commented out of `rule all` (`workflow/Snakefile`) — 0
benchmark rows exist for it in either `benchmarks_multiome.tsv` or `benchmarks_scatac.tsv`
(confirmed by direct grep). It is not characterised further here beyond noting it never ran.

`n_frag` is byte-identical between `size_manifest_multiome.tsv` and `size_manifest_scatac.tsv`
for all six clusters (verified directly), because all three rules sit upstream of the
per-cluster ARC/Neither branch at `checkpoint features_required` (briefing §2) — so every
cluster gives a true n=2 replicate on identical input, per the briefing's "shared ABC/fragment
stack" note.

**Overlap with unit 9**: `frag_to_tagAlign` calls `workflow/scripts/bedSplitSort.sh`
(unit 9's script) at `frag_to_tagAlign.smk:51-58`. Unit 9 already derived that script's own
complexity class (`O(n_frag log n_frag)`, from the sequential per-chromosome comparison-sort
loop) and its I/O round-trip through Lustre in detail
(`phase1_reports/09-bedSplitSort.md` §2, §4, §6). I do not re-derive that here. What I add:
the surrounding rule's pipe structure (which stages besides the sort are single-threaded, and
where the hard barriers sit), and — per the task brief — an explicit statement of which DAG
node Phase 2 should charge the sort seconds to (§6d).

## 2. What the code actually does

### 2a. `process_fragment_file` — the conditional that selects it

`frag_to_tagAlign.smk:64-108` is a Python `if/else` evaluated once at **workflow-parse time**,
not a runtime/wildcard branch:

```
if config["fragments_preprocessed"]:
    rule get_fragment_count:   # :66-77
else:
    rule process_fragment_file:  # :79-108
```

Both configs shipped in this repo set `fragments_preprocessed: False`
(`config/config.yaml:13`, `config/config_training.yaml:11`), so **only `process_fragment_file`
is ever a node in the DAG** for any run using these configs — `get_fragment_count` does not
exist as a rule at all once Snakemake parses the file with that config value. This is confirmed
empirically: `benchmarks_multiome.tsv` and `benchmarks_scatac.tsv` contain zero rows for
`get_fragment_count` (checked directly) and 6+6 rows for `process_fragment_file`. `config[...]`
governs whether the *fragment input is already chr-filtered upstream*; when it isn't (`False`),
`process_fragment_file` does the chr-filtering itself and also derives `fragments_filtered`,
which is what `frag_to_tagAlign`/`frag_to_norm_bigWig` actually read via
`get_processed_fragment_file` (`frag_to_tagAlign.smk:1-6`).

Shell (`frag_to_tagAlign.smk:100-108`):

```
LC_ALL=C
awk 'NR==FNR {keep[$1]; next} $1 in keep' {params.chrSizes} <(zcat {input.frag_file}) | bgzip > {output.fragments_filtered}
tabix -p bed {output.fragments_filtered}
zcat {output.fragments_filtered} | wc -l > {output.fragment_count}
```

Three **strictly sequential** statements (no shared pipe across them):
1. `zcat` (process substitution) → `awk` (hash-filter against the ~25-line chrom-sizes file,
   loaded via the classic `NR==FNR` two-file AWK idiom) → `bgzip` (no `-@`, so single-threaded)
   — one linear pass, two/three concurrently-pipelined single-threaded processes.
2. `tabix -p bed` on the now-complete file — cannot start until (1)'s `bgzip` has closed the
   file; tabix has no threading option.
3. `zcat {output.fragments_filtered} | wc -l` — **a second, fully redundant full decompression
   pass over the file just written in step (1)**, solely to get a row count. Single-threaded
   `zcat` piped to single-threaded `wc`.

`threads: 8` is declared (`:89`) but **the literal string `{threads}` never appears anywhere in
this shell block** (grepped `:100-108` directly) — no tool invoked here (`awk`, `bgzip`,
`tabix`, `zcat`, `wc`) receives a thread or parallel flag of any kind. The declared 8 CPUs are
100% decorative for this rule; nothing in the code path even attempts to use more than one.

### 2b. `frag_to_tagAlign`

Shell (`frag_to_tagAlign.smk:43-62`):

```
zcat {input.frag_file} | \
    awk ... '{print two tagAlign records per fragment}' | \
    bash {params.bedSplitSort} -i /dev/stdin -g {chrSizes} -t $SPLIT_DIR -p {threads} -s $BUFFER_SIZE -T {temp_dir} | \
    bgzip -c > {output.tagAlign_sort_file}
tabix -p bed {output.tagAlign_sort_file}
```

- `zcat` (:49) and the fragment→tagAlign `awk` (:50) are single-threaded and pipe-connected —
  the only genuine producer/consumer overlap in this rule.
- `bedSplitSort.sh` (:51-57) receives `-p {threads}` (=8), but per unit 9's derivation the
  script itself has an internal hard barrier: a single-threaded split `awk`
  (`bedSplitSort.sh:71-75`) must consume *all* of stdin before the per-chromosome sort loop
  (`bedSplitSort.sh:87-91`) can begin, and that loop calls `sort` (the only thread-aware tool in
  the whole rule) **once per chromosome, sequentially** — never more than one chromosome's sort
  running at a time, each on a small (~1/25th-scale) partition.
- `bgzip -c` (:58, no `-@`) and `tabix -p bed` (:61, no threading option) close out the rule,
  both single-threaded, both strictly after the pipe/sort work completes.
- Net: of ~6 distinct process invocations behind `threads: 8`, exactly one (`sort`, inside
  `bedSplitSort.sh`) is thread-aware, and even that one never runs more than one chromosome's
  chunk at a time.

### 2c. `frag_to_norm_bigWig`

Shell (`frag_to_tagAlign.smk:157-166`):

```
frag_count=$(<{input.fragment_count})
scale_factor=$(awk "BEGIN {print 1000000 / $frag_count}")
zcat {input.frag_file} | \
    bedtools genomecov -bg -i stdin -g {chrSizes} -scale $scale_factor | \
    sort -k1,1 -k2,2n --parallel={threads} -S $BUFFER_SIZE -T {temp_dir} > {output.bedGraph_file}
bedGraphToBigWig {output.bedGraph_file} {chrSizes} {output.bigWig_file}
```

- `zcat` (:163) — single-threaded decompression of the full fragment file.
- `bedtools genomecov -bg` (:164) — **the one categorically un-parallelizable stage.**
  `bedtools` (pinned `2.29.2`, `workflow/envs/sc_e2g.yml:8`) has never exposed a threading
  option for `genomecov` in any released version — there is no `-p`/`--threads` flag to pass
  even if one wanted to. It is a single serial streaming pass that walks the interval stream
  once, maintaining a running per-base depth accumulator, by construction inherently sequential
  (each output interval depends on the running depth state carried from the previous one).
  It sits in the middle of the pipe, between two stages that at least *could* be threaded
  (`zcat` isn't, but `sort` is).
- `sort --parallel={threads}` (:165, threads=16) — this is the one place in all three rules
  where the thread flag is correctly wired to a real, substantial job: **the whole genome's**
  bedGraph stream, unlike `frag_to_tagAlign`'s per-chromosome-partitioned sort. Confirmed on
  disk: the intermediate `ATAC_norm.bg` for `telohaec_crispri` is **10.8 GB** uncompressed
  (`/scratch/users/kaybrand/scE2G_optimize_results/igvf10_multiome/telohaec_crispri/ATAC_norm.bg`,
  `ls -la` verified), so this is not a trivially small sort.
- `bedGraphToBigWig` (:166) — a UCSC tool with no multithreading support, run as a **separate
  sequential command after the redirection to `{output.bedGraph_file}` completes**. It cannot
  overlap with anything upstream because it needs the fully-written bedGraph file, and it makes
  its own full linear pass over that same 10.8 GB file.

Despite `sort --parallel=16` being genuinely wired to a large job here (unlike the other two
rules), measured `cpu_time/wall` still tops out at **0.99** even at the largest cluster — i.e.
no more than ~1 core's aggregate work is ever realized concurrently across the whole pipeline.
The data cannot show *why* `sort`'s 16-way parallelism didn't materialize as wall-clock savings
(§10), but the code shows the two candidate explanations precisely: (a) `bedtools genomecov` is
the strictly single-core, non-parallelizable stage feeding `sort`, so `sort` can never receive
data faster than `genomecov` can produce it, and pipe backpressure caps the *pipe's* aggregate
throughput at whatever `genomecov` alone can sustain; and (b) the separate, single-core
`bedGraphToBigWig` pass at the end adds pure serial time on top, diluting whatever multi-core
benefit `sort` achieved during its portion of the pipe. `frag_to_bigWig` (unnormalized,
`:111-136`) shares this exact shape but never ran (§1) — noted, not characterised further.

## 3. Size variable

**`n_frag`** — ATAC fragments after chr filtering (`fragment_count.txt`), Tier 1,
26,231,026 – 394,847,593, ratio **15.1×** (`size_manifest_{multiome,scatac}.tsv` column 6,
identical values in both files). All three rules make one or a small fixed number of linear
passes over the fragment stream (or a stream derived from it by a constant-factor
transformation — 2 tagAlign records per fragment in `frag_to_tagAlign`), so `n_frag` is the
correct governing variable for all three; no other manifest column enters these rules' code.
Per the addendum, `n_frag_bytes` is a digit-count measurement artifact and is not used below.

## 4. Asymptotic class derived from code

| Rule | Class in `n_frag` | Justification (file:line) |
|---|---|---|
| `process_fragment_file` | `O(n_frag)` | Three sequential single-pass streams: filter (`:104`), index (`:105`, `O(n_frag)`-ish tabix build), and a **redundant second** decompress+count pass (`:107`). No sort, no join beyond an `O(1)`-sized (~25-entry) hash lookup per record. Everything here is linear; the "redundant second pass" doubles a constant, it doesn't change the class. |
| `frag_to_tagAlign` | `O(n_frag)` (own code) `+ O(n_frag log n_frag)` (via `bedSplitSort.sh`, unit 9's derivation) `= O(n_frag log n_frag)` overall | `zcat`/tagAlign-`awk`/`bgzip`/`tabix` (`:49,50,58,61`) are each one linear pass; the sort inside `bedSplitSort.sh` is a comparison sort over `Θ(n_frag)` records split into ~25 chromosome-sized chunks, `Σᵢ mᵢ log mᵢ = Θ(n_frag log n_frag)` (unit 9, §4). I adopt that term rather than re-derive it. |
| `frag_to_norm_bigWig` | `O(n_frag)` (streaming/genomecov/output) `+ O(m log m)` (the `sort`, `m` = bedGraph line count, itself `O(n_frag)`-bounded) `= O(n_frag log n_frag)` | `zcat` (:163) and `bedtools genomecov` (:164) are each one linear pass over the fragment/interval stream (genomecov's running-depth algorithm is `O(n_frag)` — no per-record work depends on data seen more than a constant window back). `sort` (:165) is a whole-genome comparison sort over the bedGraph stream, `m ≤ O(n_frag)` lines (coverage-interval merging only shrinks the count, confirmed: `ATAC_norm.bg` for `telohaec_crispri` is 10.8 GB vs. `fragments_filtered.tsv.gz` 4.27 GB compressed — same order of magnitude, not obviously smaller). `bedGraphToBigWig` (:166) is one more linear pass over the sorted bedGraph. |

## 5. Benchmark cross-check

`n_frag` genuinely spans 15.1× here — per the briefing/addendum, this is one of the few units
where the six (twelve) points can actually constrain an exponent, and I use all 12 points
(6 clusters × 2 modalities) per rule as true replicates on identical `n_frag` input (confirmed
identical across `size_manifest_multiome.tsv`/`size_manifest_scatac.tsv`).

Pooled log-log OLS fit (wall `s` vs. `n_frag`, 12 points each):

| Rule | fitted exponent | fitted coefficient | `n_frag` range ratio | wall-time range ratio |
|---|---:|---:|---:|---:|
| `process_fragment_file` | **0.969** | 4.20×10⁻⁶ | 15.05× | 15.21× |
| `frag_to_tagAlign` | **1.016** | 3.49×10⁻⁶ | 15.05× | 18.93× (12.9× if using cluster-mean-of-2; see unit 9 §5 for the noise discussion) |
| `frag_to_norm_bigWig` | **0.883** | 1.64×10⁻⁴ | 15.05× | 16.28× |

Full 12-point tables (wall `s`, `cpu_time` s, both modalities):

**`frag_to_tagAlign`** (`s`): thp1_1 120.68/119.55, thp1_2 120.82/120.91,
jurkat_pma_cd3_4hr 257.87/251.49, jurkat 338.63/337.01, k562_crispri 488.77/387.43,
telohaec_crispri 1552.20/2262.51 (multiome/scatac). `cpu_time/wall` climbs from **0.61
(thp1_1)** to **0.92 (telohaec_crispri, multiome)** with cluster size — see §6a/unit 9 §5-6 for
why (fixed per-job overhead dilutes proportionally at larger `n`).

**`process_fragment_file`** (`s`): thp1_1 65.50/64.83, thp1_2 66.05/64.64,
jurkat_pma_cd3_4hr 116.21/139.86, jurkat 189.84/188.61, k562_crispri 214.87/209.35,
telohaec_crispri 808.04/982.98. `cpu_time/wall` climbs from **0.68 (thp1_1)** to **0.94
(telohaec_crispri, multiome)**.

**`frag_to_norm_bigWig`** (`s`): thp1_1 633.56/571.57, thp1_2 443.06/581.66,
jurkat_pma_cd3_4hr 1002.22/1316.26, jurkat 1632.98/1639.29, k562_crispri 1897.51/1726.79,
telohaec_crispri 4886.69/7213.84. `cpu_time/wall` is already high even at the smallest cluster
(**0.93 at thp1_1**, vs. 0.61-0.68 for the other two rules) and reaches **0.99 (telohaec_crispri
either modality)**.

**Agreement**: all three fitted exponents (0.883–1.016) are consistent with the code's
`O(n_frag)`-to-`O(n_frag log n_frag)` prediction at this `n_frag` range — per unit 9's point
(§10, reused here), `ln(n_frag)` itself only spans 1.16× over this 15.1× range in `n_frag`, so a
pure power-law fit cannot cleanly distinguish `O(n)` from `O(n log n)`; an exponent near 1.0 is
expected either way at this scale. `frag_to_norm_bigWig`'s exponent (0.883, mildly sub-linear)
is the one visible disagreement — I attribute this the same way unit 9 attributes
`frag_to_tagAlign`'s sub-linear-looking growth: a per-job fixed cost (subprocess startup,
`bedGraphToBigWig`'s own fixed indexing overhead, Slurm scheduling latency baked into wall `s`)
that is a larger fraction of the ~450-1600 s jobs at small `n_frag` than of the ~5000-7000 s job
at `telohaec_crispri`, which drags the fitted exponent below 1 without implying the underlying
per-record work is actually sub-linear.

**Cross-modality noise**: because these three rules carry no modality-conditional code (§9),
the multiome/scATAC pair at each cluster is a true replicate on identical input. Four of six
clusters agree to within ~5-10%; `telohaec_crispri` and (for `frag_to_norm_bigWig`) `thp1_2`/
`jurkat_pma_cd3_4hr` disagree by 30-48% between the two runs (e.g. `frag_to_norm_bigWig`
telohaec_crispri: 4886.69 s vs. 7213.84 s). Per unit 9's analysis of the same pattern in
`frag_to_tagAlign`, this looks like shared-Lustre/scheduler contention on the largest, most
I/O-heavy jobs, not a modality effect — there is no `if modality` branch in this file for any
of these three rules to produce a real difference.

## 6. Concrete overhead sources — THE declared-vs-realised threading gap

| Rule | declared `threads:` | measured cpu_time/wall (telohaec_crispri, largest cluster) |
|---|---:|---:|
| `frag_to_norm_bigWig` | 16 | 0.988 (multiome) / 0.990 (scatac) |
| `frag_to_tagAlign` | 8 | 0.919 (multiome) / 0.933 (scatac) |
| `process_fragment_file` | 8 | 0.944 (multiome) / 0.980 (scatac) |

These three numbers are the ones already flagged as established; the mechanism per rule:

**a) `process_fragment_file` — `threads: 8` is wholly inert.** No tool in the shell block
(`awk`, `bgzip`, `tabix`, `zcat`, `wc`, at `:104-107`) is ever passed a thread count; the literal
substring `{threads}` does not occur anywhere in this rule's shell. The gap here isn't "a tool
ignores its thread flag" — no thread flag is ever offered to any tool. On top of that, step (3)
of the shell (`:107`, `zcat {output.fragments_filtered} | wc -l > {output.fragment_count}`) is a
second, fully redundant full decompression of the file step (1) just wrote, purely to produce a
row count — real added wall time, still single-core, with zero parallel opportunity.

**b) `frag_to_tagAlign` — one thread-aware tool, exercised on the smallest possible chunks.**
`-p {threads}` (`:55`) reaches `bedSplitSort.sh`'s `sort_chr()` (`bedSplitSort.sh:51-67`), the
only call site in this rule that passes `--parallel` to anything. But (i) that call happens
once **per chromosome, sequentially** (`bedSplitSort.sh:87-91`, a `while read` loop with no
`&`/`wait`/`xargs -P`/GNU `parallel`), so at most one chromosome's sort is ever running; and
(ii) each chromosome's file is already a ~1/25th-scale partition of the full stream, small
enough that GNU sort's internal parallel-chunk-then-merge algorithm has little external-merge
work left to parallelize. `zcat` (:49), the tagAlign-transform `awk` (:50), the split `awk`
inside `bedSplitSort.sh` (`:71-75`, a hard barrier — must consume 100% of stdin before any sort
can start), `bgzip -c` (`:58`, no `-@`) and `tabix -p bed` (`:61`, no threading option at all)
never see a thread flag. Six/seven single-threaded process invocations behind an 8-thread
reservation; one of them gets `--parallel`, and only in a form (sequential, small-chunk) that
can't realize it.

**c) `frag_to_norm_bigWig` — the thread flag is correctly wired, but the pipe's rate-limiting
stage categorically cannot use it.** `sort --parallel=16` (`:165`) is genuinely applied to a
whole-genome, multi-GB bedGraph stream (10.8 GB uncompressed at `telohaec_crispri` — the largest
single sort job of any of these three rules, unlike `frag_to_tagAlign`'s chromosome-partitioned
sorts). But it sits immediately downstream of `bedtools genomecov -bg` (`:164`), which has *no*
thread flag to ignore — `bedtools` (any version, including the `2.29.2` pin) has never exposed
multithreading for `genomecov`; it is a single serial streaming pass by construction (each
output record depends on a running per-base depth accumulator carried from the previous
record). `zcat` (:163) feeding it is also single-threaded. Downstream, `bedGraphToBigWig`
(:166) is a separate sequential single-core pass with no threading option, run only after the
sort's redirection to `{output.bedGraph_file}` is fully written — it cannot overlap with
anything. The empirical ceiling of ~0.99 cpu/wall even at the largest cluster is consistent with
`genomecov`'s single-core streaming rate capping what the pipe as a whole can deliver to `sort`,
regardless of how many threads `sort` is offered; the code cannot show whether `sort` itself
ever had a large enough backlog to use more than ~1 thread's worth of work, only that in
aggregate the *pipeline* never did (§10).

**d) Where to attribute the sort seconds (for Phase 2 — do not double-count with unit 9).**
Snakemake's `benchmark:` directive times an entire rule invocation (the whole `shell:` block as
one subprocess tree), not individual commands inside it. There is no separate DAG node, and no
separate benchmark row, for `bedSplitSort.sh` or for the bare `sort` call inside
`frag_to_norm_bigWig` — both are subprocesses that execute and exit entirely within the wall
time already reported for `frag_to_tagAlign` and `frag_to_norm_bigWig` respectively in
`benchmarks_{multiome,scatac}.tsv`. **All sort-related seconds for both rules are already fully
contained inside those two rules' measured `s`/`cpu_time` columns.** Unit 9 characterises
`bedSplitSort.sh`'s *algorithm* (why split-then-sequential-chromosome-sort); Phase 2 should not
create a third "sort" or "bedSplitSort" cost bucket alongside `frag_to_tagAlign`'s — that would
double-count time that is already inside the `frag_to_tagAlign` total.

**e) Resource-declaration finding: `runtime_hr` is a dead resource key.**
`frag_to_norm_bigWig` and `frag_to_bigWig` declare `resources: ... runtime_hr=24 ...`
(`frag_to_tagAlign.smk:151`, `:123`), whereas `frag_to_tagAlign` and `process_fragment_file`
declare `resources: ... runtime=720*2 ...` (`:37`, `:92`) — note the different key name.
Snakemake's Slurm executor plugin recognizes the resource named exactly `runtime` (minutes) to
set the Slurm job's `--time`; `runtime_hr` is not a resource name the plugin looks for, so it is
silently ignored for both `frag_to_norm_bigWig` and `frag_to_bigWig`. In practice this means
those two rules' actual Slurm time limit comes from whatever default/profile value applies when
no `runtime` resource is set, not from the `24` (hours, evidently intended) the rule author
wrote. This is a resource-declaration bug independent of the threading question, flagged here
because both affected rules are in this unit.

**f) I/O volume** (`io_in`/`io_out`, MB, `telohaec_crispri`, multiome):

| Rule | io_in | io_out |
|---|---:|---:|
| `process_fragment_file` | 430.78 | 4078.14 |
| `frag_to_tagAlign` | 24306.38 | 27388.53 |
| `frag_to_norm_bigWig` | 494.05 | 11860.66 |

`frag_to_norm_bigWig`'s `io_out` (11.86 GB) is consistent with writing both the 10.8 GB
`ATAC_norm.bg` intermediate and the 1.67 GB `ATAC_norm.bw` final output — confirmed by direct
`ls -la` on `/scratch/.../igvf10_multiome/telohaec_crispri/{ATAC_norm.bg,ATAC_norm.bw}`. That
bedGraph intermediate is declared `temp()` (`frag_to_tagAlign.smk:146`) and only survives on
disk here because Phase 0 ran with `--notemp`; in a normal run it would be deleted immediately
after the rule, but the write+read cost of producing and then re-reading 10.8 GB of uncompressed
text (once by `sort`'s output redirect, once by `bedGraphToBigWig`) is real wall time either
way. `process_fragment_file`'s low `io_in` (430.78 MB) looks implausibly small for decompressing
a multi-GB fragment file plus re-reading the filtered output for the redundant `wc -l` pass
(§6a) — I flag this as a likely measurement-capture artifact (Snakemake's benchmark `io_in` may
undercount reads performed inside short-lived child processes spawned via process substitution,
or reads served from page cache), not a claim that the redundant pass is free.

## 7. `t_est(n)`, `t_low(n)`, `t_high(n)` — one triple per rule

All formulas in seconds, `n_frag` as the raw count (e.g. `26231026`, not millions). Fit to all
12 benchmark points per rule (6 clusters × 2 modalities). Per the task brief, `t_low`/`t_high`
bound the **complexity class**, not a point estimate: `t_low` = pure `O(n_frag)` floor (smallest
observed per-record ratio across all 12 points, i.e. a rate no job undercut); `t_high` = `O(n
log n)` ceiling (largest observed per-`n log n`-unit ratio, i.e. a rate no job exceeded). They
are deliberately different functional classes, not blended.

### `process_fragment_file`

- `t_est(n_frag) = 4.1990×10⁻⁶ · n_frag^0.9692` seconds — best-guess class, `O(n_frag)`
  (least-squares power-law fit; exponent 0.969 is indistinguishable from 1.0 given the fixed
  per-job overhead noted in §5).
- `t_low(n_frag) = 2.0465×10⁻⁶ · n_frag` — `O(n_frag)` floor (smallest observed `s/n_frag`,
  `telohaec_crispri` multiome: 808.04 / 394,847,593).
- `t_high(n_frag) = 1.0192×10⁻⁷ · n_frag · log₂(n_frag)` — `O(n_frag log n_frag)` ceiling,
  included even though this rule's own code has **no sort at all** (§4), purely as a
  conservative upper bound calibrated off the largest observed `s/(n·log₂n)` ratio; there is no
  code-level reason to expect the true class exceeds `O(n_frag)` for this rule.

### `frag_to_tagAlign`

Per §1/§4, this DAG node's cost is dominated by `bedSplitSort.sh`'s `O(n_frag log n_frag)`
comparison sort (unit 9's derivation), and unit 9 already fit that exact class to these same 12
benchmark points (the whole-rule measurement — there is no finer-grained data available, §6d).
Rather than produce a second, possibly-inconsistent fit for the same DAG node from the same 12
numbers, **I adopt unit 9's calibration verbatim** (`phase1_reports/09-bedSplitSort.md` §7),
using `x = n_frag · ln(2·n_frag)` (the `2·` reflecting the tagAlign-doubling at
`frag_to_tagAlign.smk:50`):

- `t_est(n_frag) = 18.33 + 2.34×10⁻⁷ · n_frag · ln(2·n_frag)` — `O(n_frag log n_frag)`,
  least-squares fit to all 12 points; intercept absorbs subprocess/pipeline-startup overhead
  (§6b) that does not scale with `n_frag`.
- `t_low(n_frag) = 3.93×10⁻⁶ · n_frag` — `O(n_frag)` lower bound (smallest observed
  per-fragment ratio, `telohaec_crispri` multiome). This value is identical to what I
  independently computed by the same min-ratio method on the same 12 points — a useful
  cross-check that the two agents' floors agree exactly.
- `t_high(n_frag) = 18.33 + 2.87×10⁻⁷ · n_frag · ln(2·n_frag)` — `O(n_frag log n_frag)` upper
  bound, largest observed per-unit-work ratio (`k562_crispri` multiome).

### `frag_to_norm_bigWig`

- `t_est(n_frag) = 1.6353×10⁻⁴ · n_frag^0.8832` seconds — best-guess class, close to `O(n_frag)`
  with the sub-linear-looking exponent attributed to fixed per-job overhead (§5), not to
  genuinely sub-linear per-record work (the code has no mechanism for that: `genomecov` and
  `bedGraphToBigWig` are both single streaming passes).
- `t_low(n_frag) = 1.2376×10⁻⁵ · n_frag` — `O(n_frag)` floor (smallest observed `s/n_frag`,
  `telohaec_crispri` multiome: 4886.69 / 394,847,593).
- `t_high(n_frag) = 9.8005×10⁻⁷ · n_frag · log₂(n_frag)` — `O(n_frag log n_frag)` ceiling,
  reflecting the whole-genome `sort` stage (§4, §6c), calibrated off the largest observed
  `s/(n·log₂n)` ratio.

## 8. All three evaluated at mean/min/max cluster size

Mean `n_frag` across the 6 clusters = 110,823,580 (`(26,231,026+26,292,230+55,803,917+
75,342,313+86,424,398+394,847,593)/6`); min = `thp1_1` (26,231,026); max = `telohaec_crispri`
(394,847,593).

| Rule | `n_frag` | `t_low` | `t_est` | `t_high` | actual (multiome / scatac) |
|---|---|---:|---:|---:|---:|
| `process_fragment_file` | min (26.2M) | 53.7 s | 65.1 s | 65.9 s | 65.50 / 64.83 |
| | mean (110.8M) | 226.8 s | 263.2 s | 301.8 s | (no single job at mean `n`; interpolated only) |
| | max (394.8M) | 808.0 s | 901.7 s | 1149.2 s | 808.04 / 982.98 |
| `frag_to_tagAlign` | min (26.2M) | 103.1 s | 127.4 s | 152.1 s | 120.68 / 119.55 |
| | mean (110.8M) | 435.5 s | 516.6 s | 629.5 s | — |
| | max (394.8M) | 1551.6 s | 1911.2 s | 2339.9 s | 1552.20 / 2262.51 |
| `frag_to_norm_bigWig` | min (26.2M) | 324.6 s | 583.8 s | 633.6 s | 633.56 / 571.57 |
| | mean (110.8M) | 1371.6 s | 2084.6 s | 2902.5 s | — |
| | max (394.8M) | 4886.7 s | 6403.1 s | 11050.6 s | 4886.69 / 7213.84 |

(`frag_to_tagAlign` row uses unit 9's calibration, §7; the "actual" columns are the two real
replicate measurements at that cluster, not a fit.) At the max cluster, `frag_to_norm_bigWig`'s
scATAC observation (7213.84 s) sits inside `t_high` but well above `t_est` — another instance of
the largest cluster showing more run-to-run spread than the fitted class alone predicts (§5).

## 9. Modality: class or constant?

**Only the constant could conceivably differ, and the data show no clean, consistent
difference.** None of the three shell blocks contains any modality-conditional logic — no
reference to `to_generate`, `config["fragments_preprocessed"]` beyond the parse-time rule
selection in §2a (which is identical for both modalities: both `config.yaml` and
`config_training.yaml` set it to `False` regardless of modality), or any other
modality-determining value. `n_frag` is verified identical between
`size_manifest_multiome.tsv` and `size_manifest_scatac.tsv` for all six clusters. The
per-cluster multiome-vs-scATAC differences reported in §5 (mostly ≤10%, up to 48% at the
largest/most I/O-heavy clusters) match unit 9's characterisation of the same phenomenon in
`frag_to_tagAlign`: shared-filesystem/scheduler contention on a genuinely identical computation,
not a modality effect.

## 10. What the measurement cannot tell us

- **Cannot separate `O(n_frag)` from `O(n_frag log n_frag)` empirically for any of the three
  rules.** `n_frag` spans 15.1×, but `ln(n_frag)` spans only ~1.16× over that range (17.08 to
  19.79) — the same point unit 9 makes for `frag_to_tagAlign` applies identically to
  `frag_to_norm_bigWig`'s whole-genome sort. All three rules' fitted exponents (0.88-1.02) are
  compatible with either class at this scale; the class in §4 comes from reading the code, not
  from these fits.
- **Cannot confirm *why* `sort --parallel=16` in `frag_to_norm_bigWig` fails to raise cpu/wall
  above ~0.99.** §6c gives two code-supported hypotheses (genomecov-limited feed rate; a
  bedGraph stream small/uniform enough that GNU sort's internal parallel merge has little to do)
  but the benchmark data cannot distinguish them — that would need per-process (not per-rule)
  CPU accounting, which Snakemake's `benchmark:` directive does not provide (it measures the
  whole shell-block subprocess tree as one unit).
- **Cannot rule out that `io_in` for `process_fragment_file` is undercounted** (§6f) — the
  reported 430.78 MB at `telohaec_crispri` looks too small for decompressing a multi-GB fragment
  file plus a full redundant re-read; I cannot confirm from this data alone whether that reflects
  a real (surprisingly efficient) read pattern, page-cache hits not counted as I/O, or a gap in
  how Snakemake's benchmark wrapper attributes child-process I/O through a process substitution.
- **Cannot attribute the two largest clusters' 30-48% multiome/scATAC spread to a specific
  cause** beyond "contention on a shared filesystem," consistent with unit 9's finding for the
  same rule — with n=2 replicates per cluster, "this run was unlucky" and "large jobs on this
  stack are inherently more variable" cannot be distinguished.
- **Cannot measure `bedSplitSort.sh` or the bare `sort` call in `frag_to_norm_bigWig` in
  isolation from the rest of their enclosing rule** — no finer-grained benchmark than the whole
  `shell:` block exists (§6d), so all the class/overhead statements above about individual
  pipe stages are read from the code, not independently timed.
- **`frag_to_bigWig` is entirely unmeasured** (0 benchmark rows in either modality, confirmed
  directly) — its shape is identical to `frag_to_norm_bigWig` minus the `-scale` normalization,
  but nothing here characterises its actual runtime.
- **Whether nothing genuinely downstream reads `ATAC_norm.bw`'s contents**: confirmed by
  reading `workflow/Snakefile:100-102,115` and `workflow/rules/save_configs.smk:56-60` — the
  only other reference to `ATAC_norm.bw` is `save_configs.smk` writing its *path string* into an
  IGV session config file, not reading its data, and the file is only requested as a `rule all`
  target when `config["make_IGV_tracks"]` is true. This is a factual dependency observation for
  Phase 2's critical-path computation; whether that changes what should be prioritized is
  explicitly out of scope for this report.

# Unit 9 — `bedSplitSort.sh` (Tier B)

## 1. Unit identity

- Script: `workflow/scripts/bedSplitSort.sh` (91 lines)
- Called from rule `frag_to_tagAlign`, `workflow/rules/frag_to_tagAlign.smk:9-62`, at
  `frag_to_tagAlign.smk:51-58` (inside a `zcat | awk | bash bedSplitSort.sh | bgzip` pipe).
- Modality: **both**. 6 jobs per modality (12 total), one per cluster; `frag_to_tagAlign`'s
  input (`get_processed_fragment_file`, `frag_to_tagAlign.smk:1-6`) is the same fragment file
  regardless of which downstream branch (ARC vs Neither) the cluster resolves to — confirmed
  on disk: `n_frag` is byte-identical between `size_manifest_multiome.tsv` and
  `size_manifest_scatac.tsv` for all six clusters. That makes this rule part of the "shared
  ABC/fragment stack" the briefing flags as having true n=2 replicates per cluster (§4.6).
- Declares `threads: 8` (`frag_to_tagAlign.smk:34`); measured `cpu_time/wall` ranges
  **0.61 → 0.92** across the six clusters (rises with cluster size — see §5/§6), matching the
  headline 0.92 figure, which is essentially the `telohaec_crispri` value.
- Max observed wall time: 37.7 min (`telohaec_crispri`, scATAC run, 2262.51 s) — higher than
  the 25.9 min figure in the brief, which is the *multiome* run of the same cluster (1552.20 s).
  Both are real; see §5 for why they disagree by 46%.

## 2. What the code actually does

The calling rule first turns each fragment into two single-ended tagAlign records
(`frag_to_tagAlign.smk:50`: `print $1,$2,mid,"N",1000,"+"; print $1,mid+1,$3,"N",1000,"-"`), so
`bedSplitSort.sh`'s actual input stream is **2 × n_frag** lines, not n_frag.

`bedSplitSort.sh` then runs in three strictly sequential phases (bash executes each top-level
statement to completion before the next one starts — there is no `&`/`wait`, `xargs -P`, or
GNU `parallel` anywhere in the script):

1. **Split by chromosome** (`bedSplitSort.sh:71-75`). Because the caller passes
   `-i /dev/stdin`, this is a single `awk` process reading the whole 2·n_frag-line stream and
   writing each line to `$SPLITDIR/<chr>.bed` (`bedSplitSort.sh:72`). One pass, one core,
   no way to parallelize a single `awk` invocation. This is a **hard synchronization
   barrier**: nothing in phase 3 can start until this `awk` has consumed 100% of stdin,
   because bash does not proceed past line 72/74 until that command exits.
2. **Chromosome-presence warning** (`bedSplitSort.sh:78-84`). For each of ≤25 split files, a
   fresh `awk` linearly scans the 25-line sizes file for a match. This is `O(25×25)`, a fixed
   ~625-comparison cost independent of `n_frag`. Negligible.
3. **Per-chromosome sort, one chromosome at a time** (`bedSplitSort.sh:87-91`). A `while read`
   loop over the (constant, 25-line: chr1–22, X, Y, M —
   `ENCODE_rE2G/reference/GRCh38_EBV.no_alt.chrom.sizes.tsv`) sizes file, calling `sort_chr()`
   on each chromosome's split file **sequentially**. `sort_chr()` (`bedSplitSort.sh:54-67`)
   invokes `LC_ALL=C sort -k2,2n -k3,3n --parallel "$THREADS" -S "$BUFFER_SIZE" -T
   "$SORT_TMPDIR"` — a comparison sort, only ever on one chromosome's worth of data at a time.

The `sort-bed` (bedops) fast path (`bedSplitSort.sh:50-52`) is dead in production: `bedops` is
not in `workflow/envs/sc_e2g.yml` (checked — only `bedtools`, `tabix`,
`ucsc-bedgraphtobigwig`, no bedops), so `command -v sort-bed` fails and every run takes the
GNU-`sort` branch (`bedSplitSort.sh:54-68`).

Around the script (not in it, but load-bearing for the CPU/wall gap the brief asks about):
`frag_to_tagAlign.smk:49` (`zcat`), `:50` (tagAlign-transform `awk`), `:58` (`bgzip -c`, no
`-@`), `:61` (`tabix -p bed`, no threading option exists for it). All four are single-threaded
tools; `zcat`/transform-`awk` overlap with phase 1's split `awk` via the pipe, but `bgzip` and
`tabix` cannot start their real work until phase 3 has produced (and, for `tabix`, until
`bgzip` has finished writing) output.

## 3. Size variable

`n_frag` — ATAC fragments after chromosome filtering (`fragment_count.txt`), 26,231,026 –
394,847,593, ratio **15.1×** (Tier 1, genuinely variable — from `size_manifest_multiome.tsv`
column 6; scATAC's manifest carries identical values, verified). Per §2, the script's actual
working set is `2·n_frag` tagAlign records; the factor of 2 is a constant and does not change
the class.

Per the addendum, `n_frag_bytes` is a measurement artifact (digit count of
`fragment_count.txt`) and is not used anywhere below.

## 4. Asymptotic class derived from code

| Phase | file:line | Class in `n_frag` | Threaded? |
|---|---|---|---|
| Split by chromosome | `bedSplitSort.sh:71-75` (writing) via `frag_to_tagAlign.smk:49-50` (reading) | `O(n_frag)` — one linear pass over `2·n_frag` records | No — single `awk` process, always |
| Chromosome-presence check | `bedSplitSort.sh:78-84` | `O(1)` (≤25×25) | N/A |
| Per-chromosome sort loop | `bedSplitSort.sh:87-91` calling `:54-67` | `O(Σᵢ mᵢ log mᵢ)`, `Σᵢ mᵢ = 2·n_frag`, 25 sequential invocations | Only the individual `sort --parallel` call; loop itself is strictly sequential |

`sort` is a comparison sort (GNU coreutils 8.22, confirmed via `sort --version` on this
cluster); each per-chromosome call costs `O(mᵢ log mᵢ)`. Human chromosome-length shares are
skewed but fixed (chr1 ≈ 8% of the 25-chromosome total, per
`ENCODE_rE2G/reference/GRCh38_EBV.no_alt.chrom.sizes.tsv`), so `mᵢ` for the largest chromosome
is a fixed fraction of `2·n_frag` and `log(mᵢ) = log(n_frag) - O(1)` for every chromosome.
That makes `Σᵢ mᵢ log mᵢ = Θ(n_frag · log n_frag)` — the same class as sorting the whole
`2·n_frag`-record stream in one shot, just executed as 25 back-to-back smaller sorts instead
of one big one. The `-p threads`/`--parallel` flag only threads *within* one of those 25 calls
(GNU sort's internal parallel-chunk-then-merge algorithm); it does not run chromosomes
concurrently with each other.

**Overall class: `O(n_frag) + O(n_frag · log n_frag) = O(n_frag · log n_frag)`**, with the
split phase's `O(n_frag)` term as a guaranteed-single-core, guaranteed-disk-round-trip
additive cost (see §6).

## 5. Benchmark cross-check

`n_frag` and wall `s` for `frag_to_tagAlign` (whole rule; no separate benchmark exists for
`bedSplitSort.sh` alone — see §10), both modalities (true n=2 replicates on identical input,
per briefing §4.6):

| cluster | n_frag | multiome s | scATAC s | cpu_time/wall (multiome) |
|---|---:|---:|---:|---:|
| thp1_1 | 26,231,026 | 120.68 | 119.55 | 0.608 |
| thp1_2 | 26,292,230 | 120.82 | 120.91 | 0.622 |
| jurkat_pma_cd3_4hr | 55,803,917 | 257.87 | 251.49 | 0.705 |
| jurkat | 75,342,313 | 338.63 | 337.01 | 0.713 |
| k562_crispri | 86,424,398 | 488.77 | 387.43 | 0.755 |
| telohaec_crispri | 394,847,593 | 1552.20 | 2262.51 | 0.919 |

Observations:

- **`n_frag` spans 15.1×, wall time spans only 12.9×** (120.1 s mean-of-two → 1907.4 s
  mean-of-two, using the two smallest and largest cluster means). That is *sub-linear*,
  disagreeing with both `O(n_frag)` (predicts ≥15.1×) and `O(n_frag log n_frag)` (predicts
  ~17.4×, since `ln(n)` grows from 17.08 to 19.79, +16%, over the range). **This is a
  code-vs-data disagreement, and I don't smooth it.**
- The likely reason is directly visible in the CPU/wall column: it climbs monotonically from
  0.608 (thp1_1) to 0.919 (telohaec_crispri) with cluster size. `bedSplitSort.sh` plus its
  calling pipe spawns ~5 single-threaded tools per stage (`zcat`, transform-`awk`, split-`awk`,
  25× `sort`, `bgzip`, `tabix`) — subprocess start latency, `mkdir`/`trap` setup, and 25
  separate small-file opens under `$SPLITDIR` are a roughly **fixed** wall-clock cost. At
  ~120 s (thp1_1) that fixed cost is a large fraction of the total; at ~1550-2260 s
  (telohaec_crispri) it is a small fraction. A fixed overhead diluting proportionally faster
  than the size-driven term grows produces exactly this apparent-sub-linear shape — it masks,
  rather than refutes, the `n_frag log n_frag` compute term underneath.
- **`telohaec_crispri` disagrees with itself by 46% between the two modality runs**
  (1552.20 s vs 2262.51 s) on byte-identical input, while the four smaller clusters agree to
  within ~3% between runs. `k562_crispri` also disagrees by 21% (488.77 vs 387.43 s), in the
  opposite direction. Per briefing §4.6 these two runs are true replicates of the same
  computation, so this spread is measurement noise, not signal — and it is concentrated in the
  two largest, most I/O-heavy clusters. That is consistent with the split-then-sort design
  round-tripping tens of GB through `$SCRATCH` (Lustre, shared, contended) rather than through
  node-local storage (see §6): the bigger the job, the more its wall time is exposed to
  transient filesystem/scheduler contention on a shared cluster, and the *less* the six
  (twelve) points can be trusted as a clean scaling curve.
- **`-T`/`-S` verified, not assumed**: `frag_to_tagAlign.smk:47` creates
  `SPLIT_DIR=$(mktemp -d -p {resources.temp_dir} ...)` and passes it to `bedSplitSort.sh -t`;
  `resources.temp_dir = os.path.join(RESULTS_DIR, "tmp")` (`frag_to_tagAlign.smk:38-41`), and
  `RESULTS_DIR` for this run is
  `results_dir: /scratch/users/kaybrand/scE2G_optimize_results/igvf10_multiome`
  (`configs/igvf10_multiome_config.yaml:57`) — confirmed to exist on disk at that exact path.
  `/scratch` is mounted as Lustre (`mount` output: `...lustre (rw,...)`), **not**
  `$L_SCRATCH`/node-local (`/lscratch`, `xfs`, also verified via `mount`). So both the
  chromosome-split intermediate files (`bedSplitSort.sh -t "$SPLIT_DIR"`) and GNU `sort`'s own
  external-merge spill directory (`bedSplitSort.sh:64`, `-T "$SORT_TMPDIR"`) sit on a shared
  network filesystem, not local SSD. `-S "$BUFFER_SIZE"` is generous (derived from
  `determine_mem_mb`, `ENCODE_rE2G/ABC/workflow/rules/utils.smk:17-24`: roughly
  `4 × compressed_input_size_MB × 8` capped at 250 GB, divided by `threads×2` — e.g. ≈576 MB
  for thp1_1, ≈8 GB for telohaec_crispri), which likely keeps most/all individual
  chromosomes' sorts in-memory (no actual spill to `-T`) at *this* data range — I flag this as
  inference, not a confirmed fact, since GNU sort's own temp files are created and cleaned up
  within the same process and I did not catch one mid-run. If a future dataset makes even one
  chromosome's share exceed its buffer, the constant behind the `n_frag log n_frag` term would
  jump sharply, because every spill/read pass would cross Lustre instead of local disk.

## 6. Concrete overhead sources

- **Split-then-sort forces every record through disk twice, on Lustre.** The chromosome split
  (`bedSplitSort.sh:71-75`) writes all `2·n_frag` records to `$SPLITDIR/*.bed` files under
  `$SCRATCH`; the sort loop then reads them back. This is directly visible in the I/O columns:
  `telohaec_crispri` shows `io_in≈24.3 GB`, `io_out≈27.4 GB` — roughly double the ~4 GB
  on-disk compressed fragment file, consistent with writing the ~31.6 GB uncompressed
  `2·n_frag`-record intermediate to `$SPLITDIR` and reading it back, on top of the final
  compressed output. A design that piped directly into 25 concurrent sort processes (e.g.
  named pipes or per-chromosome FIFOs) would avoid this write/read round trip entirely.
- **A hard serialization barrier between phase 1 and phase 3** (§2): the split `awk` must
  finish consuming the *entire* input before the first `sort_chr` call can start, so the
  pipeline's inherent concurrency (`zcat`‖transform-`awk`‖split-`awk`‖…‖`bgzip`) only actually
  overlaps for as long as phase 1 is running; `bgzip`/`tabix` are idle until phase 3 produces
  output, and `zcat`/transform-`awk` are already done by the time phase 3 starts.
- **Declared 8 threads, realized ≈1**: `--parallel "$THREADS"` (`bedSplitSort.sh:57-58`) only
  ever touches GNU `sort`'s internal chunk-sort phase, once per chromosome, inside a strictly
  sequential 25-iteration loop (`bedSplitSort.sh:87-91`). `zcat`, both `awk` invocations,
  `bgzip -c` (no `-@`, `frag_to_tagAlign.smk:58`) and `tabix -p bed` (no threading option,
  `frag_to_tagAlign.smk:61`) never look at `threads:` at all. The measured cpu_time/wall
  (0.61–0.92, §5) is consistent with "≈1 core continuously busy, topped up by brief
  multi-core bursts from `--parallel` only on the 2-3 largest chromosomes (chr1, chr2, chrX)
  once their share of `2·n_frag` is large enough for GNU sort's internal parallel merge to
  matter" — nowhere close to the 8-core parallelism the rule requests from Slurm.
- **`sort-bed`'s ~2× speedup (script's own comment, `bedSplitSort.sh:8`) is unrealized**:
  `bedops` is absent from `workflow/envs/sc_e2g.yml`, so `command -v sort-bed`
  (`bedSplitSort.sh:50`) always fails and the GNU-`sort` branch always runs.
- **Fixed per-job overhead is a larger overhead source than the algorithm at small `n_frag`,
  and shrinks proportionally at large `n_frag`** — this is the mechanism behind §5's
  sub-linear-looking wall-time growth: ~5 single-threaded tool invocations plus 25 sequential
  `sort` process spawns plus `mkdir`/`trap`/`mktemp` bookkeeping cost roughly the same
  wall-clock regardless of `n_frag`, and dominate more of a 120 s job than a 1550-2260 s one.
- **Memory (`-S`) is sized off the fragment file's *compressed* size** (`determine_mem_mb`,
  `ENCODE_rE2G/ABC/workflow/rules/utils.smk:17-24`, ×8 for `.gz`, ×4 again, ÷16 for
  `threads×2`) rather than the actual `2·n_frag`-record stream size — a heuristic, not a
  measurement of what `sort_chr` will actually need per chromosome.

## 7. `t_est(n_frag)`, `t_low(n_frag)`, `t_high(n_frag)`

All formulas in seconds, `n_frag` as the raw fragment count (e.g. 26231026, not millions).
Calibrated on all 12 `frag_to_tagAlign` job measurements (6 clusters × 2 modalities — true
replicates on identical input, per §5/briefing §4.6), using `x = n_frag · ln(2·n_frag)` as the
regression variable for the `n log n` forms (the `2·` reflects the tagAlign-doubling at
`frag_to_tagAlign.smk:50`, per §2-3).

- **`t_low(n_frag) = 3.93×10⁻⁶ · n_frag`** — `O(n_frag)` lower bound on the class: the
  scenario where the `log(n_frag)` factor is swamped by per-record I/O/parse constants at this
  scale (consistent with the sub-linear-looking fit in §5). Constant = the *smallest* observed
  per-fragment ratio across all 12 jobs (`telohaec_crispri`, multiome run: 1552.20 s /
  394,847,593 = 3.9312×10⁻⁶ s/frag).

- **`t_high(n_frag) = 18.33 + 2.87×10⁻⁷ · n_frag · ln(2·n_frag)`** — `O(n_frag log n_frag)`
  upper bound: same class as `t_est`, but with the constant set to the *largest* observed
  per-unit-work ratio across all 12 jobs (`k562_crispri`, multiome run:
  `(488.77 − 18.33) / (86424398 · ln(172848796)) = 2.87×10⁻⁷`), which upper-bounds every one
  of the 12 observed points, including the contended `telohaec_crispri` scATAC outlier
  (2262.51 s observed vs. 2339.9 s predicted).

- **`t_est(n_frag) = 18.33 + 2.34×10⁻⁷ · n_frag · ln(2·n_frag)`** — best-guess class:
  `O(n_frag log n_frag)`, the class actually implemented by the sequential per-chromosome
  comparison-sort loop (§4), least-squares fit to all 12 points. Fixed intercept (18.33 s)
  absorbs the subprocess/pipeline-startup overhead described in §6, which does not scale with
  `n_frag`.

These are **not blended**. `t_low` is a different functional class (`O(n_frag)`) from
`t_est`/`t_high` (`O(n_frag log n_frag)`), each anchored at a real extreme observation, exactly
as instructed — averaging an `O(n)` guess with an `O(n log n)` guess would be wrong at every
scale.

## 8. All three evaluated

| n_frag | cluster | `t_low` | `t_est` | `t_high` | actual (multiome / scATAC) |
|---|---|---:|---:|---:|---:|
| 110,823,580 (mean of 6) | — | 435.5 s (7.3 min) | 516.6 s (8.6 min) | 629.5 s (10.5 min) | mean-of-cluster-means 543.9 s |
| 26,231,026 (min) | thp1_1 | 103.1 s (1.7 min) | 127.4 s (2.1 min) | 152.1 s (2.5 min) | 120.68 s / 119.55 s |
| 394,847,593 (max) | telohaec_crispri | 1551.6 s (25.9 min) | 1911.2 s (31.9 min) | 2339.9 s (39.0 min) | 1552.20 s / 2262.51 s |

At the max, the two real observations (1552.20 s, 2262.51 s) land almost exactly on `t_low`
and inside `t_high` respectively — i.e. the observed spread between the two modality runs of
the largest cluster is itself comparable to the full `t_low`–`t_high` bracket. That is the
clearest evidence in this dataset that run-to-run noise at scale (§5) is at least as large as
the uncertainty in the complexity class itself.

## 9. Modality: class or constant?

**Only the constant could conceivably differ, and even that shows no clear systematic
difference.** `bedSplitSort.sh` contains no modality-conditional logic at all — it operates
purely on the fragment stream handed to it by `frag_to_tagAlign.smk`, which is upstream of the
per-cluster ARC/Neither branch decided at `checkpoint features_required`
(briefing §2). `n_frag` is identical between `size_manifest_multiome.tsv` and
`size_manifest_scatac.tsv` for all six clusters (verified by direct comparison, §3). The
per-cluster multiome-vs-scATAC wall-time differences reported in §5 (±3% for four clusters,
+46%/−21% for the two largest) look like shared-filesystem/scheduler noise on a genuinely
identical computation, not a modality effect — there is no code path in this script (or in
`frag_to_tagAlign.smk`'s shell block) that reads `config["fragments_preprocessed"]`,
`to_generate`, or any other modality-determining value.

## 10. What the measurement cannot tell us

- **Cannot separate `O(n_frag)` from `O(n_frag log n_frag)` empirically at this range.**
  `n_frag` spans 15.1×, but `ln(n_frag)` — the only thing that would distinguish the two
  classes — spans just **1.16×** (17.08 → 19.79) over that same 15.1× range. A genuinely
  "Tier 1" variable's own logarithm behaves like a near-constant (Tier 3-style) quantity here,
  which is a subtler version of the briefing's core warning: even a 15× range in `n` is not
  enough to detect a `log n` factor, because `log` compresses range so aggressively. The class
  in §4 is derived from the code (a real, unambiguous comparison sort), not confirmed by the
  six/twelve data points.
- **Cannot attribute the sub-linear-looking wall-time growth to a single cause.** §5/§6
  proposes "fixed per-job overhead is a larger share of a small job" as the mechanism, and the
  CPU/wall trend (0.61→0.92) is consistent with it, but I have not instrumented the pipeline
  stage-by-stage (e.g. timing the split `awk` separately from the sort loop) to confirm the
  split contributes the fraction I attribute to it versus, say, `tabix`'s indexing pass.
- **Cannot confirm or rule out actual disk spilling in GNU `sort`'s external merge.** §5 notes
  the `-S` buffer is generous enough that most/all chromosomes' data may fit in memory at this
  data range, which would mean the Lustre-vs-local-disk distinction for `-T` matters little
  *today* — but this is inferred from the `determine_mem_mb` formula and file sizes, not from
  observing an actual spill. If per-chromosome data ever exceeds the buffer (larger future
  datasets, or a lower `max_memory_allocation_mb`), the Lustre location of `-T` would matter a
  lot more, and this benchmark set cannot tell us where that threshold is.
- **Cannot isolate `bedSplitSort.sh`'s own wall time from the rest of `frag_to_tagAlign`.**
  There is no separate benchmark for the script — only for the whole rule (`zcat` +
  tagAlign-transform `awk` + `bedSplitSort.sh` + `bgzip` + `tabix`). I have argued (§6) that
  the script's split+sort work is the dominant, asymptotically-important term (the other
  stages are single-pass `O(n_frag)` streaming operations with a smaller per-record constant
  than a disk-round-tripping comparison sort), but I cannot *measure* the script in isolation
  from this data.
- **The two largest clusters' multiome/scATAC disagreement (46%, 21%) is unexplained beyond
  "shared filesystem/scheduler contention."** I cannot distinguish "this specific measurement
  was unlucky" from "large jobs on this rule are inherently more variable" with n=2 replicates
  per cluster — that would need more runs, which is out of scope for characterisation.
- **`sort-bed`'s counterfactual speedup is not measured at all** — the script's own comment
  claims ~2×, but since `bedops` is absent from the conda env (§2, §6) that path has never
  executed in any of these 12 jobs, so there is no data point for it either way.

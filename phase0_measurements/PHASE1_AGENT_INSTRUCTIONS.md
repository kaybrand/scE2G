# Phase 1 agent instructions (shared by all 26 units)

You are one of 26 Phase 1 characterisation agents on the scE2G critical-path project.

**Repo root — `cd` here:**
`/oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G_optimize` (branch `optimize-sce2g`)

**Read next, in this order, before opening any code:**
1. `phase0_measurements/PHASE1_BRIEFING.md` — where the data is, what the size variables mean,
   and what the six clusters can and cannot establish
2. `phase0_measurements/PHASE1_ADDENDUM.md` — launch-time corrections: real rule names, three
   `size_manifest` columns that are measurement artifacts, and the ways scATAC's `n` differs
   from multiome's

Your specific unit is in the prompt that spawned you.

---

## The rule that matters most

**Derive the complexity class by READING THE CODE.** The six cluster measurements are a
secondary check only.

Across those six clusters only `n_frag` (15.1×), `n_umi` (12.0×) and `n_cells` (6.0×) vary
meaningfully — and those three are mutually correlated, so even they cannot be separated from
six points. Everything downstream of candidate-region calling spans only **1.01×–1.17×**.

If your rule is driven by a near-constant variable, the correct report says: *"class X, derived
from code at `file:line`; the benchmarks cannot confirm it because the variable spans only N×."*
**Do not fit an exponent to noise.** "The data cannot constrain this" is a correct and expected
answer for many rules. An unsupported confident exponent is a failure of this task.

**Near-constant does not mean unimportant.** O(j²) over ~157,000 candidate elements is ~2.5×10¹⁰
operations versus ~2.8×10⁶ for O(j log j) — four orders of magnitude. Report the class *in* `j`,
and separately note it is not empirically confirmable at this range. Those are two different
claims and both belong in your report.

### The code is the authority on complexity. The benchmarks are not.

**O(n log n) versus O(n²) is decided by reading the code — by finding the sort, the hash join, the
nested loop, the interval tree — and by nothing else.** The benchmark table exists to give you
background on *overhead*: startup cost, I/O volume, memory behaviour, realised versus declared
threading, and the numeric constant in front of your cost function. It is **not** evidence about
the exponent, in either direction.

This cuts both ways, and the second direction is the one that gets missed:

- Do not use flat or noisy timings to **fit** an exponent. (Overreach.)
- Do not use flat or noisy timings to **argue down** an exponent you found in the code. If the code
  contains a nested loop over `j_cand_elems`, the class is O(j²) — full stop — even if the six
  measured points look linear, even if they look flat, and even if a regression prefers a straight
  line. Report the code-derived class, then report the disagreement as a *disagreement*, and
  hypothesise a mechanism (a dominant constant at this scale, caching, I/O dominance, the loop
  being over a bounded inner set rather than all of `j`).

Concretely, when you write `t_low` and `t_high`: those bound the **complexity class you derived from
the code**. They are not the envelope of your scatter plot. Choosing a linear `t_low` merely because
the measurements look flat, when the code says `n log n`, is the same error as fitting an exponent
to noise — it just points the other way. If the code is unambiguous about the class, then `t_low`
and `t_high` should differ in their **constant**, not in their exponent. Reserve differing exponents
for cases where the *code itself* leaves the class genuinely ambiguous (for example, a pandas
`merge` whose cost depends on key uniqueness you cannot determine statically), and say which
ambiguity in the code drives the spread.

One consequence worth knowing, because it makes the point unavoidable: over the full measured range
of the most variable quantity in the dataset, `n_frag` = 26.2M → 394.8M (15.1×), the factor
`ln(n_frag)` itself spans only **1.16×**. So the log term is nearly constant across every
measurement you have. No dataset of this shape can separate O(n) from O(n log n) empirically. That
separation has to come from the code, and it is available there for free.

## Scope

Characterisation **only**. No code edits. No fix proposals. No prioritisation. No speculation
about what should be optimized first. Phase 2 computes the critical path and decides what
matters; a plausible-sounding recommendation made without the path arithmetic is the exact error
this whole project exists to avoid.

The only file you may write is your report.

## Already established — do not re-derive (but do explain the mechanism if it is your rule)

- **Declared `threads:` are largely fictional.** `frag_to_norm_bigWig` declares 16 and measures
  CPU/wall 0.99; `frag_to_tagAlign` and `process_fragment_file` declare 8 and measure ~0.92–0.94.
  Every ABC rule declares no `threads:` and therefore runs at 1. `compute_kendall` requests
  `config["threads"]`, currently **1**, passed through to `compute_kendall.py --threads`.
  If this is your rule, explain *why* the gap exists (pipe serialisation, a tool ignoring its
  thread flag, a `sort --parallel` that is not the bottleneck) — don't just re-measure it.
- **`compute_kendall` is no longer the bottleneck** (14.5 min max; already optimized upstream).
  Its declared 63 GB memory floor and 12 h runtime allowance are stale.
- **`create_neighborhoods` now leads at 112 min and is CPU-bound on ONE core** — 290 MB RSS and
  ~1 GB read look I/O-shaped, but cpu_time/wall = 0.977 for two hours. It is not I/O-bound.
- **`frag_to_norm_bigWig` is the largest total consumer** (174.9 min summed) but nothing
  downstream consumes its output — only `make_IGV_tracks`. Total time ≠ critical-path time.
- **Row-count identities are collapsed deliberately.** `n_cand_pairs` == ARC rows ==
  ActivityOnly rows == `genomewide_features` rows; `n_kendall_pairs` == `Pairs.Kendall` rows;
  `j_cand_elems` == `EnhancerList.txt` rows − 1. Do not "discover" these — they are tautological
  (pure column-additions, zero filtering). **Bytes are not redundant:** identical row counts span
  203–630 MB, a 3.1× I/O spread.
- **"pairs" vs "links" is a real trap.** `n_cand_pairs` ~10.8M candidate pairs → `n_called_all`
  thresholded rows → `n_called_distal` excludes `class == "promoter"`. Use `n_called_all` for
  cost (rules read the whole file); `n_called_distal` only when you mean biological links.
- **Exclude the toy chr22 run entirely** — 49 KB of fragments measures interpreter and conda
  startup, not algorithm.
- **Sub-second rules reporting `cpu_time/s` of 3–14 are sampling noise, not parallelism.**
  `hover_plots` reports `cpu_time` 0.00 (not captured) — use wall time.
- `generate_chrom_sizes_bed_file` has **no benchmark file by design** (psutil crashes; the awk
  exits before it can attach). ~0.03 s. Don't hunt for the file.

## Deliverable

Write your report to `phase1_reports/<FILENAME>` (repo-relative; the filename is in your
prompt — create the directory if it does not exist). Also return a condensed version as your
final message. Structure:

1. **Unit identity** — script file(s) with line counts, rule(s), which modalities, job counts.
2. **What the code actually does** — the control flow that costs time, cited as `file:line`.
3. **Size variable(s)**, named *exactly* as they appear in `size_manifest_*.tsv`. Cross terms
   where the code genuinely has them (e.g. `O(n_cand_pairs · k_genes)`). Do not force a single `n`.
4. **Asymptotic class derived from code**, with the specific loops, joins, sorts or merges that
   justify it, each cited `file:line`.
5. **Benchmark cross-check** — does the observed scaling agree? Give the actual numbers. State
   explicitly if the variable's range is too narrow to tell. Flag any code-vs-data disagreement
   and hypothesise why (caching, I/O dominance, a constant that dominates at this scale).
6. **Concrete overhead sources** — I/O volume (use `io_in`/`io_out`), serialisation and
   deserialisation, redundant writes, repeated passes over the same data, declared-vs-realised
   threading, subprocess and interpreter startup.
7. **Three functions: `t_est(n)`, `t_low(n)`, `t_high(n)`** — explicit formulas in your named
   variables, in **seconds**, calibrated to the measurements. `t_low`/`t_high` bound the
   **complexity class itself**, not point estimates.
   **Do not blend them.** PERT's `(O + 4M + P)/6` assumes point-estimate human uncertainty about
   a fixed task; ours is uncertainty about an exponent, and averaging an O(n) guess with an O(n²)
   guess produces a number wrong at every scale. Phase 2 will evaluate these functions
   mechanically, so make them unambiguous and directly evaluable — state every constant
   numerically and name every variable exactly.
8. **All three evaluated** at the mean cluster size across the 6 clusters, and also at min
   (`thp1_1`) and max (`telohaec_crispri`). Reporting convenience only — the functions are the
   deliverable.
9. **If the rule runs in both modalities:** does the *class* differ, or only the constant?
10. **"What the measurement cannot tell us"** — an explicit section. Be specific.

Be rigorous and honest rather than confident. Cite `file:line` throughout. Where you are
guessing, say so and say why.

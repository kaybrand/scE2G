# Phase 1 unit #3 — `create_predictions` (ABC prediction generation)

Tier A. Characterisation only — no code edits, no fix proposals, no prioritisation.

---

## 1. Unit identity

| item | value |
|---|---|
| Rule | `create_predictions`, defined at `ENCODE_rE2G/ABC/workflow/rules/predictions.smk:12` |
| Name in exported DAG | `abc_create_predictions` (ABC is included as a Snakemake module) |
| Scripts | `ENCODE_rE2G/ABC/workflow/scripts/predict.py` (316 lines) — the driver<br>`ENCODE_rE2G/ABC/workflow/scripts/predictor.py` (572 lines) — the per-chromosome kernel<br>`ENCODE_rE2G/ABC/workflow/scripts/getVariantOverlap.py` (71 lines) — see §2.6 |
| Also reads | `ENCODE_rE2G/ABC/workflow/scripts/tools.py` (`df_to_pyranges:124`, `determine_expressed_genes:101`, `run_piped_commands_safe:37`), `hic.py` (`get_powerlaw_at_distance:257`) |
| Modalities | **BOTH.** 6 jobs each = 12 jobs total |
| Inputs | `<cluster>/Neighborhoods/EnhancerList.txt` (33.0–33.3 MB, `j_cand_elems`+1 rows × 25 cols)<br>`<cluster>/Neighborhoods/GeneList.txt` (5.57–5.62 MB, `k_genes`+1 rows × 34 cols) |
| Declared outputs | `Predictions/EnhancerPredictionsAllPutative.tsv.gz` (281.5–335.3 MB, **`n_cand_pairs` rows × 29 cols**)<br>`Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz` (82.1–87.1 MB, 3.77–3.96 M rows × 29 cols) |
| Undeclared outputs | `Predictions/EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.tsv.gz` (6.15–6.24 MB) and its never-deleted `.tmp` (5.83–5.92 MB) — §2.6 |
| `threads:` | **Not declared.** `grep -n threads ENCODE_rE2G/ABC/workflow/rules/*.smk` returns nothing → runs at 1. Slurm `AllocCPUS=1` confirms. |
| `resources:` | `mem_mb = partial(determine_mem_mb, min_gb=20)` (`predictions.smk:35-36`). Inputs are 38.9 MB and uncompressed, so `max(4·38.9, 20000) = 20000 MB` (`utils.smk:17-24`). Measured peak RSS 8.08–8.82 GB → **40–44 % of the reservation.** |
| Wall time | 549.3–634.4 s (multiome), 552.3–822.6 s (scATAC) |

**Powerlaw path only.** `HiC_file` is empty for all six clusters in both
`configs/tables/igvf10_cell_clusters_multiome.tsv` and `..._scatac.tsv` (verified column 4 = `""`
for all 12 rows). `_get_run_predictions_hic_params` (`predictions.smk:3-10`) therefore returns
`--score_column powerlaw.Score` and passes **no** `--hic_file`. The written
`Predictions/parameters.predict.txt` confirms `hic_file None`, `window 5000000`,
`score_column powerlaw.Score`. Everything under `if args.hic_file:` (`predictor.py:29-60`) —
`add_hic_from_hic_file:211`, `add_hic_from_directory:280`, `qc_hic:482`,
`scale_hic_with_powerlaw:436`, `add_hic_pseudocount:465`, `fill_diagonals:109`,
`create_df_from_records:136` — **did not execute in either run.** I characterise the powerlaw
path and say nothing quantitative about the HiC branch (§10).

Two knock-on facts about the dead branch that matter for reading the output:

- `--score_column` is passed **twice** on the command line (`predictions.smk:42` gives
  `params.score_column` = `ABC.Score` from `config.yaml:47`; `predictions.smk:50` expands
  `params.hic_params` which appends `--score_column powerlaw.Score`). argparse keeps the last, so
  the effective value is `powerlaw.Score`. This only affects `test_variant_overlap`'s threshold.
- `params.flags` = `--scale_hic_using_powerlaw` (`config.yaml:41`) is parsed
  (`parameters.predict.txt` line 13 = `True`) but is read only inside
  `scale_hic_with_powerlaw` (`predictor.py:439`), which is HiC-only. It is inert here.

---

## 2. What the code actually does

### 2.1 The cost profile in one table

`predict.py` prints its own per-chromosome timings. Summing them against the benchmark wall time
decomposes the job:

| phase | code | thp1_1 (s) | share |
|---|---|---:|---:|
| conda activate + `import hicstraw/pandas/pyranges` + argparse + `write_params` | `predict.py:1-10, 172` | ~4 | 0.7 % |
| `read_csv` × 2, column subsetting | `predict.py:175-225` | ~2 | 0.4 % |
| **per-chromosome loop, all 23 chromosomes** | `predict.py:249-270` | **17.05** | **3.1 %** |
|  ├ `make_pred_table` (the pyranges/NCLS join) | `predictor.py:78-106` | 9.84 | 1.8 % |
|  └ annotate + powerlaw + 2× `compute_score` | `predictor.py:26,27,62-74` | 7.21 | 1.3 % |
| `pd.concat` + **3 × `to_csv(compression="gzip")`** + mask copies + `bedtools` pipe | `predict.py:274-298` | **~526** | **95.8 %** |
| | **total** | **549.3** | |

The `~4 s` startup is from the mtime of `parameters.predict.txt` (written at `predict.py:172`,
17:54:08) minus the job's own start line in the Slurm log (17:54:04). The 17.05 s loop total is
the sum of the script's own `Completed chromosome: … Elapsed time:` lines. Loop totals across all
12 jobs: 15.1–18.0 s, except scATAC/telohaec_crispri at 32.2 s (that whole job was contended —
§5).

**The algorithm is 3 % of this rule. Serialising ~3.8 GB of CSV text through zlib level 9 on one
core is 96 % of it.**

### 2.2 The enhancer↔gene pairing rule — the thing that sets `n_cand_pairs`

This is the deliverable the prompt flags as most important, so I state it precisely.

`make_pred_table` (`predictor.py:78-106`) does, per chromosome:

```python
enh["enh_midpoint"] = (enh["start"] + enh["end"]) / 2      # predictor.py:81
enh_pr   = df_to_pyranges(enh)                             # predictor.py:84  → [start, end)
genes_pr = df_to_pyranges(genes,
              start_col="TargetGeneTSS", end_col="TargetGeneTSS",
              start_slop=window, end_slop=window,          # predictor.py:85-92 → [TSS−W, TSS+W)
              chrom_sizes_map=chrom_sizes_map)
pred = enh_pr.join(genes_pr).df.drop([...])                # predictor.py:93-95
pred["distance"] = abs(pred["enh_midpoint"] - pred["TargetGeneTSS"])   # predictor.py:96
pred = pred.loc[pred["distance"] < window, :]              # predictor.py:97
```

`window` comes from `--window`, default **5,000,000** (`predict.py:123-128`). It is not set by the
rule (`predictions.smk:37-52` passes no `--window`) and `parameters.predict.txt:21` confirms
`window 5000000`.

So the rule is:

> **A candidate element is paired with a gene iff they are on the same chromosome and
> `|midpoint(element) − TSS(gene)| < 5,000,000 bp`.**

It is a **distance window**, not a fixed number of nearest genes. There is no `head(k)`, no
`nsmallest`, no k-NN anywhere. Mean element width is 659 bp (measured on
`Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed`), i.e. 6.6×10⁻⁵ of the window, so the
interval-overlap join at `predictor.py:93` and the midpoint distance filter at `predictor.py:97`
select essentially the same set; the filter exists "for backwards compatability" (the code's own
comment) and trims the ≤659 bp fringe.

Therefore:

```
n_cand_pairs + n_nonexpressed_pairs  =  N_out
                                     =  Σ_chromosomes Σ_{e ∈ elements(c)} |{ g ∈ genes(c) : |mid(e) − TSS(g)| < W }|
                                     =  j_cand_elems · ḡ            where ḡ = mean genes per ±5 Mb window
```

**`n_cand_pairs` is O(`j_cand_elems` · ḡ), NOT O(`j_cand_elems` · `k_genes`).** Measured ḡ (all
pairs / `j_cand_elems`) is **88.7–97.3**, against `k_genes` = 20,531 — a **211–231× reduction**.
The full cross product would be 158,548 × 20,531 = **3.254×10⁹** rows; the window gives
**1.406×10⁷**.

The mechanism is directly visible in the per-chromosome pair counts the script prints. ḡ tracks
local TSS density, exactly as a fixed 10 Mb window predicts (thp1_1, multiome):

| chrom | elements | TSSs | pairs | ḡ | TSS/Mb | 10 Mb × density |
|---|---:|---:|---:|---:|---:|---:|
| chr18 | 2,943 | 308 | 115,153 | **39.1** | 3.8 | 38 |
| chr4 | 6,956 | 784 | 318,331 | 45.8 | 4.1 | 41 |
| chr8 | 8,757 | 703 | 451,993 | 51.6 | 4.8 | 48 |
| chr10 | 6,833 | 765 | 422,347 | 61.8 | 5.7 | 57 |
| chr6 | 10,313 | 1,077 | 833,254 | 80.8 | 6.3 | 63 |
| chr16 | 5,316 | 859 | 579,524 | 109.0 | 9.5 | 95 |
| chr11 | 7,991 | 1,292 | 989,025 | 123.8 | 9.6 | 96 |
| chr17 | 6,520 | 1,197 | 1,008,532 | 154.7 | 14.4 | 144 |
| chr19 | 6,318 | 1,423 | 1,673,689 | **264.9** | 24.3 | 243 |

ḡ spans **6.8×** across chromosomes and is predicted to within ~10–20 % by
`2 · window · (TSS density)`. If the rule were "all genes on the chromosome" chr18 would give
ḡ = 308, not 39.1; if it were "all genes genome-wide" it would give 20,531. The window is
unambiguously what caps the table.

(ḡ runs slightly *above* `2·W·density` because candidate elements are accessibility peaks and are
themselves enriched in gene-dense neighbourhoods, so the element-weighted local density exceeds
the chromosome-mean density.)

### 2.3 The join is an interval tree, not a nested loop

`df_to_pyranges` (`tools.py:124-142`) returns `pr.PyRanges`, and `PyRanges.join` (pyranges
**0.0.129**, `nb_cpu=1` default, `pyranges_main.py:2249-2260`) dispatches to
`pyranges/methods/join.py:_both_indexes`, which is:

```python
it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)   # join.py:12
_self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)   # join.py:15
```

`ncls` 0.0.68 — a **Nested Containment List**. Build over the gene intervals of chromosome *c* is
O(k_c log k_c); querying with j_c element intervals is O(j_c log k_c + hits). Result assembly is
`scdf.reindex(_self_indexes)` / `ocdf.reindex(_other_indexes)` (`join.py:117-118`) — fancy
indexing, O(pairs_c · n_cols).

**There is no j² anywhere on this path.** No element×element operation exists in `predict.py` or
in the powerlaw half of `predictor.py`. The only cross terms are element×gene, and those are
bounded by the window.

### 2.4 Everything after the join is a linear pass or a hash operation

Per chromosome, in order:

| step | code | class | note |
|---|---|---|---|
| chromosome mask + `.copy()` on the full enhancer and gene frames | `predict.py:252-253` | O(`j_cand_elems` + `k_genes`) **per chromosome** | 23 full scans of the 158 k-row enhancer table = 3.6 M row copies |
| `df["chr"].apply(lambda x: chrom_sizes_map[x])`, `.apply(lambda x: max(x,0))` | `tools.py:138-139` | O(k_c) Python-level calls | ~41 k lambda invocations total; cheap but not vectorised |
| `pr.PyRanges(df)` construction | `tools.py:142` → `pyranges/methods/init.py:45` `{k: v for k, v in df.groupby(grpby_key)}` | O(rows) groupby per construction | 2 per chromosome = 46 total; source of the ~90 `FutureWarning` lines in every log |
| NCLS build + query + reindex | `predictor.py:93` | O(k_c log k_c + j_c log k_c + pairs_c · n_cols) | §2.3 |
| `distance`, `powerlaw_contact`, `powerlaw_contact_reference` | `predictor.py:96`, `452-461` → `hic.py:257-268` (`np.clip`/`np.log`/`np.exp`) | O(pairs_c) vectorised | `powerlaw_contact_reference` is **only consumed by `scale_hic_with_powerlaw` (HiC-only)** — dead compute and a dead output column here |
| `annotate_predictions`: `isSelfPromoter`, then `pred.merge(genes_subset, on="TargetGene", how="left")`, then `isSelfGenic` | `predictor.py:522-551`, merge at `534-539` | O(pairs_c + k_c) hash join | key is an **object/string** column → Python string hashing constant. **Verified many-to-one:** `symbol` in `GeneList.txt` has 20,531 distinct values in 20,531 rows, max multiplicity 1 → no row multiplication. This closes the only place the code left the class ambiguous. |
| `compute_score(..., "ABC")` and `compute_score(..., "powerlaw")` | `predictor.py:62-67` and `69-74`, body at `501-519` | 2 × O(pairs_c) | each does `np.column_stack(...).prod(axis=1)` (`:502`) and a `groupby(["TargetGene","TargetGeneTSS"]).transform("sum")` (`:507-511`) — hash groupby on an (object, float64) composite key |

**The two `compute_score` calls are identical on the powerlaw path.** `predictor.py:62-67` passes
`[activity_base_enh, powerlaw_contact]`; `predictor.py:69-74` passes the same two arrays. Verified
on disk over 400,000 rows of `thp1_1`: `ABC.Score.Numerator == powerlaw.Score.Numerator` in
400,000/400,000 rows and `ABC.Score == powerlaw.Score` in 400,000/400,000 rows. The second
groupby-transform over 14.06 M rows is 100 % redundant work, and it adds 2 duplicate columns
(18 bytes/row) to the 282 MB output.

### 2.5 The write phase — where 96 % of the time is

`predict.py:274-298`:

```python
all_putative = pd.concat(all_putative_list)                          # :274
all_putative["CellType"] = args.cellType                             # :275
all_putative.loc[all_putative.TargetGeneIsExpressed, :].to_csv(       # :279-287
    all_pred_file_expressed, sep="\t", index=False, header=True,
    compression="gzip", float_format="%.6f", na_rep="NaN")
all_putative.loc[~all_putative.TargetGeneIsExpressed, :].to_csv(...)  # :288-296
test_variant_overlap(args, all_putative)                              # :298
```

Three things make this expensive, all of them fixed by the code as written:

1. **zlib level 9.** `compression="gzip"` reaches pandas `io/common.py:761-769`, which calls
   `gzip.GzipFile(filename=..., mode=...)` with an **empty** `compression_args` dict. `GzipFile`'s
   default `compresslevel` is **9**. Nothing in `predict.py` overrides it.
2. **`float_format="%.6f"` forces the Python formatting path.** `pandas/io/formats/csvs.py`
   chunks at `_DEFAULT_CHUNKSIZE_CELLS // len(cols)` = `100000 // 29` = **3,448 rows**
   (`csvs.py:174-177`), so the 10.17 M-row table is written in **2,951 chunks**, each doing
   `df.iloc[slicer]` (a copy), `_get_values_for_csv(float_format="%.6f")` — a per-element `%`
   format over 15 float columns — and `libwriters.write_csv_rows` (`csvs.py:314-330`).
3. **The volume is large and ~20 % of it is redundant.** Measured 269.3 uncompressed bytes/row for
   the expressed table and 272.0 for the non-expressed one:

| file | rows | uncompressed | on disk (gz) | ratio |
|---|---:|---:|---:|---:|
| `EnhancerPredictionsAllPutative.tsv.gz` | 10,173,077 | 2,740 MB | 281.5 MB | 9.7× |
| `…NonExpressedGenes.tsv.gz` | 3,890,952 | 1,058 MB | 85.0 MB | 12.4× |
| `…ForVariantOverlap…tsv.gz.tmp` | 128,527 | ~34 MB | 5.88 MB | 5.8× |
| `…ForVariantOverlap…tsv.gz` (re-read + `bedtools` + re-gzip) | 128,527 | ~34 MB | 6.21 MB | 5.5× |
| | | **~3.87 GB formatted+compressed** | **378.6 MB** | |

3.80 GB of CSV text produced in ~526 s = **7.2 MB/s**, on one core, i.e. **1.28 µs per output
cell**. That number is the whole rule.

**Direct measurement of the write phase.** I ran the exact `to_csv` call from `predict.py:279-287`
on the first 500,000 rows of the real `thp1_1` output table (Slurm job 40217360, 1 core, 24 GB,
same conda env as the pipeline: `.snakemake/conda/9ec2e87051abf89ef063c76fc96d82bf_`; script at
`/scratch/users/kaybrand/p1u3/bench.py`, log at `/scratch/users/kaybrand/p1u3/bench.log`):

| operation | s / 500 k rows | s / M rows | output bytes | note |
|---|---:|---:|---:|---|
| **`to_csv(compression="gzip", float_format="%.6f")` — what the code does** | **17.57** | **35.14** | 14,272,864 | |
| ├ pandas CSV assembly + `%.6f` formatting (`to_csv` uncompressed) | 9.22 | 18.44 | 134,741,537 | **52.5 %** of the total |
| └ zlib level 9 (by subtraction) | 8.35 | 16.70 | | **47.5 %** of the total |
| external `gzip -9` on those same bytes | 8.48 | 16.96 | 14,270,756 | confirms the subtraction |
| external `gzip -6` on those same bytes | 2.68 | 5.36 | 14,632,787 | **3.16× faster; output +2.5 %** |
| `to_csv(compression={"method":"gzip","compresslevel":1})` | 10.00 | 20.00 | 17,821,270 | output +24.9 % |
| `to_csv` uncompressed, **no** `float_format` | 6.46 | 12.92 | 127,861,779 | −30 % time, −5.1 % bytes |
| `pd.concat` of 4 copies (2 M rows) | 0.16 | 0.08 | | negligible |
| boolean `.loc` copy of 2 M rows | 0.49 | 0.25 | | negligible |
| `read_csv` of 500 k rows from the gz | 1.49 | 2.98 | | |

Three things this pins down:

- **The split is ~50/50 between pandas' Python-level float formatting and zlib level 9.** Neither
  alone explains the cost; both do, additively (9.22 + 8.48 = 17.70 s ≈ the 17.57 s in-process
  figure, so pandas adds no penalty of its own for compressing inline — it is simply serial).
- **The measured 35.14 s/M agrees with the 37.50 s/M write coefficient fitted from the 10
  uncontended production jobs (§7) to within 6.3 %.** That is an independent confirmation that
  §2.1's attribution of ~96 % of the rule to the write phase is correct, arrived at from a
  completely different direction.
- **269.5 uncompressed bytes/row** in this microbenchmark against 269.3 measured by `awk` on the
  production file — the two agree to 0.07 %.

Scaling forward: 14.06 M rows across the two big tables plus ~129 k for the variant-overlap `.tmp`
at 35.14 s/M = **498.7 s**, plus concat (~1.1 s), the two mask copies (~6.8 s), the variant-overlap
masks (~1 s) and the `bedtools` pipes (~4 s) = **~512 s**, against a measured non-loop non-startup
residual of **526 s** for `thp1_1` — agreement to **2.7 %**. The write phase is accounted for.

Redundant bytes in the 29-column output (per-column mean widths measured over 200,000 rows):

| col | field | bytes/row | why it is redundant on this path |
|---:|---|---:|---|
| 7 | `activity_base_enh` | 9.07 | byte-identical to col 6 `activity_base` — `predict.py:218` puts `activity_base` in the subset and `:219` copies it again |
| 8 | `activity_base_squared_enh` | 9.33 | computed at `predict.py:220`, never used in powerlaw scoring |
| 11 | `enh_idx` | 4.53 | internal join key (`predictor.py:82`); dropped **only** inside `add_hic_from_hic_file` (`predictor.py:257-271`) and `add_hic_from_directory` (`predictor.py:412-430`), both HiC-only → leaks into the output |
| 19 | `gene_idx` | 3.72 | same, `predictor.py:83` |
| 24 | `powerlaw_contact_reference` | 9.00 | `predictor.py:459-461`; consumed only by the HiC-only `scale_hic_with_powerlaw` |
| 27 | `powerlaw.Score.Numerator` | 9.00 | byte-identical to col 25 (§2.4) |
| 28 | `powerlaw.Score` | 9.00 | byte-identical to col 26 (§2.4) |
| | **total** | **53.65 / 269.3** | **20.1 % of uncompressed bytes** (measured directly) |

Separately, `float_format="%.6f"` is applied to two float columns holding integers:
`enh_midpoint` (15.58 bytes/row, e.g. `12345.500000`) and `distance` (14.70 bytes/row, e.g.
`1234567.000000`). Roughly 14 of those ~30 bytes are trailing-zero padding.

**Memory.** `all_putative_list` is created at `predict.py:234`, appended at `:264`, concatenated at
`:274`, and **never cleared or rebound** — the 23 per-chromosome frames stay resident through all
three writes alongside the concatenated table. Plus each `.loc[mask, :]` at `:279`/`:288` is a full
copy of its slice. A single copy of the 29-column frame is ~211 B/row (6 object pointers, 3 bool,
5 int64, 15 float64). Measured peak RSS / total output rows is **574.4–579.8 bytes** across all 12
jobs — a **0.9 % spread**, and ≈2.7× the single-copy estimate, consistent with list + concat +
transient mask copy all resident.

### 2.6 `getVariantOverlap.py` — live, but its product is consumed by nothing

The prompt asks whether it runs at all. **It does.** `predict.py:9` imports
`test_variant_overlap` and `predict.py:298` calls it unconditionally — there is no `if` guard. The
file it writes exists on disk for all 12 jobs (6.15–6.24 MB), and the two `zcat | … | gzip` pipes
appear verbatim in every Slurm log.

What it costs: six boolean masks over the full 14.06 M-row table (`getVariantOverlap.py:26-31`),
`dropna` (`:33`), a `distance <= 2000000` filter (`:34-36`), a `to_csv(gzip)` of ~128 k rows
(`:38-45`), then `zcat|head -1|gzip` and `zcat|sed 1d|bedtools slop -b -150|gzip`
(`:49-63` via `tools.py:37-70`). From file mtimes on thp1_1 the `.tmp` lands at 18:03:10 and the
final file at 18:03:14, so the two subprocess pipes take **~4 s**; the masks over 14 M rows are
vectorised and cost O(1 s). Total ≈ **5–8 s of the 549 s (~1 %)**. Small.

Three facts worth recording rather than characterising further:

- `getVariantOverlap.py:30` computes `all_putative[(score_t & not_promoter) | (is_promoter & score_one)]`
  and **discards the result**, then `:31` recomputes the identical expression and assigns it. A
  duplicated (cheap, ~128 k-row) materialisation.
- `EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.tsv.gz` is referenced **nowhere
  else in the repository**. `grep -rn ForVariantOverlap` over all `.smk`/`.py`/`.R`/`.yaml`
  matches only its own definition at `getVariantOverlap.py:23`. It is not a declared rule output,
  so Snakemake does not track it and cannot clean it.
- The `.tmp` file (`:37`) is never unlinked — 5.83–5.92 MB of stale intermediate is left behind by
  every job (12 files, 70 MB across the two runs).

`getVariantOverlap.py:68-71` (the `__main__` block, which would `read_csv` the whole 282 MB
prediction file) is genuinely dead — nothing invokes the script standalone.

---

## 3. Size variables

Named exactly as in `size_manifest_{multiome,scatac}.tsv`:

| variable | role in this rule | measured range (both modalities identical) |
|---|---|---|
| **`j_cand_elems`** | rows of the enhancer input; the outer term of the pair generation | 156,420 – 158,548 (**1.014×**) |
| **`k_genes`** | rows of the gene input; NCLS is built over these | 20,531 (**exactly constant**) |
| **`n_cand_pairs`** | rows of the primary output — the dominant cost driver | 10,173,077 – 11,447,564 (**1.13×**) |
| `n_cand_pairs_bytes` | on-disk size of that output | 281,522,213 – 335,277,506 (1.19×) |
| `enh_list_bytes` | input volume | 32,939,288 – 33,323,136 (1.012×) |

Two derived quantities the manifest does not carry, both needed to state the cost honestly:

- **`N_out`** = total putative pairs generated = `n_cand_pairs` + non-expressed pairs. Read
  directly off the script's own per-chromosome log lines: **14,064,029 – 15,218,153**.
  `N_out / n_cand_pairs` = 1.329 – 1.382 (mean **1.3593**).
- **ḡ** = `N_out / j_cand_elems` = mean genes within ±5 Mb of a candidate element = **88.7 – 97.3**.

The cost is driven by `N_out`, not by `n_cand_pairs` — the rule pays full pair-generation and full
CSV-serialisation cost for the 27 % of rows that go to the non-expressed file. Normalising the 10
uncontended jobs by `N_out` is tighter (38.53–40.36 s/M, ±2.3 %) than normalising by
`n_cand_pairs` (52.27–54.70 s/M, ±4.5 %), which is the empirical signature of exactly that.

**Why `N_out/n_cand_pairs` sits near 0.73⁻¹, and why it has nothing to do with RNA.**
`Expression` is `NaN` for all 20,531 rows of `GeneList.txt` — **in both modalities.**
`determine_expressed_genes` (`tools.py:107-113`) therefore falls entirely through to
`PromoterActivityQuantile >= promoter_activity_quantile_cutoff`, default **0.30**
(`predict.py:103-108`). So `isExpressed` is decided by an ATAC promoter-activity quantile, ~70 % of
genes pass, and the pair-weighted fraction is 0.723–0.752 (slightly above 0.70 because expressed
genes sit in gene-dense regions with more nearby elements). **The ABC path never sees RNA at all.**

The cross term the code genuinely has is:

```
O( n_chrom · (j_cand_elems + k_genes)          # 23 chromosome masks, predict.py:252-253
 + Σ_c ( k_c log k_c + j_c log k_c )           # NCLS build + query, join.py:12-15
 + N_out · n_out_cols )                        # materialise, score, format, compress
```

with `n_chrom` = 23 (22 autosomes + chrX; chrY discarded at `predict.py:239-240`, and the set is
the gene∩enhancer chromosome intersection at `predict.py:238`) and `n_out_cols` = 29.

---

## 4. Asymptotic class derived from the code

**Overall: O(`N_out` · `n_out_cols`) = O(`j_cand_elems` · ḡ · 29), i.e. linear in `j_cand_elems`
with a large constant, plus an O(`j_cand_elems` log `k_genes`) query term that is dominated at
every realistic scale.**

Term by term, each with the citation that fixes it:

| term | cost | citation |
|---|---|---|
| read inputs | O(`j_cand_elems` + `k_genes`) | `predict.py:175`, `:180` |
| 23 chromosome masks + copies | O(`n_chrom` · (`j_cand_elems` + `k_genes`)) | `predict.py:252-253` |
| PyRanges construction (groupby per call, 46 calls) | O(`j_cand_elems` + `k_genes`) total | `tools.py:142` → `pyranges/methods/init.py:45` |
| chrom-size slop via `.apply` | O(`k_genes`) Python calls | `tools.py:138-139` |
| NCLS build over genes | O(Σ_c `k_c` log `k_c`) = O(`k_genes` log `k_genes`) | `pyranges/methods/join.py:12` |
| NCLS overlap query | O(Σ_c `j_c` log `k_c` + `N_out`) | `pyranges/methods/join.py:15` (`all_overlaps_both`) |
| join result materialisation | O(`N_out` · cols) | `pyranges/methods/join.py:117-118` (`reindex`) |
| distance + 2× powerlaw | O(`N_out`) vectorised | `predictor.py:96`, `452-461`; `hic.py:257-268` |
| gene-bounds merge on `TargetGene` | O(`N_out` + `k_genes`), **strictly** many-to-one | `predictor.py:534-539`; uniqueness verified §2.4 |
| 2× score + 2× groupby-transform | 2 × O(`N_out`) hash | `predictor.py:501-511`, called `:62-67` and `:69-74` |
| `pd.concat` | O(`N_out` · cols) copy | `predict.py:274` |
| **3× `to_csv(gzip level 9, float_format)`** | **O(`N_out` · cols)**, dominant constant | `predict.py:279-296`, `getVariantOverlap.py:38-45`; level 9 via `pandas/io/common.py:761-769` |
| variant-overlap masks + `bedtools` pipe | O(`N_out`) + O(small) | `getVariantOverlap.py:26-63` |

**Not O(`j_cand_elems` · `k_genes`).** The window at `predict.py:126-128` /
`predictor.py:85-97` replaces `k_genes` with ḡ ≈ 89–97. That single parameter is the difference
between 1.4×10⁷ and 3.25×10⁹ rows — a **231×** factor, and it is the reason this rule takes 9
minutes rather than a day and a half.

**Not O(`j_cand_elems`²).** There is no element×element construct on the powerlaw path. Every
grouping (`predictor.py:507`) is by target gene, i.e. over pairs, not over element×element.

**Not superlinear in `N_out`.** Every operation on the concatenated table is a vectorised pass, a
hash groupby, or a chunked serialisation. The one place key cardinality could have produced a
row-multiplying join (`predictor.py:534`) is verified many-to-one.

The class is therefore **unambiguous from the code**, so per the instructions `t_low` and `t_high`
below differ in their **constant only**, not in their exponent.

---

## 5. Benchmark cross-check

All 12 jobs (`/scratch/users/kaybrand/scE2G_optimize_results/analysis/benchmarks_*.tsv`):

| modality | cluster | `s` | `cpu_time` | cpu/s | `max_rss` MB | B/row | s per M `N_out` | s per M `n_cand_pairs` | ḡ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| multiome | thp1_1 | 549.3 | 525.8 | 0.957 | 8078 | 574.4 | 39.06 | 54.00 | 88.7 |
| scatac | thp1_1 | 556.5 | 515.0 | 0.925 | 8080 | 574.5 | 39.57 | 54.70 | 88.7 |
| multiome | thp1_2 | 555.0 | 525.8 | 0.948 | 8143 | 574.6 | 39.16 | 53.99 | 90.0 |
| scatac | thp1_2 | 552.2 | 525.9 | 0.952 | 8156 | 575.5 | 38.97 | 53.73 | 90.0 |
| multiome | jurkat_pma_cd3_4hr | 567.4 | 555.9 | 0.980 | 8491 | 576.7 | **38.53** | 52.27 | 93.5 |
| scatac | jurkat_pma_cd3_4hr | 577.1 | 556.0 | 0.963 | 8500 | 577.3 | 39.20 | 53.17 | 93.5 |
| multiome | jurkat | 586.9 | 556.3 | 0.948 | 8614 | 577.1 | 39.32 | 53.18 | 94.2 |
| scatac | jurkat | 582.0 | 549.2 | 0.944 | 8612 | 576.9 | 38.99 | 52.74 | 94.2 |
| multiome | k562_crispri | 600.9 | 585.2 | 0.974 | 8817 | 579.4 | 39.49 | 52.49 | 97.3 |
| scatac | k562_crispri | 614.2 | 557.7 | 0.908 | 8823 | 579.8 | 40.36 | 53.65 | 97.3 |
| multiome | telohaec_crispri | 634.4 | 562.3 | 0.886 | 8691 | 576.9 | *42.11* | *57.12* | 95.6 |
| scatac | telohaec_crispri | **822.6** | 731.0 | 0.889 | 8701 | 577.6 | *54.61* | *74.06* | 95.6 |

### Does the observed scaling agree with the code?

**Yes, and unusually cleanly for this project — but only because the driving variable is
`n_cand_pairs` (1.13×), which is Tier 2, not Tier 3.**

- **Linear in `N_out`:** over the 10 jobs on the 5 uncontended clusters, `s / N_out` = **38.53 –
  40.36 s per million rows**, a **±2.3 % spread** while `N_out` moves 1.082×. Predicted-linear,
  observed-linear.
- **Memory strictly linear in `N_out`:** `max_rss / N_out` = **574.4 – 579.8 B/row** across all 12
  jobs, a **0.9 % spread** — including both contended telohaec jobs, since contention costs time
  but not memory. This is the single tightest scaling law I found.
- **Single core:** `cpu_time / s` = 0.886 – 0.980, mean 0.94. No `threads:` is declared
  (`predictions.smk` has none) and Slurm confirms `AllocCPUS=1`. The 0.94 rather than 1.00 is not
  hidden parallelism; it is the ~4 s conda/import startup plus the `zcat|sed|bedtools|gzip` pipe,
  where the parent Python blocks in `communicate()` (`tools.py:62`) while four short-lived
  subprocesses run. There is no `--threads` flag to lose here: zlib level 9 inside
  `gzip.GzipFile` is serial by construction, and `PyRanges.join` defaults to `nb_cpu=1`
  (`pyranges_main.py:2257`). **Declared and realised parallelism agree at 1.**
- **Fitting an intercept fails, and that is the honest result.** Least squares of the 5-cluster
  two-run means on `N_out` gives slope 45.18 s/M with intercept **−86.4 s** (R² = 0.946). A
  negative intercept is unphysical; it is what a 1.08× lever arm plus 2 % noise produces. I do not
  read it as superlinearity. I calibrate §7 by fixing the mechanistically known constants (startup
  from mtimes, loop from the script's own timers) and fitting only the per-row write coefficient,
  which lands at 37.10 – 38.18 s/M — **±1.4 %**.

### Replicate noise, and the one disagreement

This rule sits upstream of the ARC/Kendall branch point, so the two runs are true n=2 replicates.
`GeneList.txt` md5 is identical between modalities for all six clusters; `n_cand_pairs`,
`j_cand_elems`, ḡ and all output byte counts match; and
`EnhancerPredictionsAllPutative.tsv.gz` is **bit-identical** between the two runs for 5 of 6
clusters (`checksums_*.tsv`).

| cluster | multiome `s` | scatac `s` | ratio |
|---|---:|---:|---:|
| jurkat | 586.9 | 582.0 | 0.992 |
| thp1_2 | 555.0 | 552.2 | 0.995 |
| thp1_1 | 549.3 | 556.5 | 1.013 |
| jurkat_pma_cd3_4hr | 567.4 | 577.1 | 1.017 |
| k562_crispri | 600.9 | 614.2 | 1.022 |
| **telohaec_crispri** | 634.4 | **822.6** | **1.297** |

Five clusters replicate to within **2.2 %**. telohaec_crispri differs by **29.7 %** — and its
scATAC job's *self-reported* per-chromosome loop time is also 2.1× the multiome one (32.19 s vs
16.96 s) at bit-comparable input, which proves the discrepancy is node/Lustre contention and not
data. This matches the briefing's stated noise floor exactly. **The only apparent code-vs-data
disagreement — telohaec being 30 % slow — is measurement noise, and the largest single job in the
table (822.6 s) is not a data point about scaling.** Excluding it, `s/N_out` is flat.

Small caveat on the replicate claim, recorded for accuracy: `EnhancerList.txt` is *not*
bit-identical between the two runs (quantile tie-breaking upstream in `create_neighborhoods`).
For thp1_1 and jurkat the columns `predict.py` actually reads are unaffected —
`activity_base` differs in **0 / 158,548** and **0 / 158,397** rows respectively, hence the
bit-identical outputs. For telohaec_crispri `activity_base` differs in ~0.85 % of rows at the 6th
decimal, which is why its output differs by 452 bytes out of 316 MB (1.4×10⁻⁶) and its md5
differs. Row counts are identical for all six. For cost purposes these are replicates.

### What the code says that the data cannot see

- `j_cand_elems` spans **1.014×** and `k_genes` is **exactly constant**. Nothing in this dataset
  can confirm any exponent in either. The O(`j_cand_elems` · ḡ) claim and the "not
  O(`j·k_genes`)" claim rest entirely on `predict.py:126-128` + `predictor.py:85-97`, cross-checked
  by the per-chromosome ḡ-vs-density table in §2.2 — which is a *within-run* check across 23
  chromosomes spanning 6.8× in ḡ, not a cross-cluster fit.
- `ln(k_genes)` is constant, so the NCLS log factor is invisible by construction.

---

## 6. Concrete overhead sources

Ordered by measured cost.

1. **zlib level 9 + pandas' Python float-formatting path, on one core: ~526 s of 549 s (96 %).**
   3.80 GB of CSV text at 7.2 MB/s. Fixed by `compression="gzip"` (`predict.py:283`) reaching
   `gzip.GzipFile` with default `compresslevel=9` (`pandas/io/common.py:761-769`) and by
   `float_format="%.6f"` (`predict.py:285`) forcing per-element Python `%` formatting over 15 float
   columns in 2,951 chunks of 3,448 rows (`pandas/io/formats/csvs.py:174-177, 314-330`).
   **Measured split (§2.5): 52.5 % pandas formatting/CSV assembly, 47.5 % zlib level 9.** For
   reference at fixed input bytes, level 9 costs 3.16× the CPU of level 6 and yields 2.5 % smaller
   output; dropping `float_format` cuts the formatting half by 30 %. Recording the cost structure,
   not proposing a change — Phase 2 decides what matters.
2. **~20.1 % of the written bytes are redundant on this path** (measured, §2.5 table): two exactly
   duplicated column pairs, two internal join keys that escape because their `drop` lives inside
   HiC-only functions (`predictor.py:257-271`, `412-430`), and two HiC-only quantities computed
   unconditionally.
3. **One fully redundant 14 M-row hash groupby.** `compute_score` is called twice with identical
   arguments on the powerlaw path (`predictor.py:62-67` vs `69-74`); verified bit-identical
   outputs. Part of the 7.2 s post-join loop cost.
4. **~2× peak memory from a list that is never released.** `all_putative_list` (`predict.py:234`)
   stays referenced through `pd.concat` (`:274`) and all three writes. 574–580 B/row measured
   versus ~211 B/row for a single copy. Peak 8.08–8.82 GB against a 20 GB reservation
   (`predictions.smk:36`, `min_gb=20`) — **40–44 % utilised**; Slurm's cgroup `MaxRSS` is
   8.63–9.45 GB (it also counts the nested Snakemake interpreter).
5. **A ~34 MB round-trip and a leaked temp file.** `getVariantOverlap.py:38-45` writes a gzipped
   `.tmp`, then `:49-63` reads it back twice through `zcat` and re-compresses it twice
   (header, then body through `bedtools slop`), buffering the whole result in the parent's memory
   via `communicate()` (`tools.py:62`) before writing. ~4 s. The `.tmp` is never unlinked.
6. **The product of step 5 is consumed by nothing.** `EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.tsv.gz`
   has no reader anywhere in the repository, and is not a declared output so Snakemake neither
   tracks nor cleans it (§2.6).
7. **23 full scans of the enhancer table.** `predict.py:252-253` re-masks and `.copy()`-s the
   whole 158 k-row frame once per chromosome. Small in absolute terms (part of the 17 s loop) but
   it is a repeated pass over the same data, and it is O(`n_chrom` · `j_cand_elems`).
8. **Nested-Snakemake wrapper: 26–46 s per job (5–7 %).** Slurm `Elapsed` (9:39–11:20) exceeds the
   `benchmark:` `s` (9:09–10:34) by that much: the profile submits a Slurm job that starts a second
   Snakemake, which builds its own DAG, activates conda, and afterwards runs "Storing output in
   storage" (visible in every log, and it is why the two declared outputs' mtimes are *later* than
   the job's own final log line while the undeclared variant-overlap file's mtime is not).
   `benchmark:` does not include this.
9. **~90 pandas `FutureWarning` lines per log** from `pyranges/methods/init.py:45` (two
   `PyRanges` constructions × 23 chromosomes). Cosmetic, but each one is a `groupby` on the
   Chromosome column at construction time.
10. **`io_in`/`io_out` cross-check.** `io_out` = 340.4–412.0 MB against a true on-disk
    378.6–430.9 MB, so for *this* rule it happens to land within ~10 % (it misses the
    subprocess-written variant-overlap file). `io_in` = 37.1–138.4 MB against a true 38.9 MB
    input; the 122.9/138.4 MB readings on the two telohaec jobs are spurious. Consistent with the
    briefing: I use on-disk sizes and `*_bytes`, not these columns.

---

## 7. `t_est`, `t_low`, `t_high`

Variables named exactly as in `size_manifest_*.tsv`. All formulas in **seconds**. Each is a plain
arithmetic expression in `n_cand_pairs` and `j_cand_elems` — directly evaluable.

The class is unambiguous from the code (§4), so **all three have the same functional form and
differ only in their constants.** They are not a blend and must not be averaged.

```
t_est(n_cand_pairs, j_cand_elems)  =  6.5  +  1.60e-5 * j_cand_elems  +  5.253e-5 * n_cand_pairs

t_low(n_cand_pairs, j_cand_elems)  =  4.0  +  1.00e-5 * j_cand_elems  +  4.550e-5 * n_cand_pairs

t_high(n_cand_pairs, j_cand_elems) = 10.0  +  3.00e-5 * j_cand_elems  +  7.600e-5 * n_cand_pairs
```

### How `t_est` is built (each constant is mechanistic, not fitted to the exponent)

| constant | value | origin |
|---|---|---|
| fixed | 6.5 s | conda activate + `import hicstraw/pandas/pyranges` + argparse + `write_params` ≈ 4 s (mtime of `parameters.predict.txt` minus job start), + ~2.5 s for the two `read_csv` calls |
| per `j_cand_elems` | 1.60e-5 s | `read_csv` of the 33 MB `EnhancerList.txt` plus the 23 chromosome masks at `predict.py:252-253` — 2.52 s at `j` = 157,635 |
| per `N_out`, loop | 1.15e-6 s | the script's own timers: 16.9 s at `N_out` = 14.69 M (NCLS join, powerlaw, merge, 2× `compute_score`) |
| per `N_out`, write | 3.750e-5 s | residual after subtracting the above from the 5 uncontended clusters' two-run means: 37.10, 37.28, 37.40, 37.52, 38.18 s/M → mean **37.50**, spread **±1.4 %** |
| `N_out` → `n_cand_pairs` | ×1.3593 | measured `N_out / n_cand_pairs` mean over the six clusters (range 1.329–1.382) |

Combined per-`n_cand_pairs` coefficient: `(1.15e-6 + 3.750e-5) × 1.3593 = 5.253e-5`.

Residuals of `t_est` against the two-run means for the five uncontended clusters: thp1_1 −1.7 %,
thp1_2 −0.8 %, jurkat_pma +1.2 %, jurkat +0.7 %, k562 +0.4 %. **Max error 1.7 %.**

### What `t_low` and `t_high` bound

Same class; the constant is uncertain because of, in order of measured magnitude:

1. **Node and Lustre contention — the dominant term.** Measured directly: 0.992–1.022 replicate
   ratio on five clusters, **1.297** on telohaec_crispri. Across all 12 jobs `s/N_out` spans
   38.53 → 54.61 s/M = **1.42×**. `t_high`'s per-pair coefficient is set from the worst observed
   rate (74.06 s per M `n_cand_pairs`, scATAC/telohaec) rounded up; `t_low`'s from the best
   observed (52.27 s/M, multiome/jurkat_pma) shaved ~13 % for a quiet node with faster storage.
2. **The expressed fraction, which sets `N_out / n_cand_pairs`.** 1.329–1.382 here, driven purely
   by `promoter_activity_quantile_cutoff = 0.30` with all-NaN `Expression` (§3). On a cluster where
   real RNA populated `GeneList.Expression`, or where the cutoff moved, this multiplier could
   plausibly reach ~2.0 — which is inside `t_high` at fixed `n_cand_pairs`.
3. **ḡ, the mean genes per ±5 Mb window (88.7–97.3, 1.10×).** For a fixed reference and window
   this is stable, but it is the lever that would move if the TSS annotation, the window, or the
   spatial distribution of candidate elements changed. It enters only through `N_out`.
4. **Mean output row width (269.3 B/row here).** Fixed by the 29 columns and `%.6f`; it would move
   if the column set changed.

Both bounds preserve the exponents. I deliberately did **not** give `t_low` a lower exponent
because the measurements look flat: they look flat *because the class is linear*, and the code says
so at `predictor.py:78-106` and `predict.py:279-296`.

---

## 8. Evaluated (reporting convenience only — the functions are the deliverable)

Cluster sizes: `j_cand_elems` and `n_cand_pairs` from `size_manifest_*.tsv` (identical in both
modalities).

| | `j_cand_elems` | `n_cand_pairs` | `t_low` | `t_est` | `t_high` | measured (multiome / scATAC) |
|---|---:|---:|---:|---:|---:|---|
| **mean of 6 clusters** | 157,635 | 10,816,106 | **497.7 s** (8.29 min) | **577.2 s** (9.62 min) | **836.8 s** (13.94 min) | mean of 12 jobs = 599.9 s |
| **min — `thp1_1`** | 158,548 | 10,173,077 | **468.5 s** (7.81 min) | **543.4 s** (9.06 min) | **787.9 s** (13.13 min) | 549.3 / 556.5 |
| **max — `telohaec_crispri`** | 157,551 | 11,106,687 | **510.9 s** (8.52 min) | **592.5 s** (9.87 min) | **858.8 s** (14.31 min) | 634.4 / **822.6** |

`t_est` is 1.7 % low at `thp1_1` and 6.6 % low against telohaec's multiome measurement (and 28 %
low against its contended scATAC one, which `t_high` covers).

Note that "min" and "max" here are by cluster *depth* (`n_frag`), per the instructions. By this
rule's actual driving variable the ordering differs: the largest `n_cand_pairs` is
**`k562_crispri`** at 11,447,564, where `t_est` = **610.3 s** against a measured 600.9 / 614.2 —
a **+0.4 %** residual. `n_frag` spans 15.1× across these clusters and is irrelevant to this rule;
`n_cand_pairs` spans 1.13× and is what governs it.

Per-job resource prediction, if Phase 2 wants it: **`max_rss` ≈ 577 bytes × `N_out`
≈ 784 bytes × `n_cand_pairs`** (0.9 % spread over 12 jobs; 8.0–8.8 GB observed against a 20 GB
reservation).

---

## 9. Both modalities: does the class differ, or only the constant?

**Neither. The class *and* the constant are the same, because the rule receives the same data.**

`create_predictions` runs upstream of `checkpoint features_required`, so the ARC/Kendall branch has
not yet diverged. Concretely:

- `j_cand_elems`, `k_genes`, `n_cand_pairs`, `n_cand_pairs_bytes`, `enh_list_bytes` and ḡ are
  **numerically identical** row-for-row between `size_manifest_multiome.tsv` and
  `size_manifest_scatac.tsv` for all six clusters.
- `GeneList.txt` md5 matches between modalities for all six clusters
  (`checksums_multiome.tsv` / `checksums_scatac.tsv`).
- `EnhancerPredictionsAllPutative.tsv.gz` is bit-identical between modalities for 5 of 6 clusters
  and differs by 452 bytes in 316 MB (1.4×10⁻⁶) for telohaec_crispri (§5).
- `Expression` is all-NaN in `GeneList.txt` in **both** modalities, so even the notionally
  RNA-dependent `isExpressed` split (`tools.py:107-113`) is computed identically from ATAC
  promoter-activity quantiles. The ABC path never sees RNA.

So the 12 jobs are **6 clusters × 2 replicates**, not two modality arms. No modality term belongs
in §7, and the measured `s/N_out` distributions overlap completely (multiome 38.53–39.49 excluding
telohaec; scATAC 38.97–40.36 excluding telohaec).

---

## 10. What the measurement cannot tell us

Explicitly and specifically:

1. **Nothing about any exponent in `j_cand_elems` or `k_genes`.** `j_cand_elems` spans **1.014×**
   (156,420–158,548) and `k_genes` is **exactly constant** at 20,531. The O(`j_cand_elems` · ḡ)
   class, and the crucial negative claim that it is *not* O(`j_cand_elems` · `k_genes`), come
   entirely from `predict.py:126-128` and `predictor.py:85-97`. My strongest empirical support is
   the within-run per-chromosome table in §2.2 — 23 points spanning 6.8× in ḡ, agreeing with
   `2·window·density` — which tests the *mechanism*, not the cluster-level exponent.
2. **Nothing about the NCLS log factor.** `ln(k_genes)` is constant by construction. O(`j log k`)
   and O(`j`) are indistinguishable here and always will be on a fixed genome annotation.
3. **`window` was never varied.** It is 5,000,000 in all 12 jobs. That `N_out` scales ~linearly in
   `window` (while `window` ≪ chromosome length) is a code-derived prediction, untested. It is
   also the single largest lever on this rule's cost and the measurements say nothing about it.
4. **The formatting-vs-compression split *is* now measured, but only at one point.** §2.5 puts it
   at 52.5 % / 47.5 % on 500,000 real rows of the `thp1_1` table (Slurm job 40217360). Two limits:
   (a) those are the first 500 k rows, i.e. chr1 and expressed-only, whereas the non-expressed
   table has more `NaN`s and compresses 12.4× rather than 9.7× — its per-row cost is probably
   slightly lower and I did not measure it; (b) the microbenchmark ran on one node at one moment,
   so its 35.14 s/M carries the same contention exposure as any single job here. The agreement
   with the production-fitted 37.50 s/M to within 6.3 % is reassuring but is not a replicate.
5. **The expressed fraction was never varied.** `N_out / n_cand_pairs` ∈ [1.329, 1.382] across all
   six clusters, because `Expression` is uniformly NaN and the split is a fixed 0.30 promoter-
   activity quantile. Any dataset with real gene expression in `GeneList.txt` breaks that
   multiplier, and `t_est` written in `n_cand_pairs` inherits the error. Phase 2 should prefer the
   `N_out` form if it ever has `N_out` available.
6. **The HiC branch produced zero measurements.** All six clusters have an empty `HiC_file`;
   `parameters.predict.txt` records `hic_file None`; neither powerlaw feature table contains
   `contactFrequency`. I can say the branch *exists* and where — `predictor.py:29-60` dispatching to
   `add_hic_from_hic_file:211` or `add_hic_from_directory:280`, plus `qc_hic:482`,
   `fill_diagonals:109`, `create_df_from_records:136` — and I can say the rule's own comment
   anticipates it being expensive (`predictions.smk:36`: *"Use 100GB if using average HiC"*). I
   **cannot** give it a class or a constant. Two structural features are visible by reading and
   would need their own unit if HiC is ever enabled: `add_hic_from_hic_file` re-`merge`s the entire
   per-chromosome `pred` table once per 8,000-row matrix block inside a loop
   (`predictor.py:231-251`), and `fill_diagonals` (`predictor.py:122-133`) is a Python-level
   per-row loop with a nested bin search. Neither is characterised here, and neither should be
   assumed cheap or expensive from this report.
7. **`telohaec_crispri` cannot be used as a scaling point.** Its two replicates differ by 29.7 %
   (634.4 vs 822.6 s) on bit-comparable input, and even the script's own internal loop timer
   differs 2.1× (16.96 vs 32.19 s). Anything below ~30 % at that cluster is contention.
   `n_frag` = 394.8 M makes it the "largest" cluster in the briefing's Tier 1 sense, but `n_frag`
   does not enter this rule at all.
8. **`io_in`/`io_out` are not used.** For this rule `io_out` happens to land within ~10 % of the
   true on-disk volume, but `io_in` reports 37–138 MB for a fixed 38.9 MB input, so I treat both as
   unusable per the briefing and use on-disk sizes and the `*_bytes` manifest columns instead.
9. **`benchmark:` excludes 26–46 s of per-job nested-Snakemake wrapper** (Slurm `Elapsed` minus
   `s`, 5–7 %). `t_est` in §7 is calibrated to `benchmark: s`, so Phase 2 should add that wrapper
   term separately if it is modelling wall-clock occupancy rather than in-script time.
10. **Six clusters, one reference genome, one gene annotation, one `nStrongestPeaks` cap.** ḡ,
    `k_genes`, `n_chrom` = 23 and the 29-column output width are all properties of that single
    configuration, and every constant in §7 is conditioned on it.

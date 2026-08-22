# Phase 1 unit #1 (Tier A) — `compute_kendall`

Author: Phase 1 characterisation agent, 2026-08-21.
Scope: characterisation only. No code edits, no fix proposals, no prioritisation.

---

## 0. Headline

The central question posed to this unit was: **is `fast_kendall_sc`'s per-pair Kendall τ an
O(c log c) merge-sort or an O(c²) naive double loop, where c = `n_cells`?**

**Answer: neither.** The compiled path this pipeline actually takes is asymptotically *better*
than the classical O(c log c). It is

> **O( m_π · log(q_γ · m_π) ) per (gene γ, peak π) pair**, where
> `m_π` = number of cells in which peak π is accessible (nnz of that ATAC column) and
> `q_γ` = number of cells in which gene γ is detected (nnz of that RNA column).

Both `m_π` and `q_γ` are ≪ `n_cells`: measured pair-weighted mean `m_π` is **2.8 %–7.2 % of
`n_cells`**. The algorithm never touches a zero cell. It is a closed-form rank-gap accumulation
over the peak's accessible cells only, cited at `src/lib.rs:312-360`.

I obtained the actual Rust source (see §1.1) rather than inferring it, and separately confirmed it
against symbols and instantiated monomorphisations in the installed `.so`.

Two secondary findings that matter more than the exponent at current scale:

1. At the largest cluster the Rust kernel is **43.6 %** of wall time and the single line
   `writeMM(data.RNA, rna_mtx_path)` (`compute_kendall.R:60`) — ASCII MatrixMarket serialisation
   of the normalised RNA matrix for the `system2` handoff — is **29.6 %**. At the smallest cluster
   the kernel is only **9.9 %**, and fixed interpreter/annotation startup is **30 %**.
2. The R→Python `system2` boundary round-trips **0.82 GB (thp1_1) to 3.96 GB (telohaec_crispri)**
   of plain ASCII through `/tmp`, and that volume is charged against the job's Slurm memory
   cgroup (§6.2). It is a first-class cost, not a rounding error.

---

## 1. Unit identity

| item | value |
|---|---|
| Rule | `compute_kendall`, `workflow/rules/compute_kendall.smk:4` |
| Entry point | `workflow/scripts/feature_computation/compute_kendall.R` — **298 lines** (R; `script:` directive at `compute_kendall.smk:45-46`) |
| Numeric CLI | `workflow/scripts/feature_computation/compute_kendall.py` — **72 lines** (invoked via `system2`, `compute_kendall.R:92-104`) |
| Python library | `fast_kendall_sc==1.0.0` (`workflow/envs/sc_e2g.yml:30`) → `fast_kendall_sc/metrics.py` — **206 lines** |
| Rust core | `_fast_kendall_sc` cdylib, `src/lib.rs` — **440 lines** (pyo3 0.22.6, numpy 0.22.1, rayon 1.12.0) |
| Modality | **MULTIOME ONLY.** 0 rows in `benchmarks_scatac.tsv`; 6 rows in `benchmarks_multiome.tsv`. |
| Jobs | **6** (one per cluster) |
| Declared resources | `mem_mb = determine_mem_mb(min_gb=63)` (`compute_kendall.smk:39`), `runtime = attempt*12*60` min (`:40`), `threads = config["threads"]` = **1** (`:41-42`, `configs/igvf10_multiome_config.yaml:70`) |
| Inputs | `Kendall/Pairs.tsv.gz`, `Kendall/atac_matrix.rds`, the cluster's 10x RNA directory |
| Outputs | `Kendall/Pairs.Kendall.tsv.gz`, `umi_count.txt`, `cell_count.txt`, `Kendall/gene_expression_metrics.tsv.gz` |

### 1.1 How I got the Rust source, and why you can trust it

`fast_kendall_sc` ships as a prebuilt `manylinux2014` wheel; the installed package contains only
`__init__.py`, `metrics.py` and a 914 KB stripped-ish `.so`. PyPI, however, also publishes an
**sdist** for 1.0.0. I downloaded it to `$SCRATCH/fks_src/`:

```
fast_kendall_sc-1.0.0.tar.gz
sha256 2f0e56aae7a3b082b403f5d89fceeb345563c467e5ab4380a8f1120cc7e54844
→ src/lib.rs (440 lines), fast_kendall_sc/metrics.py, Cargo.toml, Cargo.lock
```

Provenance check: `metrics.py` and `__init__.py` in the sdist are **byte-identical** (md5
`80f9c0e1…`, `920a9efb…`) to the copies installed in the run's conda prefix
`…/scE2G/.snakemake/conda/c1ae1ceebb37236cbea110033b350574_/lib/python3.11/site-packages/fast_kendall_sc/`.

Independent confirmation from the installed binary itself (`readelf -sW`, `strings`) — the symbol
table names the exact private items in `src/lib.rs` and the monomorphisations that pin the sorts:

- `__fast_kendall_sc::score_peak`, `…::GeneSortInfo::build`, `…::__batch_kendall_tau`,
  `…::__batch_kendall_tau_sparse`, `…::_kendall_concordance_diff`
- `core::slice::sort::stable::driftsort_main::<(usize,f64), …GeneSortInfo::build{closure}>`
  → the `explicit.sort_by` at `src/lib.rs:234` is a **stable comparison sort over nnz_gene
  elements**, not over `n_cells`
- `core::slice::sort::unstable::ipnsort::<(usize,usize), …sort_unstable_by_key… GeneSortInfo::build>`
  → the `by_cell_id.sort_unstable_by_key` at `src/lib.rs:275`
- `core::slice::sort::unstable::ipnsort::<usize, <usize as PartialOrd>::lt>` → the bare
  `ranks.sort_unstable()` at `src/lib.rs:319`
- `rayon::iter::plumbing::bridge_producer_consumer<rayon::range::IterProducer<usize>, …
  MapConsumer<CollectConsumer<Vec<Vec<f64>>>>, …__batch_kendall_tau_sparse{closure}>`
  → parallelism is a `(0..n_genes).into_par_iter().map(…).collect()` (`src/lib.rs:408-429`)
- `pyo3` `PyList_Append` / `PyFloat_FromDouble` present → the `PyResult<Vec<f64>>` return
  (`src/lib.rs:387`) marshals as a Python **list**, not an ndarray (empirically confirmed, §6.4)
- `Cargo.lock`/SBOM: no external sort or Fenwick-tree crate. Deps are pyo3, numpy, ndarray,
  rayon/crossbeam only.

So: the class below is derived from source, and the source is verified to be the source of the
binary that ran.

---

## 2. What the code actually does

### 2.1 R stage (`compute_kendall.R`)

| line(s) | operation | cost order |
|---|---|---|
| `:188` | `fread(Pairs.tsv.gz)` — 6 columns, 9.7–11.3 M rows. `PairName` is unique per row, so this interns ~10⁷ CHARSXPs into R's global string cache | Θ(`n_kendall_pairs`) |
| `:191` | `readRDS(atac_matrix.rds)` — peaks × cells `dgCMatrix` of raw fragment counts | Θ(nnz_atac) |
| `:192-193` | `BinarizeCounts()` — full second copy of the matrix | Θ(nnz_atac) |
| `:196-207` | RNA load. For this dataset `rna_matrix_file` is a **10x directory**, so `Read10X(…, gene.column=1)` at `:204` (61,217 features × `n_cells`, gene *symbols* as rownames) | Θ(nnz_rna) |
| `:209` | subset RNA columns to `colnames(matrix.atac)` — no-op here, cell sets already match | Θ(nnz_rna) |
| `:216` | `as(…, "CsparseMatrix")` — **no-op on this path** (`Read10X` already returns dgCMatrix). The 2.8 h → 1 s comment at `:211-215` refers to the h5ad/dgRMatrix branch, which this dataset does not use |
| `:219-223` | `sum()` / `ncol()` → `umi_count.txt`, `cell_count.txt`. **This is where the manifest's `n_umi` and `n_cells` are produced** | Θ(nnz_rna) |
| `:226` | `NormalizeData()` — log-normalisation, third RNA copy, sparsity pattern preserved | Θ(nnz_rna) |
| `:229-232` | `rowMeans`, `rowSums`, and `matrix.rna_count > 0` (an extra lgCMatrix temporary) | Θ(nnz_rna) |
| `:132-168` | `map_gene_names()` — `fread` the 53 MB gzipped GENCODE v43 GTF (3.4 M lines), filter to `type == "gene"`, then `unlist(lapply(gene_ref$attributes, extract_attributes, …))` **twice** over ~62 k rows (`:136-137`). Then a `left_join` against the TSS500 BED and a `group_by/filter(n()==1)` de-duplication | **fixed** (annotation-sized, cluster-independent) |
| `:246-250` | `fwrite(gene_expression_metrics.tsv.gz)` — 20,530 rows, tiny | Θ(`k_genes`) |
| `:48-49` | pair filter: `TargetGene %in% rownames(RNA) & PeakName %in% rownames(ATAC)`. **Measured: drops 0 rows** on both instrumented clusters (`pairs_in == pairs_scored`), consistent with `row_identities_ok = yes` | Θ(`n_kendall_pairs`) |
| `:60-70` | **the handoff write** — `writeMM` ×2, `writeLines` ×2, `fwrite` of the (gene, peak) pair list, all into `tempfile("kendall_")` | Θ(nnz_rna + nnz_atac + `n_kendall_pairs`) |
| `:261-262` | `rm(); gc()` — deliberate, so R and the Python child do not both hold the matrices (documented at `:30-39`, `:252-254`). It works: R drops to ~3.0 GB before `system2` |
| `:92-104` | `system2("python", …)` — fork/exec of a fresh CPython |
| `:110` | `as.numeric(readLines(results_path))` — 10⁷-element character vector, then `strtod`. The literal `"NA"` lines emitted by `compute_kendall.py:68` produce the `NAs introduced by coercion` warning seen in every job log. Measured NA fraction: **3.60 %** (thp1_1), **1.32 %** (telohaec_crispri) |
| `:273` | `pairs.E2G[, (gex_cols) := df.exp_filt[pairs.E2G$TargetGene, gex_cols]]` — `data.frame` row-indexing by a 10⁷-long **character** vector. Triggers the `Invalid .internal.selfref` shallow copy of the whole 10⁷-row table (visible in every job log) | Θ(`n_kendall_pairs`) |
| `:287-292` | column subset + `setorder(chr, start, end, TargetGene)` — data.table radix forder | Θ(`n_kendall_pairs`) |
| `:294-298` | `fwrite(*.tsv.gz)` — 10 columns × 10⁷ rows, single-threaded deflate | Θ(`n_kendall_pairs`) |

### 2.2 Python stage (`compute_kendall.py`)

`:36` sets `RAYON_NUM_THREADS` via `os.environ.setdefault` **before** importing the extension
(`:40`) — correct ordering, since rayon sizes its global pool once on first use. Note `setdefault`:
an inherited `RAYON_NUM_THREADS` would silently win over `--threads`.

`:42-43` `scipy.io.mmread(...).T.tocsc()` for both matrices. `mmread` returns COO ⇒
**`sp.issparse()` is True for both**, which selects the sparse Rust kernel at `metrics.py:144-147`.
This is the branch decision that fixes the complexity class.

`:55-59` a **pure-Python loop over ~10⁷ pair lines**, two dict lookups and two `list.append` per
line. `:61-62` converts to `int64` ndarrays. `:64` calls `batch_kendall_tau`. `:66-68` writes one
decimal per line with an f-string and a per-scalar `np.isnan`.

### 2.3 `fast_kendall_sc.metrics.batch_kendall_tau` (`metrics.py:90-206`)

- `:137` `np.argsort(gene_indices, kind="stable")` — full sort of `n_kendall_pairs` int64 keys to
  group pairs by gene (CSR-style).
- `:141` `np.unique(gene_sorted, return_counts=True)` — `np.unique` sorts unconditionally, so this
  is a **second full sort of an array that `:137-138` already sorted**.
- `:144` the sparse/dense branch. Sparse ⇒ `_run_sparse_batch` (`:168-206`).
- `:173-181` `csc_matrix(...)` + `eliminate_zeros()` for both; `:191` computes the ATAC tie term
  `n_y` from `np.diff(indptr)` (with an in-code note that `sum(axis=0)` was 5 orders of magnitude
  slower).
- `:198-204` `.astype(np.float64)` / `.astype(np.uintp)` on `data`, `indices`, `indptr` for both
  matrices. scipy stores `indices` as **int32**; the Rust signature needs `usize`, so every index
  array is duplicated at 8 bytes.
- `:163-164` scatters the flat result back into caller pair order.

### 2.4 The Rust kernel — the actual numerics

**Dispatch:** `_batch_kendall_tau_sparse` (`src/lib.rs:374-433`).

`:394-398` borrows the numpy buffers as slices — the full matrices are **not copied** into Rust
(the in-code rationale is at `:388-393`).

`:407-429` `py.allow_threads(|| (0..n_genes).into_par_iter().map(|g| …).collect())` — rayon,
parallel **over distinct genes**, one task per gene.

**Per gene** — `GeneSortInfo::build` (`:223-278`):
- `:230-232` gather the gene's explicit `(cell_id, value)` entries → `Vec` of length `nnz_gene`
- `:234` `explicit.sort_by(descending by value)` — **stable comparison sort, O(nnz_gene log nnz_gene)**
- `:240-262` one linear pass to record tie-group start/size per rank and accumulate `n_x`
- `:265-268` the implicit-zero block's `n_x` contribution in **closed form**, no scan
- `:270-275` build `by_cell_id` (cell_id → rank) and `sort_unstable_by_key(cell_id)` —
  **O(nnz_gene log nnz_gene)**
- `:280-288` `rank_of(cell)` = `binary_search_by_key` over `by_cell_id`; a miss means "implicit
  zero" and its rank is derived arithmetically from the insertion point — **O(log nnz_gene)**,
  and crucially **no `n_cells`-sized array is ever allocated** (documented `:189-207`)

**Per pair** — `score_peak` (`:312-360`), called once per (gene, peak) at `:419`:
- `:315-318` map the peak's `m_π` accessible cell ids to ranks → `m_π` binary searches
- `:319` `ranks.sort_unstable()` — **O(m_π log m_π)**
- `:325-330` one linear pass computing concordant/discordant from **rank gaps**:
  `discordant += r − k`, `concordant += k·(r − prev_rank − 1)` — the zero cells' aggregate
  contribution is captured by the gaps, so they are never visited
- `:336-357` the same accumulation restricted to each RNA tie-group; the inner `while` at `:343`
  advances monotonically and `i = j`, so the whole block is **O(m_π)** amortised
- `:420` denominator `sqrt((n_0 − n_x)(n_0 − n_y[peak]))`, NaN if zero

**Where `n_cells²` appears:** only as the scalar `n_0 = n_cells·(n_cells−1)/2` (`:405`) and in
`n_y` (`metrics.py:193`). Closed-form arithmetic, **O(1)**. There is no O(`n_cells`²) loop
anywhere in the crate.

**The dead sibling paths** (present in the binary, not taken here):
- `_batch_kendall_tau` (`:80-187`) — the *dense* kernel. `:116-121` allocates and sorts an
  `n_cells`-long index array per gene, and `:154-170` sweeps **all `n_cells` rows per pair** ⇒
  O(`n_kendall_pairs` · `n_cells`). Selected only if neither input is `scipy.sparse`.
- `_kendall_concordance_diff` (`:16-65`) — one gene per call, O(`n_cells`) per column; used by
  `kendall_tau_score` / `pairwise_kendall_tau`, which this pipeline never calls.

This is the one genuine code-level branch that changes the *class*. For this workflow it is
resolved (sparse), because `compute_kendall.R:60-61` writes MatrixMarket and
`compute_kendall.py:42-43` reads it back as COO. See §7.4.

---

## 3. Size variables

Manifest columns (`size_manifest_multiome.tsv`) that actually govern this rule:

| manifest variable | role in this rule | range (6 clusters) | ratio |
|---|---|---|---|
| `n_kendall_pairs` | rows of `Pairs.tsv.gz`; number of pairs scored (filter is a measured no-op); drives every R/Python glue pass | 9,694,978 – 11,317,834 | **1.17×** |
| `n_cells` | cells in both matrices; sets `n_0`, and (with per-cell depth) sets `q_γ` | 2,593 – 15,555 | **6.0×** |
| `n_frag` | best manifest proxy for nnz of the binarised ATAC matrix and hence for `m_π` | 26.2 M – 394.8 M | **15.1×** |
| `n_umi` | weaker proxy for nnz of the RNA matrix (1.5× spread in the ratio vs 1.38× for `n_cells`) | 14.5 M – 174.2 M | 12.0× |
| `k_genes` (20,532) | the gene universe the RNA matrix is collapsed onto; sets the *number* of rayon tasks | constant | 1.00× |
| `j_cand_elems` | ≈ the peak universe (unique `PeakName` measured 151,054 – 156,399) | 156,420 – 158,548 | 1.01× |

**The cross terms the code genuinely has.** The kernel's work is a sum over pairs of the peak's
nnz, so the honest cost variable is a *product*, not any single manifest column:

```
W_ATAC ≡ Σ_{(γ,π) ∈ pairs} m_π          (pair-weighted total accessible cells)
W_RNA  ≡ Σ_{γ ∈ distinct genes} q_γ     (≈ nnz of the gene-filtered RNA matrix)
```

I measured these directly for all six clusters (job 40217748, `$SCRATCH/kendall_probe/probe.R`,
independently reproduced by the instrumented replays):

| cluster | `n_cells` | `n_frag` | nnz_atac | nnz_rna(raw) | uniq peaks | uniq genes | mean `m_π` (pair-wt) | **W_ATAC** | mean `q_γ` | **W_RNA** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| thp1_1 | 2,676 | 26.23 M | 13.36 M | 7.157 M | 151,054 | 14,358 | 102.5 | 0.994e9 | 410.0 | 5.887 M |
| thp1_2 | 2,593 | 26.29 M | 13.09 M | 6.995 M | 151,555 | 14,347 | 100.1 | 0.992e9 | 401.3 | 5.757 M |
| jurkat_pma_cd3_4hr | 5,124 | 55.80 M | 21.35 M | 15.499 M | 155,224 | 14,354 | 156.1 | 1.680e9 | 942.3 | 13.525 M |
| k562_crispri | 5,109 | 86.42 M | 37.71 M | 18.380 M | 153,663 | 14,379 | 280.9 | 3.179e9 | 1064.4 | 15.305 M |
| jurkat | 8,239 | 75.34 M | 31.28 M | 25.921 M | 154,356 | 14,376 | 231.1 | 2.493e9 | 1587.3 | 22.820 M |
| telohaec_crispri | 15,555 | 394.85 M | 161.66 M | 57.211 M | 156,399 | 14,370 | 1113.5 | **12.401e9** | 3451.0 | 49.591 M |

**`W_ATAC` spans 12.5× — the widest-spanning cost quantity available to this unit**, wider than
`n_cells` (6.0×) and second only to `n_frag` (15.1×), with which it is nearly proportional.
This is what makes the kernel's exponent empirically testable here (§5.2) — unusual for this
project.

Empirical bridges to manifest variables (needed so Phase 2 can evaluate the functions from the
manifest alone):

| bridge | constant | measured spread |
|---|---|---|
| `m̄_π ≈ κ_m · n_frag` | κ_m = 3.3e-6 | 2.80e-6 – 3.91e-6 (1.40×) |
| `W_ATAC ≈ κ_m · n_kendall_pairs · n_frag` | 3.3e-6 | same |
| nnz_atac ≈ κ_A · `n_frag` | κ_A = 0.44 | 0.383 – 0.509 (1.33×) |
| nnz_rna(raw) ≈ κ_R · `n_cells` | κ_R = 3,070 | 2,674 – 3,678 (1.38×) |
| nnz_rna(gene-filtered) ≈ 0.90 · nnz_rna(raw) | 0.90 | 0.884, 0.916 (2 measured) |
| `q̄_γ ≈ 0.182 · n_cells` | 0.182 | ±20 % |
| distinct genes in pairs | 14,360 | 14,347 – 14,379 (**1.002×**) |

Note the last row: the number of rayon tasks is effectively a constant ~14,360, not a variable.

---

## 4. Asymptotic class derived from code

### 4.1 The Rust kernel (the numerics)

```
T_kernel = O(  Σ_γ  q_γ log q_γ                          [src/lib.rs:234, :275]
             + Σ_{(γ,π)} m_π · (log q_γ + log m_π)       [src/lib.rs:315-318, :319]
             + Σ_{(γ,π)} m_π                             [src/lib.rs:325-330, :339-357]  )
```

Justifying constructs, each cited:

| construct | file:line | cost |
|---|---|---|
| stable `sort_by` on the gene's explicit `(cell, value)` entries | `src/lib.rs:234` | O(q_γ log q_γ) per gene |
| `sort_unstable_by_key` on `by_cell_id` | `src/lib.rs:275` | O(q_γ log q_γ) per gene |
| tie-group scan | `src/lib.rs:244-262` | O(q_γ) per gene |
| implicit-zero tie term, closed form (no loop) | `src/lib.rs:265-268` | O(1) |
| `binary_search_by_key` in `rank_of`, once per accessible cell | `src/lib.rs:281`, called at `:317` | O(m_π log q_γ) per pair |
| `ranks.sort_unstable()` | `src/lib.rs:319` | O(m_π log m_π) per pair |
| rank-gap concordance pass | `src/lib.rs:325-330` | O(m_π) per pair |
| tie-group-local pass, monotone inner cursor | `src/lib.rs:339-357` | O(m_π) per pair |
| `n_0 = n_cells(n_cells−1)/2` | `src/lib.rs:405` | O(1), **not a loop** |
| rayon fan-out over genes | `src/lib.rs:408-409` | ~14,360 independent tasks |

In manifest variables, since q_γ, m_π ≤ `n_cells`:

```
T_kernel = O( W_RNA · log n_cells  +  W_ATAC · log n_cells )
         = O( n_kendall_pairs · n_frag · log n_cells )        [via W_ATAC ≈ 3.3e-6·P·n_frag]
```

`W_ATAC` dominates `W_RNA` by 169× (thp1_1) to 250× (telohaec_crispri), so the per-gene
`GeneSortInfo::build` term is numerically negligible even though it is the only place a
`log`-linear sort over the *gene* axis occurs.

**Explicitly, for the record:** this is *not* O(`n_cells` log `n_cells`) per pair (the classical
merge-sort τ) and *not* O(`n_cells`²) per pair (naive). It is O(m_π log m_π) with
m_π ≈ 0.028–0.072 · `n_cells`, i.e. a factor ~14–36 cheaper than the classical form at these
accessibility rates, and it is the sparsity of the ATAC column — not any property of `n_cells` —
that buys that.

### 4.2 The R + Python glue

```
T_glue = Θ( n_kendall_pairs · log n_kendall_pairs )   [two np sorts, one data.table forder]
       + Θ( nnz_rna )                                  [Read10X, NormalizeData, writeMM, mmread]
       + Θ( nnz_atac )                                 [readRDS, BinarizeCounts, writeMM, mmread]
       + Θ( 1 )                                        [R library load, GTF parse, py imports]
```

Superlinear terms, all with near-constant log factors (log₂ 10⁷ ≈ 23.3, spanning 1.007× across
the six clusters — unfittable, as the briefing predicts):
`np.argsort` (`metrics.py:137`), `np.unique`'s internal sort (`metrics.py:141`),
`setorder` (`compute_kendall.R:292`).

### 4.3 Which term dominates, measured

Instrumented replays (§5.1) give the split. Shares of wall time — a **selection of the largest
groups, not an exhaustive partition** (the ~4–5 % remainder is the other Θ(nnz_rna) stages:
`Read10X`, `NormalizeData`, `mmread(RNA)`, `gc`):

| | thp1_1 (min) | telohaec_crispri (max) |
|---|---:|---:|
| Rust kernel — O(W_ATAC · log) | 9.9 % | **43.6 %** |
| `writeMM(RNA)` — Θ(nnz_rna) | 12.3 % | **29.6 %** |
| Θ(`n_kendall_pairs`) glue (R + Python) | 38.8 % | 12.4 % |
| Fixed startup (R libs + GTF + py imports) | **30.0 %** | 4.7 % |
| remaining Θ(nnz_atac) | 3.9 % | 5.1 % |

The rule's *shape* changes qualitatively across the measured range: at thp1_1 it is a
startup-plus-row-count job; at telohaec_crispri it is a kernel-plus-serialisation job.

---

## 5. Benchmark cross-check

### 5.1 Measurements used

Real Phase 0 rows (`benchmarks_multiome.tsv`, rule `compute_kendall`):

| cluster | s | max_rss MB | cpu_time | cpu/wall | `seff` MaxRSS |
|---|---:|---:|---:|---:|---:|
| thp1_1 | 219.99 | 5,414 | 162.02 | 0.736 | 6.36 GB |
| thp1_2 | 219.35 | 5,388 | 153.50 | 0.700 | — |
| jurkat_pma_cd3_4hr | 278.03 | 5,939 | 235.71 | 0.848 | — |
| k562_crispri | 338.01 | 6,495 | 300.23 | 0.888 | — |
| jurkat | 449.01 | 6,484 | 409.62 | 0.912 | 8.11 GB |
| telohaec_crispri | 871.70 | 10,444 | 816.47 | 0.937 | 14.74 GB |

Because the six benchmark points alone cannot decompose this rule, I built two **instrumented
replays** — verbatim copies of `compute_kendall.R` / `.py` with only the snakemake inputs
literalised and per-stage timers added, run at the same `--cpus-per-task=1` / `--mem=63GB`
(jobs 40217957, 40217960; sources and logs in `$SCRATCH/kendall_probe/`). Fidelity:

| cluster | real `s` | replay total | Δ |
|---|---:|---:|---:|
| thp1_1 | 219.99 | 209.96 | −4.6 % |
| telohaec_crispri | 871.70 | **871.65** | **−0.006 %** |

Stage table (seconds):

| stage | thp1_1 | telohaec |
|---|---:|---:|
| R: `library()` (Signac/Seurat/anndata/tidyverse) | 27.62 | 16.35 |
| R: `fread(Pairs.tsv.gz)` | 17.41 | 23.56 |
| R: `readRDS(atac_matrix.rds)` | 3.57 | 10.47 |
| R: `BinarizeCounts` | 2.09 | 2.91 |
| R: `Read10X` | 3.86 | 21.79 |
| R: subset + `NormalizeData` + gene metrics | 0.79 | 5.56 |
| R: `map_gene_names` (GTF fread 16.93/11.57 + attrs 4.03/3.85 + key 3.12/1.34) | 24.08 | 16.76 |
| R: **`writeMM(RNA)`** | **25.80** | **258.09** |
| R: `writeMM(ATAC)` | 1.91 | 25.38 |
| R: `writeLines` + `fwrite(pairs.tsv)` | 0.33 | 0.39 |
| R: pair filter + `rm`/`gc` | 2.41 | 3.65 |
| Py: interpreter + numpy/scipy/ext import | 10.50 | 5.23 |
| Py: `mmread` ×2 + `.T.tocsc()` ×2 | 2.15 | 10.17 |
| Py: pair-line loop + `np.array` | 2.75 | 4.14 |
| Py: `argsort` + `np.unique` group-by-gene | 0.79 | 1.94 |
| Py: `csc_matrix`/`eliminate_zeros`/`n_y`/dtype casts | 0.17 | 0.75 |
| **Rust `_batch_kendall_tau_sparse`** | **20.89** | **379.84** |
| Py: PyList→ndarray + scatter | 0.43 | 0.52 |
| Py: write `kendall_results.txt` | 10.67 | 16.39 |
| subprocess spawn residual | 0.82 | 3.05 |
| R: `readLines` + `as.numeric` | 8.63 | 11.62 |
| R: gex join `df[<10⁷ char>, ]` | 11.24 | 10.56 |
| R: column subset + `setorder` | 4.12 | 1.28 |
| R: `fwrite(*.tsv.gz)` | 24.59 | 36.70 |
| R: `unlink(tmp_dir)` | 0.08 | 0.31 |

### 5.2 Does the observed kernel scaling agree with the code? **Yes, to 2.3 %.**

This is the rare case in this project where the data *can* discriminate, because `W_ATAC` spans
12.5× between the two instrumented clusters.

Observed kernel ratio telohaec / thp1_1 = 379.84 / 20.89 = **18.18×**.

| candidate class | predicted ratio | verdict |
|---|---:|---|
| O(`n_kendall_pairs` · `n_cells`) — the *dense* kernel (`src/lib.rs:154`) | 1.149 × 5.813 = **6.68×** | rejected (2.7× low) |
| O(`n_kendall_pairs` · `n_cells` · log `n_cells`) — classical merge-sort τ | 6.68 × 1.222 = **8.16×** | rejected (2.2× low) |
| O(W_ATAC), no log | **12.48×** | rejected (1.46× low) |
| **O(W_ATAC · (log₂ q̄ + log₂ m̄))** — `src/lib.rs:312-360` | 12.48 × 1.424 = **17.77×** | **accepted (−2.3 %)** |
| O(`n_kendall_pairs` · `n_cells`²) | 1.149 × 33.8 = 38.8× | rejected wildly |

Fitting the single constant: κ = t_kernel / (W_ATAC · L), L = log₂ q̄ + log₂ m̄:

- thp1_1: 20.89 / (0.994e9 × 15.359) = **1.369e-9 s**
- telohaec: 379.84 / (12.401e9 × 21.874) = **1.401e-9 s**

Two points 12.5× apart in `W_ATAC` agreeing on the constant to 2.3 % is strong confirmation of
the code-derived form, *including* the log factor. I still flag the obvious limitation: n = 2
instrumented points, and the four intermediate clusters' kernel times are interpolated from κ,
not measured.

### 5.3 Whole-rule cross-check, and the one disagreement

Calibrating each term's rate constant from the stage table (not by fitting the six wall times)
and then predicting all six:

| cluster | observed s | `t_est` | error |
|---|---:|---:|---:|
| thp1_1 | 220.0 | 204.0 | −7.3 % |
| thp1_2 | 219.4 | 205.1 | −6.5 % |
| jurkat_pma_cd3_4hr | 278.0 | 284.0 | +2.1 % |
| k562_crispri | 338.0 | 326.5 | −3.4 % |
| **jurkat** | **449.0** | **353.9** | **−21.2 %** |
| telohaec_crispri | 871.7 | 884.6 | +1.5 % |

**Disagreement, flagged not smoothed:** `jurkat` runs 21 % longer than the structural model, and
this is not a kernel-exponent effect — its `cpu_time` (409.6 s) alone exceeds the whole
prediction. Candidate mechanisms, none confirmable from one measurement:

1. **Cold-cache startup variance.** The two replays measured the "fixed" term at **63 s**
   (thp1_1, cold conda prefix on Lustre) vs **41 s** (telohaec, warm) — a 22 s / 1.5× spread on a
   term I model as constant. That covers ~5 of the 95 s gap.
2. **`writeMM` rate variance.** Measured 4.08 µs/nnz (thp1_1) vs 4.93 µs/nnz (telohaec), a 1.21×
   spread. At jurkat's 23.3 M RNA nonzeros, the same spread is ±20 s.
3. **Different scheduling epoch.** jurkat's job id is 40079572, ~160 k ids after the other five
   (39918185–39949395): it was re-run in a later batch, on a differently loaded node. Every
   Phase 0 row is a single measurement, so node heterogeneity is unmodelled.

None of these individually explains 95 s; together they plausibly do. I am not adjusting the
model to absorb it — I widened `t_high` instead (§7).

### 5.4 Where the benchmarks cannot help

`n_kendall_pairs` spans **1.17×** and `j_cand_elems` **1.01×**. The Θ(`n_kendall_pairs`) glue —
which is 38.8 % of the smallest job — is therefore *not* empirically separable from the fixed
term by these six points. Its rate constant (8.4–9.7 µs/pair) comes from the two instrumented
runs, not from the benchmark table. Likewise the O(P log P) sorts at `metrics.py:137,141` cannot
be distinguished from O(P): log₂ P spans 1.007× here.

---

## 6. Concrete overhead sources

### 6.1 The R→Python `system2` boundary — what is serialised, in what format, at what cost

This is the single biggest structural overhead in the rule. `prepare_kendall_inputs`
(`compute_kendall.R:40-81`) writes **five plain-text files** into `tempfile("kendall_")`;
`compute_kendall.py` reads all five and writes a sixth back. Measured exactly (`file.info`):

| file | written at | format | thp1_1 bytes | telohaec bytes | bytes/record |
|---|---|---|---:|---:|---:|
| `rna_matrix.mtx` | `.R:60` `writeMM` | ASCII MatrixMarket, `i j <double>` | 180,564,564 | 1,519,394,971 | 28.5–29.0 /nnz |
| `atac_matrix.mtx` | `.R:61` `writeMM` | ASCII MatrixMarket, `i j 1` | 143,747,476 | 1,865,520,107 | 10.8–11.5 /nnz |
| `rna_genes.txt` | `.R:62` | one name per line | 137,631 | 137,631 | — |
| `atac_peaks.txt` | `.R:63` | one name per line | 3,610,180 | 3,731,755 | — |
| `pairs.tsv` | `.R:64-70` `fwrite` | TSV, (TargetGene, PeakName), no header | 294,291,701 | 337,601,645 | 30.3 /pair |
| `kendall_results.txt` | `.py:66-68` | one decimal per line, or `NA` | 201,165,554 | 233,821,068 | 20.8–21.0 /pair |
| **total round-tripped** | | | **823.5 MB** | **3,960.2 MB** | |

Cost of that boundary, in seconds:

| | thp1_1 | telohaec |
|---|---:|---:|
| R writes (`writeMM` ×2 + `fwrite` + `writeLines`) | 28.04 | 283.86 |
| Python reads (`mmread` ×2 + `.T.tocsc()` ×2 + pair loop + `np.array`) | 4.90 | 14.31 |
| Python writes result text | 10.67 | 16.39 |
| R reads it back (`readLines` + `as.numeric`) | 8.63 | 11.62 |
| fork/exec + fresh CPython + numpy/scipy/ext import | 11.32 | 8.28 |
| **boundary total** | **63.6 (30.3 %)** | **334.5 (38.4 %)** |

Three specific observations:

- **The serialisation is grossly asymmetric.** R's `writeMM` formats the RNA matrix at
  **4.1–4.9 µs per nonzero**; scipy's `mmread` parses the *same file* at **0.065–0.23 µs per
  nonzero** — a **21–75× asymmetry**. Writing 1.52 GB of ASCII costs 258 s; reading it costs 3.4 s.
  The handoff cost is essentially *all* on the R side, and specifically in double→decimal
  formatting: the ATAC matrix (same code path, but values are all the literal `1`) writes at
  0.14–0.16 µs/nnz, i.e. **30× faster per nonzero** than the RNA matrix.
- **The ATAC values are carried end-to-end and never read.** `src/lib.rs:315` consumes only
  `atac_indices`; `metrics.py:180-193` documents that only column *membership* means "accessible".
  Yet the all-ones value vector is materialised five times: R `dgCMatrix@x` (1.29 GB of `1.0`
  doubles at telohaec), the `BinarizeCounts` copy (another 1.29 GB), the third ASCII field of
  1.87 GB of `.mtx` text, scipy's COO `data`, and scipy's CSC `data`.
- **Where it lands.** Every job log shows `resources: tmpdir=/tmp`, and R's `tempfile()` honours
  `TMPDIR=/tmp`. Verified on a compute node (job 40222013): `/tmp` is a per-user bind mount of
  **node-local NVMe** (`/dev/nvme1n1p1[/kaybrand]`, xfs, 3.5 TB), *not* tmpfs and *not* Lustre.
  So the 0.82–3.96 GB round trip is local-SSD I/O — cheap in latency terms, which is why the
  Python read side is so fast, and why the cost is concentrated in R's ASCII *formatting* rather
  than in the storage.
- **`disk_mb` is under-declared.** The job log shows `disk_mb=1000` / `disk=1 GB` while the actual
  temp footprint is up to 3.96 GB. Nothing enforces it on Sherlock (no `--tmp` is passed) and
  there is 3.5 TB available, so it is inert bookkeeping — but it is wrong by 4×.

### 6.2 Memory: what is resident, in what layout, is it sparse

**It is sparse throughout. No dense `n_cells × n_genes` or `n_cells × n_peaks` array is ever
materialised**, on either side of the boundary. Confirmed in code:
`src/lib.rs:189-207` (GeneSortInfo is sized to `nnz_gene`, never `n_cells`);
`src/lib.rs:394-398` (the Rust side borrows the numpy CSC buffers rather than copying);
`metrics.py:100-110`. The only dense allocations in the crate are in the *unused* dense kernel
(`src/lib.rs:116` `order`, `:126` `new_group`, both `n_cells`-long per gene).

Resident set, by owner, at telohaec_crispri (peak, from R's own `gc()` plus the layout arithmetic):

| item | layout | bytes |
|---|---|---:|
| `pairs.E2G` data.table, 11.1 M × 6 | 4 character + 2 integer columns → ~0.4 GB of pointers | ~0.4 GB |
| `PairName` string pool | 11.1 M **unique** CHARSXPs interned in R's global string cache; never released | ~1.1 GB |
| `matrix.atac_count` | dgCMatrix, 161.7 M nnz: `i` int32 0.65 GB + `x` float64 1.29 GB | 1.94 GB |
| `matrix.atac` (`BinarizeCounts` copy, `.R:192`) | second full dgCMatrix | 1.94 GB |
| `matrix.rna_count` | 57.2 M nnz → 0.69 GB | 0.69 GB |
| `matrix.rna` (`NormalizeData` copy) | 0.69 GB | 0.69 GB |
| `matrix.rna_count > 0` lgCMatrix temporary (`.R:230`) | 0.46 GB | 0.46 GB |
| `matrix.rna_filt` (52.4 M nnz) | 0.63 GB | 0.63 GB |
| **R total** | | **≈ 8.4 GB** |

R's own accounting agrees: peak `Ncells` 1,905.7 MB + `Vcells` 7,084.7 MB ≈ **9.0 GB**, against
`max_rss` 10,444 MB. After the deliberate `rm(); gc()` at `.R:261-262` R falls to 958.6 +
2,031.0 = **3.0 GB** for the duration of the subprocess — the mitigation documented at
`.R:30-39` demonstrably works.

Python side, peak (telohaec): the ASCII→COO→CSC→`uintp`-cast chain, not the algorithm —
ATAC COO (161.7 M × 16 B = 2.6 GB) coexisting with its CSC form during `.tocsc()`, then
`indices.astype(np.uintp)` doubling 0.65 GB → 1.29 GB; RNA equivalently ~1.1 GB; the two Python
`list`s of ~11.1 M boxed ints from `compute_kendall.py:53-59` (~0.8 GB, *still referenced* through
the kernel call); and the result marshalling (§6.4, ~0.5 GB). ≈ 6.5 GB, consistent with
10.4 GB total minus R's retained 3.0 GB.

**The handoff bytes are also charged to the memory cgroup.** `seff` MaxRSS exceeds Snakemake's
psutil `max_rss` by 0.95 GB (thp1_1), 1.63 GB (jurkat) and 4.30 GB (telohaec) — against measured
/estimated temp-file volumes of 0.82 GB, ~1.59 GB and 3.96 GB. Three-for-three within ~15 %.
`/tmp` is node-local NVMe, not tmpfs (verified, §6.1), so the mechanism is **page cache for the
handoff files being accounted to the job's memory cgroup** rather than anonymous tmpfs pages.
Either way the practical consequence is the same: **the ASCII handoff consumes `mem_mb` on top of
RSS**, and it does so in proportion to nnz_rna + nnz_atac — i.e. to `n_cells` and `n_frag`, not
to `n_kendall_pairs`.

For completeness on the declared floor: `determine_mem_mb(min_gb=63)`
(`ENCODE_rE2G/ABC/workflow/rules/utils.smk:17-24`) computes `max(4 × 8 × input_size_mb, 63000)`;
even at telohaec's ~750 MB of declared inputs that is 24,000 MB, so **the 63 GB floor binds on all
six clusters**. Observed cgroup peak 6.4–14.7 GB ⇒ 10–24 % efficiency. (Not re-litigated; the
briefing already establishes the floor is stale. Recorded because the *shape* — RSS plus handoff
bytes — is what a future sizing would need.)

### 6.3 Declared vs realised threading

- `compute_kendall.smk:41-42` requests `threads: config["threads"]`;
  `configs/igvf10_multiome_config.yaml:70` sets it to **1** (with an in-repo comment at `:68-69`
  that it is left at scE2G's default deliberately for parity).
- Path: `snakemake@threads` → `cores` (`.R:183`) → `--threads 1` (`.R:102`) →
  `os.environ.setdefault("RAYON_NUM_THREADS", "1")` (`.py:36`, before the extension import at
  `:40`) → rayon global pool = 1 worker.
- **Realised parallelism at threads=1 is exactly 1.** Measured cpu_time/wall is 0.70–0.94 (mean
  0.83); `/usr/bin/time` on the replays gives 83 % (thp1_1) and 97 % (telohaec). The shortfall
  below 1.0 is **I/O wait, not missing parallelism** — the small clusters wait proportionally
  longer on the Lustre reads and `/tmp` writes because their compute term is smaller. The trend
  is monotone in cluster size (0.700 → 0.937), which is the signature of a fixed I/O cost against
  growing compute.
- What the code *would* do with more threads: `src/lib.rs:408-409` (sparse) and `:105-106`
  (dense) fan out over **distinct genes** — measured 14,347–14,379 tasks with 675–775 pairs each,
  so there is no task-count ceiling below any Sherlock node width. Per-gene work is skewed
  (pairs/gene and `m_π` both vary widely), but rayon's adaptive range splitting plus work stealing
  is the standard remedy. `py.allow_threads` (`:104`, `:407`) releases the GIL for the whole
  parallel region, so nothing serialises on the interpreter.
- **But only the kernel is inside the parallel region.** Everything else — all of the R script,
  and the Python import / `mmread` / pair-loop / result-write stages — is strictly serial. The
  kernel's measured share of wall is 9.9 % (thp1_1) to 43.6 % (telohaec_crispri), so Amdahl caps
  the achievable whole-rule speedup at **1.11× (thp1_1) to 1.77× (telohaec_crispri)** no matter
  how many cores are given. Stated as characterisation; **no recommendation is made** — Phase 2
  owns that call.
- Incidental: both replays reported `cpus=2` from a `--cpus-per-task=1` request (Sherlock
  hyperthread pairing), so `data.table`/OpenMP may see 2. Rayon does not, being pinned by the
  env var.

### 6.4 Redundant passes over the same data

Seven representations of the same P-length τ vector, between the Rust `f64` and the output file:

1. `Vec<Vec<f64>>` — one `Vec` per gene, ~14,360 allocations (`src/lib.rs:407`)
2. `.concat()` — a second full flat copy (`src/lib.rs:432`)
3. `PyResult<Vec<f64>>` marshals through pyo3 as a **Python `list` of P `float` objects**
   (`src/lib.rs:387`). Verified empirically: `type(raw) == list`, `len == 11,136,922`. At CPython
   3.11 that is ~32 B/element ⇒ **~356 MB** for telohaec, transient but coexisting with (1) and (2).
4. `np.asarray(...)` back to a contiguous float64 array (`metrics.py:145`)
5. `results[order] = flat_results` — an O(P) random-access scatter (`metrics.py:164`)
6. `kendall_results.txt` — P lines of decimal text (`compute_kendall.py:66-68`), 201–234 MB
7. R `character` vector of P strings, then `as.numeric` (`compute_kendall.R:110`)

Other repeated work:

- `np.argsort(gene_indices, kind="stable")` at `metrics.py:137` sorts P keys; `np.unique(...)` at
  `:141` then sorts **the already-sorted result again** (numpy's `_unique1d` sorts
  unconditionally). Two full P-length sorts where the second's input is known-sorted.
  Measured jointly: 0.79 s / 1.94 s — small absolutely, listed for completeness.
- `sp.csc_matrix(rna_matrix, dtype=np.float64)` at `metrics.py:173` on a matrix that is *already*
  CSC float64, followed by `.data.astype(np.float64)` at `:198` — two copy opportunities on the
  same array.
- `.astype(np.uintp)` on `indices`/`indptr` for both matrices (`metrics.py:199-204`): scipy stores
  int32, the Rust ABI needs `usize`, so every index array is duplicated at 8 B. +1.29 GB (ATAC)
  and +0.42 GB (RNA) at telohaec.
- `BinarizeCounts` (`compute_kendall.R:192`) duplicates the ATAC matrix; `NormalizeData` (`:226`)
  duplicates the RNA matrix; `matrix.rna_count > 0` (`:230`) makes a third RNA-shaped temporary.
- `pairs.E2G[, (gex_cols) := …]` (`.R:273`) triggers data.table's `Invalid .internal.selfref`
  shallow copy of the entire 10⁷-row table — visible as a warning in all six job logs.
- The pairs list traverses **eight** representations: `Pairs.tsv.gz` → R data.table → filtered
  copy → `pairs.tsv` ASCII → Python `str` lines → two Python `list[int]` → two `int64` ndarrays →
  two sorted `uintp` copies.

### 6.5 Fixed startup

| component | thp1_1 | telohaec |
|---|---:|---:|
| R `library(Signac, Seurat, anndata, data.table, Matrix, dplyr, tibble, tools)` | 27.62 | 16.35 |
| `map_gene_names` GTF work (53 MB gz, 3.4 M lines, two `lapply` passes over ~62 k attribute strings) | 24.08 | 16.76 |
| CPython + numpy + scipy + `_fast_kendall_sc` import | 10.50 | 5.23 |
| fork/exec residual | 0.82 | 3.05 |
| **total** | **63.0 (30.0 %)** | **41.4 (4.7 %)** |

Wholly independent of cluster size; the 1.5× spread between the two runs is filesystem cache
state on the conda prefix (Lustre). Note this is *inside* the benchmarked `s` — it is the script's
own startup, not Slurm latency. Slurm/Snakemake overhead is additional and outside `s`: the job
logs show ~3 s between the outer `rule compute_kendall` and the inner `localrule compute_kendall`
(the SLURM executor re-invokes Snakemake inside the allocation, so each job pays a second DAG
build and a conda activation).

### 6.6 On `io_in` / `io_out`

Per instruction I did not rely on them; I used `file.info()` on the actual temp files and on-disk
sizes. For the record, they happen to be *consistent* for this rule: `io_out` 3,220 MB (thp1_1)
and 6,429 MB (telohaec) against `ru_oublock`-derived 3.53 GB and 6.83 GB from my replays; `io_in`
133 MB and 705 MB against 162 MB and 752 MB. Both within ~10 %. I mention this only so the next
reader does not mistake agreement for a contradiction of the project-wide verdict — this rule's
I/O is large, uncached and block-level, which is exactly the regime where those columns behave.

---

## 7. `t_est`, `t_low`, `t_high`

### 7.1 Form B — evaluable from `size_manifest_multiome.tsv` alone (**the deliverable**)

Variables, all exactly as named in the manifest: `n_kendall_pairs`, `n_cells`, `n_frag`.
Seconds. `log2` is base-2 logarithm.

```
LB = log2(0.182 * n_cells) + log2(3.3e-6 * n_frag)

t_est(n_kendall_pairs, n_cells, n_frag) =
      55.0
    + 9.06e-6  * n_kendall_pairs
    + 1.463e-2 * n_cells
    + 1.54e-7  * n_frag
    + 4.59e-15 * n_kendall_pairs * n_frag * LB

t_low(n_kendall_pairs, n_cells, n_frag) =
      38.0
    + 7.00e-6  * n_kendall_pairs
    + 1.030e-2 * n_cells
    + 1.10e-7  * n_frag
    + 4.13e-15 * n_kendall_pairs * n_frag * LB

t_high(n_kendall_pairs, n_cells, n_frag) =
      80.0
    + 14.0e-6  * n_kendall_pairs
    + 2.260e-2 * n_cells
    + 2.42e-7  * n_frag
    + 6.11e-15 * n_kendall_pairs * n_frag * LB
```

Term meanings: constant = R library load + GTF parse + Python imports (§6.5);
`n_kendall_pairs` term = all Θ(P) glue passes in R and Python (§4.2);
`n_cells` term = Θ(nnz_rna) work, dominated by `writeMM(RNA)` (via nnz_rna ≈ 2,760·`n_cells`);
`n_frag` term = Θ(nnz_atac) work (via nnz_atac ≈ 0.44·`n_frag`);
product term = the Rust kernel, `O(W_ATAC · log)` (via `W_ATAC` ≈ 3.3e-6·P·`n_frag`).

**`t_low` and `t_high` differ from `t_est` only in their constants — the exponents and the cross
term are identical.** That is deliberate and follows the instructions: the code leaves the class
unambiguous (§4.1, single sparse code path, every loop identified and cited), so the uncertainty
is in the numeric constants, not in the exponent. The spread is set from measured variation, not
invented: startup 41–63 s (1.5×), `writeMM` 4.08–4.93 µs/nnz (1.21×), glue 8.4–9.7 µs/pair
(1.16×), kernel κ 1.369–1.401e-9 (1.02×), plus enough headroom on the upper side to bracket the
unexplained +21 % `jurkat` residual (§5.3). `t_low` is below and `t_high` above **all six**
observations.

### 7.2 Form A — tighter, if Phase 2 can supply the measured intermediates

Uses `W_ATAC`, nnz_rna(filtered), nnz_atac, `q̄`, `m̄` from the §3 table (I measured all six
clusters, so they are available today; they are *not* manifest columns):

```
L = log2(q_bar) + log2(m_bar)

t_est = 55.0
      + 9.06e-6  * n_kendall_pairs
      + 5.3e-6   * nnz_rna_filtered
      + 0.35e-6  * nnz_atac
      + 1.39e-9  * W_ATAC * L
```

Residuals: −7.8 %, −7.1 %, −1.5 %, +0.2 %, **−21.7 % (jurkat)**, −1.1 %. Mean |error| 6.6 %
excluding jurkat: 3.5 %.

### 7.3 Assumptions a future dataset could break

- `distinct genes in pairs ≈ 14,360` and `distinct peaks ≈ 154,000` are treated as constants.
  Both are set by `k_genes` and the `nStrongestPeaks` cap and are constant here to 1.002× and
  1.035×. A different gene annotation or peak cap invalidates the `q̄`/`m̄` bridges.
- The bridges `m̄ ≈ 3.3e-6·n_frag` and nnz_rna ≈ 3,070·`n_cells` are per-cell-depth relations.
  They hold to ±20 % across these six clusters but are properties of the assay, not the code.
- `max_cell_count: 20000` (`configs/igvf10_multiome_config.yaml:63`) is not binding here (max
  15,555). Above it, `generate_atac_matrix.R` subsamples, and `n_cells` — and hence `m̄`, `q̄` and
  the whole kernel term — **saturates**. All three functions become optimistic ceilings past that
  point, in the sense that they would keep growing while the real cost plateaus.

### 7.4 The dense-path contingency (not used in `t_est`)

If the RNA or ATAC matrix ever reached `batch_kendall_tau` as a dense array, `metrics.py:144`
would select `_batch_kendall_tau` and the kernel term would become

```
t_kernel_dense ≈ c_d * n_kendall_pairs * n_cells   +   c_s * 14360 * n_cells * log2(n_cells)
```

i.e. O(P·`n_cells`) from the full sweep at `src/lib.rs:154-170`. At telohaec that inner work is
11.14e6 × 15,555 = 1.73e11 element-visits versus the sparse path's 1.24e10 — **14× more**. I did
not calibrate `c_d` (that path never runs). Recorded because it is the only place the *class*
is genuinely branch-dependent, and the branch is one line.

---

## 8. All three functions evaluated

Form B, seconds. "obs" is the Phase 0 `s` column.

| point | `n_kendall_pairs` | `n_cells` | `n_frag` | `t_low` | **`t_est`** | `t_high` | obs |
|---|---:|---:|---:|---:|---:|---:|---:|
| **min — thp1_1** | 9,694,978 | 2,676 | 26,231,026 | 152.4 | **204.0** | 306.4 | 220.0 |
| **mean of 6 clusters** | 10,601,259 | 6,549 | 110,823,580 | 282.8 | **365.0** | 537.7 | 396.0 (mean) |
| **max — telohaec_crispri** | 11,136,922 | 15,555 | 394,847,593 | 715.8 | **884.6** | 1,269.1 | 871.7 |

Term decomposition of `t_est` (seconds), which is where the interesting story is:

| point | const | Θ(P) glue | Θ(nnz_rna) | Θ(nnz_atac) | **kernel** |
|---|---:|---:|---:|---:|---:|
| thp1_1 | 55.0 | 87.8 | 39.2 | 4.0 | **18.0** |
| mean | 55.0 | 96.0 | 95.8 | 17.1 | **100.8** |
| telohaec_crispri | 55.0 | 100.9 | 227.6 | 60.8 | **439.9** |

The kernel goes from 8.8 % to 49.7 % of `t_est` across the measured range while the Θ(P) glue term
barely moves (87.8 → 100.9 s, because `n_kendall_pairs` spans only 1.17×). Any extrapolation of
this rule is an extrapolation in `n_frag` and `n_cells`, essentially not at all in
`n_kendall_pairs`.

For reference, the other three clusters: `t_est` = 205.1 (thp1_2), 284.0 (jurkat_pma_cd3_4hr),
326.5 (k562_crispri), 353.9 (jurkat).

---

## 9. Modality

**The rule does not exist on the scATAC path.** It is one of the four Multiome-only rules
(`make_kendall_pairs → generate_atac_matrix → compute_kendall → arc_e2g`), gated by
`checkpoint features_required` resolving to `to_generate == "ARC"`. Confirmed: 0 rows in
`benchmarks_scatac.tsv`, and `n_kendall_pairs` / `n_cells` / `n_umi` / `kendall_scored_bytes` are
all empty in `size_manifest_scatac.tsv`. There is therefore no cross-modality class or constant
comparison to make. Six jobs total, all multiome.

---

## 10. What the measurement cannot tell us

1. **The four un-instrumented clusters' internal split.** I timed stages only for `thp1_1` and
   `telohaec_crispri`. jurkat_pma_cd3_4hr, k562_crispri, jurkat and thp1_2 have their stage
   shares *inferred* from the two calibrated rate constants. In particular the `jurkat` +21 %
   residual (§5.3) is **unattributed**: I cannot say whether it was startup, `writeMM`, node
   speed, or something I have not modelled, and one measurement cannot distinguish them.
2. **The kernel constant rests on n = 2.** κ = 1.369e-9 vs 1.401e-9 is a beautiful agreement
   across a 12.5× range in `W_ATAC`, but it is two points. A third instrumented cluster (k562 is
   the informative one: high `m̄` at low `n_cells`) would be the cheapest way to break the
   `n_frag`/`n_cells` collinearity that the briefing warns about.
3. **`n_frag`, `n_umi` and `n_cells` remain mutually confounded.** My model attributes the
   Θ(nnz_rna) term to `n_cells` and the Θ(nnz_atac) and kernel terms to `n_frag`, on the basis of
   which measured intermediate each *stage* consumes (from the instrumented runs) — not on the
   basis of the six wall times, which cannot separate them. If a future dataset has an unusual
   depth-per-cell, the attribution, and hence the extrapolation, will be wrong even though the
   class is right.
4. **`n_kendall_pairs` spans 1.17× and `j_cand_elems` 1.01×.** The Θ(`n_kendall_pairs`) glue term
   — 38.8 % of the smallest job — is not empirically separable from the constant term by these
   six points, and the O(P log P) sorts at `metrics.py:137,141` are indistinguishable from O(P)
   (log₂ P spans 1.007×). Those coefficients come from the stage table only.
5. **The claim that `m_π` scales with `n_frag` is empirical, not structural.** The *code's* cost
   is O(nnz of the ATAC column); that nnz happens to track sequencing depth in these six
   datasets (ratio 2.80e-6 – 3.91e-6, a 1.40× spread). Nothing in the code guarantees it. A
   deeply-sequenced but low-complexity library, or a change to `BinarizeCounts` thresholds, would
   move `m̄` independently of `n_frag`.
6. **Behaviour above `threads: 1` is entirely unmeasured.** I characterised what the code *would*
   do (rayon over ~14,360 gene tasks, GIL released, Amdahl ceiling 1.11×–1.77×), all read from
   `src/lib.rs:407-429` and the measured serial/parallel split. No multi-thread run exists in this
   dataset, so scaling efficiency, per-gene load imbalance and any memory multiplication from
   `T` concurrent `GeneSortInfo` allocations are predictions, not observations.
7. **The cgroup-memory attribution of the handoff is inferred, not directly instrumented.** The
   arithmetic is 3-for-3 (seff − psutil `max_rss` ≈ temp-file volume) and `/tmp` is confirmed
   node-local NVMe rather than tmpfs, so page-cache accounting is the only mechanism left standing
   — but I did not read `memory.stat` inside a live `compute_kendall` job to prove it. If Phase 2
   ever resizes `mem_mb` on the strength of RSS alone, this is the term that would bite.
8. **The `max_cell_count: 20000` saturation point is untested.** No cluster reaches it, so the
   kernel term's growth has never been observed to plateau.
9. **Single measurements, no replicates.** Unlike the shared ABC/fragment stack, this rule exists
   only in the multiome run, so every one of the six points is n = 1. Node heterogeneity, page-cache
   state and Lustre contention are unmodelled and — from the 1.5× spread I measured on the
   supposedly-constant startup term alone — are worth ±10 % on any single number here.

---

## Appendix: artifacts produced (all under `$SCRATCH`, nothing in the repo)

| path | what |
|---|---|
| `/scratch/users/kaybrand/fks_src/fast_kendall_sc-1.0.0/` | the verified Rust source (`src/lib.rs`) |
| `/scratch/users/kaybrand/kendall_probe/probe.R`, `probe.40217748.log` | per-cluster shape probe: nnz_atac, per-peak nnz distribution, `W_ATAC`, `W_RNA` |
| `/scratch/users/kaybrand/kendall_probe/probe_kendall.{R,py}` | instrumented verbatim replays |
| `/scratch/users/kaybrand/kendall_probe/replay_thp1_1.40217957.log`, `replay_telohaec_crispri.40217960.log` | stage timings, handoff byte counts, `PYTYPE`/`PYWORK` lines |
| `/scratch/users/kaybrand/kendall_probe/model2.{awk,dat}` | the six-cluster validation of §5.3/§7 |

Per-peak `m_π` distribution, for anyone modelling the kernel's tail (accessible cells per peak):

| cluster | min | q25 | median | q75 | q90 | q99 | max |
|---|---:|---:|---:|---:|---:|---:|---:|
| thp1_1 | 3 | 21 | 39 | 98 | 268 | 828 | 2,518 |
| telohaec_crispri | 45 | 226 | 400 | 1,085 | 3,087 | 8,244 | 12,883 |

The distribution is heavily right-skewed (q99/median ≈ 21), which is why per-gene rayon tasks are
load-imbalanced and why `m̄` (pair-weighted) sits ~1.14× above the unweighted per-peak mean in
every cluster.

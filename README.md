# scE2G

A computational pipeline for predicting genome-wide enhancer-gene regulatory links from single-cell ATAC-seq or paired ATAC and RNA-seq (multiome) data.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2024.11.23.624931v1-red.svg)](https://www.biorxiv.org/content/10.1101/2024.11.23.624931v1)
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/EngreitzLab/scE2G/tree/main.svg?style=shield)](https://dl.circleci.com/status-badge/redirect/gh/EngreitzLab/scE2G/tree/main)

<hr>

## Overview

**Input:** Single-cell ATAC-seq or paired ATAC and RNA-seq (multiome) data per cell cluster

**Output:** Genome-wide enhancer-gene regulatory link predictions per cell cluster

### Pipeline components

1. **ABC model predictions** - Compute ABC model predictions for each cell cluster
2. **E2G feature generation** - Generate element-gene features from ABC predictions  
3. **Correlation analysis** *(multiome only)* - Compute Kendall correlation and/or ARC-E2G score for each cell cluster
4. **Feature integration** - Combine components 2 & 3 to construct feature file for predictive model
5. **Model training** *(optional)* - Train predictive model using CRISPR-validated E-G pairs from K562 dataset
6. **Prediction** - Apply trained model to assign scores to each element-gene pair

<hr>

## System requirements

### Hardware requirements

The scE2G pipeline necessitates a standard computer furnished with ample RAM to facilitate the operations as defined by a user. 
For optimal performance, we suggest using a computer equipped with 32+ GB RAM.

### Software requirements

The scE2G pipeline is compatible with Linux. It has been tested successfully on the following systems:
- Linux: Red Hat Enterprise Linux 8.10 (Ootpa)
- Linux: CentOS Linux 7 (Core)

The software dependencies and versions on which the software has been tested are listed in 
- `workflow/envs/run_snakemake.yml` 
- `workflow/envs/sc_e2g.yml files`

<hr>

## Installation

### Clone repository

```bash
# Quick clone with submodules
git clone --recurse-submodules --shallow-submodules --depth 1 https://github.com/EngreitzLab/scE2G.git

# Initialize and update nested submodules
cd scE2G
git submodule update --init --recursive
```

### Set up environment

We highly recommend using the provided conda environment for compatibility:

```bash
# Configure conda for flexible channel packaging
conda config --set channel_priority flexible

# Create and activate environment
conda env create -f workflow/envs/run_snakemake.yml
conda activate run_snakemake
```

<hr>

## Usage

### Prerequisites

Before running scE2G, perform clustering to define cell clusters through standard single-cell analysis (e.g., Seurat & Signac).

### Input data requirements

See example data in `resources/example_chr22_multiome_cluster/` folder.

#### 1. Pseudobulk fragment files

- **Format:** One file per cell cluster with corresponding `*.tbi` index files
- **Requirements:**
  - Sorted by coordinates, compressed using `bgzip`, and indexed with `tabix`
  - 5 columns (no header) corresponding to: `chr`, `start`, `end`, `cell_name`, `read_count`
  - Cell names must match RNA count matrix column names
  - Must contain exact set of cells represented in RNA matrix
 
**Preparation:**
```bash
# Sort if needed
sort -k1,1 -k2,2n atac_fragments.unsorted.tsv > atac_fragments.tsv

# Compress and index
bgzip atac_fragments.tsv
tabix -p bed atac_fragments.tsv.gz
```

#### 2. RNA count matrix *(multiome only)*

- **Format:** Gene × cell matrix for each cell cluster
- **Requirements:**
  - Use unnormalized (raw) counts
  - No duplicated gene names
  - Supported formats:
    - `.csv.gz`
    - `.h5ad` or `.h5` (may require matching `anndata` version)
    - Sparse matrix directory (`matrix.mtx.gz`, `genes.tsv`/`features.tsv.gz`, `barcodes.tsv.gz`)

**Gene mapping:** By default, genes are mapped via Ensembl ID using GENCODE v43 annotations. Modify `gene_annotation` in `config/config.yaml` for different versions (e.g., GENCODE v32 for CellRanger data).

<details>
<summary><h4>Advanced file format options</h4></summary>

**Pre-processed fragment files**

If your fragment files are already properly sorted and filtered to main chromosomes, you can skip preprocessing steps by setting `fragments_preprocessed: True` in your config file. Use this option only if you are certain that your files:
1. Are sorted with `sort -k1,1 -k2,2n` (chromosome then numerical position)
2. Only contain fragments on chromosomes present in the chromosome sizes reference file

This option allows you to skip the sorting and filtering steps of the pipeline, which can be very resource intensive for large fragment files.

**Cell filtering configurations**

The default pipeline settings assume the each cell cluster has a corresponding fragment file and RNA matrix that contain the exact same cells. If you instead have an RNA matrix containing cells from many clusters, you can avoid making cluster-specific matrices by setting `RNA_matrix_filtered: False` in your config file, and using the same RNA matrix for all clusters in your `cell_cluster_config`. The pipeline will then use the *intersection of cells contained in the cluster-specific fragment file and combined RNA matrix* to compute the Kendall correlation.
Please note:
1. You still must provide cluser-specific fragment files
2. The RNA matrix must meet the formatting requirements indicated above
3. The memory requirements to load a very large RNA matrix may exceed the default estimations in the pipeline.

</details>

### Configuration (for generating predictions)

1. **Main config** - Edit `config/config.yaml`:
   - Set `results_dir` path
   
2. **Cell clusters** - Edit `config/config_cell_clusters.tsv`:
   - Specify cluster name, ATAC fragment file path, RNA matrix path
   - For ATAC-only analysis: leave RNA matrix path empty but include the column
   
3. **Model selection** - Specify path to model directory, or a comma-separated list of multiple models. Current supported models estimate contact using a power law function of genomic distance:
   - For multiome predictions: `models/multiome_powerlaw_v3`
   - For ATAC-only predictions: `models/scATAC_powerlaw_v3`

### Running the pipeline

#### Option 1: Using Singularity container (recommended)

This approach uses a pre-built container from SyLabs Cloud, so you only need to create the small `run_snakemake` environment locally.

1. Create and activate the environment:
```bash
conda env create -f workflow/envs/run_snakemake.yml
conda activate run_snakemake
```

2. Run the pipeline:
```bash
snakemake --configfile config/config.yaml -j1 --use-conda --use-singularity
```

#### Option 2: Using conda only

This approach builds all conda environments locally, which can take significant time on first run.

```bash
snakemake -j1 --use-conda --configfile config/config.yaml
```

> **Note:** First run may take time to build conda environments, usually around 30-40 minutes according to CircleCI tests. If it exceeds 1 hour, ensure you're using mamba and have sufficient memory.

### Output

#### Key outputs
- **Tabular predictions**
  - `{results_dir}/{cell_cluster}/{model_name}/scE2G_predictions.tsv.gz`: All putative enhancer-gene predictions for a cell cluster
  - `{results_dir}/{cell_cluster}/{model_name}/scE2G_predictions_threshold{model_threshold}.tsv.gz`: Thresholded predictions containing enhancer-gene pairs that pass score threshold and other filtering steps
  - Key score column in these files is `E2G.Score.qnorm` (quantile-normalized scE2G score)
- **Genome-browser files** (produced if `make_IGV_tracks: True` in your config file)
  - `{IGV_dir}/{cell_cluster}/ATAC_norm.bw`: bigWig file with read-depth normalized pseudobulk ATAC signal 
  - `{IGV_dir}/{cell_cluster}/scE2G_predictions_threshold{model_threshold}.bedpe`: bedpe file with filtered enhnancer-gene predictions
- **QC report**
  - `{results_dir}/qc_plots/predictions_qc_report.html`: Report summarizing properties of predictions in comparison to reference values

<details>
<summary><h4>Full output structure</h4></summary>
  
```
{results_dir}/                                         # Main results directory
├── {cell_cluster}/                                      # Outputs for each cell cluster
│   ├── ActivityOnly_features.tsv.gz                       # All element-gene pairs with activity-based features
│   ├── ActivityOnly_plus_external_features.tsv.gz         # All element-gene pairs with activity-based and other features
│   ├── ARC/                                               # ARC-E2G results and intermediate files (multiome only)
│   ├── external_features_config.tsv                       # Configuration for external features
│   ├── feature_table.tsv                                  # Feature table reflecting all models designated for the cell cluster
│   ├── genomewide_features.tsv.gz                         # Genome-wide element-gene pairs with features, formatted for model application
│   ├── Kendall/                                           # Kendall correlation results and intermediate results (multiome only)
│   ├── {model_name}/                                      # scE2G model-specific results (e.g., multiome_powerlaw_v3)
│   │   ├── scE2G_predictions.tsv.gz                         # All predictions with scores
│   │   ├── scE2G_predictions_threshold{threshold}.tsv.gz    # Filtered predictions
│   │   ├── scE2G_predictions_threshold{threshold}_stats.tsv # Properties of filtered predictions
│   │   ├── scE2G_element_list.tsv.gz                        # List of candidate elements and associated features
│   │   └── scE2G_gene_list.tsv.gz                           # List of genes in reference file and associated features
│   ├── Neighborhoods/                                    # Results from "Neighborhoods" step of ABC
│   ├── new_features/                                     # Additional computed features
│   ├── Peaks/                                            # Results from "Peaks" step of ABC
│   ├── Predictions/                                      # Results from "Predictions" step of ABC
│   ├── processed_genes_file.bed                          # Processed gene annotations
│   ├── tagAlign/                                         # Results from converting fragment to tagAlign file
│   └── to_generate.txt                                   # Pipeline generation tracking file
└── qc_plots/                                           # Quality control outputs across clusters
  └── predictions_qc_report.html                          # Report of prediction properties compared to reference
  └── [PDFs of prediction property plots]

{IGV_dir (= results_dir if not defined}/             # Genome-browser results directory (only if make_IGV_tracks is True)
├── {cell_cluster}/                                    # Outputs for each cell cluster
│   ├── ATAC.bw                                          # bigWig with unnormalized pseudobulk ATAC signal
│   ├── ATAC_norm.bw                                     # bigWig with read-depth normalized pseudobulk ATAC signal
└── └── {model_name}/                                    # scE2G model results (e.g., multiome_powerlaw_v3)
       └── scE2G_predictions_threshold{threshold}.bedpe   # bedpe file corresponding to filtered predictions
```
  
</details>

<hr>

## Model training

scE2G predicts enhancer-gene links from genomic features.  It learns from CRISPR perturbation data how to use these features to understand the regulatory landscape.  This section explains how to train new versions of the scE2G model using CRISPR data in additional cell types or new features.

### How it works

For each row in the cell clusters configuration table:
1. scE2G takes ATAC fragments files and RNA count matrices, generates candidate E2G links, and computes genome-wide features (number of intralink TSS or elements, normalized ATAC at promoter, ARC E2G score, etc.) for each candidate link, just as it would when you run the model.
2. These genome-wide features are overlapped with the CRISPRi links for the relevant cell type.  The resulting coordinate universe is defined by the CRISPRi links file.  Where a CRISPR perturbation does not map to exactly one candidate E2G link, features are merged or filled in accordance with the feature_table specified in the model configuration table.

For each row in the model configuration table:

3. Resulting CRISPRi links annotated with genome-wide features for each cell type are appended together and sorted into one combined CRISPR dataset.  scE2G Multiome is trained on this dataset, defining weights for each of the features<sup>*</sup> specified in the feature_table.

<sup>\*</sup> *scE2G Multiome's default features are numTSSEnhGene, normalizedATAC_prom, numNearbyEnhancers, ubiqExpresed, numCandidateEnhGene, & ARC.E2G.Score.  scE2G scATAC uses ABC.Score instead of ARC.E2G.Score.*

### Training Configuration

The crispr_cell_types in config_training.yaml dictates the CRISPR data and therefore the CRISPR cell types.  The model_config specifies which ATAC(+RNA) datasets to use to produce the genomewide features and which CRISPRi datafile has the relevant cell types. The dataset config specifies which ATAC(+RNA) datasets match which crispr cell types.

#### 1 Edit **`config/config_training.yaml`**:

Example:
```
### INPUT
model_config: config/model_config.tsv
cell_clusters: config/cell_cluster_config.tsv

### OUTPUT
results_dir: results/out_folder_name

### RESOURCES
crispr_dataset: # crispr_datasets named in the model_config must match one of these keys
  training_multi: path to CRISPRi links in multiple cell types (datafile)
  training_K562: path to CRISPRi links in K562 alone (datafile)
crispr_cell_types: # list cell clusters available in CRISPRi datasets; list items must match cell types in CRISPRi datafile
  training_multi: ["K562", "WTC11"]
  training_K562: ["K562"]
```

**You may use CRISPR datasets available in [CRISPR Comparison GitHub](https://github.com/EngreitzLab/CRISPR_comparison/blob/main/resources/crispr_data).** Be sure to filter your data to only include cell types described in your cell cluster configuration table.

<details>
<summary> More details on scE2G configurations: </summary>

##### OPTIONS

Only set to benchmark_performance True if all cell clusters are K562
```
benchmark_performance: False
```

If True, create ATAC bw and prediction bedpe files.
```
make_IGV_tracks: False
```

Only set fragments_preprocessed to True if you are sure that your fragment files (1) are sorted with sort -k1,1 -k2,2n and (2) only contain fragments on chromosomes in the chromosome sizes reference file (default if not specified: False)
```
fragments_preprocessed: False
```

Set RNA_matrix_filtered to True if the RNA matrix contains the exact set of cells as the ATAC fragment file, and False if it contains more cells (default if not specified: True)
```
RNA_matrix_filtered: True
```

If the cell count in a cell type exceeds the max_cell_count, randomly extract max_cell_count cells.
```
max_cell_count: 20000
```

Maximum memory that can be allocated for submitted jobs (Megabytes)
```
max_memory_allocation_mb: 250000
```

Number of threads for parallel computing of Kendall correlation; must be less than or equal to snakemake's global '-j' parameter
```
threads: 1
```

##### REFERENCE FILES
Should be compatible with the RNA matrix(ces)
```
gene_annotations: "resources/genome_annotations/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz"
```
Gene list with promoter bounds
```
gene_TSS500: "resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed"
```
Gene list with gene body bounds
```
genes: "resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.bed"
```

</details>

#### 2 Make the model configuration table
**`model_config`** columns:
  - `model`: Model name
  - `dataset`: Cell cluster dataset unique identifier(s) (use commas without spaces to separate multiple cell types)
  - `ABC_directory`: ABC results directory
  - `crispr_dataset`: unique identifier indicating the CRISPRi perturbation dataset to use for training; corresponds to a key under crispr_dataset in the training configuration yaml
  - `feature_table`: Path to a configuration table of model features
  - `polynomial`: Use polynomial features? (Note: models with polynomial features cannot be directly used in Apply model workflow)
  - `override_params`: Are there model training parameters you would like to change from the default logistic regression settings specified in `config/config_training.yaml`?

override_params: *See [this example](https://pastebin.com/zt1868R3) `model_config` for how to specify override parameters. If there are no override_params, leave the column blank but still include the header.*

Example of model configuration:
| model | dataset | ABC_directory | crispr_dataset | feature_table | polynomial | override_params |
|-------|---------|---------------|----------------|---------------|------------|-----------------|
| scE2G_train_on_multi | K562_Xu_et_al,WTC11 | | training_multi | resources/feature_tables/multiome_arc_n6.tsv | FALSE | |
| scE2G_train_on_K562 | K562_Xu_et_al | | training_K562 | resources/feature_tables/multiome_arc_n6.tsv | FALSE | |

##### `feature_table` columns: *Feature tables must be specified for each model*
  - feature (name in final table)
  - input_col (name in ABC output)
  - second_input (multiplied by input_col if provided)
  - aggregate_function (how to combine feature values when a CRISPR element overlaps more than one ABC element)
  - fill_value (how to replace NAs)
  - nice_name (used when plotting)

Example of feature configuration:
| feature | input_col | second_input | aggregate_function | fill_value | nice_name |
|---------|-----------|--------------|-------------------|------------|-----------|
| numTSSEnhGene | numTSSEnhGene | NA | max | 0 | # TSSs between E and P |
| normalizedATAC_prom | normalized_atac_prom | NA | mean | 0 | ATAC signal at P |
| numNearbyEnhancers | numNearbyEnhancers | NA | max | 0 | # peaks within 5Kb of E |
| ubiqExpressed | is_ubiquitous_uniform | NA | max | 0 | Ubiquitous expression |
| numCandidateEnhGene | numCandidateEnhGene | NA | max | 0 | # peaks between E and P |
| ARC.E2G.Score | ARC.E2G.Score | NA | mean | 0 | ARC-E2G score |

#### 3 Make the cell cluster configuration table
**`cell_cluster_config`** columns:
If `ABC_directory` not specified, must contain required ABC biosample parameters
  - `cluster`: Each "dataset" in `model_config` must match a unique cell `cluster` row identifier here
  - `crispr_cell_type`: Name of the cell cluster in the CRISPRi links table
  - `rna_matrix_file`: Path to RNA count matrix for this cell cluster
  - `atac_frag_file`: Path to ATAC fragment file for this cell cluster
  - `model_dir`: scE2G models to train, *ex: models/multiome_powerlaw_v3,models/scATAC_powerlaw_v3*
  <details>
  <summary> click to expand required ABC biosample parameters </summary>

  - `HiC_file` *ex: https://s3.us-central-1.wasabisys.com/aiden-encode-hic-mirror/bifocals_iter2/tissues.hic*
  - `HiC_type` *ex: hic*
  - `HiC_resolution` *ex: 5000*
  </details>

Example:
| cluster | crispr_cell_type | rna_matrix_file | atac_frag_file | HiC_file | HiC_type | HiC_resolution | model_dir |
|---------|------------------|-----------------|----------------|----------|----------|----------------|-----------|
| WTC11 | WTC11 | path/to/wtc11_rna_count_matrix.csv.gz | path/to/wtc11_atac_fragments.tsv.gz | path/to/hic | hic | 5000 | models/multiome_powerlaw_v3,models/scATAC_powerlaw_v3 |
| K562_Xu_et_al | K562 | path/to/k562_rna_count_matrix.csv.gz | path/to/k562_atac_fragments.tsv.gz | path/to/hic | hic | 5000 | models/multiome_powerlaw_v3,models/scATAC_powerlaw_v3 |

### Assembling model directory to use in application workflow

Each model directory must contain:
1. `model.pkl`
2. `feature_table.tsv`
3. `score_threshold_{score_threshold}`, where `score_threshold` is a value from 0–1 (e.g., `0.177`)
4. `tpm_threshold_{tpm_threshold}`, where `tpm_threshold` is any non-negative value  (use 0 for ATAC-only models)
5. `qnorm_reference.tsv.gz` (single column with header `E2G.Score` containing raw scores)

### Running training workflow

With Singularity (recommended):
```bash
snakemake -s workflow/Snakefile_training --configfile config/config_training.yaml -j1 --use-conda --use-singularity
```

Or with conda only:
```bash
snakemake -s workflow/Snakefile_training --configfile config/config_training.yaml -j1 --use-conda
```

Output appears at the path to `results_dir` specified in `config_training.yaml`.

<hr>

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Citation

If you use scE2G in your research, please cite our preprint:
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2024.11.23.624931v1)

### Support

For questions and issues, please use the [GitHub Issues](https://github.com/EngreitzLab/scE2G/issues) page.

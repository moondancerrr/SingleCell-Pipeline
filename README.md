# Quality Control, Doublet Detection, Normalization, Feature Selection for Single-Cell Data, Clustering and comprehensive cell-type annotation.

This project provides a Snakemake workflow for performing quality control, normalization, and feature selection on single-cell RNA sequencing data. The workflow includes steps for filtering low-quality cells, correcting ambient RNA contamination, detecting doublets, normalizing data, and selecting highly variable features. It leverages Python, R, and several bioinformatics tools, including [SoupX](https://github.com/constantAmateur/SoupX) and [scDblFinder](https://github.com/plger/scDblFinder) for contamination correction and doublet detection. It also provides comprehensive cell-type annotation through both manual and automated approaches. Automated cell-type annotation using [CellTypist].

## Requirements

### Software

- **Python** (version 3.8+)
- **Snakemake** (version 6+)
- **R** (version 4.0+)

### Python Packages
- `scanpy`
- `seaborn`
- `matplotlib`
- `scipy`
- `anndata2ri`
- `rpy2`
- `python-igraph`
- `leidenalg`
- `scarches`
- `celltypist`

You can install these packages with:
```bash
pip install scanpy seaborn matplotlib scipy anndata2ri rpy2 python-igraph leidenalg scarches celltypist
```

### Adding Package Installation to Snakemake Workflow

If you want to ensure that these packages are installed within the Snakemake workflow, you can add a `conda` environment configuration file to install dependencies automatically. Here’s how:

1. **Create a Conda environment file (e.g., `environment.yaml`)**:

    ```yaml
    name: scRNAseq_qc
    channels:
      - conda-forge
      - bioconda
      - defaults
    dependencies:
      - python=3.8
      - boost=1.85
      - scanpy
      - seaborn
      - matplotlib
      - scipy
      - anndata2ri
      - rpy2
      - celltypist
      - scarches
      - python-igraph
      - leidenalg
      - r-base
      - bioconductor-scater
      - bioconductor-scdbfinder
      - bioconductor-soupx
      - bioconductor-scry
    ```

2. **Activate the Conda environment in Snakemake** by adding this line to your `Snakefile`:

    ```python
    # At the top of your Snakefile
    conda: "environment.yaml"
    ```

This ensures that the environment is set up automatically with all required packages when running the Snakemake workflow.

## Usage

1. Clone the repository:

- git clone https://github.com/moondancerrr/SingleCell-Pipeline.git
- cd SingleCell-Pipeline

2. Run the workflow

- snakemake --cores <number_of_cores>

## Output

The workflow generates:

- s4d8_quality_control.h5ad: Processed AnnData file with quality control and doublet detection results.
- s4d8_normalization.h5ad: AnnData file after applying normalization methods.
- s4d8_feature_selection.h5ad: AnnData file after feature selection with highly deviant genes.
- figures/: Directory containing .png files that visualize key steps of the analysis:
- total_counts_log1p_norm.png: Histogram of total counts and shifted log1p normalization.- 
- scran_normalization.png: Histogram of log1p with Scran estimated size factors.
- analytic_pearson_residuals.png: Histogram of analytic Pearson residuals.
- dispersion_vs_mean_highly_deviant.png: Dispersion vs mean plot for highly deviant genes.

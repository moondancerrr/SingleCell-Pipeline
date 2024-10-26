import argparse
import anndata2ri
import logging
import numpy as np
import os
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import sys
from scipy.stats import median_abs_deviation

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro

parser = argparse.ArgumentParser(description="Quality Control Script")
parser.add_argument("--input", required=True, help="Input filtered_feature_bc_matrix file")
parser.add_argument("--raw", required=True, help="Input raw_feature_bc_matrix file")
parser.add_argument("--output", required=True, help="Output path for processed AnnData file")
parser.add_argument("--figure_dir", required=True, help="Directory to save figures")
args = parser.parse_args()

# Load filtered_feature_bc_matrix
adata = sc.read_10x_h5(filename=args.input)

# Load raw_feature_bc_matrix as SoupX needs it
adata_raw = sc.read_10x_h5(filename=args.raw)

# Create the figure directory if it doesn't exist
os.makedirs(args.figure_dir, exist_ok=True)

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

adata.var_names_make_unique()

# First step: calculate the QC covariates or metric

# Define mitochondrial, ribosomal and hemoglobin genes
# Mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# Ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# Hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

# Calculate the respective QC metrics with scanpy
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)

# Plot the three QC covariates n_genes_by_counts, total_counts and pct_counts_mt per sample to assess how well the respective cells were captured
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
sns.distplot(adata.obs["total_counts"], bins=100, kde=False, ax=ax[0])
ax[0].set_title("Total Counts - Before")

sc.pl.violin(adata, "pct_counts_mt", ax=ax[1], show=False)
ax[1].set_title("pct_counts_mt - Before")
fig.suptitle("Initial QC Metrics - Before Filtering")
plt.savefig("figures/qc_metrics_before_filtering.png")
plt.close()

# Define a function that takes a metric, i.e. a column in .obs and the number of MADs (nmad) that is still permissive within the filtering strategy
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

# Apply this function to the log1p_total_counts, log1p_n_genes_by_counts and pct_counts_in_top_20_genes QC covariates each with a threshold of 5 MADs
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)

# Cells with a percentage of mitochondrial counts exceeding 8 % are filtered out
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)

# Filter our AnnData object based on these two additional columns (outliers)
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

# Plot after filtering
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
sns.distplot(adata.obs["total_counts"], bins=100, kde=False, ax=ax[0])
ax[0].set_title("Total Counts - After")
sc.pl.violin(adata, "pct_counts_mt", ax=ax[1], show=False)
ax[1].set_title("pct_counts_mt - After")
fig.suptitle("QC Metrics - After Filtering")
plt.savefig("figures/qc_metrics_after_filtering.png")
plt.close()
#p1 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# Load the required python and R packages needed for running SoupX - Method for removal of ambient mRNA ￼
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

# Attempt to load SoupX; if not available, install it via BiocManager
try:
    ro.r('library(SoupX)')
except Exception:  # Catch any exception if SoupX is not installed
    ro.r('if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")')
    ro.r('BiocManager::install("SoupX")')
    ro.r('library(SoupX)')

# Copy of our AnnData object, normalize and log1p transform it
adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp)
sc.pp.log1p(adata_pp)

# Compute the principle components of the data to obtain a lower dimensional representation
sc.pp.pca(adata_pp)
# Generate a neighbourhood graph of the data
sc.pp.neighbors(adata_pp)
# Run leiden clustering on the KNN-graph
sc.tl.leiden(adata_pp, key_added="soupx_groups")

# Preprocess variables for SoupX
soupx_groups = adata_pp.obs["soupx_groups"]

# Delete the copy of our AnnData object as we generated a vector of clusters which can be used in soupX
del adata_pp

# Save the cell names, gene names and the data matrix of the filtered cellranger output￼
cells = adata.obs_names
genes = adata.var_names
# Transpose .X as SoupX required a matrix of shape features x barcodes
data = adata.X.T

adata_raw.var_names_make_unique()
data_tod = adata_raw.X.T

del adata_raw

# Load necessary R packages
ro.r('''
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!requireNamespace("SoupX", quietly = TRUE)) BiocManager::install("SoupX")
library(Matrix)
library(SoupX)
''')

# Pass Python data to R environment
ro.globalenv['data'] = data
ro.globalenv['data_tod'] = data_tod
ro.globalenv['genes'] = genes
ro.globalenv['cells'] = cells
ro.globalenv['soupx_groups'] = soupx_groups

# Construct a so called SoupChannel from the table of droplets and the table of cells. Next, we add metadata to the SoupChannel object which can be any metadata in the form of a data.frame
# Execute the code as a block in R
# As a first step, SoupX calculates the profile of the soup. It estimates the ambient mRNA expression profile from empty droplets as given by the unfiltered Cellranger matrix. Next, SoupX estimates the cell specific contamination fraction. Lastly, it corrects the expression matrix according to the ambient mRNA expression profile and the estimated contamination.
ro.r('''
# Specify row and column names of data
rownames(data) <- genes
colnames(data) <- cells

# Ensure correct sparse format for table of counts and table of droplets
data <- as(data, "sparseMatrix")
data_tod <- as(data_tod, "sparseMatrix")

# Generate SoupChannel Object for SoupX 
sc <- SoupChannel(data_tod, data, calcSoupProfile = FALSE)

# Add extra meta data to the SoupChannel object
soupProf <- data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
sc <- setSoupProfile(sc, soupProf)

# Set cluster information in SoupChannel
sc <- setClusters(sc, soupx_groups)

# Estimate contamination fraction
sc <- autoEstCont(sc, doPlot=FALSE)

# Infer corrected table of counts and round to integer
out <- adjustCounts(sc, roundToInt = TRUE)
''')

# Retrieve the output in Python
out = ro.globalenv['out']

adata.layers["counts"] = adata.X
adata.layers["soupX_counts"] = out.T
adata.X = adata.layers["soupX_counts"]

# Min 20 cells - filters out 0 count genes
sc.pp.filter_genes(adata, min_cells=20)

# Doublet detection step
# Load necessary R packages
ro.r('''
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("scater", quietly = TRUE)) install.packages("scater")
if (!requireNamespace("scDblFinder", quietly = TRUE)) BiocManager::install("scDblFinder")
if (!requireNamespace("BiocParallel", quietly = TRUE)) BiocManager::install("BiocParallel")
library(Seurat)
library(scater)
library(scDblFinder)
library(BiocParallel)
''')

# Transfer `adata.X.T` to R as a matrix
data_mat = adata.X.T
ro.globalenv['data_mat'] = data_mat

# Run scDblFinder in R and retrieve results
ro.r('''
set.seed(123)
sce <- scDblFinder(
    SingleCellExperiment(
        list(counts=data_mat)
    )
)
# The final doublet score (the higher the more likely that the cell is a doublet)
doublet_score <- sce$scDblFinder.score
# The classification (doublet or singlet)
doublet_class <- sce$scDblFinder.class
''')

# Transfer `doublet_score` and `doublet_class` back to Python
doublet_score = np.array(ro.globalenv['doublet_score'])
doublet_class = np.array(ro.globalenv['doublet_class'])

# Add results to `adata.obs` in Python
adata.obs["scDblFinder_score"] = doublet_score
adata.obs["scDblFinder_class"] = doublet_class

# Display doublet class counts
print(adata.obs["scDblFinder_class"].value_counts())

# Save the modified AnnData object
adata.write("s4d8_quality_control.h5ad")

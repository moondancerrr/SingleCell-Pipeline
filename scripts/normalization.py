import scanpy as sc
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import anndata2ri
import logging
from scipy.sparse import issparse

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro

# Set up argument parsing
parser = argparse.ArgumentParser(description="Normalization on single-cell data.")
parser.add_argument("--input", required=True, help="Path to input AnnData file")
parser.add_argument("--output_adata", required=True, help="Path to save the output AnnData file")
parser.add_argument("--output_dir", required=True, help="Directory to save output figures")
args = parser.parse_args()

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    # color_map="YlGnBu",
    frameon=False,
)

rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

adata = sc.read(
    filename=args.input
)
# Setting the inplace parameter to False as we want to explore three different normalization techniques 
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# Log1p transform
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

# Save histogram for total counts and shifted log normalization
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total Counts")
sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Shifted Logarithm")
plt.savefig(os.path.join(args.output_dir, "total_counts_log1p_norm.png"))
plt.close()

from scipy.sparse import csr_matrix, issparse
ro.r('''
library(scran)
library(BiocParallel)
''')

# Preliminary clustering for differentiated normalisation
adata_pp = adata.copy()
sc.pp.normalize_total(adata_pp)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="groups")

data_mat = adata_pp.X.T
# Convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
if issparse(data_mat):
    if data_mat.nnz > 2**31 - 1:
        data_mat = data_mat.tocoo()
    else:
        data_mat = data_mat.tocsc()
ro.globalenv["data_mat"] = data_mat
ro.globalenv["input_groups"] = adata_pp.obs["groups"]

ro.r('''
size_factors = sizeFactors(
    computeSumFactors(
        SingleCellExperiment(
            list(counts=data_mat)), 
            clusters = input_groups,
            min.mean = 0.1,
            BPPARAM = MulticoreParam()
    )
)
''')
# Compute the size factors based on the groups of cells calculated before
size_factors = np.array(ro.globalenv['size_factors'])
# Save size_factors in .obs and normalize the data and subsequently apply a log1p transformation
adata.obs["size_factors"] = size_factors
scran = adata.X / adata.obs["size_factors"].values[:, None]
adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))

# Plot and save histogram for Scran normalization
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total Counts")
sns.histplot(adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Log1p with Scran Estimated Size Factors")
plt.savefig(os.path.join(args.output_dir, "scran_normalization.png"))
plt.close()

# Analytic Pearson residuals normalization
analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])

# Plot and save histogram for Analytic Pearson residuals
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total Counts")
sns.histplot(adata.layers["analytic_pearson_residuals"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Analytic Pearson Residuals")
plt.savefig(os.path.join(args.output_dir, "analytic_pearson_residuals.png"))
plt.close()

# Save the updated AnnData object
adata.write("s4d8_normalization.h5ad")






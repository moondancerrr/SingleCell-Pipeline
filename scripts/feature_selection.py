import argparse
import scanpy as sc
import anndata2ri
import logging
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro

# Set up argument parsing
parser = argparse.ArgumentParser(description="Feature selection on single-cell data.")
parser.add_argument("--input", required=True, help="input AnnData file")
args = parser.parse_args()

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

ro.r('''
library(scry)
''')

# Load the already normalized dataset
adata = sc.read(
    filename=args.input
)

# Save the AnnData object in our R environment
ro.globalenv["adata"] = adata
# Call feature selection with deviance on the non-normalized counts matrix￼
ro.r('''
sce = devianceFeatureSelection(adata, assay="X")
''')
# Export the binomial deviance values as a vector
binomial_deviance = ro.r("rowData(sce)$binomial_deviance").T
# Sort the vector an select the top 4,000 highly deviant genes￼
idx = binomial_deviance.argsort()[-4000:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True
# Save them as an additional column in .var as ‘highly_deviant’
adata.var["highly_deviant"] = mask
# Save the computed binomial deviance in case we want to sub-select a different number of highly variable genes afterwards
adata.var["binomial_deviance"] = binomial_deviance

sc.pp.highly_variable_genes(adata, layer="scran_normalization")
# Plotting dispersion versus mean for the genes and color by ‘highly_deviant’
plt.figure(figsize=(8, 6))
ax = sns.scatterplot(data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5)
ax.set_xlim(None, 1.5)
ax.set_ylim(None, 3)
ax.set_title("Dispersion vs Mean for Highly Deviant Genes")
plt.savefig("figures/dispersion_vs_mean_highly_deviant.png")
plt.close()
adata.write("s4d8_feature_selection.h5ad")


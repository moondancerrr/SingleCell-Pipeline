import warning

# Ignore deprecation warnings specifically
warnings.filterwarnings("ignore", category=DeprecationWarning)

# Import additional necessary libraries for data processing and visualization
import numba
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning

# Ignore deprecation warnings related to numba
warnings.simplefilter("ignore", category=NumbaDeprecationWarning)

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import os
import re
from scipy.sparse import csr_matrix
import seaborn as sns
import matplotlib.pyplot as plt
matplotlib.use("Agg")  # Use non-GUI backend to prevent tkinter issues
import celltypist
from celltypist import models
import scarches as sca
import urllib.request

# Suppress performance warnings from pandas
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

# Set default figure parameters for Scanpy plots￼
sc.set_figure_params(figsize=(5, 5))
￼
# Argument parsing
parser = argparse.ArgumentParser(description="Single-cell data processing and visualization")
parser.add_argument("--input", required=True, help="Input AnnData file (h5ad format)")
parser.add_argument("--output_dir", required=True, help="Directory to save output plots")
args = parser.parse_args()

# Load data
adata = sc.read(
    filename="s4d8_clustered.h5ad")

# Set the output directory and create it if it doesn't exist
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

### Manual Annotation ###
# Define marker genes for each cell type
marker_genes = {
     "Progenitor Cells": ["NES", "PROM1"],
     "Astrocytes": ["GFAP", "AQP4"],
     "Oligodendrocytes": ["MBP", "MOG"],
     "Microglia": ["TMEM119", "P2RY12"],
     "Neurons": ["RBFOX3", "MAP2"],
     "Pericytes": ["PDGFRB", "RGS5"],
     "Ependymal Cells": ["S100B", "FOXJ1"],
     "Mesenchymal Stem Cells": ["CD44", "ENG"],
     "Oligodendrocyte Precursor Cells (OPCs)": ["PDGFRA", "CSPG4"],
     "Fibroblasts": ["COL1A1", "FAP"],
     "Mast Cells": ["TPSAB1", "KIT"],
     "Dendritic Cells": ["ITGAX", "HLA-DRA"],
     "Schwann Cells": ["MPZ", "S100B"],
     "Erythrocytes": ["HBB", "HBA1"],
}

# Check if markers are present in the data and store them in a dictionary
marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found

# Set up layers and normalize data using precomputed normalization layers
adata.layers["counts"] = adata.X
adata.X = adata.layers["scran_normalization"]
adata.var["highly_variable"] = adata.var["highly_deviant"]
￼
# Run PCA on the data for dimensionality reduction
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
￼
# Compute the neighborhood graph for clustering
sc.pp.neighbors(adata)
￼
# Compute UMAP for visualization in 2D space
sc.tl.umap(adata)
￼
# Define a list of brain cell types for visualization purposes
brain_cells = [
      "Progenitor Cells",
      "Astrocytes",
      "Oligodendrocytes",
      "Microglia",
      "Neurons",
      "Pericytes",
      "Ependymal Cells",
      "Mesenchymal Stem Cells",
     "Oligodendrocyte Precursor Cells (OPCs)",
     "Fibroblasts",
     "Mast Cells",
     "Dendritic Cells",
     "Schwann Cells",
     "Erythrocytes"
     ]

# Function to clean up filenames by replacing invalid characters
def sanitize_filename(name):
    return re.sub(r'[^A-Za-z0-9_]+', '_', name)

# Loop through each brain cell type to generate UMAP plots for each marker
for ct in brain_cells:
    gene_list = marker_genes_in_data.get(ct)
    if gene_list:
        print(f"{ct.upper()}:")
        
        # Generate UMAP plot without saving, then save manually
        sc.pl.umap(
            adata,
            color=gene_list,
            vmin=0,
            vmax="p99",
            sort_order=False,
            frameon=False,
            cmap="Reds",
            show=False  # Disable display
        )
        
        # Sanitize the cell type name for the filename
        sanitized_ct = sanitize_filename(ct)
        plot_filename = os.path.join(output_dir, f"umap_{sanitized_ct}.png")
        
        # Save the plot manually with plt.savefig
        plt.savefig(plot_filename, bbox_inches="tight")
        plt.close()
        
        print(f"Plot saved: {plot_filename}\n\n\n")
    else:
        print(f"Warning: No genes found for cell type {ct}")

# Perform clustering with a specified resolution
sc.tl.leiden(adata, resolution=2, key_added="leiden_2")
￼
# Filter marker genes to include only those present in data
brain_cells_markers = {
    ct: [m for m in ct_markers if m in adata.var.index]
    for ct, ct_markers in marker_genes.items()
    if ct in brain_cells
    }

# Define the filename for the dot plot
dotplot_filename = os.path.join(output_dir, "dotplot_leiden_2_B_markers.png")

# Generate and save the dot plot
sc.pl.dotplot(
    adata,
    groupby="leiden_2",
    var_names=brain_cells_markers,
    standard_scale="var",  # Standard scale: normalize each gene to range from 0 to 1
    save=True, 
    show=False  # Prevent interactive display
)

# Save the figure manually using plt.savefig
plt.savefig(dotplot_filename, bbox_inches="tight")
plt.close()

# Manual annotation mapping clusters to cell types according to what we see in the dot plot
cl_annotation = {
     "2": "Denderitic Cells",
     "6": "Ependymal Cells",
     "4": "Astrocytes",
     "5": "Mesenchymal Cells"
     }

adata.obs["manual_celltype_annotation"] = adata.obs.leiden_2.map(cl_annotation)

# Save UMAP plot with manual annotations
umap_filename = os.path.join(output_dir, "umap_manual_celltype_annotation.png")
sc.pl.umap(
    adata,
    color=["manual_celltype_annotation"],
    show=False,  # Prevent interactive display
    save=True 
)
# Save the figure manually using plt.savefig
plt.savefig(umap_filename, bbox_inches="tight")

# Identify differentially expressed genes in clusters
sc.tl.rank_genes_groups(
    adata, groupby="leiden_2", method="wilcoxon", key_added="dea_leiden_2"
)

# Filter differential expression results
sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="dea_leiden_2",
    key_added="dea_leiden_2_filtered",
)

# Define the filename for the rank genes groups dot plot
rank_genes_dotplot_filename = os.path.join(output_dir, "rank_genes_groups_dotplot_leiden_2.png")

# Generate and save the rank genes groups dot plot
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden_2",
    standard_scale="var",  # Standard scale: normalize each gene to range from 0 to 1
    n_genes=5,
    key="dea_leiden_2_filtered",
    show=False,  # Prevent interactive display
    save=True 
    )

# Save the figure manually using plt.savefig
plt.savefig(rank_genes_dotplot_filename, bbox_inches="tight")

# Define the filename 
microglia_markers_filename = os.path.join(output_dir, "microglia_markers_plot.png")

# Generate and save UMAP plot of microglia markers
sc.pl.umap(
     adata,
     color=["CTSS", "TYROBP"], vmax="p99",
     legend_loc="on data",
     frameon=False,
     cmap="Reds",
     how=False,  # Prevent interactive display
     save=True
)

# Save the figure manually using plt.savefig
plt.savefig(microglia_markers_filename, bbox_inches="tight")

# Update annotations and save final annotated dataset
cl_annotation["3"] = "microglia"
adata.obs["manual_celltype_annotation"] = adata.obs.leiden_2.map(cl_annotation)

# Final save
adata.write("annotation_adata_out.h5ad")

### Automatic Annotation ###
# Copy data and set raw counts for celltypist model
adata_celltypist = adata.copy()  # make a copy of our adata

# Normalize and log-transform data for celltypist
adata_celltypist.X = adata.layers["counts"]  # set adata.X to raw counts
sc.pp.normalize_per_cell(
    adata_celltypist, counts_per_cell_after=10**4
)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform
# Make .X dense instead of sparse, for compatibility with celltypist:
adata_celltypist.X = adata_celltypist.X.toarray()

# Download the celltypist models for Developing_Human_Brain cells
models.download_models(
     force_update=True, model=["Developing_Human_Brain.pkl"]
     )

# Load model
model_high = models.Model.load(model="Developing_Human_Brain.pkl")

# Run model
predictions_high = celltypist.annotate(
    adata_celltypist, model=model_high, majority_voting=True
)

# Transform the predictions to adata to get the full output
predictions_high_adata = predictions_high.to_adata()
￼
# Copy the results to our original AnnData object
adata.obs["celltypist_cell_label_coarse"] = predictions_high_adata.obs.loc[
    adata.obs.index, "majority_voting"
]
adata.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs.loc[
    adata.obs.index, "conf_score"
]

# Plot
celltypist_umap = os.path.join(output_dir, "celltypist_umap.png")

sc.pl.umap(
    adata,
    color=["celltypist_cell_label_coarse", "celltypist_conf_score_coarse"],
    frameon=False,
    sort_order=False,
    wspace=1,
    save=True,
    show=False  # Prevent interactive display
)

# Save the figure manually using plt.savefig
plt.savefig(celltypist_umap, bbox_inches="tight")
plt.close()

# Dendrogram partly reflects prior knowledge on cell type relation
celltypist_dendrogram = os.path.join(output_dir, "celltypist_dendrogram.png")

# Generate dendrogram for cell relationships in celltypist annotations
sc.pl.dendrogram(adata,
      groupby="celltypist_cell_label_fine",
      save=True,
      show=False  # Prevent interactive display
)

# Save the figure manually using plt.savefig
plt.savefig(celltypist_umap, bbox_inches="tight")
plt.close()

celltypist_manual = os.path.join(output_dir, "celltypist_manual.png")

# Save UMAP plot comparing manual and celltypist annotations
sc.pl.umap(
    adata,
    color=["manual_celltype_annotation"],
    frameon=False,
    save=True,
    show=False  # Prevent interactive display
)

# Save the figure manually using plt.savefig
plt.savefig(celltypist_manual, bbox_inches="tight")
plt.close()

# Breakdown of cluster 3 in terms of fine celltypist labels
pd.crosstab(adata.obs.leiden_2, adata.obs.celltypist_cell_label_fine).loc["3", :].sort_values(ascending=False)

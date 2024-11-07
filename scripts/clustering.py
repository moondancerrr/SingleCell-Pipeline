import argparse
import matplotlib.pyplot as plt
import scanpy as sc

# Set up argument parsing
parser = argparse.ArgumentParser(description="Feature selection on single-cell data.")
parser.add_argument("--input", required=True, help="input AnnData file")
parser.add_argument("--output_dir", required=True, help="Directory to save output figures")
parser.add_argument("--output_adata", required=True, help="Output path for the clustered AnnData file")
args = parser.parse_args()

# Configure Scanpy settings for output verbosity and figure parameters
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80,
      facecolor="white",
      frameon=False
      )

# Load data
adata = sc.read(
      filename=args.input
      )

# Calculate the KNN graph on a lower-dimensional gene expression representation￼
sc.pp.neighbors(adata, n_pcs=30)

# Embed the cells into a UMAP embedding
sc.tl.umap(adata)

# Call the Leiden algorithm for clustering cells
sc.tl.leiden(adata)

# Inspect the impact of different resolutions on the clustering result
sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

# Set up a figure for visualizing the UMAP embedding with different clustering resolutions
plt.figure(figsize=(8, 6))
sc.pl.umap(
    adata,
    color=["leiden_res0_25", "leiden_res0_5", "leiden_res1"],
    legend_loc="on data",
    show=False  # Disable interactive display
    )

# Save the UMAP plot with different clustering resolutions to the specified directory
output_path = os.path.join(args.output_dir, "leiden_clustering.png")
plt.savefig(output_path)
plt.close()

# Save the clustered AnnData object to the specified output path
adata.write(args.output_adata)


rule all:
    input:
        "s4d8_quality_control.h5ad",
        directory("figures")

# Define the URLs for the input files
URL_FILTERED = "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Brain_Tumor_3p_LT/Brain_Tumor_3p_LT_filtered_feature_bc_matrix.h5"
URL_RAW = "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Brain_Tumor_3p_LT/Brain_Tumor_3p_LT_raw_feature_bc_matrix.h5"

# Define paths to the downloaded files
rule download_filtered:
    output:
        "data/Brain_Tumor_3p_LT_filtered_feature_bc_matrix.h5"
    shell:
        """
        curl -o {output} -L {URL_FILTERED} || wget -O {output} {URL_FILTERED}
        """

rule download_raw:
    output:
        "data/Brain_Tumor_3p_LT_raw_feature_bc_matrix.h5"
    shell:
        """
        curl -o {output} -L {URL_RAW} || wget -O {output} {URL_RAW}
        """

# Main QC rule that depends on downloaded files
rule QC:
    input:
        input_file="data/Brain_Tumor_3p_LT_filtered_feature_bc_matrix.h5",
        raw_file="data/Brain_Tumor_3p_LT_raw_feature_bc_matrix.h5"
    output:
        processed_adata="s4d8_quality_control.h5ad",
        figures_dir=directory("figures")
    shell:
        """
        python scripts/quality_control.py --input {input.input_file} --raw {input.raw_file} --output {output.processed_adata} --figure_dir {output.figures_dir}
        """

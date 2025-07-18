"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
This script generates a boxplot and stripplot of normalized gene expression
for a specified gene across sample groups.

USAGE:
    python script_name.py <norm_counts_file> <design_file> <gene_symbol>

ARGUMENTS:
    norm_counts_file : CSV file with normalized gene expression values (genes as rows, samples as columns).
    design_file      : CSV file with study design metadata, must include 'runID' and the grouping column.
    gene_symbol      : Gene name to plot (must exist in norm_counts_file index).
    
NOTES:
    - By default, the grouping variable used is 'antiCCP'. You can modify this in the function call.
    - The plot is saved as a PNG in the 'plots' directory.
    - In the provided code, file paths and gene symbol are hardcoded for testing. Uncomment sys.argv lines to enable CLI usage.

REQUIREMENTS:
    - pandas
    - seaborn
    - matplotlib
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

def plot_gene_expression(norm_counts_file, design_file, gene, group_col='antiCCP', output_dir="plots"):
    norm_counts = pd.read_csv(norm_counts_file, index_col=0)
    design = pd.read_csv(design_file)

    if gene not in norm_counts.index:
        print(f"[ERROR] Gene '{gene}' not found in normalized counts.")
        return

    merged = pd.DataFrame({
        'expression': norm_counts.loc[gene],
        'group': design.set_index('runID').loc[norm_counts.columns][group_col]
    }).reset_index()

    sns.boxplot(data=merged, x='group', y='expression')
    sns.stripplot(data=merged, x='group', y='expression', color='black', alpha=0.5)
    plt.title(f"Expression of {gene}")
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f"{output_dir}/{gene}_expression.png")
    plt.close()

if __name__ == "__main__":
    #norm_counts_file = sys.argv[1]
    #design_file = sys.argv[2]
    #gene = sys.argv[3]
    #******** Adjust the file path ****************
    norm_counts_file = "/THE_PATH/results/counts/norm_counts_with_symbols.csv"
    design_file = "/THE_PATH/meta/study_design.csv"
    gene = "PLCH2"
    plot_gene_expression(norm_counts_file, design_file, gene)


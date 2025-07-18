"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
This script generates a heatmap of the top differentially expressed genes (DEGs) using normalized count data.

Usage:
    python heatmap_plot.py <norm_counts_file> <design_file> <top_genes_file>

Arguments:
    norm_counts_file : CSV file of normalized gene expression counts (gene symbols as row indices).
    design_file      : CSV file containing the experimental design with sample metadata (must include 'runID' and the grouping column).
    top_genes_file   : CSV file with a column named 'Gene' listing top DEGs to be plotted.

Output:
    A hierarchical clustered heatmap will be saved as a PNG file to the path specified in the `output_path` argument (default: "plots/top_degs_heatmap.png").

Example:
    python heatmap_plot.py norm_counts_with_symbols.csv study_design.csv top_genes.csv

Note:
    - Ensure all input files are properly formatted.
    - Edit the 'group_col' parameter in the `plot_heatmap` function if you want to color samples by a different metadata field.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import sys
import os

def plot_heatmap(norm_counts_file, design_file, top_genes_file, group_col='antiCCP', output_path="plots/top_degs_heatmap.png"):
    norm_counts = pd.read_csv(norm_counts_file, index_col=0)  # now gene symbols
    design = pd.read_csv(design_file)
    top_genes = pd.read_csv(top_genes_file)['Gene'].tolist()

    # Filter to genes present in normalized data
    top_genes = [g for g in top_genes if g in norm_counts.index]
    if not top_genes:
        print("[WARNING] No top genes found in norm_counts_with_symbols.")
        return

    # Subset and scale
    subset = norm_counts.loc[top_genes]
    scaled = StandardScaler().fit_transform(subset.T)
    scaled_df = pd.DataFrame(scaled, index=norm_counts.columns, columns=top_genes)

    # Add group colors from metadata
    group_map = design.set_index('runID')[group_col].to_dict()
    sample_groups = [group_map.get(sample, 'Unknown') for sample in scaled_df.index]
    unique_groups = sorted(set(sample_groups))
    palette = sns.color_palette("Set2", n_colors=len(unique_groups))
    lut = dict(zip(unique_groups, palette))
    col_colors = pd.Series(sample_groups, index=scaled_df.index).map(lut)

    # Plot clustermap
    g = sns.clustermap(
        scaled_df,
        cmap="vlag",
        row_colors=col_colors,
        figsize=(12, 8),
        col_cluster=True,
        yticklabels=True,
        xticklabels=True
    )

    # Add legend
    for label in unique_groups:
        g.ax_col_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(loc="center", ncol=len(unique_groups), title=group_col)

    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path)
    plt.close()


if __name__ == "__main__":
    #norm_counts_file = sys.argv[1]
    #design_file = sys.argv[2]
    #top_genes_file = sys.argv[3]
    #***************Adjust the right path **********************
    norm_counts_file = "RIGHT_PATH/one_w_results/counts/norm_counts_with_symbols.csv"
    design_file = "RIGHT_PATH/meta/study_design.csv"
    top_genes_file = "RIGHT_PATH/one_w_results/DEGs/top_genes.csv"
    plot_heatmap(norm_counts_file, design_file, top_genes_file)


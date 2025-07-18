"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################

Small RNA Heatmap Plotter
-------------------------

This script generates a clustered heatmap of the top differentially expressed small RNA sequences
based on normalized expression counts (CPM). It is useful for visualizing expression patterns 
across samples grouped by experimental conditions (e.g., antiCCP status).

Usage:
- Prepare the following CSV files:
    1. Normalized counts file: a matrix with sequences as rows and samples as columns.
    2. Design file: metadata with sample information including a 'runID' column and a grouping column (e.g., 'antiCCP').
    3. Top sequences file: a list of the most differentially expressed sequences, with a column named 'TopSequences'.

- You can run this script directly or call the `plot_small_rna_heatmap()` function in another program.

Default file paths are set in the `__main__` block, but they can be modified to fit your directory structure.
The resulting heatmap will be saved in: plots/top_small_rna_heatmap.png

Dependencies:
- pandas
- seaborn
- matplotlib
- scikit-learn
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import os

def plot_small_rna_heatmap(norm_counts_file, design_file, top_seqs_file, group_col='antiCCP', output_path="plots/top_small_rna_heatmap.png"):
    """
    Plots a heatmap for top differentially expressed small RNA sequences.
    
    Parameters:
    - norm_counts_file: Path to normalized CPM matrix (sequences x samples)
    - design_file: Path to sample metadata (must contain 'runID' and group_col)
    - top_seqs_file: Path to CSV with top small RNA sequences (column name: 'TopSequences')
    - group_col: Column name for sample grouping (e.g., treatment or condition)
    - output_path: File path to save the heatmap
    """
    norm_counts = pd.read_csv(norm_counts_file, index_col=0)
    design = pd.read_csv(design_file)
    top_seqs = pd.read_csv(top_seqs_file)['TopSequences'].tolist()

    # Ensure top sequences are in the matrix
    top_seqs = [s for s in top_seqs if s in norm_counts.index]
    if not top_seqs:
        print("[WARNING] No top sequences found in normalized counts.")
        return

    # Subset and scale
    subset = norm_counts.loc[top_seqs]
    scaled = StandardScaler().fit_transform(subset.T)
    scaled_df = pd.DataFrame(scaled, index=subset.columns, columns=top_seqs)

    # Group annotations
    group_map = design.set_index('runID')[group_col].to_dict()
    sample_groups = [group_map.get(sample, 'Unknown') for sample in scaled_df.index]
    unique_groups = sorted(set(sample_groups))
    palette = sns.color_palette("Set2", n_colors=len(unique_groups))
    lut = dict(zip(unique_groups, palette))
    col_colors = pd.Series(sample_groups, index=scaled_df.index).map(lut)

    # Clustered heatmap
    g = sns.clustermap(
        scaled_df,
        cmap="vlag",
        row_colors=col_colors,
        figsize=(12, 8),
        col_cluster=True,
        yticklabels=True,
        xticklabels=True
    )

    # Add legend for sample groups
    for label in unique_groups:
        g.ax_col_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(loc="center", ncol=len(unique_groups), title=group_col)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path)
    plt.close()

if __name__ == "__main__":
    #Replace RIGHT_PATH
    norm_counts_file = "RIGHT_PATH/results/counts/norm_counts.csv"
    design_file = "RIGHT_PATH/meta/study_design.csv"
    top_seqs_file = "RIGHT_PATH/results/DEGs/top_sequences.csv"
    plot_small_rna_heatmap(norm_counts_file, design_file, top_seqs_file)


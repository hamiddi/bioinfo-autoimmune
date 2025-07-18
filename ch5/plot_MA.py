"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
MA Plot Generator for Differential Expression Analysis

This script creates an MA plot (log2 fold change vs. mean expression) from two input CSV files:
1. A file containing differential expression results (with columns: 'Gene' or 'GeneSymbol', 'log2FC', and 'adj-p').
2. A file containing normalized gene expression counts (with genes as row indices and samples as columns).

Usage:
------
You can run the script from the command line as:
    python ma_plot_script.py <degs_file.csv> <normalized_counts.csv>

Alternatively, you can modify the file paths in the __main__ block for manual execution.

Arguments:
----------
- degs_file         : CSV file with differential expression results.
- norm_counts_file  : CSV file with normalized counts per gene.
- output_path       : Path to save the MA plot image (default: "plots/ma_plot.png").
- pval_thresh       : Adjusted p-value threshold to consider significance (default: 0.05).
- fc_thresh         : Fold change threshold to consider up/down regulation (default: 1.0).
- show_plot         : Whether to display the plot in an interactive window (default: True).

Output:
-------
- A PNG image of the MA plot showing upregulated, downregulated, and non-significant genes.
- Top 10 most significant genes are labeled on the plot.

Dependencies:
-------------
- pandas
- numpy
- matplotlib
- seaborn
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def ma_plot(degs_file, norm_counts_file, output_path="plots/ma_plot.png",
            pval_thresh=0.05, fc_thresh=1.0, show_plot=True):
    degs = pd.read_csv(degs_file).dropna(subset=['adj-p'])
    norm_counts = pd.read_csv(norm_counts_file, index_col=0)

    # Average duplicates if present
    norm_counts = norm_counts.groupby(norm_counts.index).mean()

    # Average expression per gene (log2)
    avg_expr = norm_counts.mean(axis=1)
    avg_expr_log2 = np.log2(avg_expr + 1e-6)

    id_field = 'Gene'
    if degs[id_field].iloc[0] not in avg_expr_log2.index and 'GeneSymbol' in degs.columns:
        id_field = 'GeneSymbol'

    degs['A'] = degs[id_field].map(avg_expr_log2)
    df = degs.dropna(subset=['A']).copy()

    def label(row):
        if row['adj-p'] < pval_thresh:
            if row['log2FC'] > fc_thresh:
                return 'Upregulated'
            elif row['log2FC'] < -fc_thresh:
                return 'Downregulated'
        return 'Not significant'

    df['status'] = df.apply(label, axis=1)

    top_genes = df.nsmallest(10, 'adj-p').copy()
    top_genes['label'] = top_genes.get('GeneSymbol', top_genes['Gene'])
    df['label'] = ""
    df.loc[top_genes.index, 'label'] = top_genes['label']

    # Plot
    plt.figure(figsize=(12, 8), dpi=300)
    sns.scatterplot(data=df, x='A', y='log2FC', hue='status',
                    palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not significant': 'grey'},
                    alpha=0.7, edgecolor=None)

    plt.axhline(fc_thresh, linestyle='--', color='black')
    plt.axhline(-fc_thresh, linestyle='--', color='black')
    plt.axhline(0, linestyle='-', color='black')

    for _, row in top_genes.iterrows():
        plt.text(row['A'], row['log2FC'], row['label'], fontsize=10, color='black', ha='right')

    plt.xlabel("Average Expression (log2 CPM)", fontsize=16)
    plt.ylabel("log2 Fold Change", fontsize=16)
    plt.title("MA Plot with Gene Labels", fontsize=20)
    plt.legend(title="Regulation", fontsize=12, title_fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path)
    if show_plot:
        plt.show()
    plt.close()

if __name__ == "__main__":
    #degs_file = sys.argv[1]           # e.g. results/DEGs/diff_expr_with_symbols.csv
    #norm_counts_file = sys.argv[2]    # e.g. results/counts/norm_counts_with_symbols.csv
    #Adjust the right path
    degs_file = "RIGHT_PATH/one_w_results/DEGs/diff_expr_with_symbols_log2fc.csv"
    norm_counts_file = "RIGHT_PATH/one_w_results/counts/norm_counts_with_symbols.csv"
    ma_plot(degs_file, norm_counts_file)


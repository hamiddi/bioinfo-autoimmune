"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################

Volcano Plot Generator for Differential Gene Expression Results

This script reads a CSV file containing differential gene expression data and generates a volcano plot
highlighting significantly upregulated and downregulated genes based on specified thresholds.

Usage:
    python volcano_plot.py <path_to_DEGs_csv_file>

Arguments:
    <path_to_DEGs_csv_file>: Path to the CSV file containing differential expression results.
                             The file must include columns: 'log2FC', 'adj-p', and either 'GeneSymbol' or 'Gene'.

Optional Parameters (can be modified inside the script or the function call):
    - output_path: File path to save the plot (default: plots/volcano_plot.png)
    - pval_thresh: Adjusted p-value threshold for significance (default: 0.05)
    - fc_thresh: Log2 fold change threshold for biological significance (default: 1.0)
    - show_plot: Whether to display the plot interactively (default: True)

Example:
    python volcano_plot.py ../../../one_w_results/DEGs/diff_expr_with_symbols_log2fc.csv

Note:
    The script will automatically create the output directory if it does not exist.
    The top 10 most significant genes are labeled on the plot.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def plot_volcano(degs_file, output_path="plots/volcano_plot.png", pval_thresh=0.05, fc_thresh=1.0, show_plot=True):
    df = pd.read_csv(degs_file).dropna(subset=['adj-p'])
    df['-log10p'] = -np.log10(df['adj-p'] + 1e-10)

    # Define status
    def label(row):
        if row['adj-p'] < pval_thresh:
            if row['log2FC'] > fc_thresh:
                return 'Upregulated'
            elif row['log2FC'] < -fc_thresh:
                return 'Downregulated'
        return 'Not significant'

    df['status'] = df.apply(label, axis=1)

    # Label top genes
    top_genes = df.nsmallest(10, 'adj-p').copy()
    top_genes['label'] = top_genes.get('GeneSymbol', top_genes['Gene'])
    df['label'] = ""
    df.loc[top_genes.index, 'label'] = top_genes['label']

    # Plot
    plt.figure(figsize=(12, 8), dpi=300)
    sns.scatterplot(data=df, x='log2FC', y='-log10p', hue='status',
                    palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not significant': 'grey'},
                    alpha=0.7, edgecolor=None)

    plt.axhline(-np.log10(pval_thresh), linestyle='--', color='black')
    plt.axvline(fc_thresh, linestyle='--', color='black')
    plt.axvline(-fc_thresh, linestyle='--', color='black')

    for _, row in top_genes.iterrows():
        plt.text(row['log2FC'], row['-log10p'], row['label'], fontsize=10, color='black', ha='right')

    plt.xlabel("log2 Fold Change", fontsize=16)
    plt.ylabel("-log10 Adjusted P-value", fontsize=16)
    plt.title("Volcano Plot with Top Gene Labels", fontsize=20)
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
    #degs_file = sys.argv[1]
    #Replace RIGHT_PATH with the path  
    degs_file = "RIGHT_PATH/one_w_results/DEGs/diff_expr_with_symbols_log2fc.csv"
    plot_volcano(degs_file)


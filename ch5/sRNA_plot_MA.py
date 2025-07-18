"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################

MA Plot Generator for Small RNA-Seq Data
----------------------------------------

This Python script generates an MA plot for small RNA-Seq data based on differential expression results
and normalized count data. The plot highlights upregulated and downregulated sequences and labels the
top 10 most significant sequences based on adjusted p-values.

USAGE:
    Run the script directly, or import the `ma_plot_small_rna` function in another Python module.

    By default, the script expects:
        - Differential expression results CSV file at:    ../../../results/DEGs/diff_expr_log2fc.csv
        - Normalized counts CSV file at:                  ../../../results/counts/norm_counts.csv

    These files must include:
        - A 'Sequence' column in the DEGs file (for small RNA identifiers)
        - Normalized counts with sequence IDs as the row index and samples as columns

    OUTPUT:
        - A high-resolution MA plot will be saved to: plots/ma_plot_small_rna.png
        - Optionally displays the plot if `show_plot=True`

FUNCTION PARAMETERS:
    ma_plot_small_rna(degs_file, norm_counts_file,
                      output_path="plots/ma_plot_small_rna.png",
                      pval_thresh=0.05, fc_thresh=1.0, show_plot=True)

    You can adjust the p-value threshold, fold change threshold, and output path.

REQUIREMENTS:
    pandas, numpy, matplotlib, seaborn
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def ma_plot_small_rna(degs_file, norm_counts_file, output_path="plots/ma_plot_small_rna.png",
                      pval_thresh=0.05, fc_thresh=1.0, show_plot=True):
    """
    Creates an MA plot for small RNA-Seq data using sequence-level expression.
    """
    degs = pd.read_csv(degs_file).dropna(subset=['adj-p'])
    norm_counts = pd.read_csv(norm_counts_file, index_col=0)

    # Average duplicates if any (some sequences may be identical)
    norm_counts = norm_counts.groupby(norm_counts.index).mean()

    # Compute average expression per sequence
    avg_expr = norm_counts.mean(axis=1)
    avg_expr_log2 = np.log2(avg_expr + 1e-6)

    id_field = 'Sequence'
    if id_field not in degs.columns:
        raise ValueError("The differential expression file must contain a 'Sequence' column for small RNA IDs.")

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

    top_seqs = df.nsmallest(10, 'adj-p').copy()
    top_seqs['label'] = top_seqs['Sequence']
    df['label'] = ""
    df.loc[top_seqs.index, 'label'] = top_seqs['label']

    # Plot MA
    plt.figure(figsize=(12, 8), dpi=300)
    sns.scatterplot(data=df, x='A', y='log2FC', hue='status',
                    palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not significant': 'grey'},
                    alpha=0.7, edgecolor=None)

    plt.axhline(fc_thresh, linestyle='--', color='black')
    plt.axhline(-fc_thresh, linestyle='--', color='black')
    plt.axhline(0, linestyle='-', color='black')

    for _, row in top_seqs.iterrows():
        plt.text(row['A'], row['log2FC'], row['label'], fontsize=10, color='black', ha='right')

    plt.xlabel("Average Expression (log2 CPM)", fontsize=16)
    plt.ylabel("log2 Fold Change", fontsize=16)
    plt.title("MA Plot of Small RNA Sequences", fontsize=20)
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
    #Replace RIGHT_PATH
    degs_file = "RIGHT_PATH/results/DEGs/diff_expr_log2fc.csv"
    norm_counts_file = "RIGHT_PATH/results/counts/norm_counts.csv"
    ma_plot_small_rna(degs_file, norm_counts_file)


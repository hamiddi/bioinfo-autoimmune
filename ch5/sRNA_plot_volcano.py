"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
This script generates a volcano plot for small RNA-Seq differential expression results.

USAGE:
- Make sure your input CSV file (degs_file) contains at least the following columns:
    'Sequence'   : Identifier of the small RNA sequence
    'log2FC'     : Log2 fold change of expression
    'adj-p'      : Adjusted p-value

- The default file path for input is:
    '../../../results/DEGs/diff_expr_log2fc.csv'
  Modify it as needed or call the function directly with a different path.

- The script will output a volcano plot to:
    'plots/volcano_plot_small_rna.png'
  You can change this path by passing a different `output_path` to the function.

- Thresholds for significance and labeling can be adjusted via:
    - `pval_thresh`: Default = 0.05
    - `fc_thresh`: Default = 1.0

- If run directly (as main), it will execute with the default input and output paths.

REQUIREMENTS:
- pandas
- numpy
- matplotlib
- seaborn
- Python 3.6+

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_volcano_small_rna(degs_file, output_path="plots/volcano_plot_small_rna.png", pval_thresh=0.05, fc_thresh=1.0, show_plot=True):
    """
    Generates a volcano plot for small RNA-Seq differential expression results.

    Parameters:
    - degs_file: Path to differential expression results (must include 'Sequence', 'log2FC', 'adj-p')
    - output_path: Output file path for the plot
    - pval_thresh: Significance threshold for adjusted p-values
    - fc_thresh: Fold-change threshold for labeling
    - show_plot: Whether to display the plot interactively
    """
    df = pd.read_csv(degs_file).dropna(subset=['adj-p'])
    df['-log10p'] = -np.log10(df['adj-p'] + 1e-10)

    # Label differential status
    def label(row):
        if row['adj-p'] < pval_thresh:
            if row['log2FC'] > fc_thresh:
                return 'Upregulated'
            elif row['log2FC'] < -fc_thresh:
                return 'Downregulated'
        return 'Not significant'

    df['status'] = df.apply(label, axis=1)

    # Label top sequences
    top_seqs = df.nsmallest(10, 'adj-p').copy()
    top_seqs['label'] = top_seqs['Sequence']
    df['label'] = ""
    df.loc[top_seqs.index, 'label'] = top_seqs['label']

    # Plot
    plt.figure(figsize=(12, 8), dpi=300)
    sns.scatterplot(data=df, x='log2FC', y='-log10p', hue='status',
                    palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not significant': 'grey'},
                    alpha=0.7, edgecolor=None)

    plt.axhline(-np.log10(pval_thresh), linestyle='--', color='black')
    plt.axvline(fc_thresh, linestyle='--', color='black')
    plt.axvline(-fc_thresh, linestyle='--', color='black')

    for _, row in top_seqs.iterrows():
        plt.text(row['log2FC'], row['-log10p'], row['label'], fontsize=10, color='black', ha='right')

    plt.xlabel("log2 Fold Change", fontsize=16)
    plt.ylabel("-log10 Adjusted P-value", fontsize=16)
    plt.title("Volcano Plot of Small RNA Sequences", fontsize=20)
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
    plot_volcano_small_rna(degs_file)


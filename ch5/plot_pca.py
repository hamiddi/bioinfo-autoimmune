"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################

This script performs Principal Component Analysis (PCA) on RNA-Seq gene expression data 
and generates a 2D PCA scatter plot, where each point represents a sample colored by a 
specified metadata group (e.g., condition or treatment).

Usage:
1. Ensure you have two CSV files:
   - One containing normalized gene expression data with genes as rows and samples as columns.
   - One containing metadata with sample IDs as the index and at least one column for grouping.

2. Set the paths to your expression and metadata files in the `expression_file` and 
   `metadata_file` variables under the `__main__` section.

3. Update the `group_column` parameter in the `plot_pca` function call to match the column 
   in the metadata you want to use for coloring the samples (e.g., 'condition', 'treatment', 'antiCCP').

4. Run the script:
       python pca_plot_script.py

Dependencies:
- pandas
- matplotlib
- seaborn
- scikit-learn

The PCA plot will be saved to the specified output file path ('plots/pca_plot.png' by default) 
and displayed on screen.
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns

def plot_pca(expression_data: pd.DataFrame, metadata: pd.DataFrame, group_column: str, n_components=2, output_file="plots/pca_plot.png"):
    """
    Generates and saves a PCA plot from RNA-Seq expression data with sample labels.

    Parameters:
    - expression_data: DataFrame of gene expression with samples as columns and genes as rows.
    - metadata: DataFrame with sample metadata. Must include a column matching the sample names.
    - group_column: Column in metadata used for coloring (e.g., treatment group or condition).
    - n_components: Number of PCA components (default: 2).
    - output_file: File path to save the plot (default: 'pca_plot.png').
    """
    
    # Transpose expression data to have samples as rows
    expression_data_t = expression_data.T

    # Standardize the data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(expression_data_t)

    # PCA
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(X_scaled)

    # Create PCA DataFrame
    pca_df = pd.DataFrame(pcs, columns=[f'PC{i+1}' for i in range(n_components)])
    pca_df[group_column] = metadata.loc[expression_data_t.index, group_column].values
    pca_df["Sample"] = expression_data_t.index

    # Plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=pca_df,
        x='PC1', y='PC2',
        hue=group_column,
        palette='Set2',
        s=100, edgecolor='k'
    )

    # Add sample name labels with smaller font
    for _, row in pca_df.iterrows():
        plt.text(
            row['PC1'], row['PC2'],
            row['Sample'],
            fontsize=8,
            ha='right',
            va='bottom'
        )

    plt.title('PCA Plot of RNA-Seq Samples', fontsize=18)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)', fontsize=14)
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(title=group_column, title_fontsize=13, fontsize=12)
    plt.grid(True)
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_file, dpi=300)
    plt.show()
    plt.close()



if __name__ == "__main__":
    # Load example expression data and metadata
    #Replace RIGHT_PATH with the path
    expression_file = "RIGHT_PATH/one_w_results/counts/norm_counts_with_symbols.csv"  # Rows = genes, Columns = samples
    metadata_file = "RIGHT_PATH/meta/study_design.csv"      # Columns = sample info (e.g., condition)

    # Read data
    expression_data = pd.read_csv(expression_file, index_col=0)
    metadata = pd.read_csv(metadata_file, index_col=0)

    # Assume 'condition' column in metadata
    plot_pca(expression_data, metadata, group_column='antiCCP')


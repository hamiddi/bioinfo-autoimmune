"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
PCA Plot Generator for Small RNA-Seq Data
-----------------------------------------
This script generates a PCA (Principal Component Analysis) scatter plot using normalized small RNA-Seq expression data 
(counts at the sequence level) and associated metadata. It is useful for visualizing sample clustering based on biological 
or clinical groupings.

Usage:
1. Make sure the following input files exist:
   - A CSV file containing normalized expression counts with sequences as rows and sample IDs as columns.
   - A CSV metadata file with sample IDs as rows and metadata fields (e.g., condition, treatment) as columns.

2. Update the file paths in the `__main__` section:
   - `expression_file`: Path to the normalized expression data CSV.
   - `metadata_file`: Path to the sample metadata CSV.

3. Set the `group_column` to the name of the column in the metadata file that contains group labels (e.g., 'antiCCP').

4. Run the script:
   ```bash
   python this_script.py
"""
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns

def plot_small_rna_pca(expression_data: pd.DataFrame, metadata: pd.DataFrame, group_column: str, n_components=2, output_file="plots/pca_plot_small_rna.png"):
    """
    Generates a PCA plot for small RNA-Seq expression data using sequence-level counts.

    Parameters:
    - expression_data: DataFrame of expression with sequences as rows and samples as columns.
    - metadata: DataFrame with sample metadata.
    - group_column: Column in metadata for group labeling (e.g., condition).
    - n_components: Number of PCA components.
    - output_file: Where to save the PCA image.
    """
    # Transpose expression data: rows = samples, columns = sequences
    expression_data_t = expression_data.T

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(expression_data_t)

    # Perform PCA
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(X_scaled)

    # Create DataFrame for PCA results
    pca_df = pd.DataFrame(pcs, columns=[f'PC{i+1}' for i in range(n_components)])
    pca_df[group_column] = metadata.loc[expression_data_t.index, group_column].values
    pca_df["Sample"] = expression_data_t.index

    # Plot PCA
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=pca_df,
        x='PC1', y='PC2',
        hue=group_column,
        palette='Set2',
        s=100, edgecolor='k'
    )

    for _, row in pca_df.iterrows():
        plt.text(row['PC1'], row['PC2'], row['Sample'], fontsize=8, ha='right', va='bottom')

    plt.title('PCA Plot of Small RNA-Seq Samples', fontsize=18)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)', fontsize=14)
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(title=group_column, title_fontsize=13, fontsize=12)
    plt.grid(True)
    plt.tight_layout()

    # Save
    plt.savefig(output_file, dpi=300)
    plt.show()
    plt.close()

if __name__ == "__main__":
    # File paths
    #Replace RIGHT_PATH
    expression_file = "RIGHT_PATH/results/counts/norm_counts.csv"  # Sequence-level normalized expression
    metadata_file = "RIGHT_PATH/meta/study_design.csv"

    # Load
    expression_data = pd.read_csv(expression_file, index_col=0)
    metadata = pd.read_csv(metadata_file, index_col=0)

    # Plot using 'antiCCP' or any other metadata column
    plot_small_rna_pca(expression_data, metadata, group_column='antiCCP')


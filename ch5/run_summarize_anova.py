"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
This script performs differential expression analysis using one-way ANOVA
based on a selected grouping variable (e.g., Ethnic group, pathotype, etc.).
It compares gene expression levels across groups and outputs significant results.

INSTRUCTIONS:
1. Prepare a CPM-normalized expression matrix in tab-delimited format
   and place it at the path specified in `count_file`.
   Rows should represent genes and columns should be sample IDs.

2. Prepare a study design CSV file containing metadata for each sample.
   It should include a column corresponding to the grouping variable (e.g., "Ethnic").

3. Update the following parameters in the script as needed:
   - `group_column`: The name of the metadata column for grouping.
   - `control_group`: The reference group for comparisons.
   - `min_samples`: Minimum number of samples required per group.
   - `use_posthoc`: Set to True to perform posthoc pairwise Tukey HSD tests.

4. Run the script using:
       python script_name.py

5. Output will be saved to the path specified in `output_file`.

Ensure that the custom function `summarize_anova_results` is implemented
and accessible via import from `summarize_anova.py`.
"""

import pandas as pd
from summarize_anova import summarize_anova_results  # Function we defined earlier

def main():
    # === File paths ===
    count_file = "results/counts/gene_counts_simulated.txt"     # Your CPM-normalized expression matrix
    design_file = "meta/study_design.csv"             # Your detailed metadata
    output_file = "results/DEGs/anova_results_ethnic.csv"  # Where to save the result

    # === Parameters ===
    group_column = "Ethnic"            # You can also use 'pathotype', 'antiCCP', etc.
    control_group = "White"            # Compare all other ethnic groups to 'White'
    min_samples = 2
    use_posthoc = False                # Set to True if you want pairwise Tukey HSD

    # === Load input data ===
    print("Reading normalized expression matrix...")
    norm_counts = pd.read_csv(count_file, sep='\t', index_col=0)

    print("Reading study design...")
    design = pd.read_csv(design_file)

    # === Run the analysis ===
    print(f"Running differential expression grouped by {group_column}...")
    results = summarize_anova_results(
        norm_counts=norm_counts,
        design=design,
        group_col=group_column,
        control_group=control_group,
        min_samples_per_group=min_samples,
        posthoc=use_posthoc
    )

    # === Save output ===
    print(f"Saving results to {output_file}")
    results.to_csv(output_file, index=False)
    print("Done!")

if __name__ == "__main__":
    main()


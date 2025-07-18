"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
This script performs gene-level one-way ANOVA analysis using normalized gene expression data and metadata.

USAGE:
- Call the `summarize_anova_results` function with the following inputs:
    - `norm_counts`: A pandas DataFrame of normalized gene expression (genes x samples)
    - `design`: A pandas DataFrame of sample metadata that includes a 'runID' column and a grouping variable (e.g., 'Ethnic')
    - Optional parameters:
        - `group_col`: Name of the column in `design` to group by (default = 'Ethnic')
        - `control_group`: A reference group for computing log2 fold-change and Cohen’s d (default = None)
        - `min_samples_per_group`: Minimum required samples per group to include in the analysis (default = 2)
        - `posthoc`: Whether to perform Tukey HSD post hoc tests (default = False)

OUTPUT:
- Returns a DataFrame with the results of one-way ANOVA per gene, including:
    - p-value and FDR-adjusted p-value
    - Eta-squared (effect size)
    - log2 fold change and Cohen’s d (if control group is provided or inferred)
    - Optional: Tukey HSD post hoc test results if `posthoc=True`

NOTE:
- Ensure `runID` values in `design` match the column names in `norm_counts`.
- Designed for exploratory differential expression analysis in transcriptomic studies.
"""

import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from itertools import combinations

def summarize_anova_results(norm_counts, design, group_col='Ethnic', control_group=None, min_samples_per_group=2, posthoc=False):
    """
    Performs one-way ANOVA for each gene and calculates:
    - log2FC (if control group is specified)
    - Cohen's d
    - Eta-squared (η²)
    - Optional Tukey HSD post hoc testing

    Parameters:
        norm_counts (DataFrame): Normalized expression data (genes x samples)
        design (DataFrame): Metadata with sample info
        group_col (str): Grouping variable
        control_group (str): Optional control group for log2FC
        min_samples_per_group (int): Min number of samples required in a group
        posthoc (bool): Whether to perform Tukey HSD pairwise comparisons

    Returns:
        results_df (DataFrame): Summary of differential expression stats
    """
    design['runID'] = design['runID'].astype(str).str.strip()
    design[group_col] = design[group_col].astype(str).str.strip()

    results = []

    for gene in norm_counts.index:
        data = []
        group_data = {}
        sample_sizes = {}

        for group in design[group_col].unique():
            sample_ids = design[design[group_col] == group]['runID']
            sample_ids = [s for s in sample_ids if s in norm_counts.columns]
            if len(sample_ids) >= min_samples_per_group:
                expr = norm_counts.loc[gene, sample_ids]
                data.append(expr)
                group_data[group] = expr
                sample_sizes[group] = len(expr)

        if len(group_data) < 2:
            continue

        # One-way ANOVA
        stat, pval = f_oneway(*data)
        total_mean = pd.concat(group_data.values()).mean()
        ss_between = sum(len(g) * (g.mean() - total_mean)**2 for g in group_data.values())
        ss_total = sum(((val - total_mean) ** 2).sum() for val in group_data.values())
        eta_sq = ss_between / ss_total if ss_total != 0 else np.nan

        # Post hoc comparisons (Tukey HSD)
        if posthoc:
            df_long = pd.DataFrame({
                'Expression': np.concatenate([group_data[g] for g in group_data]),
                'Group': np.concatenate([[g] * len(group_data[g]) for g in group_data])
            })
            tukey = pairwise_tukeyhsd(df_long['Expression'], df_long['Group'], alpha=0.05)
            tukey_df = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])
            tukey_df.insert(0, 'Gene', gene)
            results.extend(tukey_df.to_dict('records'))
        else:
            # Control-based or top-two log2FC + effect size
            if control_group and control_group in group_data:
                for g in group_data:
                    if g == control_group:
                        continue
                    mean_g = group_data[g].mean()
                    mean_c = group_data[control_group].mean()
                    std_pooled = np.sqrt((group_data[g].var() + group_data[control_group].var()) / 2)
                    cohen_d = (mean_g - mean_c) / (std_pooled + 1e-6)
                    log2fc = np.log2((mean_g + 1e-6) / (mean_c + 1e-6))
                    results.append({
                        'Gene': gene,
                        'Group': g,
                        'Baseline': control_group,
                        'log2FC': log2fc,
                        'Cohen_d': cohen_d,
                        'Eta_squared': eta_sq,
                        'p-value': pval
                    })
            else:
                top2 = sorted(group_data.items(), key=lambda x: x[1].mean(), reverse=True)[:2]
                log2fc = np.log2((top2[0][1].mean() + 1e-6) / (top2[1][1].mean() + 1e-6))
                std_pooled = np.sqrt((top2[0][1].var() + top2[1][1].var()) / 2)
                cohen_d = (top2[0][1].mean() - top2[1][1].mean()) / (std_pooled + 1e-6)
                results.append({
                    'Gene': gene,
                    'Group': top2[0][0],
                    'Baseline': top2[1][0],
                    'log2FC': log2fc,
                    'Cohen_d': cohen_d,
                    'Eta_squared': eta_sq,
                    'p-value': pval
                })

    # Final output
    results_df = pd.DataFrame(results)
    if 'p-value' in results_df.columns:
        results_df['adj-p'] = multipletests(results_df['p-value'], method='fdr_bh')[1]

    return results_df


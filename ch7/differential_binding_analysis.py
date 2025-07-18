"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script performs differential binding analysis for ChIP-Seq peak regions across experimental conditions.
It processes BAM files and peak calls (in narrowPeak format), quantifies read counts over merged peak regions, and 
applies statistical testing to identify significantly differentially bound regions between groups (e.g., control vs. treated). 
The output is a CSV file with peak-level statistics, including log2 fold change, p-values, and FDR-adjusted q-values.

Required Python Packages:
- os
- subprocess
- pandas
- numpy (imported at runtime for backward compatibility)
- pybedtools
- statsmodels
- scipy

External Tools Required:
- samtools (for BAM indexing)

Input Files:
- Peak files: `results/peaks/*.narrowPeak` (MACS3 output)
- BAM files: `results/aligned/*.bam` (aligned reads)
- Metadata CSV file: `meta/metadata.csv` with columns `runID` and `condition` (must include values like `control` and `treated`)

Output Files:
- Merged peak regions (used internally as BEDTool object)
- Count matrix (in memory)
- Differential binding results: `differential_binding_results.csv` (contains peak IDs, log2 fold change, p-values, and q-values)

Workflow Summary:
1. Index BAM files using samtools (if not already indexed)
2. Merge all peak regions from input files
3. Count overlapping reads per region per BAM file using `pybedtools`
4. Construct a count matrix with peaks as rows and samples as columns
5. Group samples using metadata and perform two-sample t-tests for each peak
6. Adjust p-values using Benjamini-Hochberg FDR correction
7. Save differential binding results as a CSV file

Note:
- HOMER-style peak annotations and motif discovery can be integrated downstream.
- Ensure samtools and pybedtools dependencies are met and properly configured.
"""

#pip pybedtools statsmodels scipy
import os
import subprocess
import pandas as pd
import pybedtools
from glob import glob
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind
import subprocess
def index_bam_files(bam_files):
    for bam in bam_files:
        bai_file = bam + ".bai"
        if not os.path.exists(bai_file):
            print(f"Indexing BAM file: {bam}")
            subprocess.run(["samtools", "index", bam], check=True)
        else:
            print(f"Index already exists: {bai_file}")

def load_peaks(peak_files):
    all_peaks = []
    for f in peak_files:
        bed = pybedtools.BedTool(f)
        all_peaks.extend(bed)
    merged = pybedtools.BedTool(all_peaks).sort().merge()
    return merged

def count_reads(bam_file, merged_peaks):
    counts = []
    bam = pybedtools.BedTool(bam_file)
    for idx, interval in enumerate(merged_peaks):
        region_str = f"{interval.chrom}\t{interval.start}\t{interval.end}"
        region = pybedtools.BedTool(region_str, from_string=True)
        try:
            coverage = region.coverage(bam).to_dataframe()

            if coverage.empty:
                print(f"[Warning] Empty coverage at region {region_str}")
                counts.append(0)
            elif coverage.shape[1] <= 6:
                print(f"[Warning] Unexpected column count ({coverage.shape[1]}) "
                     f"at region "
                     f"{region_str}")
                counts.append(0)
            else:
                count = int(coverage.iloc[0, 6])
                counts.append(count)
        except Exception as e:
            print(f"[Error] Failed coverage at region {region_str}: {e}")
            counts.append(0)
    return counts

def create_count_matrix(bam_files, merged_peaks):
    matrix = {}
    for bam in bam_files:
        sample_name = os.path.basename(bam).split('.')[0]
        matrix[sample_name] = count_reads(bam, merged_peaks)
    df = pd.DataFrame(matrix)
    df.index = [f"peak_{i+1}" for i in range(len(merged_peaks))]
    return df

def load_metadata(metadata_file):
    df = pd.read_csv(metadata_file)
    return dict(zip(df['runID'].astype(str), df['condition']))

def perform_ttest(count_df, group_labels):
    group1_samples = [s for s in count_df.columns
                  if group_labels.get(s, '') == 'control']
    group2_samples = [s for s in count_df.columns
                  if group_labels.get(s, '') == 'treated']
    
    group1 = count_df[group1_samples]
    group2 = count_df[group2_samples]

    pvals = []
    lfc = []
    for i in range(count_df.shape[0]):
        stat, p = ttest_ind(group1.iloc[i], group2.iloc[i], equal_var=False)
        pvals.append(p)
        log_fc = (group2.iloc[i].mean() + 1) / (group1.iloc[i].mean() + 1)
        lfc.append(pd.np.log2(log_fc))  # Will replace np below

    qvals = multipletests(pvals, method='fdr_bh')[1]
    results = pd.DataFrame({
        'peak_id': count_df.index,
        'log2FoldChange': lfc,
        'pvalue': pvals,
        'qvalue': qvals
    })
    return results

def save_results(results, output_file='differential_binding_results.csv'):
    results.to_csv(output_file, index=False)

def main():
    peak_files = glob("results/peaks/*.narrowPeak")
    bam_files = glob("results/aligned/*.bam")
    metadata_file = "meta/metadata.csv"

    index_bam_files(bam_files)  # Ensure indexing before anything else
    merged_peaks = load_peaks(peak_files)
    count_matrix = create_count_matrix(bam_files, merged_peaks)

    group_labels = load_metadata(metadata_file)
    result = perform_ttest(count_matrix, group_labels)
    save_results(result)

if __name__ == "__main__":
    import numpy as np
    pd.np = np  # For backward compatibility in log2
    main()


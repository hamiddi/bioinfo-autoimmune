"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
Small RNA-Seq Multi-Group Differential Expression Pipeline
----------------------------------------------------------

This Python pipeline processes small RNA-Seq FASTQ files for multi-group differential expression 
analysis using normalized counts and log2 fold-change. It includes adapter trimming, alignment 
with Bowtie1, unique sequence counting, normalization, and ANOVA-based statistical testing 
between groups (e.g., ethnicities).

Requirements:
-------------
- Input FASTQ files named as: data/raw/<runID>.fastq.gz
- Study design CSV file: meta/study_design.csv
    Must include columns 'runID' and a group label (e.g., 'Ethnic')

Tools used (must be installed and in PATH):
-------------------------------------------
- fastp
- bowtie
- samtools
- Python packages: pandas, numpy, scipy, statsmodels, pysam

How to Run:
-----------
1. Prepare the input FASTQ files and the study design CSV.
2. Make sure the Bowtie1 index for miRBase is built and set correctly in `BOWTIE_INDEX`.
3. Run the script from the command line:
       python <this_script>.py
4. Results will be saved to:
   - Trimmed reads: data/processed/
   - Alignment files: data/processed/
   - Counts and normalized counts: results/counts/
   - Differential expression results: results/DEGs/

Output:
-------
- Differentially expressed small RNAs between groups with log2FC and adjusted p-values.
- Top 20 sequences by significance: results/DEGs/top_sequences.csv

"""

# Small RNA-Seq pipeline with log2FC for multi-group comparison
import os
import subprocess
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests

THREADS = 8
BOWTIE_INDEX = "data/ref/mirbase_bowtie"
STUDY_DESIGN_FILE = "meta/study_design.csv"

# Trim reads for small RNA using typical 3' adapter
def trim_small_rna_reads(fq_in, fq_out, sample_id):
    os.makedirs("data/processed", exist_ok=True)
    os.makedirs("results/fastp", exist_ok=True)
    cmd = [
        "fastp", "-i", fq_in, "-o", fq_out,
        "--adapter_sequence=TGGAATTCTCGGGTGCCAAGG",
        "-h", f"results/fastp/{sample_id}.html",
        "-j", f"results/fastp/{sample_id}.json"
    ]
    subprocess.run(cmd, check=True)

# Map with bowtie1 (short reads)
def map_small_rna_reads(fq_trimmed, sample_id):
    sam_file = f"data/processed/{sample_id}.sam"
    cmd = ["bowtie", "-v", "1", "-p", str(THREADS), "-S", BOWTIE_INDEX, fq_trimmed, sam_file]
    subprocess.run(cmd, check=True)
    return sam_file

# Convert to BAM and sort

def convert_and_sort_sam(sam_file):
    bam_file = sam_file.replace(".sam", ".bam")
    sorted_bam = bam_file.replace(".bam", ".sorted.bam")
    subprocess.run(["samtools", "view", "-bS", sam_file, "-o", bam_file], check=True)
    subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam], check=True)
    subprocess.run(["samtools", "index", sorted_bam], check=True)
    return sorted_bam

# Count unique read sequences
def count_unique_sequences(bam_files, output_file):
    from collections import Counter
    import pysam

    count_dict = {}
    for bam in bam_files:
        sample = os.path.basename(bam).split(".")[0]
        samfile = pysam.AlignmentFile(bam, "rb")
        seq_counter = Counter()
        for read in samfile.fetch():
            if not read.is_unmapped:
                seq = read.query_sequence
                seq_counter[seq] += 1
        samfile.close()
        count_dict[sample] = seq_counter

    all_seqs = sorted(set(seq for d in count_dict.values() for seq in d))
    df = pd.DataFrame(index=all_seqs)
    for sample, counter in count_dict.items():
        df[sample] = [counter.get(seq, 0) for seq in all_seqs]

    df.to_csv(output_file)
    return df

# Normalize and filter low-expression sequences
def normalize_and_filter(counts, design_file, min_cpm=1, min_samples=3):
    design = pd.read_csv(design_file)
    counts = counts.loc[:, counts.columns.intersection(design['runID'])]
    lib_sizes = counts.sum(axis=0)
    cpm = counts.divide(lib_sizes, axis=1) * 1e6
    keep = (cpm > min_cpm).sum(axis=1) >= min_samples
    filtered = counts[keep]
    norm_counts = filtered.divide(filtered.sum(axis=0), axis=1) * 1e6
    return norm_counts, design

# Multi-group ANOVA with log2FC vs control or top2 groups
def differential_expression(norm_counts, design, group_col='Ethnic', control_group='White', min_samples_per_group=2):
    design['runID'] = design['runID'].astype(str).str.strip()
    design[group_col] = design[group_col].astype(str).str.strip()
    groups = design[group_col].dropna().unique()

    results = []
    for seq in norm_counts.index:
        seq_data = []
        group_means = {}
        group_exprs = {}
        for g in groups:
            sample_ids = design[design[group_col] == g]['runID']
            sample_ids = [s for s in sample_ids if s in norm_counts.columns]
            if len(sample_ids) >= min_samples_per_group:
                expr = norm_counts.loc[seq, sample_ids]
                seq_data.append(expr)
                group_exprs[g] = expr
                group_means[g] = expr.mean()

        if len(seq_data) >= 2:
            stat, pval = f_oneway(*seq_data)
            if control_group and control_group in group_exprs:
                for g in group_exprs:
                    if g == control_group:
                        continue
                    log2fc = np.log2((group_means[g] + 1e-6) / (group_means[control_group] + 1e-6))
                    results.append((seq, g, control_group, log2fc, pval))
            else:
                top2 = sorted(group_means.items(), key=lambda x: x[1], reverse=True)[:2]
                log2fc = np.log2((top2[0][1] + 1e-6) / (top2[1][1] + 1e-6))
                results.append((seq, top2[0][0], top2[1][0], log2fc, pval))

    df = pd.DataFrame(results, columns=['Sequence', 'Group', 'Baseline', 'log2FC', 'p-value'])
    df['adj-p'] = multipletests(df['p-value'], method='fdr_bh')[1]
    return df.sort_values(['Sequence', 'adj-p'])

# Main execution

def run_small_rna_log2fc_pipeline():
    os.makedirs("results/DEGs", exist_ok=True)
    os.makedirs("results/counts", exist_ok=True)
    os.makedirs("results/fastp", exist_ok=True)

    study_design = pd.read_csv(STUDY_DESIGN_FILE)
    bam_files = []

    for _, row in study_design.iterrows():
        sid = row['runID']
        fq_in = f"data/raw/{sid}.fastq.gz"
        fq_out = f"data/processed/{sid}.trim.fastq.gz"
        trim_small_rna_reads(fq_in, fq_out, sid)
        sam = map_small_rna_reads(fq_out, sid)
        sorted_bam = convert_and_sort_sam(sam)
        bam_files.append(sorted_bam)

    counts_df = count_unique_sequences(bam_files, "results/counts/seq_counts.csv")
    norm_counts, design = normalize_and_filter(counts_df, STUDY_DESIGN_FILE)
    norm_counts.to_csv("results/counts/norm_counts.csv")
    degs = differential_expression(norm_counts, design, group_col='Ethnic', control_group='White')
    degs.to_csv("results/DEGs/diff_expr_log2fc.csv", index=False)
    top = degs.head(20)['Sequence'].tolist()
    pd.DataFrame(top, columns=['TopSequences']).to_csv("results/DEGs/top_sequences.csv", index=False)

if __name__ == "__main__":
    run_small_rna_log2fc_pipeline()


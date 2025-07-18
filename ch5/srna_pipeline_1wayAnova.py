"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################

Small RNA-Seq ANOVA Pipeline
----------------------------

This script implements a complete pipeline for small RNA sequencing analysis,
including preprocessing, mapping, quantification, normalization, and 
differential expression analysis using one-way ANOVA.

**Usage Instructions:**
1. Prepare your input files:
   - Raw FASTQ files must be placed in `data/raw/` and named as <runID>.fastq.gz.
   - A metadata file `meta/study_design.csv` is required with a column named 'runID'
     matching the FASTQ file names and a grouping column (e.g., 'Ethnic').

2. Set the Bowtie index path in the `BOWTIE_INDEX` variable.
   - The index should be built from mature microRNA sequences (e.g., from miRBase).

3. Ensure required tools are installed and available in your environment:
   - fastp, bowtie, samtools, and required Python packages (pandas, numpy, pysam, scipy, statsmodels).

4. Run the script from the command line:
   python script_name.py

**Output:**
- Trimmed reads in `data/processed/`
- Mapping results and BAM files in `data/processed/`
- Quality reports in `results/fastp/`
- Raw and normalized count tables in `results/counts/`
- Differential expression results (p-value, adjusted p-value, log2FC) in `results/DEGs/`
"""

# Small RNA-Seq ANOVA Pipeline
import os
import subprocess
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests

THREADS = 8
BOWTIE_INDEX = "data/ref/mirbase_bowtie"
STUDY_DESIGN_FILE = "meta/study_design.csv"

# Adapter trimming for small RNA reads
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

# Mapping with Bowtie1
def map_small_rna_reads(fq_trimmed, sample_id):
    output_sam = f"data/processed/{sample_id}.sam"
    cmd = ["bowtie", "-v", "1", "-p", str(THREADS), "-S", BOWTIE_INDEX, fq_trimmed, output_sam]
    subprocess.run(cmd, check=True)
    return output_sam

# Convert and sort SAM to BAM
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

# One-way ANOVA across groups
def differential_expression(norm_counts, design, group_col='Ethnic', control_group=None, min_samples_per_group=2):
    design['runID'] = design['runID'].astype(str).str.strip()
    design[group_col] = design[group_col].astype(str).str.strip()
    groups = design[group_col].dropna().unique()

    results = []
    for seq in norm_counts.index:
        seq_data = []
        group_means = {}
        for g in groups:
            sample_ids = design[design[group_col] == g]['runID']
            sample_ids = [s for s in sample_ids if s in norm_counts.columns]
            if len(sample_ids) >= min_samples_per_group:
                expr = norm_counts.loc[seq, sample_ids]
                seq_data.append(expr)
                group_means[g] = expr.mean()
        if len(seq_data) >= 2:
            stat, pval = f_oneway(*seq_data)
            if control_group and control_group in group_means:
                top_group = max(group_means.items(), key=lambda x: x[1])[0]
                if top_group == control_group:
                    top_group = sorted(group_means.items(), key=lambda x: x[1], reverse=True)[1][0]
                log2fc = np.log2((group_means[top_group] + 1e-6) / (group_means[control_group] + 1e-6))
            else:
                top2 = sorted(group_means.items(), key=lambda x: x[1], reverse=True)[:2]
                log2fc = np.log2((top2[0][1] + 1e-6) / (top2[1][1] + 1e-6))
            results.append((seq, pval, log2fc))

    df = pd.DataFrame(results, columns=['Sequence', 'p-value', 'log2FC'])
    df['adj-p'] = multipletests(df['p-value'], method='fdr_bh')[1]
    return df.sort_values('adj-p')

# Main runner
def run_small_rna_anova_pipeline():
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
    degs = differential_expression(norm_counts, design, group_col='Ethnic')
    degs.to_csv("results/DEGs/diff_expr_anova.csv", index=False)
    top = degs.head(20)['Sequence'].tolist()
    pd.DataFrame(top, columns=['TopSequences']).to_csv("results/DEGs/top_sequences.csv", index=False)

if __name__ == "__main__":
    run_small_rna_anova_pipeline()


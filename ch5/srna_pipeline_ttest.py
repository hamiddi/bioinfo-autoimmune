"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
Small RNA-Seq Analysis Pipeline
-------------------------------
This Python script performs a complete small RNA-Seq data analysis workflow using raw FASTQ files and a study design table. 
It is adapted from a standard mRNA-Seq pipeline and tailored for microRNA and other small RNA quantification.

Pipeline Steps:
1. Adapter trimming using `fastp` with a standard small RNA 3' adapter.
2. Alignment to reference using `bowtie1` (no splicing).
3. SAM to sorted BAM conversion and indexing using `samtools`.
4. Counting unique small RNA sequences from aligned reads.
5. Normalization and filtering based on counts per million (CPM).
6. Differential expression analysis using t-tests and FDR correction.

Usage Instructions:
- Place your raw FASTQ files in `data/raw/` with filenames matching the `runID` column in the study design CSV.
- Ensure `meta/study_design.csv` includes a column named `runID` and a group column (e.g., "antiCCP").
- Reference files must include:
    - Bowtie index located at `data/ref/mirbase_bowtie.*`
    - 3' adapter sequence file at `data/ref/illumina_3prime_adapter.fa` (optional here; fastp uses inline adapter)
- Required tools: `fastp`, `bowtie`, `samtools`, `pysam`, and Python libraries (`pandas`, `numpy`, `scikit-learn`, `statsmodels`).
- Run the script with: `python <script_name>.py`

Output:
- Trimmed FASTQ files: `data/processed/*.trim.fastq.gz`
- Alignment BAM files: `data/processed/*.sorted.bam`
- Count matrix: `results/counts/seq_counts.csv`
- Normalized counts: `results/counts/norm_counts.csv`
- Differential expression results: `results/DEGs/diff_expr_smallRNA.csv`
- Top sequences: `results/DEGs/top_seqs.csv`
"""

# Small RNA-Seq Pipeline (Adapted from mRNA-Seq version)
import os
import subprocess
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# Configuration
THREADS = 8
BOWTIE_INDEX = "data/ref/mirbase_bowtie"
FASTA_ADAPTER = "data/ref/illumina_3prime_adapter.fa"
STUDY_DESIGN_FILE = "meta/study_design.csv"

# Step 1: Adapter Trimming (for small RNAs)
def trim_small_rna_reads(fq_in, fq_out, sample_id):
    os.makedirs("data/processed", exist_ok=True)
    os.makedirs("results/fastp", exist_ok=True)
    cmd = [
        "fastp", "-i", fq_in, "-o", fq_out,
        "--adapter_sequence=TGGAATTCTCGGGTGCCAAGG",  # Typical 3' adapter for small RNA
        "-h", f"results/fastp/{sample_id}.html",
        "-j", f"results/fastp/{sample_id}.json"
    ]
    subprocess.run(cmd, check=True)

# Step 2: Mapping with Bowtie1 (no splicing)
def map_small_rna_reads(fq_trimmed, sample_id):
    output_sam = f"data/processed/{sample_id}.sam"
    cmd = [
        "bowtie", "-v", "1", "-p", str(THREADS), "-S", BOWTIE_INDEX, fq_trimmed, output_sam
    ]
    subprocess.run(cmd, check=True)
    return output_sam

# Step 3: Convert SAM to BAM and sort

def convert_and_sort_sam(sam_file):
    bam_file = sam_file.replace(".sam", ".bam")
    sorted_bam = bam_file.replace(".bam", ".sorted.bam")
    subprocess.run(["samtools", "view", "-bS", sam_file, "-o", bam_file], check=True)
    subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam], check=True)
    subprocess.run(["samtools", "index", sorted_bam], check=True)
    return sorted_bam

# Step 4: Count unique read sequences (collapsed reads)
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

    # Create DataFrame
    all_seqs = sorted(set(seq for d in count_dict.values() for seq in d))
    df = pd.DataFrame(index=all_seqs)
    for sample, counter in count_dict.items():
        df[sample] = [counter.get(seq, 0) for seq in all_seqs]

    df.to_csv(output_file)
    return df

# Step 5: Normalize & Filter
def normalize_and_filter(counts, design_file, min_cpm=1, min_samples=3):
    design = pd.read_csv(design_file)
    counts = counts.loc[:, counts.columns.intersection(design['runID'])]
    lib_sizes = counts.sum(axis=0)
    cpm = counts.divide(lib_sizes, axis=1) * 1e6
    keep = (cpm > min_cpm).sum(axis=1) >= min_samples
    filtered = counts[keep]
    norm_counts = filtered.divide(filtered.sum(axis=0), axis=1) * 1e6
    return norm_counts, design

# Step 6: Differential Expression
def differential_expression(norm_counts, design, group_col='antiCCP'):
    group1_ids = design[design[group_col].str.lower() == 'positive']['runID']
    group2_ids = design[design[group_col].str.lower() == 'negative']['runID']

    group1 = [s for s in group1_ids if s in norm_counts.columns]
    group2 = [s for s in group2_ids if s in norm_counts.columns]

    results = []
    for seq in norm_counts.index:
        expr1 = norm_counts.loc[seq, group1]
        expr2 = norm_counts.loc[seq, group2]

        if expr1.var() == 0 and expr2.var() == 0:
            results.append((seq, np.nan, np.nan))
            continue

        stat, pval = ttest_ind(expr1, expr2, equal_var=False, nan_policy='omit')
        log2fc = np.log2((expr1.mean() + 1e-6) / (expr2.mean() + 1e-6))
        results.append((seq, log2fc, pval))

    df = pd.DataFrame(results, columns=['Sequence', 'log2FC', 'p-value'])
    df.dropna(subset=['p-value'], inplace=True)
    df['adj-p'] = multipletests(df['p-value'], method='fdr_bh')[1]
    return df.sort_values('adj-p')

# Main Pipeline

def run_small_rna_pipeline():
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

    degs = differential_expression(norm_counts, design)
    degs.to_csv("results/DEGs/diff_expr_smallRNA.csv", index=False)
    top_seqs = degs.head(20)['Sequence'].tolist()
    pd.DataFrame(top_seqs, columns=['TopSequences']).to_csv("results/DEGs/top_seqs.csv", index=False)

if __name__ == "__main__":
    run_small_rna_pipeline()


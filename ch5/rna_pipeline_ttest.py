"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################

RNA-Seq Analysis Pipeline with Gene Symbol Annotation and Differential Expression

This Python script implements a complete RNA-Seq data analysis workflow, including:
- Reference genome and annotation file download
- Trimming and quality control of raw FASTQ files using fastp
- Alignment with STAR aligner
- BAM indexing with SAMtools
- Gene-level quantification with featureCounts
- Low-expression filtering and CPM normalization
- Differential expression analysis between two groups using t-test
- Gene symbol annotation and result export

USAGE:
1. Place the paired-end raw FASTQ files in the `data/raw/` directory, named as:
   <sample_id>_1.fastq.gz and <sample_id>_2.fastq.gz
2. Ensure a metadata file `meta/study_design.csv` exists with a 'runID' column matching sample IDs and a grouping column (e.g., 'antiCCP').
3. Run the script:
   - Set `simulation_mode=False` for full pipeline execution (download, trimming, alignment, counting, etc.)
   - Set `simulation_mode=True` to skip alignment and use precomputed counts (e.g., for debugging or simulated data)

OUTPUT:
- Gene count matrix: results/counts/gene_counts.txt
- CPM-normalized expression matrix: results/counts/norm_counts.csv
- Differential expression results: results/DEGs/diff_expr_with_symbols.csv
- Top 20 differentially expressed genes: results/DEGs/top_genes.csv
- Top 100 DEGs with gene symbols: results/DEGs/top100_genes.csv

DEPENDENCIES:
External tools required (must be in PATH):
- fastp
- STAR
- samtools
- featureCounts

Python packages required:
- pandas, numpy, matplotlib, seaborn, sklearn, statsmodels, scipy, requests

IMPORTANT:
- The first time you run the script, the reference genome and GTF file will be downloaded and indexed.
- You can comment/uncomment sections inside the `run_pipeline()` function depending on which stages to run.
"""

#include gene symbols and additional count file
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import requests
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests

# Configuration
THREADS = 8
REFERENCE_URL = "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
REFERENCE_FA = "data/ref/hg38.fa"
STAR_INDEX = "data/ref/star_index"
GTF_FILE = "data/ref/genes.gtf"
STUDY_DESIGN_FILE = "meta/study_design.csv"

def trim_reads(fq1, fq2, out1, out2, sample_id):
    os.makedirs("data/processed", exist_ok=True)
    os.makedirs("results/fastp", exist_ok=True)
    cmd = [
        "fastp", "-i", fq1, "-I", fq2, "-o", out1, "-O", out2,
        "-h", f"results/fastp/{sample_id}.html",
        "-j", f"results/fastp/{sample_id}.json"
    ]
    subprocess.run(cmd, check=True)

def prepare_reference():
    os.makedirs("data/ref", exist_ok=True)
    if not os.path.exists(REFERENCE_FA):
        command = (
                 "wget -O - "
                 f"{REFERENCE_URL} "
                 "--no-check-certificate "
                 "| gunzip -c "
                 f"> {REFERENCE_FA}"
              )
        subprocess.run(
           command,
           shell=True,
           check=True
        )
    if not os.path.exists(GTF_FILE):
        gtf_url = (
         "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/"
         "Homo_sapiens.GRCh38.110.gtf.gz"
        )
        command = (
          "wget -O data/ref/genes.gtf.gz "
          f"{gtf_url} "
          "--no-check-certificate"
        )

        subprocess.run(
          command,
          shell=True,
          check=True
        )
        subprocess.run("gunzip -f data/ref/genes.gtf.gz", shell=True, check=True)
    if not os.path.exists(os.path.join(STAR_INDEX, "genomeParameters.txt")):
        os.makedirs(STAR_INDEX, exist_ok=True)
        cmd = [
            "STAR", "--runThreadN", str(THREADS),
            "--runMode", "genomeGenerate",
            "--genomeDir", STAR_INDEX,
            "--genomeFastaFiles", REFERENCE_FA,
            "--sjdbGTFfile", GTF_FILE,
            "--sjdbOverhang", "100"
        ]
        subprocess.run(cmd, check=True)

def map_reads(fq1, fq2, sample_id):
    output_prefix = f"data/processed/{sample_id}_"
    cmd = [
        "STAR", "--runThreadN", str(THREADS),
        "--genomeDir", STAR_INDEX,
        "--readFilesIn", fq1, fq2,
        "--readFilesCommand", "zcat",
        "--outFileNamePrefix", output_prefix,
        "--outSAMtype", "BAM", "SortedByCoordinate"
    ]
    subprocess.run(cmd, check=True)

def index_bam(bam_file):
    subprocess.run(["samtools", "index", bam_file])

def count_reads(bam_files, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    cmd = [
        "featureCounts", "-T", str(THREADS), "-p",
        "-a", GTF_FILE, "-o", output_file
    ] + bam_files
    subprocess.run(cmd, check=True)

def filter_low_expression_genes(counts, min_cpm=1, min_samples=3):
    """
    Filters out genes with low expression across samples.

    Parameters:
        counts (pd.DataFrame): Raw count matrix (genes x samples).
        min_cpm (float): Minimum CPM threshold.
        min_samples (int): Minimum number of samples where gene must exceed CPM.

    Returns:
        pd.DataFrame: Filtered count matrix.
    """
    lib_sizes = counts.sum(axis=0)
    cpm = counts.divide(lib_sizes, axis=1) * 1e6
    keep = (cpm > min_cpm).sum(axis=1) >= min_samples
    filtered = counts[keep]
    print(f"[INFO] Filtered from {counts.shape[0]} to {filtered.shape[0]} genes based on CPM > {min_cpm} in ≥ {min_samples} samples.")
    return filtered

def normalize_counts(count_matrix_path, design_file):
    counts = pd.read_csv(count_matrix_path, sep='\t', comment='#', index_col=0)
    counts = counts.iloc[:, 5:]
    counts.columns = [col.split('/')[-1].replace('_Aligned.sortedByCoord.out.bam', '') for col in counts.columns]
    study_design = pd.read_csv(design_file)
    available_samples = [s for s in study_design['runID'] if s in counts.columns]
    counts = counts[available_samples]
    study_design = study_design[study_design['runID'].isin(available_samples)]
    #norm_counts = counts.div(counts.sum(axis=0), axis=1) * 1e6
    #return norm_counts, study_design
    counts = filter_low_expression_genes(counts)
    norm_counts = counts.div(counts.sum(axis=0), axis=1) * 1e6
    return norm_counts, study_design

def normalize_counts_sim(count_matrix_path, design_file):
    #no need to skip the first 5 lines created by featureCount
    """
    Normalizes a simulated RNA-Seq count matrix using counts per million (CPM),
    and aligns it with the study design.

    Parameters:
        count_matrix_path (str): Path to the raw simulated count matrix (genes x samples).
        design_file (str): Path to the study design CSV with a 'runID' column.

    Returns:
        norm_counts (pd.DataFrame): CPM-normalized expression matrix.
        study_design (pd.DataFrame): Filtered design matching available samples.
    """
    # Load data
    counts = pd.read_csv(count_matrix_path, sep='\t', index_col=0)
    study_design = pd.read_csv(design_file)

    # Standardize sample ID formatting
    counts.columns = counts.columns.astype(str).str.strip()
    study_design['runID'] = study_design['runID'].astype(str).str.strip()

    # Match samples present in both count matrix and study design
    available_samples = [s for s in study_design['runID'] if s in counts.columns]
    counts = counts[available_samples]
    study_design = study_design[study_design['runID'].isin(available_samples)]

    # Normalize using CPM (Counts Per Million)
    counts = filter_low_expression_genes(counts)
    norm_counts = counts.div(counts.sum(axis=0), axis=1) * 1e6
    return norm_counts, study_design

def load_gene_symbols(gtf_file):
    symbols = {}
    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] != 'gene':
                continue
            info = fields[8]
            gene_id = None
            gene_name = None
            for entry in info.strip().split(';'):
                entry = entry.strip()
                if entry.startswith("gene_id"):
                    gene_id = entry.split(' ')[1].replace('"', '')
                elif entry.startswith("gene_name"):
                    gene_name = entry.split(' ')[1].replace('"', '')
            if gene_id and gene_name:
                symbols[gene_id] = gene_name
    return symbols


def differential_expression(norm_counts, design, group_col='antiCCP'):
    '''
    Perform differential expression analysis between two groups.

    Parameters:
        norm_counts (pd.DataFrame): Gene expression matrix (genes x samples).
        design (pd.DataFrame): Metadata with at least 'runID' and grouping column.
        group_col (str): Column name in design that defines the group (e.g., 'antiCCP').

    Returns:
        pd.DataFrame: Table with Gene, log2 Fold Change, p-value, and adjusted p-value.
    '''
    group1_ids = design[design[group_col].str.lower() == 'positive']['runID']
    group2_ids = design[design[group_col].str.lower() == 'negative']['runID']

    group1 = [s for s in group1_ids if s in norm_counts.columns]
    group2 = [s for s in group2_ids if s in norm_counts.columns]

    if not group1 or not group2:
        raise ValueError("One of the groups has no samples in the expression matrix.")

    results = []
    for gene in norm_counts.index:
        expr1 = norm_counts.loc[gene, group1]
        expr2 = norm_counts.loc[gene, group2]

        if expr1.var() == 0 and expr2.var() == 0:
            # Both groups have no variation → skip
            results.append((gene, np.nan, np.nan))
            continue

        stat, pval = ttest_ind(expr1, expr2, equal_var=False, nan_policy='omit')
        log2fc = np.log2((expr1.mean() + 1e-6) / (expr2.mean() + 1e-6))
        results.append((gene, log2fc, pval))

    df = pd.DataFrame(results, columns=['Gene', 'log2FC', 'p-value'])

    # Drop rows where p-value is missing (e.g. due to t-test failure)
    missing = df['p-value'].isna().sum()
    if missing > 0:
        print(f"[WARNING] {missing} genes had NaN p-values and were excluded from adjustment.")
        df = df.dropna(subset=['p-value'])

    # Adjust p-values
    df['adj-p'] = multipletests(df['p-value'], method='fdr_bh')[1]

    return df.sort_values('adj-p')


def export_counts_with_symbols(count_matrix_path, gene_symbol_map, output_path):
    counts = pd.read_csv(count_matrix_path, sep='\t', comment='#', index_col=0)
    counts = counts.iloc[:, 5:]

    # Clean sample column names (e.g., ERR9539375)
    cleaned_columns = [col.split('/')[-1].replace('_Aligned.sortedByCoord.out.bam', '') for col in counts.columns]
    counts.columns = cleaned_columns

    # Add gene symbols
    counts.reset_index(inplace=True)
    counts.rename(columns={'Geneid': 'accession'}, inplace=True)
    counts['gene_symbol'] = counts['accession'].map(gene_symbol_map)

    # Reorder columns
    ordered_cols = ['accession', 'gene_symbol'] + [col for col in counts.columns if col not in ['accession', 'gene_symbol']]
    counts = counts[ordered_cols]

    # Save full file (including genes without symbols)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    counts.to_csv(output_path, index=False)

    # Filter to only rows with valid gene symbols
    filtered = counts.dropna(subset=['gene_symbol'])
    filtered_path = output_path.replace('.csv', '_with_symbols_only.csv')
    filtered.to_csv(filtered_path, index=False)
    print(f"[INFO] Saved counts with all genes: {output_path}")
    print(f"[INFO] Saved counts with gene symbols only: {filtered_path}")


def run_pipeline(simulation_mode=False):
    os.makedirs("results/DEGs", exist_ok=True)
    os.makedirs("results/counts", exist_ok=True)
    os.makedirs("results/fastp", exist_ok=True)

    gene_symbols = load_gene_symbols(GTF_FILE)

    if simulation_mode:
        print("[SIMULATION] Skipping alignment and preprocessing.")
        study_design = pd.read_csv(STUDY_DESIGN_FILE)
        norm_counts, design = normalize_counts_sim("results/counts/gene_counts_simulated.txt", STUDY_DESIGN_FILE)
    else:
        #prepare_reference()
        study_design = pd.read_csv(STUDY_DESIGN_FILE)
        ###     
        for _, row in study_design.iterrows():
            sid = row['runID']
            fq1 = f"data/raw/{sid}_1.fastq.gz"
            fq2 = f"data/raw/{sid}_2.fastq.gz"
            out1 = f"data/processed/{sid}_1.trim.fastq.gz"
            out2 = f"data/processed/{sid}_2.trim.fastq.gz"
            trim_reads(fq1, fq2, out1, out2, sid)
            map_reads(out1, out2, sid)
            index_bam(f"data/processed/{sid}_Aligned.sortedByCoord.out.bam")
        ###
        #bam_files = [f"data/processed/{sid}_Aligned.sortedByCoord.out.bam" for sid in study_design['runID']]
        bam_files = []
        for sid in study_design['runID']:
           bam_path = f"data/processed/{sid}_Aligned.sortedByCoord.out.bam"
           if os.path.exists(bam_path):
              bam_files.append(bam_path)
           else:
              print(f"[WARNING] Missing BAM file for {sid}, skipping.")
        count_reads(bam_files, "results/counts/gene_counts.txt")
        export_counts_with_symbols("results/counts/gene_counts.txt", gene_symbols, 
                                  "results/counts/gene_counts_with_symbols.csv")
        norm_counts, design = normalize_counts("results/counts/gene_counts.txt", STUDY_DESIGN_FILE)
        norm_counts.to_csv("results/counts/norm_counts.csv")
        # Create version with gene symbols
        norm_counts_with_symbols = norm_counts.copy()
        norm_counts_with_symbols['gene_symbol'] = norm_counts_with_symbols.index.map(gene_symbols)
        norm_counts_with_symbols = norm_counts_with_symbols.dropna(subset=['gene_symbol'])
        norm_counts_with_symbols = norm_counts_with_symbols.set_index('gene_symbol')
        # Save filtered version
        norm_counts_with_symbols.to_csv("results/counts/norm_counts_with_symbols.csv")

    #degs = differential_expression(norm_counts, design)
    degs = differential_expression(norm_counts, design)

    degs['GeneSymbol'] = degs['Gene'].map(gene_symbols)
    degs.to_csv("results/DEGs/diff_expr_with_symbols.csv", index=False)

    #top_genes = degs.head(20)['Gene'].tolist()
    # Select top 20 genes that have gene symbols
    top_degs = degs.dropna(subset=['GeneSymbol']).head(20)
    top_genes = top_degs['GeneSymbol'].tolist()
    # Save to CSV for downstream plotting
    pd.DataFrame(top_genes, columns=["Gene"]).to_csv("results/DEGs/top_genes.csv", index=False)
    top_100_genes = degs.sort_values('adj-p').dropna(subset=['GeneSymbol']).head(100)
    top_100_genes['GeneSymbol'].to_csv("results/DEGs/top100_genes.csv", index=False)

    #significant_genes = degs[degs['adj-p'] < 0.1]['Gene'].tolist()
    #significant_genes.to_csv("results/DEGs/diff_significant_genes.csv", index=False)
if __name__ == "__main__":
    run_pipeline(simulation_mode=False)


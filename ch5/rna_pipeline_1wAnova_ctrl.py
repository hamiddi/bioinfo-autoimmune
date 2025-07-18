"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
RNA-Seq Differential Expression Pipeline with Gene Symbol Annotation
This script performs an end-to-end RNA-Seq data analysis pipeline including:
- Reference genome preparation
- Trimming of raw paired-end FASTQ reads
- Read alignment with STAR
- Gene-level quantification using featureCounts
- Low-expression gene filtering
- CPM normalization
- Differential expression analysis using ANOVA (multi-group comparison)
- Export of results with gene symbols
- Visualization-ready output files

Usage Instructions:
1. Place raw paired-end FASTQ files in `data/raw/`, named as `<sampleID>_1.fastq.gz` and `<sampleID>_2.fastq.gz`.
2. Ensure your study design file is located at `meta/study_design.csv` with a column named `runID` and a grouping column (e.g., 'Ethnic').
3. Set `simulation_mode=True` in `run_pipeline()` to test the pipeline using precomputed counts.
4. Set `simulation_mode=False` to run the full pipeline on real FASTQ data (requires STAR, samtools, fastp, and featureCounts installed and in PATH).
5. Output directories are created automatically:
   - `results/fastp`: Quality reports from fastp
   - `data/processed`: Aligned and trimmed files
   - `results/counts`: Count matrices and normalized results
   - `results/DEGs`: Differentially expressed genes and top gene lists

Dependencies:
- Python packages: pandas, matplotlib, seaborn, numpy, sklearn, statsmodels, scipy, requests
- External tools: fastp, STAR, samtools, featureCounts (Subread)

Make sure required tools are installed and accessible from the command line before running.
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
        subprocess.run(f"wget -O - {REFERENCE_URL} --no-check-certificate | gunzip -c > {REFERENCE_FA}",
                      shell=True, check=True)
    if not os.path.exists(GTF_FILE):
        gtf_url = "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
        subprocess.run(f"wget -O data/ref/genes.gtf.gz {gtf_url} --no-check-certificate", shell=True, check=True)
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
    #norm_counts = counts.div(counts.sum(axis=0), axis=1) * 1e6
    #return norm_counts, study_design
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


# Modified Differential Expression for multi-group ANOVA with log2FC between top two groups
def differential_expression(norm_counts, design, group_col='Ethnic', control_group="White", min_samples_per_group=2):
    """
    Perform one-way ANOVA for DE across multiple groups.
    If control_group is provided, computes log2FC for all other groups vs. control.

    Parameters:
        norm_counts (pd.DataFrame): Gene expression matrix (genes x samples).
        design (pd.DataFrame): Metadata with 'runID' and group column.
        group_col (str): Column in design defining groups.
        control_group (str): Control group to compare others against.
        min_samples_per_group (int): Minimum samples per group to be included.

    Returns:
        pd.DataFrame: Gene-wise p-values and log2FCs vs control (if specified).
    """
    design['runID'] = design['runID'].astype(str).str.strip()
    design[group_col] = design[group_col].astype(str).str.strip()

    groups = design[group_col].dropna().unique()
    if len(groups) < 2:
        raise ValueError("Not enough groups for DE.")

    results = []

    for gene in norm_counts.index:
        gene_data = []
        group_means = {}
        group_exprs = {}

        for g in groups:
            sample_ids = design[design[group_col] == g]['runID']
            sample_ids = [s for s in sample_ids if s in norm_counts.columns]
            if len(sample_ids) >= min_samples_per_group:
                expr = norm_counts.loc[gene, sample_ids]
                gene_data.append(expr)
                group_exprs[g] = expr
                group_means[g] = expr.mean()

        if len(gene_data) >= 2:
            stat, pval = f_oneway(*gene_data)

            if control_group and control_group in group_exprs:
                for g in group_exprs:
                    if g == control_group:
                        continue
                    # Compute log2FC: group vs control
                    mean_control = group_means[control_group] + 1e-6
                    mean_group = group_means[g] + 1e-6
                    log2fc = np.log2(mean_group / mean_control)
                    results.append((gene, g, control_group, log2fc, pval))
            else:
                # Default: just compare top two
                top2 = sorted(group_means.items(), key=lambda x: x[1], reverse=True)[:2]
                log2fc = np.log2((top2[0][1] + 1e-6) / (top2[1][1] + 1e-6))
                results.append((gene, top2[0][0], top2[1][0], log2fc, pval))

    if not results:
        raise ValueError("No genes passed filtering criteria.")

    df = pd.DataFrame(results, columns=['Gene', 'Group', 'Baseline', 'log2FC', 'p-value'])
    df['adj-p'] = multipletests(df['p-value'], method='fdr_bh')[1]
    return df.sort_values(['Gene', 'adj-p'])


# In run_pipeline(), update the call to differential_expression:
# degs = differential_expression(norm_counts, design)
# Change to:
# degs = differential_expression(norm_counts, design, group_col='Ethnic')

def export_counts_with_symbols(count_matrix_path, gene_symbol_map, output_path):
    """
    Export raw count matrix with added gene symbols.

    Parameters:
        count_matrix_path (str): Path to raw featureCounts output.
        gene_symbol_map (dict): Ensembl ID to gene symbol mapping.
        output_path (str): Destination CSV file path.
    """
    counts = pd.read_csv(count_matrix_path, sep='\t', comment='#', index_col=0)
    counts = counts.iloc[:, 5:]  # drop annotation columns
    counts.reset_index(inplace=True)
    counts.rename(columns={'Geneid': 'accession'}, inplace=True)
    counts['gene_symbol'] = counts['accession'].map(gene_symbol_map)
    ordered_cols = ['accession', 'gene_symbol'] + [col for col in counts.columns if col not in ['accession', 'gene_symbol']]
    counts = counts[ordered_cols]
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    counts.to_csv(output_path, index=False)
    print(f"[INFO] Saved count matrix with gene symbols to {output_path}")

def run_pipeline(simulation_mode=True):
    os.makedirs("results/DEGs", exist_ok=True)
    os.makedirs("results/plots", exist_ok=True)
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
        bam_files = [f"data/processed/{sid}_Aligned.sortedByCoord.out.bam" for sid in study_design['runID']]
        count_reads(bam_files, "results/counts/gene_counts.txt")
        export_counts_with_symbols("results/counts/gene_counts.txt", gene_symbols, 
                                  "results/counts/gene_counts_with_symbols.csv")
        norm_counts, design = normalize_counts("results/counts/gene_counts.txt", STUDY_DESIGN_FILE)

        # Create version with gene symbols
        norm_counts_with_symbols = norm_counts.copy()
        norm_counts_with_symbols['gene_symbol'] = norm_counts_with_symbols.index.map(gene_symbols)
        norm_counts_with_symbols = norm_counts_with_symbols.dropna(subset=['gene_symbol'])
        norm_counts_with_symbols = norm_counts_with_symbols.set_index('gene_symbol')
        # Save filtered version
        norm_counts_with_symbols.to_csv("results/counts/norm_counts_with_symbols.csv")



    #degs = differential_expression(norm_counts, design)
    degs = differential_expression(norm_counts, design, group_col='Ethnic')
    degs['GeneSymbol'] = degs['Gene'].map(gene_symbols)
    degs.to_csv("results/DEGs/diff_expr_with_symbols_log2fc.csv", index=False)

    top_degs = degs.dropna(subset=['GeneSymbol']).head(20)
    top_genes = top_degs['GeneSymbol'].tolist()
    # Save to CSV for downstream plotting
    pd.DataFrame(top_genes, columns=["Gene"]).to_csv("results/DEGs/top_genes.csv", index=False)
    top_100_genes = degs.sort_values('adj-p').dropna(subset=['GeneSymbol']).head(100)
    top_100_genes['GeneSymbol'].to_csv("results/DEGs/top100_genes.csv", index=False)

    #significant_genes = degs[degs['adj-p'] < 0.1]['GeneSymbol'].tolist()
    #significant_genes.to_csv("results/DEGs/diff_significant_genes.csv", index=False)




if __name__ == "__main__":
    run_pipeline(simulation_mode=False)


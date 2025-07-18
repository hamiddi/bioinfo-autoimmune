"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
RNA-Seq Analysis Pipeline with Two-Way ANOVA and Functional Enrichment
-----------------------------------------------------------------------

This Python program performs end-to-end RNA-Seq data analysis, including read trimming, genome alignment,
read counting, normalization, two-way ANOVA for differential expression, visualization, and functional enrichment.

USAGE:
------
1. **Simulation Mode (default)**:
   - To generate synthetic RNA-Seq data and perform analysis for testing and demonstration purposes:
     Run the script as-is:
         python script_name.py

2. **Real Data Mode**:
   - Set `simulation_mode=False` in the `run_pipeline()` call at the end of the script.
   - Ensure you have the following in place:
     • Raw FASTQ files in `data/raw/` named as `<runID>_1.fastq.gz` and `<runID>_2.fastq.gz`
     • A study design file at `meta/study_design.csv` with a column named `runID` and metadata columns (`Ethnic`, `antiCCP`)
     • Required tools installed and accessible: `fastp`, `STAR`, `samtools`, `featureCounts`
   - The script will:
     • Download the reference genome and GTF file from Ensembl if not already present
     • Generate STAR index
     • Preprocess reads and align them
     • Count reads per gene
     • Normalize expression values (CPM)
     • Perform two-way ANOVA
     • Generate volcano and heatmap plots
     • Identify significant genes and run GO/KEGG/Reactome enrichment analysis via g:Profiler

OUTPUT:
-------
Results will be saved in the `results/` directory, including:
 - Trimmed read quality reports
 - Normalized expression matrices
 - Differential expression results
 - Volcano and heatmap plots
 - Functional enrichment summary

DEPENDENCIES:
-------------
Requires the following Python packages:
    pandas, matplotlib, seaborn, numpy, scikit-learn, statsmodels, scipy, requests

To install all dependencies:
    pip install -r requirements.txt  (create a requirements file with needed packages)

"""

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

import statsmodels.api as sm
from statsmodels.formula.api import ols
from simulate_2wAnovaData import *

# Configuration
#Replace with the right PATH
THREADS = 8
REFERENCE_URL = "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
REFERENCE_FA = "../data/ref/hg38.fa"
STAR_INDEX = "../data/ref/star_index"
GTF_FILE = "../data/ref/genes.gtf"
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

def normalize_counts(count_matrix_path, design_file):
    counts = pd.read_csv(count_matrix_path, sep='\t', comment='#', index_col=0)
    counts = counts.iloc[:, 5:]
    counts.columns = [col.split('/')[-1].replace('_Aligned.sortedByCoord.out.bam', '') for col in counts.columns]
    study_design = pd.read_csv(design_file)
    available_samples = [s for s in study_design['runID'] if s in counts.columns]
    counts = counts[available_samples]
    study_design = study_design[study_design['runID'].isin(available_samples)]
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


def two_way_anova(norm_counts, design, factor1='Ethnic', factor2='antiCCP', min_samples=2):
    design = design.copy()
    design['runID'] = design['runID'].astype(str).str.strip()
    design[factor1] = design[factor1].astype(str).str.strip()
    design[factor2] = design[factor2].astype(str).str.strip()

    results = []

    for gene in norm_counts.index:
        gene_expr = norm_counts.loc[gene]
        expr_df = pd.DataFrame({
            'Expression': gene_expr,
            'runID': gene_expr.index
        })

        merged = expr_df.merge(design, on='runID')
        merged = merged.dropna(subset=['Expression', factor1, factor2])

        if merged.shape[0] < min_samples * 2:
            continue

        try:
            model = ols(f'Expression ~ C({factor1}) + C({factor2}) + C({factor1}):C({factor2})', data=merged).fit()
            anova_table = sm.stats.anova_lm(model, typ=2)

            p1 = anova_table.loc[f'C({factor1})', 'PR(>F)']
            p2 = anova_table.loc[f'C({factor2})', 'PR(>F)']
            p_int = anova_table.loc[f'C({factor1}):C({factor2})', 'PR(>F)']

            results.append((gene, p1, p2, p_int))
        except Exception as e:
            continue

    df = pd.DataFrame(results, columns=['Gene', f'{factor1}_p', f'{factor2}_p', 'Interaction_p'])
    df[f'{factor1}_adjp'] = multipletests(df[f'{factor1}_p'], method='fdr_bh')[1]
    df[f'{factor2}_adjp'] = multipletests(df[f'{factor2}_p'], method='fdr_bh')[1]
    df['Interaction_adjp'] = multipletests(df['Interaction_p'], method='fdr_bh')[1]

    return df.sort_values('Interaction_adjp')


def plot_volcano(degs, p_col='Interaction_adjp'):
    if p_col not in degs.columns:
        print(f"[WARNING] Column '{p_col}' not found in DE results. Skipping volcano plot.")
        return

    degs = degs.copy()
    degs['log_p'] = -np.log10(degs[p_col] + 1e-10)
    degs['GeneIndex'] = range(len(degs))

    plt.figure(figsize=(10,6))
    sns.scatterplot(data=degs, x='GeneIndex', y='log_p', s=20)
    plt.xlabel("Gene Index")
    plt.ylabel(f"-log10({p_col})")
    plt.title("Volcano Plot (Two-way ANOVA)")
    plt.tight_layout()
    plt.savefig("results/plots/volcano.png")
    plt.close()
###
def plot_heatmap(norm_counts, design, top_genes):
    subset = norm_counts.loc[top_genes]
    subset = StandardScaler().fit_transform(subset.T)
    sns.clustermap(pd.DataFrame(subset, index=norm_counts.columns, columns=top_genes), cmap="vlag")
    plt.title("Top Differentially Expressed Genes")
    plt.savefig("results/plots/top_degs_heatmap.png")
    plt.close()
###

def plot_heatmap(norm_counts, design, top_genes):
    from matplotlib.colors import to_hex
    from sklearn.preprocessing import StandardScaler
    import matplotlib.patches as mpatches

    # Subset top genes and normalize
    subset = norm_counts.loc[top_genes]
    data_scaled = StandardScaler().fit_transform(subset.T)
    data_df = pd.DataFrame(data_scaled, index=subset.columns, columns=top_genes)

    # Prepare group labels
    meta = design.set_index('runID').loc[data_df.index][['Ethnic', 'antiCCP']]
    meta['Group'] = meta['Ethnic'] + "_" + meta['antiCCP']

    # Assign unique colors to each group
    unique_groups = meta['Group'].unique()
    palette = sns.color_palette("tab10", len(unique_groups))
    lut = dict(zip(unique_groups, map(to_hex, palette)))
    col_colors = meta['Group'].map(lut)

    # Clustermap with column group annotations
    g = sns.clustermap(
        data_df.T,  # genes as rows
        cmap="vlag",
        col_colors=col_colors,
        xticklabels=False,
        yticklabels=True,
        figsize=(10, 8)
    )

    # Add gene labels
    g.ax_heatmap.set_yticklabels(g.data2d.index, rotation=0, fontsize=8)
    g.ax_heatmap.set_title("Top Differentially Expressed Genes", fontsize=12, pad=20)
    #g.ax_heatmap.set_title("Top Differentially Expressed Genes", fontsize=12)

    # Create legend patches for group colors
    legend_patches = [mpatches.Patch(color=lut[group], label=group) for group in unique_groups]

    # Place the legend outside the heatmap
    g.ax_col_dendrogram.legend(
        handles=legend_patches,
        title="Sample Groups",
        bbox_to_anchor=(1, 1),  # Position outside top-right
        loc='upper left',
        fontsize=8,
        title_fontsize=9,
        frameon=False
    )

    # Save the heatmap
    plt.savefig("results/plots/top_degs_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()


def plot_gene_expression(norm_counts, design, gene):
    merged = pd.DataFrame({
        'expression': norm_counts.loc[gene],
        'group': design.set_index('runID').loc[norm_counts.columns]['antiCCP']
    }).reset_index()
    sns.boxplot(data=merged, x='group', y='expression')
    sns.stripplot(data=merged, x='group', y='expression', color='black', alpha=0.5)
    plt.title(f"Expression of {gene}")
    plt.savefig(f"results/plots/{gene}_expression.png")
    plt.close()

def functional_enrichment(gene_list):
    if not gene_list:
        print("[INFO] No significant genes for enrichment.")
        return
    print("[INFO] Querying g:Profiler...")
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    headers = {'Content-Type': 'application/json'}
    payload = {
        "organism": "hsapiens",
        "query": gene_list,
        "sources": ["GO:BP", "KEGG", "REAC"]
    }
    response = requests.post(url, headers=headers, json=payload)
    if response.ok:
        res = response.json()
        terms = [(r['name'], r['p_value']) for r in res['result'][:10]]
        print("[INFO] Top Enriched Terms:")
        for name, p in terms:
            print(f" - {name} (p={p:.2e})")
    else:
        print("[ERROR] Failed to fetch enrichment results.")

def export_counts_with_symbols(count_matrix_path, gene_symbol_map, output_path):
    counts = pd.read_csv(count_matrix_path, sep='\t', comment='#', index_col=0)
    counts = counts.iloc[:, 5:]
    counts.reset_index(inplace=True)
    counts.rename(columns={'Geneid': 'accession'}, inplace=True)
    counts['gene_symbol'] = counts['accession'].map(gene_symbol_map)
    ordered_cols = ['accession', 'gene_symbol'] + [col for col in counts.columns if col not in ['accession', 'gene_symbol']]
    counts = counts[ordered_cols]
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    counts.to_csv(output_path, index=False)

def run_pipeline(simulation_mode=False):
    os.makedirs("results/DEGs", exist_ok=True)
    os.makedirs("results/plots", exist_ok=True)
    os.makedirs("results/counts", exist_ok=True)
    os.makedirs("results/fastp", exist_ok=True)

    gene_symbols = load_gene_symbols(GTF_FILE)
    """
    if simulation_mode:
        print("[SIMULATION] Skipping alignment and preprocessing.")
        study_design = pd.read_csv(STUDY_DESIGN_FILE)
        norm_counts, design = normalize_counts_sim("results/counts/gene_counts_simulated.txt", STUDY_DESIGN_FILE)
        """
    if simulation_mode:
        print("[SIMULATION] Generating synthetic data for two-way ANOVA...")
        norm_counts, design = simulate_data(num_genes=1000, samples_per_group=3)
        print(f"[INFO] Simulated {norm_counts.shape[0]} genes across {norm_counts.shape[1]} samples.")

    else:
        prepare_reference()
        study_design = pd.read_csv(STUDY_DESIGN_FILE)

        for _, row in study_design.iterrows():
            sid = row['runID']
            fq1 = f"data/raw/{sid}_1.fastq.gz"
            fq2 = f"data/raw/{sid}_2.fastq.gz"
            out1 = f"data/processed/{sid}_1.trim.fastq.gz"
            out2 = f"data/processed/{sid}_2.trim.fastq.gz"
            trim_reads(fq1, fq2, out1, out2, sid)
            map_reads(out1, out2, sid)
            index_bam(f"data/processed/{sid}_Aligned.sortedByCoord.out.bam")
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

    degs = two_way_anova(norm_counts, design, factor1='Ethnic', factor2='antiCCP')

    degs['GeneSymbol'] = degs['Gene'].map(gene_symbols)
    degs.to_csv("results/DEGs/diff_expr_with_symbols.csv", index=False)

    #plot_volcano(degs)
    plot_volcano(degs, p_col='Interaction_adjp')

    top_genes = degs.head(20)['Gene'].tolist()
    plot_heatmap(norm_counts, design, top_genes)
    if top_genes:
        plot_gene_expression(norm_counts, design, top_genes[0])

    significant_genes = degs[degs['Interaction_adjp'] < 0.1]['Gene'].tolist()
    #significant_genes = degs[degs['Ethnic_adjp'] < 0.1]['Gene'].tolist()
    #significant_genes = degs[degs['antiCCP_adjp'] < 0.1]['Gene'].tolist()

    functional_enrichment(significant_genes)

if __name__ == "__main__":
    run_pipeline(simulation_mode=True)


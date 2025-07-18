"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This script implements a full **assembly-based shotgun metagenomics analysis pipeline** for functional profiling, gene prediction, 
differential analysis, and genome binning using Python and external bioinformatics tools.

The pipeline includes the following stages:
1. Quality control using FastQC
2. De novo assembly of paired-end reads using MEGAHIT
3. Contig binning with MetaBAT2
4. Gene prediction from contigs with Prodigal
5. Functional annotation using DIAMOND with a UniProt database
6. Functional abundance matrix construction
7. Visualization of pathway-level data and PCA
8. Statistical testing for differential function abundance
9. Pathway reconstruction with MinPath
10. BAM generation and binning preparation with BWA and Samtools

Required Software (external dependencies):
- `fastqc`, `megahit`, `prodigal`, `diamond`, `bwa`, `samtools`, `metabat2`, `MinPath.py`

Required Python Packages:
- pandas, matplotlib, seaborn, sklearn, scipy, skbio, glob, os, subprocess

Paths to Provide/Verify:
- `DATABASE = "uniprot.dmnd"`: Path to a DIAMOND-formatted UniProt database
- Make sure `MinPath.py` is installed and accessible in your environment

Input Files:
- Paired-end FASTQ reads: located under `data/raw/`, named like `sample_1.fastq.gz` and `sample_2.fastq.gz`
- Metadata file: `meta/metadata.csv` with at least `runID`, `condition`, and `gender` columns

Output Files and Folders (created in `results_assembly/`):
- `qc/`: Quality control reports
- `assembly/`: MEGAHIT-assembled contigs (`final.contigs.fa`)
- `bins/`: BAM alignments and MetaBAT bins
- `genes/`: Predicted proteins in `.faa` format
- `diamond/`: Functional annotation TSVs from DIAMOND
- `minpath/`: MinPath reconstructed pathway outputs
- `functional_abundance_matrix.tsv`: Final matrix of function counts per sample
- `function_heatmap.png`: Top 25 most abundant functions heatmap
- `pca_functions.png`: PCA plot of functional profiles
- `differential_functions.csv`: Statistical test results between case/control groups

Execution:
Ensure all software and reference databases are accessible, then run:
```bash
python assembly_metagenomics_pipeline.py

Note:
BWA is used for mapping reads to contigs to assist with coverage estimation during binning.
MinPath step extracts KEGG Orthologs (KOs) from DIAMOND annotations for minimal pathway reconstruction.
This pipeline is ideal for reconstructing metabolic potential and microbial functions from metagenomic data without relying on taxonomic pre-annotation.
"""


import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
from sklearn.decomposition import PCA
from scipy.stats import ttest_ind
from skbio.diversity import alpha_diversity

# --- CONFIGURATION ---
RAW_DIR = "data/raw"
META_FILE = "meta/metadata.csv"
OUTPUT_DIR = "results_assembly"
THREADS = "8"
DATABASE = "uniprot.dmnd"  # Diamond-formatted protein db
PRODIGAL_TYPE = "meta"

# --- SETUP ---
def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def load_metadata(meta_file):
    df = pd.read_csv(meta_file)
    df.columns = df.columns.str.strip()
    df['runID'] = df['runID'].astype(str)
    return df

def get_sample_pairs(raw_dir):
    files = glob(os.path.join(raw_dir, "*.fastq.gz"))
    samples = {}
    for f in files:
        base = os.path.basename(f)
        sid = base.split("_")[0].split(".")[0]
        if sid not in samples:
            samples[sid] = [None, None]
        if "_1" in base or "_R1" in base:
            samples[sid][0] = f
        elif "_2" in base or "_R2" in base:
            samples[sid][1] = f
    return samples

# --- PIPELINE FUNCTIONS ---
def run_fastqc(read1, read2, outdir):
    cmd = f"fastqc -t {THREADS} -o {outdir} {read1} {read2}"
    subprocess.run(cmd, shell=True, check=True)
"""
def run_megahit(read1, read2, sample_id, outdir):
    sample_out = os.path.join(outdir, sample_id)
    cmd = f"megahit -1 {read1} -2 {read2} -o {sample_out} -t {THREADS}"
    subprocess.run(cmd, shell=True, check=True)
    return os.path.join(sample_out, "final.contigs.fa")
"""
def run_megahit(read1, read2, sample_id, outdir):
    sample_out = os.path.join(outdir, sample_id)
    cmd = f"megahit -1 {read1} -2 {read2} -o {sample_out} -t {THREADS}"
    subprocess.run(cmd, shell=True, check=True)
    return os.path.join(sample_out, "final.contigs.fa")

def run_prodigal(contigs, sample_id, outdir):
    genes_faa = os.path.join(outdir, f"{sample_id}.faa")
    cmd = f"prodigal -i {contigs} -a {genes_faa} -p {PRODIGAL_TYPE} -q"
    subprocess.run(cmd, shell=True, check=True)
    return genes_faa

def run_diamond(faa, sample_id, outdir):
    diamond_out = os.path.join(outdir, f"{sample_id}_diamond.tsv")
    cmd = (
     f"diamond blastp "
     f"-d {DATABASE} "
     f"-q {faa} "
     f"-o {diamond_out} "
     f"-f 6 "
     f"-k 1 "
     f"-p {THREADS}"
    )
    subprocess.run(cmd, shell=True, check=True)
    return diamond_out

def count_function_hits(tsv_file):
    df = pd.read_csv(tsv_file, sep="\t", header=None)
    df.columns = [
     "query", "target", "pident", "length", "mismatch", "gapopen",
     "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    return df["target"].value_counts().to_dict()

def build_abundance_matrix(hit_dicts):
    all_functions = set()
    for hits in hit_dicts.values():
        all_functions.update(hits.keys())
    all_functions = sorted(all_functions)

    rows = []
    for sample, hits in hit_dicts.items():
        row = [hits.get(func, 0) for func in all_functions]
        rows.append(row)
    df = pd.DataFrame(rows, index=hit_dicts.keys(), columns=all_functions).T
    return df

# --- VISUALIZATION & ANALYSIS ---
def plot_function_heatmap(df, metadata, outdir):
    top = df.sum(axis=1).nlargest(25).index
    heat = (
     df.loc[top]
     .T
     .reset_index()
     .merge(metadata, left_on="index", right_on="runID")
     .set_index("runID")
    )
    sns.heatmap(heat[top], cmap="mako")
    plt.title("Top 25 Functions")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "function_heatmap.png"))

def plot_pca(df, metadata, outdir):
    comp = PCA(n_components=2).fit_transform(df.T)
    pc_df = pd.DataFrame(comp, columns=["PC1", "PC2"])
    pc_df["runID"] = df.columns
    pc_df = pc_df.merge(metadata, on="runID")
    sns.scatterplot(data=pc_df, x="PC1", y="PC2", hue="condition", style="gender")
    plt.title("PCA of Functional Profiles")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "pca_functions.png"))

def analyze_differential(df, metadata, outdir):
    condition_map = metadata.set_index("runID")["condition"].to_dict()
    ctrl = [c for c in df.columns if condition_map.get(c, "").lower() == "control"]
    case = [c for c in df.columns if condition_map.get(c, "").lower() == "ms"]
    results = []
    for func in df.index:
        stat, pval = ttest_ind(
                     df.loc[func, ctrl],
                     df.loc[func, case],
                     equal_var=False
        )
        results.append((func, pval))
    df = pd.DataFrame(results, columns=["function", "p_value"])
    df = df.sort_values("p_value")
    output_path = os.path.join(outdir, "differential_functions.csv")
    df.to_csv(output_path, index=False)

## ------------------
def index_contigs(contigs):
    cmd = f"samtools faidx {contigs}"
    subprocess.run(cmd, shell=True, check=True)

def run_bwa_mapping(contigs, read1, read2, bam_out):
    index_cmd = f"bwa index {contigs}"
    aln_cmd = (
      f"bwa mem -t {THREADS} {contigs} {read1} {read2} | "
      "samtools view -bS - | "
      f"samtools sort -o {bam_out}"
    )
    index_bam = f"samtools index {bam_out}"
    subprocess.run(index_cmd, shell=True, check=True)
    subprocess.run(aln_cmd, shell=True, check=True)
    subprocess.run(index_bam, shell=True, check=True)

def run_metabat(contigs, bam_file, bin_dir):
    cmd = (
      f"metabat2 "
      f"-i {contigs} "
      f"-a <(jgi_summarize_bam_contig_depths "
      f"--outputDepth depth.txt "
      f"{bam_file}) "
      f"-o {os.path.join(bin_dir, 'bin')}"
    )
    subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)
# --------- Pathway Reconstruction with MinPath ---------------
def run_minpath(annotated_file, sample_id, outdir):
    input_file = os.path.join(outdir, f"{sample_id}_ko.list")
    with open(annotated_file) as f_in, open(input_file, "w") as f_out:
        for line in f_in:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                f_out.write(f"{parts[1]}\n")
    minpath_out = os.path.join(outdir, f"{sample_id}_minpath.out")
    cmd = f"MinPath.py -any {input_file} -map KO -report {minpath_out}"
    subprocess.run(cmd, shell=True, check=True)
    return minpath_out


# --- MAIN CONTROLLER ---
def main():
    ensure_dir(OUTPUT_DIR)
    metadata = load_metadata(META_FILE)
    samples = get_sample_pairs(RAW_DIR)

    gene_hits = {}

    for sid, (r1, r2) in samples.items():
        if not r1 or not r2:
            continue

        print(f"\nrocessing sample: {sid}")

        # 1. Quality Control
        qc_dir = os.path.join(OUTPUT_DIR, "qc")
        ensure_dir(qc_dir)
        run_fastqc(r1, r2, qc_dir)

        # 2. Assembly
        asm_dir = os.path.join(OUTPUT_DIR, "assembly")
        ensure_dir(asm_dir)
        contigs = run_megahit(r1, r2, sid, asm_dir)

        # 3. Contig Indexing & Binning
        bin_dir = os.path.join(OUTPUT_DIR, "bins", sid)
        ensure_dir(bin_dir)
        bam_file = os.path.join(bin_dir, f"{sid}.bam")
        index_contigs(contigs)
        run_bwa_mapping(contigs, r1, r2, bam_file)
        run_metabat(contigs, bam_file, bin_dir)

        # 4. Gene Prediction
        gene_dir = os.path.join(OUTPUT_DIR, "genes")
        ensure_dir(gene_dir)
        faa = run_prodigal(contigs, sid, gene_dir)

        # 5. Functional Annotation with DIAMOND
        diam_dir = os.path.join(OUTPUT_DIR, "diamond")
        ensure_dir(diam_dir)
        diamond_out = run_diamond(faa, sid, diam_dir)

        # 6. Pathway Reconstruction with MinPath
        minpath_dir = os.path.join(OUTPUT_DIR, "minpath")
        ensure_dir(minpath_dir)
        run_minpath(diamond_out, sid, minpath_dir)

        # 7. Count functions from DIAMOND results
        gene_hits[sid] = count_function_hits(diamond_out)

    # 8. Build Abundance Matrix
    abundance_df = build_abundance_matrix(gene_hits)
    abundance_df.to_csv(os.path.join(OUTPUT_DIR, "functional_abundance_matrix.tsv"), sep="\t")

    # 9. Visualization and Statistical Analysis
    plot_function_heatmap(abundance_df, metadata, OUTPUT_DIR)
    plot_pca(abundance_df, metadata, OUTPUT_DIR)
    analyze_differential(abundance_df, metadata, OUTPUT_DIR)

    print("\nPipeline complete. All results are saved in:", OUTPUT_DIR)


if __name__ == "__main__":
    main()


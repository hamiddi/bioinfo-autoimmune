"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script implements a comprehensive shotgun metagenomics pipeline combining taxonomic and functional 
profiling with visualization and differential analysis. It uses Kaiju for taxonomic classification, HUMAnN for 
pathway-level functional profiling, and various Python libraries to generate plots and statistical outputs, including PCA, 
Shannon diversity, differential abundance (with and without gender stratification), and a final HTML summary report.

Required Python Packages:
- os
- re
- glob
- pandas
- numpy
- matplotlib
- seaborn
- plotly.express
- jinja2
- scipy.stats
- sklearn.decomposition
- skbio.diversity

External Tools Required:
- kaiju
- humann
- Optional tools like gzip and bowtie2 for upstream preprocessing

Paths to Provide:
- `KAIJU_DB` (replace `/provide_kaijudb_path/kaijudb/refseq/refseq/kaiju_db_refseq.fmi`)
- `NODES_DMP` (replace `/provide_kaijudb_path/kaijudb/kaijudb/refseq/nodes.dmp`)

Input Files:
- Paired-end FASTQ files: `data/raw/{sample}_1.fastq.gz` and `{sample}_2.fastq.gz`
- Metadata file: `meta/metadata.csv` (must include `runID`, `condition`, and `gender`)
- JASPAR MEME motif file (optional for extensions): `motifs/jaspar.meme`

Output Files (saved to `results/`):
- Taxonomic abundance plot: `taxonomic_abundance.png`
- Shannon diversity plot: `shannon_diversity.png`
- PCA interactive plot: `pca_taxa.html`
- Differential abundance tables: `differential_abundance.csv`, `gender_stratified_diff.csv`
- Merged Kaiju table: in memory
- HUMAnN merged pathways: `humann/merged_pathabundance.tsv`
- Functional pathway heatmap: `top_pathways_heatmap.png`
- Final HTML report: `report.html`

Workflow Summary:
1. Identify and pair input FASTQ files
2. Run Kaiju for taxonomic profiling
3. Merge and analyze taxonomic data (diversity, PCA, differential abundance)
4. Run HUMAnN for pathway-level functional profiling
5. Visualize top pathways
6. Generate interactive and static visual summaries in an HTML report

Usage:
Make sure all required tools and paths are correctly configured, then run:
```bash
python shotgun_class_pipeline.py

Note:
You must update KAIJU_DB and NODES_DMP with valid paths to your Kaiju database.
HUMAnN may require pre-indexed databases (e.g., ChocoPhlAn, UniRef) and Bowtie2/Samtools support.
This pipeline assumes paired-end reads; it can be modified for single-end if needed.
"""



from sklearn.decomposition import PCA
import plotly.express as px
from jinja2 import Environment, FileSystemLoader
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
from scipy.stats import ttest_ind
from skbio.diversity import alpha_diversity

# -------- CONFIGURATION --------
RAW_DIR = "data/raw"
META_FILE = "meta/metadata.csv"
KAIJU_DB = "/biodata/users/cse815/cse815/software/progs/kaijudb/refseq/refseq/kaiju_db_refseq.fmi"
NODES_DMP = "/biodata/users/cse815/cse815/software/progs/kaijudb/refseq/nodes.dmp"
OUTPUT_DIR = "results"
HUMANN_OUTPUT = os.path.join(OUTPUT_DIR, "humann")
THREADS = "16"

# -------- UTILITY FUNCTIONS --------
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

def load_metadata(meta_file):
    df = pd.read_csv(meta_file)
    df.columns = df.columns.str.strip()
    df['runID'] = df['runID'].astype(str)
    return df

# -------- TAXONOMIC PROFILING --------
"""
def run_kaiju(read1, read2, output_prefix, db_path):
    for i, read in enumerate([read1, read2], 1):
        out_file = f"{output_prefix}_R{i}.out"
        cmd = [
            "kaiju",
            "-t", NODES_DMP,
            "-f", db_path,
            "-i", read,
            "-o", out_file,
            "-z", THREADS
        ]
        print(f"Running Kaiju on {read}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print("Kaiju failed:", result.stderr)
            raise RuntimeError("Kaiju failed")
        output_files.append(out_file)
    return output_files
"""
def run_kaiju(read1, read2, output_prefix, db_path):
    output_files = []  # <<< YOU MUST ADD THIS

    for i, read in enumerate([read1, read2], 1):
        out_file = f"{output_prefix}_R{i}.out"
        cmd = [
            "kaiju",
            "-t", NODES_DMP,
            "-f", db_path,
            "-i", read,
            "-o", out_file,
            "-z", THREADS
        ]
        print(f"Running Kaiju on {read}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print("Kaiju failed:", result.stderr)
            raise RuntimeError("Kaiju failed")
        
        output_files.append(out_file)   # Now it will work correctly
    return output_files


def merge_kaiju_outputs(kaiju_files):
    dfs = []
    for f in kaiju_files:
        try:
            df = pd.read_csv(f, sep="\t", header=None)
            if df.shape[1] == 3:  # Correct number of columns for your Kaiju output
                df.columns = ["status", "read_id", "taxon_id"]
                sample_id = os.path.basename(f).split("_kaiju")[0]
                df["sample"] = sample_id
                dfs.append(df)
            else:
                print(f"Skipping malformed Kaiju file (wrong number of columns): {f}")
        except Exception as e:
            print(f"Skipping malformed Kaiju file ({f}): {e}")

    if not dfs:
        raise ValueError("No Kaiju output files could be merged.")

    full_df = pd.concat(dfs)

    # Only use classified reads (C), not unclassified (U)
    classified_df = full_df[full_df["status"] == "C"]

    # Count number of classified reads per taxon_id per sample
    merged = classified_df.groupby(["sample", "taxon_id"]).size().reset_index(name="count")

    # Pivot to have taxon IDs as rows, samples as columns
    pivot = merged.pivot(index="taxon_id", columns="sample", values="count").fillna(0)

    return pivot

# -------- FUNCTIONAL PROFILING --------
def run_humann(input_fastq, sample_id, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    cmd = [
        "humann", "--input", input_fastq,
        "--output", output_dir,
        "--output-basename", sample_id,
        "--threads", THREADS
    ]
    subprocess.run(cmd, check=True)

def merge_humann_tables(output_dir):
    cmd = [
        "humann_join_tables",
        "--input", output_dir,
        "--output", os.path.join(output_dir, "merged_pathabundance.tsv"),
        "--file_name", "pathabundance"
    ]
    subprocess.run(cmd, check=True)

# -------- DIVERSITY AND DIFFERENTIAL ANALYSIS --------
def plot_taxonomic_abundance(count_df, metadata, output_dir):
    df = count_df.T.reset_index().melt(id_vars="index")
    df.columns = ["sample", "taxon", "count"]
    df = df.merge(metadata, left_on="sample", right_on="runID")
    plt.figure(figsize=(14, 6))
    sns.barplot(data=df, x="taxon", y="count", hue="condition", errorbar=None)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "taxonomic_abundance.png"))

def plot_shannon_diversity(count_df, metadata, output_dir):
    counts = count_df.T
    counts.index.name = "sample"
    shannon = alpha_diversity("shannon", counts.values, ids=counts.index)
    df = pd.DataFrame({"sample": shannon.index, "shannon": shannon.values})
    df = df.merge(metadata, left_on="sample", right_on="runID")
    plt.figure(figsize=(6, 4))
    sns.boxplot(data=df, x="condition", y="shannon", hue="gender")
    plt.title("Shannon Diversity by Condition and Gender")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "shannon_diversity.png"))

def analyze_differential_abundance(count_df, metadata, output_dir):
    group_map = metadata.set_index("runID")["condition"].to_dict()
    control = [c for c in count_df.columns if group_map.get(c, '').lower() == "control"]
    treated = [c for c in count_df.columns if group_map.get(c, '').lower() == "ms"]
    result = []
    for taxon in count_df.index:
        ctrl_vals = count_df.loc[taxon, control]
        treat_vals = count_df.loc[taxon, treated]
        stat, pval = ttest_ind(ctrl_vals, treat_vals, equal_var=False)
        result.append((taxon, pval))
    df_result = pd.DataFrame(result, columns=["taxon", "p_value"]).sort_values("p_value")
    df_result.to_csv(os.path.join(output_dir, "differential_abundance.csv"), index=False)

# ------Gender-Stratified Differential Taxonomic Analysis (Optional) ----
def analyze_gender_stratified_differential(count_df, metadata, output_dir):
    results = []
    for gender in metadata["gender"].unique():
        meta_sub = metadata[metadata["gender"] == gender]
        group_map = meta_sub.set_index("runID")["condition"].to_dict()
        ctrl = [c for c in count_df.columns if group_map.get(c, "").lower() == "control"]
        treat = [c for c in count_df.columns if group_map.get(c, "").lower() == "ms"]
        for taxon in count_df.index:
            ctrl_vals = count_df.loc[taxon, ctrl]
            treat_vals = count_df.loc[taxon, treat]
            if len(ctrl_vals) >= 2 and len(treat_vals) >= 2:
                stat, pval = ttest_ind(ctrl_vals, treat_vals, equal_var=False)
                results.append((taxon, gender, pval))
    df = pd.DataFrame(results, columns=["taxon", "gender", "p_value"])
    df.to_csv(os.path.join(output_dir, "gender_stratified_diff.csv"), index=False)

# ----------------- PCA visualization of taxonomic profiles ----------------
def plot_pca(count_df, metadata, output_dir):
    df = count_df.T
    df = df.loc[metadata["runID"]]
    pca = PCA(n_components=2)
    components = pca.fit_transform(df)
    pc_df = pd.DataFrame(components, columns=["PC1", "PC2"])
    pc_df["runID"] = df.index
    pc_df = pc_df.merge(metadata, on="runID")
    fig = px.scatter(
        pc_df, x="PC1", y="PC2", color="condition", symbol="gender",
        hover_data=["runID"], title="PCA of Taxonomic Profiles"
    )
    fig.write_html(os.path.join(output_dir, "pca_taxa.html"))

# --------------- Heatmap of the top pathways from HUMAnN ------------
def plot_top_pathway_heatmap(humann_path, metadata, output_dir, top_n=20):
    df = pd.read_csv(humann_path, sep="\t", index_col=0)
    df = df.loc[df.index.str.startswith("UNMAPPED") == False]
    df = df.drop(columns=["UNMAPPED"], errors='ignore')
    top_pathways = df.sum(axis=1).nlargest(top_n).index
    top_df = df.loc[top_pathways]
    top_df.columns.name = "runID"
    top_df = top_df.T.reset_index().merge(metadata, on="runID").set_index("runID")
    heatmap_data = top_df[top_pathways]
    plt.figure(figsize=(12, 8))
    sns.heatmap(heatmap_data.T, cmap="viridis", cbar_kws={"label": "Abundance"})
    plt.title("Top Functional Pathways (HUMAnN)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "top_pathways_heatmap.png"))


def generate_html_report(output_dir):
    env = Environment(loader=FileSystemLoader("."))
    template_str = """
    <html>
    <head><title>Shotgun Metagenomics Report</title></head>
    <body>
    <h1>Shotgun Metagenomics Analysis Report</h1>
    <ul>
        <li><a href="taxonomic_abundance.png">Taxonomic Abundance</a></li>
        <li><a href="shannon_diversity.png">Shannon Diversity Plot</a></li>
        <li><a href="pca_taxa.html">PCA of Taxonomic Profiles</a></li>
        <li><a href="top_pathways_heatmap.png">Top Pathway Heatmap</a></li>
        <li><a href="differential_abundance.csv">Differential Abundance Table</a></li>
        <li><a href="humann/merged_pathabundance.tsv">Merged HUMAnN Pathways</a></li>
    </ul>
    </body>
    </html>
    """
    with open(os.path.join(output_dir, "report.html"), "w") as f:
        f.write(template_str)



# -------- MAIN PIPELINE --------
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    metadata = load_metadata(META_FILE)
    samples = get_sample_pairs(RAW_DIR)

    kaiju_list_file = os.path.join(OUTPUT_DIR, "kaiju_outputs.txt")
    kaiju_outputs = []


    
    # === STEP 1: Run Kaiju if needed ===
    if not os.path.exists(kaiju_list_file):
        print("🛠 Running Kaiju classification...")
        for sample_id, (r1, r2) in samples.items():
            if not r1 or not r2:
                continue
            print(f"Running Kaiju on sample {sample_id}")
            output_prefix = os.path.join(OUTPUT_DIR, f"{sample_id}_kaiju.out")
            output_files = run_kaiju(r1, r2, output_prefix, KAIJU_DB)
            kaiju_outputs.extend(output_files)

        # Save list of Kaiju output files
        with open(kaiju_list_file, "w") as f:
            for file in kaiju_outputs:
                f.write(file + "\n")
    else:
        print("Loading existing Kaiju outputs...")
        with open(kaiju_list_file, "r") as f:
            kaiju_outputs = [line.strip() for line in f.readlines()]
    
    # === STEP 2: Merge and Analyze Kaiju outputs ===
    kaiju_table = merge_kaiju_outputs(kaiju_outputs)

    plot_taxonomic_abundance(kaiju_table, metadata, OUTPUT_DIR)
    plot_shannon_diversity(kaiju_table, metadata, OUTPUT_DIR)
    plot_pca(kaiju_table, metadata, OUTPUT_DIR)
    analyze_differential_abundance(kaiju_table, metadata, OUTPUT_DIR)
    analyze_gender_stratified_differential(kaiju_table, metadata, OUTPUT_DIR)

    # === STEP 3: HUMAnN Functional Analysis ===
    for sample_id, (r1, r2) in samples.items():
        if not r1:
            continue
        run_humann(r1, sample_id, HUMANN_OUTPUT)

    merge_humann_tables(HUMANN_OUTPUT)
    merged_path = os.path.join(HUMANN_OUTPUT, "merged_pathabundance.tsv")
    plot_top_pathway_heatmap(merged_path, metadata, OUTPUT_DIR)

    generate_html_report(OUTPUT_DIR)

    print("Pipeline finished.")

if __name__ == "__main__":
    main()


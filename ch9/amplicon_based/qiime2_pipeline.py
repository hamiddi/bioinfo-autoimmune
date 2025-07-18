"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script automates a complete amplicon-based microbiome analysis pipeline using QIIME 2. 
It performs sequence import, quality filtering, denoising, taxonomy classification, phylogenetic tree 
construction, diversity analysis, and differential abundance testing using ANCOM. The script is structured 
as modular functions and processes QIIME 2 `.qza` artifacts and `.qzv` visualizations, making it suitable for 
16S rRNA or other marker gene datasets.

Required Python Package:
- os
- subprocess

External Tools Required:
- QIIME 2 (must be installed in your environment)

Input Files:
- Manifest file: `meta/manifest.csv` (for paired-end FASTQ input)
- Metadata file: `meta/metadata.tsv` (sample metadata for diversity and ANCOM)
- Pre-trained classifier: `classifier/silva-138-99-nb-classifier.qza` (or compatible)

Output Files (key artifacts and visualizations):
- Imported demultiplexed sequences: `data/paired-end-demux.qza`
- DADA2 outputs: `results/table.qza`, `results/rep-seqs.qza`, `results/denoising-stats.qza`
- Visual summaries: `.qzv` files for quality, tables, sequences, and taxonomy
- Taxonomy assignments: `results/taxonomy.qza`
- Phylogenetic tree: `results/rooted-tree.qza`
- Core diversity metrics: `results/core-metrics-results/`
- ANCOM differential abundance results: `results/ancom-condition.qzv`

Workflow Summary:
1. Import paired-end FASTQ files using a QIIME 2 manifest
2. Summarize demultiplexed reads to assess quality
3. Denoise reads with DADA2 (truncates to 250 bp)
4. Visualize the feature table, sequences, and denoising stats
5. Assign taxonomy using a pre-trained Naive Bayes classifier (e.g., SILVA)
6. Generate taxonomy bar plots and metadata tables
7. Build a phylogenetic tree using MAFFT and FastTree
8. Compute alpha and beta diversity metrics
9. Visualize group significance and beta diversity using Emperor plots
10. Perform differential abundance testing using ANCOM

Note:
- Ensure QIIME 2 is activated in your shell environment before running this script.
- Modify truncation lengths and metadata column names as needed for your dataset.
- Pre-trained classifiers must match the target region and primer set used in your study.
"""

import os
import subprocess


def run_command(command, description=""):
    """Run a shell command and print status."""
    print(f"\n Running: {description}\n{command}")
    process = subprocess.run(command, shell=True)
    if process.returncode != 0:
        raise RuntimeError(f"Command failed: {description}")

def import_sequences():
    """Import sequences using manifest format."""
    cmd = (
        "qiime tools import "
        "--type 'SampleData[PairedEndSequencesWithQuality]' "
        "--input-path meta/manifest.csv "
        "--output-path data/paired-end-demux.qza "
        "--input-format PairedEndFastqManifestPhred33V2"
    )
    run_command(cmd, "Importing paired-end sequences via manifest")


def summarize_sequences():
    """Summarize demux data to assess quality."""
    cmd = (
        "qiime demux summarize "
        "--i-data data/paired-end-demux.qza "
        "--o-visualization results/demux-summary.qzv"
    )
    run_command(cmd, "Summarizing demultiplexed sequences")


def denoise_dada2():
    """Run DADA2 for quality control, denoising, and chimera removal."""
    cmd = (
        "qiime dada2 denoise-paired "
        "--i-demultiplexed-seqs data/paired-end-demux.qza "
        "--p-trim-left-f 0 "
        "--p-trim-left-r 0 "
        "--p-trunc-len-f 250 "
        "--p-trunc-len-r 250 "
        "--o-table results/table.qza "
        "--o-representative-sequences results/rep-seqs.qza "
        "--o-denoising-stats results/denoising-stats.qza"
    )
    run_command(cmd, "Denoising with DADA2")


def visualize_feature_table():
    """Visualize feature table and representative sequences."""
    cmds = [
        (
            "qiime feature-table summarize "
            "--i-table results/table.qza "
            f"--o-visualization results/table.qzv "
            "--m-sample-metadata-file meta/metadata.tsv"
        ),
        (
            "qiime feature-table tabulate-seqs "
            "--i-data results/rep-seqs.qza "
            "--o-visualization results/rep-seqs.qzv"
        ),
        (
            "qiime metadata tabulate "
            "--m-input-file results/denoising-stats.qza "
            "--o-visualization results/denoising-stats.qzv"
        ),
    ]
    for cmd in cmds:
        run_command(cmd, "Visualizing feature table and stats")


def assign_taxonomy():
    """Assign taxonomy using a pre-trained classifier (e.g., SILVA or Greengenes)."""
    classifier_path = "classifier/silva-138-99-nb-classifier.qza"  # Ensure this file exists
    if not os.path.exists(classifier_path):
        raise FileNotFoundError(f"{classifier_path} not found. Download a compatible classifier.")

    cmd = (
        "qiime feature-classifier classify-sklearn "
        f"--i-classifier {classifier_path} "
        "--i-reads results/rep-seqs.qza "
        "--o-classification results/taxonomy.qza"
    )
    run_command(cmd, "Assigning taxonomy")


def visualize_taxonomy():
    """Visualize taxonomy assignments."""
    cmds = [
        (
            "qiime metadata tabulate "
            "--m-input-file results/taxonomy.qza "
            "--o-visualization results/taxonomy.qzv"
        ),
        (
            "qiime taxa barplot "
            "--i-table results/table.qza "
            "--i-taxonomy results/taxonomy.qza "
            "--m-metadata-file meta/metadata.tsv "
            "--o-visualization results/taxa-bar-plots.qzv"
        ),
    ]
    for cmd in cmds:
        run_command(cmd, "Generating taxonomy visualizations")

def generate_phylogeny():
    """Generates phylogenetic tree required for beta diversity."""
    cmds = [
        (
            "qiime phylogeny align-to-tree-mafft-fasttree "
            "--i-sequences results/rep-seqs.qza "
            "--o-alignment results/aligned-rep-seqs.qza "
            "--o-masked-alignment results/masked-alignment.qza "
            "--o-tree results/unrooted-tree.qza "
            "--o-rooted-tree results/rooted-tree.qza"
        )
    ]
    for cmd in cmds:
        run_command(cmd, "Generating phylogenetic tree")


def compute_diversity_metrics(sampling_depth=10000):
    """Run core diversity metrics for alpha and beta diversity."""
    cmd = (
        "qiime diversity core-metrics-phylogenetic "
        "--i-phylogeny results/rooted-tree.qza "
        "--i-table results/table.qza "
        f"--p-sampling-depth {sampling_depth} "
        "--m-metadata-file meta/metadata.tsv "
        "--output-dir results/core-metrics-results"
    )
    run_command(cmd, "Computing diversity metrics")


def visualize_diversity():
    """Group significance tests and visualizations for diversity."""
    cmds = [
        (
            "qiime diversity alpha-group-significance "
            "--i-alpha-diversity results/core-metrics-results/faith_pd_vector.qza "
            "--m-metadata-file meta/metadata.tsv "
            "--o-visualization results/core-metrics-results/faith-pd-group-significance.qzv"
        ),
        (
            "qiime diversity alpha-group-significance "
            "--i-alpha-diversity results/core-metrics-results/evenness_vector.qza "
            "--m-metadata-file meta/metadata.tsv "
            "--o-visualization results/core-metrics-results/evenness-group-significance.qzv"
        ),
        (
            "qiime diversity beta-group-significance "
            "--i-distance-matrix results/core-metrics-results/unweighted_unifrac_distance_matrix.qza "
            "--m-metadata-file meta/metadata.tsv "
            "--m-metadata-column condition "
            "--o-visualization results/core-metrics-results/unweighted-unifrac-condition.qzv "
            "--p-pairwise"
        ),
        (
            "qiime emperor plot "
            "--i-pcoa results/core-metrics-results/unweighted_unifrac_pcoa_results.qza "
            "--m-metadata-file meta/metadata.tsv "
            "--o-visualization results/core-metrics-results/unweighted-unifrac-emperor.qzv"
        )
    ]
    for cmd in cmds:
        run_command(cmd, "Visualizing alpha/beta diversity")


def differential_abundance():
    """Perform differential abundance analysis using ANCOM."""
    cmds = [
        (
            "qiime composition add-pseudocount "
            "--i-table results/table.qza "
            "--o-composition-table results/composition-table.qza"
        ),
        (
            "qiime composition ancom "
            "--i-table results/composition-table.qza "
            "--m-metadata-file meta/metadata.tsv "
            "--m-metadata-column condition "
            "--o-visualization results/ancom-condition.qzv"
        )
    ]
    for cmd in cmds:
        run_command(cmd, "Running differential abundance with ANCOM")



def main():
    os.makedirs("results", exist_ok=True)
    import_sequences()
    summarize_sequences()
    denoise_dada2()
    visualize_feature_table()
    assign_taxonomy()
    visualize_taxonomy()
    generate_phylogeny()
    compute_diversity_metrics()
    visualize_diversity()
    differential_abundance()
    print("\nAmplicon-based QIIME 2 analysis complete with diversity and differential analysis!")


if __name__ == "__main__":
    main()


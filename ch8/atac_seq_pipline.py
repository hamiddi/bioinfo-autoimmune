"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script implements a complete single-end ATAC-Seq data analysis pipeline. It automates 
the steps of reference genome setup, read trimming, alignment using Bowtie2, BAM filtering, and peak calling using MACS2. 
The pipeline is metadata-driven, reading sample IDs and conditions from a CSV file, and organizes output into structured directories.

Required Python Packages:
- os
- subprocess
- pandas

External Tools Required:
- wget
- gunzip
- fastp
- samtools
- bowtie2
- MACS2

Input Files:
- Metadata CSV: `meta/metadata.csv` with columns including `SRARUNID` and `condition`
- Raw FASTQ files: `data/raw/{SRARUNID}.fastq.gz` (single-end reads)

Output Files:
- Trimmed FASTQ files: `data/trimmed/{SRARUNID}_trimmed.fastq.gz`
- Aligned BAM files: `data/aligned/{SRARUNID}.bam` and filtered versions
- BAM index files: `*.bam.bai`
- Peak files: `data/peaks/{SRARUNID}_{condition}/*.narrowPeak`
- Reference genome: downloaded and indexed in `reference/` (hg38)

Workflow Summary:
1. Download and index the hg38 reference genome (Bowtie2 index and samtools FASTA index)
2. Trim adapters and low-quality bases using `fastp`
3. Align trimmed reads to the genome using `bowtie2` and sort with `samtools`
4. Filter out mitochondrial reads and generate BAM index
5. Call peaks using `macs2` for each sample

Note:
- Designed for single-end sequencing reads.
- Modify the trimming and alignment steps for paired-end data if needed.
- Ensure all required external tools are installed and in the system's PATH.
"""

import os
import subprocess
import pandas as pd

# Load metadata
metadata = pd.read_csv("meta/metadata.csv", sep=',')

# Define directories
RAW_DIR = "data/raw"
TRIMMED_DIR = "data/trimmed"
ALIGNED_DIR = "data/aligned"
PEAKS_DIR = "data/peaks"

for d in [TRIMMED_DIR, ALIGNED_DIR, PEAKS_DIR]:
    os.makedirs(d, exist_ok=True)


# Step 1: Download and index reference genome
def setup_reference_genome():
    ref_dir = "reference"
    os.makedirs(ref_dir, exist_ok=True)

    ref_url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    ref_gz = os.path.join(ref_dir, "hg38.fa.gz")
    ref_fa = os.path.join(ref_dir, "hg38.fa")
    bowtie2_index_base = os.path.join(ref_dir, "hg38")

    if not os.path.exists(ref_fa):
        print("Downloading hg38 reference genome...")
        subprocess.run(f"wget -O {ref_gz} {ref_url}", shell=True, check=True)
        subprocess.run(f"gunzip -f {ref_gz}", shell=True, check=True)

    # Indexing with samtools
    if not os.path.exists(ref_fa + ".fai"):
        print("Indexing with samtools...")
        subprocess.run(f"samtools faidx {ref_fa}", shell=True, check=True)

    # Indexing with bowtie2
    if not all(os.path.exists(f"{bowtie2_index_base}.{i}.bt2") for i in range(1, 5)):
        print("Indexing with bowtie2...")
        subprocess.run(f"bowtie2-build {ref_fa} {bowtie2_index_base}", shell=True, check=True)

    return bowtie2_index_base


# Step 2: Trimming
def trim_reads(sample_id):
    input_file = os.path.join(RAW_DIR, f"{sample_id}.fastq.gz")
    output_file = os.path.join(TRIMMED_DIR, f"{sample_id}_trimmed.fastq.gz")
    cmd = f"fastp -i {input_file} -o {output_file} --detect_adapter_for_pe -w 4"
    subprocess.run(cmd, shell=True, check=True)

# Step 3: Alignment
def align_reads(sample_id, genome_index):
    input_file = os.path.join(TRIMMED_DIR, f"{sample_id}_trimmed.fastq.gz")
    output_bam = os.path.join(ALIGNED_DIR, f"{sample_id}.bam")
    cmd = f"bowtie2 -x {genome_index} -U {input_file} | samtools view -Sb - | samtools sort -o {output_bam}"
    subprocess.run(cmd, shell=True, check=True)


# Step 4: Filter BAM (optional: remove mitochondrial, low quality, duplicates)
def filter_bam(sample_id):
    bam_file = os.path.join(ALIGNED_DIR, f"{sample_id}.bam")
    filtered_bam = os.path.join(ALIGNED_DIR, f"{sample_id}_filtered.bam")
    cmd = f"samtools view -h {bam_file} | grep -v 'chrM' | samtools view -Sb - | samtools sort -o {filtered_bam}"
    subprocess.run(cmd, shell=True, check=True)
    subprocess.run(f"samtools index {filtered_bam}", shell=True, check=True)

# Step 5: Peak calling
def call_peaks(sample_id, condition):
    bam_file = os.path.join(ALIGNED_DIR, f"{sample_id}_filtered.bam")
    out_dir = os.path.join(PEAKS_DIR, f"{sample_id}_{condition}")
    os.makedirs(out_dir, exist_ok=True)
    cmd = f"macs2 callpeak -t {bam_file} -f BAM -g hs -n {sample_id} --outdir {out_dir} --nomodel --shift -100 --extsize 200 -q 0.01"
    subprocess.run(cmd, shell=True, check=True)

# Pipeline runner
def run_pipeline():
    genome_index = setup_reference_genome()
    for idx, row in metadata.iterrows():
        sample_id = row['SRARUNID']
        condition = row['condition']
        print(f"\nProcessing sample {sample_id} [{condition}]")
        trim_reads(sample_id)
        align_reads(sample_id, genome_index)
        filter_bam(sample_id)
        call_peaks(sample_id, condition)

if __name__ == "__main__":
    run_pipeline()


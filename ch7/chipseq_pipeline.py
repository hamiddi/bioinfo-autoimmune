"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script implements a ChIP-Seq data analysis pipeline for identifying DNA-protein interaction sites. 
It automates key steps including downloading and indexing the reference genome, trimming raw reads, aligning 
paired-end sequences using BWA, converting and sorting SAM/BAM files with SAMtools, and calling peaks with MACS3. 
The script is designed to process multiple samples based on metadata provided in a CSV file.

Required Python Packages:
- os
- subprocess
- pandas
- glob

External Tools Required:
- wget
- gunzip
- BWA
- SAMtools
- Trimmomatic
- MACS3

Input Files:
- Metadata CSV file: `meta/metadata.csv` (must include `runID` and `condition` columns)
- Paired-end FASTQ files: stored in `data/raw/`, named as `{runID}_1.fastq.gz` and `{runID}_2.fastq.gz`
- Reference genome: downloaded from UCSC (hg38.fa.gz) and indexed for BWA alignment

Output Files:
- Trimmed FASTQ files: `results/trimmed/{runID}_1.trimmed.fastq.gz`, `{runID}_2.trimmed.fastq.gz`
- Aligned BAM files: `results/aligned/{runID}.bam`
- Peak files: `results/peaks/{runID}_peaks.narrowPeak`
- BWA reference index: generated in the `reference/` directory

Note:
This pipeline is tailored for ChIP-Seq experiments and assumes paired-end sequencing data. 
MACS3 is configured for broad and sharp peak detection using BAMPE mode with optional control-based adjustments.
Ensure all external tools are installed and accessible in the system environment.
"""


import os
import subprocess
import pandas as pd
from glob import glob


data_dir = "data/raw"
meta_file = "meta/metadata.csv"
out_dir = "results"
reference_dir = "reference"
genome_fasta = os.path.join(reference_dir, "hg38.fa")
bwa_index_prefix = os.path.join(reference_dir, "hg38")

os.makedirs(out_dir, exist_ok=True)
os.makedirs(reference_dir, exist_ok=True)

def download_and_index_reference():
    if all(os.path.exists(f"{bwa_index_prefix}.{ext}") for ext in ["amb", "ann", "bwt", "pac", "sa"]):
        print("Reference genome already indexed for BWA.")
        return

    print("Downloading and indexing reference genome for BWA...")
    url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    compressed_fasta = genome_fasta + ".gz"
    try:
        subprocess.run(["wget", "-O", compressed_fasta, url], check=True)
        subprocess.run(["gunzip", compressed_fasta], check=True)
        subprocess.run(["bwa", "index", "-p", bwa_index_prefix, genome_fasta], check=True)
    except subprocess.CalledProcessError as e:
        print("Error downloading or indexing genome:", e)
        raise

def read_metadata(meta_path):
    metadata = pd.read_csv(meta_path)
    return metadata

def find_fastq_pairs(run_id):
    r1 = glob(os.path.join(data_dir, f"{run_id}*_1.fastq.gz"))
    r2 = glob(os.path.join(data_dir, f"{run_id}*_2.fastq.gz"))
    if r1 and r2:
        return r1[0], r2[0]
    else:
        raise FileNotFoundError(f"FASTQ pairs not found for {run_id}")

def trim_reads(run_id, r1, r2):
    trimmed_dir = os.path.join(out_dir, "trimmed")
    os.makedirs(trimmed_dir, exist_ok=True)
    trimmed_r1 = os.path.join(trimmed_dir, f"{run_id}_1.trimmed.fastq.gz")
    trimmed_r2 = os.path.join(trimmed_dir, f"{run_id}_2.trimmed.fastq.gz")
    cmd = [
        "trimmomatic", "PE", "-phred33",
        r1, r2,
        trimmed_r1, "/dev/null",
        trimmed_r2, "/dev/null",
        "SLIDINGWINDOW:4:20", "MINLEN:36"
    ]
    subprocess.run(cmd, check=True)
    return trimmed_r1, trimmed_r2

def align_reads(run_id, trimmed_r1, trimmed_r2):
    aligned_dir = os.path.join(out_dir, "aligned")
    os.makedirs(aligned_dir, exist_ok=True)
    sam_output = os.path.join(aligned_dir, f"{run_id}.sam")
    bam_output = os.path.join(aligned_dir, f"{run_id}.bam")
    cmd_align = [
        "bwa", "mem", bwa_index_prefix,
        trimmed_r1, trimmed_r2
    ]
    with open(sam_output, "w") as samfile:
        subprocess.run(cmd_align, stdout=samfile, check=True)
    cmd_view = ["samtools", "view", "-bS", sam_output]
    cmd_sort = ["samtools", "sort", "-o", bam_output, "-"]
    p1 = subprocess.Popen(cmd_view, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd_sort, stdin=p1.stdout)
    p1.stdout.close()
    p2.communicate()
    os.remove(sam_output)
    return bam_output

def call_peaks(run_id, bam_file, condition):
    peaks_dir = os.path.join(out_dir, "peaks")
    os.makedirs(peaks_dir, exist_ok=True)
    peak_output = os.path.join(peaks_dir, f"{run_id}_peaks.narrowPeak")
    cmd = [
        "macs3", "callpeak",
        "-t", bam_file,
        "-n", run_id,
        "--outdir", peaks_dir,
        "-f", "BAMPE",
        "-g", "hs",
        "--keep-dup", "all",
        "-q", "0.01",
        "--nomodel"
    ]
    if condition == "control":
        cmd += ["--nolambda"]
    subprocess.run(cmd, check=True)
    return peak_output

def run_pipeline():
    download_and_index_reference()
    metadata = read_metadata(meta_file)
    for _, row in metadata.iterrows():
        run_id = row['runID']
        condition = row['condition']
        print(f"Processing {run_id} ({condition})...")
        r1, r2 = find_fastq_pairs(run_id)
        trimmed_r1, trimmed_r2 = trim_reads(run_id, r1, r2)
        bam_file = align_reads(run_id, trimmed_r1, trimmed_r2)
        peak_file = call_peaks(run_id, bam_file, condition)
        print(f"Finished {run_id}, peaks at: {peak_file}")

if __name__ == "__main__":
    run_pipeline()


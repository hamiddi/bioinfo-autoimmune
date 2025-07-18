"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script implements a pipeline for transcription factor footprinting and differential binding 
analysis using the TOBIAS suite. It performs signal correction, footprint scoring, and motif-based 
differential footprint detection across experimental groups from ATAC-Seq or ChIP-Seq data.

The script assumes peak regions are already called and merged into a BED file, and BAM alignments are 
filtered and indexed. TOBIAS tools (`ATACorrect`, `FootprintScores`, and `BINDetect`) are then used 
to compute binding scores and compare accessibility between sample groups defined in a metadata file.

Required Python Packages:
- os
- subprocess
- pandas
- shutil
- warnings

External Tools Required:
- TOBIAS (https://github.com/loosolab/TOBIAS)
- samtools (for BAM indexing, not used directly here but required for preprocessing)

Input Files:
- Reference genome FASTA: `reference/hg38.fa`
- Merged peak regions BED: `data/peaks/merged_peaks.bed`
- Filtered BAM files: `data/aligned/{sample_id}_filtered.bam`
- JASPAR or other motif file in MEME format: `motifs/jaspar.meme`
- Metadata CSV file: `meta/metadata.csv`, with sample IDs (`SRARUNID`) and grouping columns (`condition`, `CD38`, etc.)

Output Files:
- TOBIAS-corrected signal files: `data/tobias/{sample_id}_filtered_corrected.bw`
- Footprint score tracks: `data/tobias/{sample_id}_footprints.bw`
- Differential binding results (BINDetect): 
  - `data/tobias/differential_condition/`
  - `data/tobias/differential_CD38/`

Workflow Summary:
1. Correct ATAC-Seq signals for Tn5 bias using `TOBIAS ATACorrect`
2. Calculate footprint scores across peaks with `TOBIAS FootprintScores`
3. Compare footprinting signal across experimental conditions using `TOBIAS BINDetect`
4. Organize results by group identifiers (e.g., `condition`, `CD38`)

Note:
- Ensure TOBIAS and its dependencies (e.g., bedtools, samtools) are installed and in your system path.
- This script currently only executes the BINDetect phase (step 3). To enable ATACorrect and FootprintScores,
  uncomment the corresponding block in `run_pipeline()`.
"""

import os
import subprocess
import pandas as pd
import shutil
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

REFERENCE_FASTA = "reference/hg38.fa"
PEAKS_BED = "data/peaks/merged_peaks.bed"
BAM_DIR = "data/aligned"
OUTPUT_DIR = "data/tobias"
MOTIFS_FILE = "motifs/jaspar.meme"
METADATA_FILE = "meta/metadata.csv"

os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_metadata(metadata_path):
    df = pd.read_csv(metadata_path, sep=",")
    return df

def run_tobias_atacorrect(sample_id):
    bam_file = os.path.join(BAM_DIR, f"{sample_id}_filtered.bam")

    # TOBIAS will create this subdirectory and put outputs in it
    out_basename = os.path.join(OUTPUT_DIR, f"{sample_id}_filtered_corrected")

    # Final expected flat output location
    final_output_bw = os.path.join(OUTPUT_DIR, f"{sample_id}_filtered_corrected.bw")

    # If previous output is a directory, remove it
    if os.path.isdir(final_output_bw):
        print(f"Warning: {final_output_bw} exists as directory, removing.")
        shutil.rmtree(final_output_bw)
    elif os.path.isfile(final_output_bw):
        os.remove(final_output_bw)

    # Run TOBIAS ATACorrect
    cmd = [
        "TOBIAS", "ATACorrect",
        "--bam", bam_file,
        "--genome", REFERENCE_FASTA,
        "--peaks", PEAKS_BED,
        "--out", out_basename  # Do not add .bw
    ]
    subprocess.run(cmd, check=True)

    # Actual output path TOBIAS will create
    produced_bw = os.path.join(out_basename, f"{sample_id}_filtered_corrected.bw")

    # Move it to flat OUTPUT_DIR
    if os.path.isfile(produced_bw):
        shutil.move(produced_bw, final_output_bw)
    else:
        raise FileNotFoundError(f"TOBIAS did not create expected output: {produced_bw}")

def run_tobias_footprintscore(sample_id):
    corrected_bw = os.path.join(OUTPUT_DIR, f"{sample_id}_filtered_corrected.bw")
    out_fp = os.path.join(OUTPUT_DIR, f"{sample_id}_footprints.bw")

    # Optional safeguard: check if input .bw exists before running
    if not os.path.isfile(corrected_bw):
        raise FileNotFoundError(f"Expected corrected bigWig not found: {corrected_bw}")

    cmd = [
        "TOBIAS", "FootprintScores",
        "--signal", corrected_bw,
        "--regions", PEAKS_BED,
        "--output", out_fp
    ]
    subprocess.run(cmd, check=True)

def run_differential_binding(df, group_by="condition"):
    grouped = df.groupby(group_by)
    signal_files = []
    labels = []
    """
    # if condition in the metadata.csv: HC1, HC2, SLE1, SLE2, ...
    for name, group in grouped:
        label = str(name)
        for sample_id in group["SRARUNID"]:
            footprint_bw = os.path.join(OUTPUT_DIR, f"{sample_id}_footprints.bw")
            if os.path.isfile(footprint_bw):
                signal_files.append(footprint_bw)
                labels.append(label)
    """
    for name, group in grouped:
      for idx, sample_id in enumerate(group["SRARUNID"]):
        label = f"{name}{idx+1}"  #this will fix condition to e.g., HC1, HC2, SLE1, ...
        footprint_bw = os.path.join(OUTPUT_DIR, f"{sample_id}_footprints.bw")
        if os.path.isfile(footprint_bw):
            signal_files.append(footprint_bw)
            labels.append(label)

    out_dir = os.path.join(OUTPUT_DIR, f"differential_{group_by}")
    os.makedirs(out_dir, exist_ok=True)

    # 🛠 Correct order: --signals → --labels → [others]
    cmd = [
        "TOBIAS", "BINDetect",
        "--signals"
    ] + signal_files + [
        "--cond-names"
    ] + labels + [
        "--motifs", MOTIFS_FILE,
        "--peaks", PEAKS_BED,
        "--genome", REFERENCE_FASTA,
        "--outdir", out_dir
    ]
    print(cmd)
    subprocess.run(cmd, check=True)
    print(f"BINDetect completed for group: {group_by}")

def run_pipeline():
    df = load_metadata(METADATA_FILE)
    """
    for sample_id in df["SRARUNID"]:
        print(f"\nRunning TOBIAS for sample: {sample_id}")
        run_tobias_atacorrect(sample_id)
        run_tobias_footprintscore(sample_id)
    """
    run_differential_binding(df, group_by="condition")
    run_differential_binding(df, group_by="CD38")

if __name__ == "__main__":
    run_pipeline()


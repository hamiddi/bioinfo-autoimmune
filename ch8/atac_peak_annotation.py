"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script automates genomic peak annotation using HOMER's `annotatePeaks.pl`. It processes MACS2 peak 
files for each sample based on metadata, annotates genomic regions (e.g., promoters, exons, introns), and 
saves the results into a structured output directory.

The script reads sample identifiers and conditions from a CSV file, constructs the path to each sample's peak 
file, and generates annotated text files using HOMER for the specified genome (default: hg38).

Required Python Packages:
- os
- subprocess
- pandas

External Tools Required:
- HOMER (annotatePeaks.pl must be installed and configured in the system PATH)

Input Files:
- Metadata CSV file: `meta/metadata.csv` with columns `SRARUNID` and `condition`
- Peak files: `data/peaks/{SRARUNID}_{condition}/{SRARUNID}_peaks.narrowPeak` (MACS2 output)

Output Files:
- Annotated peak files: `data/annotated/{SRARUNID}_annotated.txt`

Note:
- Update the `GENOME` variable if using a genome other than `hg38` (e.g., `mm10`).
- Ensure HOMER is installed and the appropriate genome is configured via `configureHomer.pl`.
"""



import os
import subprocess
import pandas as pd

# Load metadata
metadata = pd.read_csv("meta/metadata.csv", sep=",")

# Directories
PEAKS_DIR = "data/peaks"
ANNOTATION_DIR = "data/annotated"
GENOME = "hg38"  # Change this if you use another genome

os.makedirs(ANNOTATION_DIR, exist_ok=True)

# Function to annotate peaks using HOMER
def annotate_peaks(sample_id, condition):
    peak_dir = os.path.join(PEAKS_DIR, f"{sample_id}_{condition}")
    peak_file = os.path.join(peak_dir, f"{sample_id}_peaks.narrowPeak")

    if not os.path.exists(peak_file):
        print(f"Peak file not found for {sample_id}, skipping...")
        return

    output_file = os.path.join(ANNOTATION_DIR, f"{sample_id}_annotated.txt")

    cmd = f"annotatePeaks.pl {peak_file} {GENOME} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"Annotated peaks for {sample_id} saved to {output_file}")

# Runner function
def run_annotation():
    for idx, row in metadata.iterrows():
        sample_id = row["SRARUNID"]
        condition = row["condition"]
        annotate_peaks(sample_id, condition)

if __name__ == "__main__":
    run_annotation()


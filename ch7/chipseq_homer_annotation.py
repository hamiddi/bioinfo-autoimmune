"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script performs genomic annotation of MACS3 peak files using HOMER's `annotatePeaks.pl` tool. 
It supports both `.narrowPeak` and `.bed` input formats by converting narrowPeak files to BED when necessary. 
The script is command-line driven, allowing users to specify the input file, genome build (e.g., hg38), and 
output file name. It automates peak-to-gene annotation for ChIP-Seq or ATAC-Seq data, providing biological 
context to enriched regions such as promoter proximity, gene names, and functional regions.

Required Python Packages:
- os
- subprocess
- argparse
- pandas

External Tools Required:
- HOMER (annotatePeaks.pl must be installed and accessible in the system PATH)

Input Files:
- MACS3 `.narrowPeak` file or standard `.bed` file containing peak coordinates
- Genome reference tag supported by HOMER (e.g., `hg38`, `mm10`)

Output Files:
- Annotated text file (default: `annotated_peaks.txt`) with HOMER-based functional peak annotations

Usage Example:
```bash
python annotate_peaks_homer.py -i sample_peaks.narrowPeak -g hg38 -o annotated_peaks.txt

Note:
Ensure HOMER is properly installed and configured with the target genome before running this script.
"""


import os
import subprocess
import argparse
import pandas as pd

def run_command(command):
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True)

def convert_narrowpeak_to_bed(narrowpeak_file, bed_file):
    print(f"Converting {narrowpeak_file} to BED format...")
    df = pd.read_csv(narrowpeak_file, sep="\t", header=None)
    bed = df.iloc[:, [0, 1, 2]]
    bed.to_csv(bed_file, sep="\t", header=False, index=False)

def run_homer_annotation(bed_file, genome, output_file):
    print(f"Running HOMER annotation for {bed_file} on genome {genome}...")
    command = f"annotatePeaks.pl {bed_file} {genome} > {output_file}"
    run_command(command)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Annotate MACS3 peaks using HOMER."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input MACS3 narrowPeak or BED file."
    )
    parser.add_argument(
        "-g",
        "--genome",
        default="hg38",
        help="Genome version (e.g., hg38, mm10)."
    )
    parser.add_argument(
        "-o",
        "--output",
        default="annotated_peaks.txt",
        help="Output annotated file."
    )
    return parser.parse_args()

def main():
    args = parse_arguments()

    input_file = args.input
    genome = args.genome
    output_file = args.output
    ext = os.path.splitext(input_file)[-1]

    # Convert if narrowPeak
    if ext == ".narrowPeak":
        bed_file = input_file.replace(".narrowPeak", ".bed")
        convert_narrowpeak_to_bed(input_file, bed_file)
    else:
        bed_file = input_file

    run_homer_annotation(bed_file, genome, output_file)
    print(f"Annotation completed. Output saved to {output_file}")

if __name__ == "__main__":
    main()


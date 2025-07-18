"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script automates de novo motif discovery from MACS3 peak files using HOMER's `findMotifsGenome.pl`. 
It reformats MACS3 `.narrowPeak` output by centering each peak region around the summit and resizing to a fixed 
window (default 200 bp), then passes the resulting BED file to HOMER for motif enrichment analysis.

The pipeline is suitable for ChIP-Seq or ATAC-Seq datasets and simplifies HOMER-based analysis into a single command.

Required Python Packages:
- os
- subprocess

External Tools Required:
- HOMER (findMotifsGenome.pl must be installed and configured in the system PATH)

Input Files:
- MACS3 narrowPeak file (e.g., `results/peaks/SRR26147696_peaks.narrowPeak`)

Output Files:
- Formatted BED file for HOMER input: `formatted_peaks.bed` (written to specified output directory)
- HOMER motif discovery results: written to the output directory (e.g., `results/motifs_SRR26147696/`)

Usage:
This script can be executed directly. Modify the example call under `if __name__ == "__main__":` to specify:
- Input peak file
- Reference genome (e.g., `hg38`, `mm10`)
- Output directory
- Peak window size around summits (e.g., 200 bp)
- Motif length (e.g., 8 for 8-mers)

Example:
```bash
python discover_motifs.py

Note:
Ensure HOMER and the appropriate genome package are installed and configured before running this script.
"""

import os
import subprocess

def create_output_dir(output_dir):
    """Create directory to store motif results."""
    os.makedirs(output_dir, exist_ok=True)
    print(f"[INFO] Output directory created at: {output_dir}")

def prepare_peak_file(input_peak_file, formatted_peak_file, peak_size=200):
    """
    Format MACS3 peak file to have fixed-size regions centered at summits.
    This uses awk to reformat each line for HOMER compatibility.
    """
    with open(formatted_peak_file, 'w') as fout, open(input_peak_file, 'r') as fin:
        for line in fin:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            summit = start + ((end - start) // 2)
            fixed_start = max(0, summit - peak_size // 2)
            fixed_end = fixed_start + peak_size
            fout.write(f"{chrom}\t{fixed_start}\t{fixed_end}\n")
    print(f"[INFO] Peaks formatted and written to: {formatted_peak_file}")

def run_homer_find_motifs(formatted_peak_file, genome, output_dir, motif_length=8):
    """
    Run HOMER to discover motifs from formatted peaks.
    """
    cmd = [
        "findMotifsGenome.pl",
        formatted_peak_file,
        genome,
        output_dir,
        "-len", str(motif_length),
        "-size", "given"
    ]
    print(f"[INFO] Running HOMER motif discovery...")
    try:
        subprocess.run(cmd, check=True)
        print(f"[INFO] HOMER motif discovery completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] HOMER failed: {e}")

def discover_motifs_from_macs3(input_peak_file, genome='hg38', output_dir='motif_results', peak_size=200, motif_length=8):
    """
    Full pipeline to discover motifs from MACS3 output using HOMER.
    """
    create_output_dir(output_dir)
    formatted_peak_file = os.path.join(output_dir, "formatted_peaks.bed")
    prepare_peak_file(input_peak_file, formatted_peak_file, peak_size)
    run_homer_find_motifs(formatted_peak_file, genome, output_dir, motif_length)

# Example usage
if __name__ == "__main__":
    discover_motifs_from_macs3(
        input_peak_file="results/peaks/SRR26147696_peaks.narrowPeak",
        genome="hg38",
        output_dir="results/motifs_SRR26147696",
        peak_size=200,
        motif_length=8
    )


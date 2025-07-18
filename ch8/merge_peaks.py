"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script merges all narrowPeak files from ChIP-Seq or ATAC-Seq experiments into a single non-redundant 
set of genomic intervals. It recursively scans the `data/peaks/` directory for MACS2 narrowPeak outputs, extracts 
chromosome, start, and end coordinates, and uses `pybedtools` to merge overlapping or adjacent peaks.

This merged BED file is essential for uniform peak-based quantification, motif discovery, and footprinting analysis.

Required Python Packages:
- os
- glob
- pybedtools

Input:
- Peak files: `data/peaks/*/*_peaks.narrowPeak` (one per sample/condition)

Output:
- Merged BED file: `data/peaks/merged_peaks.bed`

Usage:
```bash
python merge_peaks.py

Notes:
Ensure that pybedtools is installed and BEDTools is accessible in your system environment.
The script will raise an error if no peak files are found, so make sure the input directory structure is correct.
"""

import os
import glob
import pybedtools

PEAKS_DIR = "data/peaks"
OUTPUT_FILE = os.path.join(PEAKS_DIR, "merged_peaks.bed")

def get_peak_files(peak_dir, pattern="*_peaks.narrowPeak"):
    # Recursively find narrowPeak files in sample subfolders
    return sorted(glob.glob(os.path.join(peak_dir, "*", pattern)))

def merge_peak_files(peak_files, output_file):
    print(f"Found {len(peak_files)} peak files to merge.")
    
    if len(peak_files) == 0:
        raise FileNotFoundError("No peak files found. Check the folder structure and pattern.")

    # Combine all peaks (chr, start, end)
    combined = pybedtools.BedTool(peak_files[0]).cut([0, 1, 2])
    for peak_file in peak_files[1:]:
        bed = pybedtools.BedTool(peak_file).cut([0, 1, 2])
        combined = combined.cat(bed, postmerge=False)

    # Sort and merge overlapping/adjacent intervals
    merged = combined.sort().merge()
    merged.saveas(output_file)
    print(f"Merged peaks written to: {output_file}")

if __name__ == "__main__":
    peak_files = get_peak_files(PEAKS_DIR)
    merge_peak_files(peak_files, OUTPUT_FILE)


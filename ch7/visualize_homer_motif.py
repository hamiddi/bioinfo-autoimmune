"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script visualizes position weight matrices (PWMs) of transcription factor binding motifs 
discovered by HOMER. It parses a HOMER motif file (usually located in `homerResults/motif*.motif`), 
extracts the consensus sequence, optional gene name, and base frequency matrix, then generates a 
stacked bar chart showing the relative frequency of A, C, G, and T at each position in the motif.

Required Python Packages:
- sys
- numpy
- re
- matplotlib.pyplot

Input:
- HOMER motif file (text file containing motif header and PWM, e.g., `motif1.motif`)

Output:
- Interactive bar chart showing the motif's base composition across positions

Usage:
```bash
python visualize_homer_motif.py motif1.motif

Notes:

This script assumes the motif file follows HOMER's format: a header line beginning with > followed
by multiple lines of four numeric values representing base frequencies for A, C, G, and T.

If a gene name is embedded in the motif header, it is extracted and shown in the plot title.
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import re

def parse_motif_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    header_line = ""
    matrix = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header_line = line[1:].strip()
            continue
        values = list(map(float, line.split()))
        if len(values) == 4:
            matrix.append(values)

    pwm = np.array(matrix)
    
    # Parse consensus and gene name if found
    parts = header_line.split()
    consensus = parts[0]
    gene_name = ""
    if len(parts) > 1:
        match = re.search(r'([A-Za-z0-9_-]+)\(', parts[1])
        if match:
            gene_name = match.group(1)

    return consensus, gene_name, pwm

def plot_pwm(consensus, gene_name, pwm):
    positions = np.arange(pwm.shape[0])
    bases = ['A', 'C', 'G', 'T']
    colors = ['green', 'blue', 'orange', 'red']
    bar_width = 0.2

    plt.figure(figsize=(12, 6))
    for i in range(4):
        plt.bar(positions + i * bar_width, pwm[:, i], width=bar_width, color=colors[i], label=bases[i])

    plt.xticks(positions + bar_width * 1.5, [f"Pos {i+1}" for i in range(pwm.shape[0])], fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Motif Position", fontsize=14)
    plt.ylabel("Base Frequency", fontsize=14)

    title = f"Motif: {consensus}"
    if gene_name:
        title += f" ({gene_name})"
    plt.title(title, fontsize=16)

    plt.legend(title="Base", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

def main():
    if len(sys.argv) != 2:
        print("Usage: python visualize_homer_motif.py <motif_file>")
        sys.exit(1)

    motif_file = sys.argv[1]
    consensus, gene_name, pwm = parse_motif_file(motif_file)
    plot_pwm(consensus, gene_name, pwm)

if __name__ == "__main__":
    main()


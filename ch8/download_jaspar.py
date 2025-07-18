"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script downloads the latest non-redundant transcription factor position frequency matrices (PFMs) 
from the JASPAR 2024 CORE database in MEME format. These motifs are widely used in transcription factor 
binding site prediction and footprinting analyses (e.g., with TOBIAS, FIMO, or HOMER).

The script retrieves the file from the official JASPAR FTP server and saves it locally under a `motifs/` directory.

Required Python Packages:
- os
- urllib.request

Input:
- Remote URL: JASPAR2024 CORE non-redundant PFMs in MEME format

Output:
- Local file: `motifs/jaspar.meme`

Usage:
Run the script to download and store the motif file for downstream analysis.
```bash
python download_jaspar_meme.py

Note:
Ensure internet access is available when running this script. The downloaded file can be used with tools such as
TOBIAS (BINDetect), HOMER (findMotifsGenome.pl), or MEME Suite applications for motif enrichment analysis.
"""

import os
import urllib.request

def download_jaspar_meme():
    motifs_dir = "motifs"
    url = "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_meme.txt"
    os.makedirs(motifs_dir, exist_ok=True)
    output_path = os.path.join(motifs_dir, "jaspar.meme")

    print(f"Downloading JASPAR MEME motifs from:\n{url}")
    urllib.request.urlretrieve(url, output_path)
    print(f"Saved motifs as: {output_path}")

if __name__ == "__main__":
    download_jaspar_meme()


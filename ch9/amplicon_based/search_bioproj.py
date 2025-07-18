"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script automates the discovery of amplicon-based BioProjects related to Crohn’s disease in NCBI, 
identifies associated SRA runs, and checks whether they are paired-end. It uses Biopython's Entrez module for 
NCBI queries and the Requests library to retrieve run information from the SRA database.

The script prints summaries of BioProjects, fetches linked SRA accessions, and reports which ones include 
paired-end sequencing—information essential for downstream microbiome or metagenomic analyses.

Required Python Packages:
- time
- requests
- Bio.Entrez (from Biopython)

Input:
- NCBI Entrez search term (currently: `'amplicon AND Crohn "bioproject"[Filter]'`)

Output:
- Console output listing:
  - BioProject ID, title, and description
  - Links to paired-end SRA runs (if found)

Usage:
```bash
python find_paired_end_sra_crohn.py

Notes:
Replace Entrez.email with your own email address to comply with NCBI's access policy.
The script includes a delay (time.sleep(0.4)) to avoid exceeding NCBI’s API rate limits.
You may increase retmax to retrieve more BioProjects or modify the query for other conditions or assays.
"""

import time
import requests
from Bio import Entrez

Entrez.email = "your.email@example.com"  # Replace with your email

# Search parameters
search_term = 'amplicon AND Crohn  "bioproject"[Filter]'

def search_bioprojects(term, retmax=10):
    handle = Entrez.esearch(db="bioproject", term=term, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def get_bioproject_summary(bioproject_id):
    handle = Entrez.esummary(db="bioproject", id=bioproject_id)
    summary = Entrez.read(handle)
    handle.close()
    return summary[0]

def fetch_sra_links(bioproject_id):
    handle = Entrez.elink(dbfrom="bioproject", db="sra", id=bioproject_id)
    record = Entrez.read(handle)
    handle.close()
    sra_ids = []
    for linkset in record:
        for link in linkset.get("LinkSetDb", []):
            for link_info in link["Link"]:
                sra_ids.append(link_info["Id"])
    return sra_ids

def check_for_paired_end(sra_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "sra",
        "id": sra_id,
        "rettype": "runinfo",
        "retmode": "text"
    }
    response = requests.get(url, params=params)
    if response.status_code == 200:
        lines = response.text.strip().split("\n")
        header = lines[0].split(",")
        if "LibraryLayout" in header:
            layout_index = header.index("LibraryLayout")
            for line in lines[1:]:
                parts = line.split(",")
                if len(parts) > layout_index and parts[layout_index].lower() == "paired":
                    return True
    return False

def main():
    bioproject_ids = search_bioprojects(search_term, retmax=20)
    for bp_id in bioproject_ids:
        summary = get_bioproject_summary(bp_id)
        print(f"\nBioProject ID: {bp_id}")
        print(f"Title: {summary.get('Title')}")
        print(f"Description: {summary.get('Description')}")

        sra_ids = fetch_sra_links(bp_id)
        paired_end_runs = []
        for sra_id in sra_ids:
            if check_for_paired_end(sra_id):
                paired_end_runs.append(sra_id)
            time.sleep(0.4)  # To avoid hitting API limits

        if paired_end_runs:
            print("Found paired-end SRA runs:")
            for run in paired_end_runs:
                print(f" - https://www.ncbi.nlm.nih.gov/sra/?term={run}")
        else:
            print("No paired-end SRA runs found for this BioProject.")

if __name__ == "__main__":
    main()


"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
##########################################################
This script retrieves a list of NCBI Gene IDs related to specific autoimmune.

Usage:
- Replace 'your_email@example.com' with your own email address.
- Call the function `fetch_gene_ids(disease, output_file)` with:
    - `disease`: the name of the disease to search for (e.g., "multiple sclerosis")
    - `output_file`: the name of the file where the retrieved gene IDs will be saved

Example:
    fetch_gene_ids("multiple sclerosis", "MS_gene_ids.txt")
    fetch_gene_ids("Systemic Lupus Erythematosus", "SLE_gene_ids.txt")
"""
from Bio import Entrez
def fetch_gene_ids(disease, output_file):
    Entrez.email = "your_email@example.com"
    query = f"{disease}[Title] AND human[ORGN] AND alive[prop]"
    handle = Entrez.esearch(db="gene", term=query, retmax=100)
    record = Entrez.read(handle)
    handle.close()
    gene_ids = record["IdList"]
    with open(output_file, "w") as f:
        for gene_id in gene_ids:
            f.write(f"{gene_id}\n")
fetch_gene_ids("multiple sclerosis", "MS_gene_ids.txt")
fetch_gene_ids("Systemic Lupus Erythematosus", "SLE_gene_ids.txt")


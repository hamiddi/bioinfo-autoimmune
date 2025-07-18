"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This script uses Biopython's Entrez module to search the NCBI BioProject database for publicly available 
ATAC-seq datasets related to systemic lupus in humans. It filters results to include only BioProjects 
with SRA data (bioproject_sra[filter]) and exclude controlled-access datasets (bioproject_gap[filter]).

Required Python Package:
- Biopython (Entrez module)

Before Use:
- Replace 'your.email@example.com' with a valid email address. NCBI requires this for API access.

Input:
- Query string: `(ATAC-seq AND human[ORGN] AND systemic lupus) AND bioproject_sra[filter] NOT bioproject_gap[filter]`

Output:
- Printed list of BioProject accessions (e.g., PRJNAxxxxxx) that match the query criteria

Note:
- This script is limited to retrieving up to 1000 BioProjects. Modify `retmax` to increase/decrease results.
- Network access is required. Ensure firewall or proxy settings do not block Entrez queries.
"""

from Bio import Entrez

# Set your email (NCBI requires this)
Entrez.email = "your.email@example.com"

# Define the search query
query = (
    '(ATAC-seq AND human[ORGN] AND systemic lupus) '
    'AND bioproject_sra[filter] '
    'NOT bioproject_gap[filter]'
)
# Search the BioProject database
handle = Entrez.esearch(db="bioproject", term=query, retmax=1000)
record = Entrez.read(handle)
handle.close()

# Retrieve BioProject IDs
bioproject_ids = record["IdList"]

# Fetch summary to get accessions
if bioproject_ids:
    handle = Entrez.esummary(db="bioproject", id=",".join(bioproject_ids))
    records = Entrez.read(handle)
    handle.close()

    print("BioProject Accessions:")
    for project in records["DocumentSummarySet"]["DocumentSummary"]:
        print(project["Project_Acc"])
else:
    print("No BioProjects found for the query.")


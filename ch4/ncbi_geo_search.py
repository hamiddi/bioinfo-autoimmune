"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autimmune Diseases
############################################################
This script searches the NCBI GEO (Gene Expression Omnibus) database for RNA-Seq datasets 
related to autoimmune diseases and retrieves corresponding dataset details and associated 
SRA (Sequence Read Archive) IDs.

Usage Instructions:
1. Make sure you have Biopython installed: pip install biopython
2. Replace 'your_email@example.com' with your actual email address to comply with NCBI's API policy.
3. Run the script in a Python environment with internet access.
4. The script will:
   - Search GEO for relevant datasets based on the query: "RNA-Seq AND autoimmune disease"
   - Save the retrieved GEO dataset IDs to a file named 'geo_ids.txt'
   - Fetch metadata (title, summary, accession, type) and SRA IDs for each dataset
   - Save this detailed information to 'geo_details.txt'

Note:
- The script respects NCBI's rate limits by including a short delay between API requests.
- You can modify the search query or retmax value to suit different research needs.
"""
from Bio import Entrez
import time
# Set email for NCBI access
Entrez.email = "your_email@example.com"
# Function to search GEO
def search_geo(query, db="gds", retmax=100):
    handle = Entrez.esearch(db=db, term=query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]
# Function to fetch GEO dataset details
def fetch_geo_details(geo_id):
    handle = Entrez.esummary(db="gds", id=geo_id)
    record = Entrez.read(handle)
    handle.close()
    return record[0] if record else None
# Function to find SRA ID from GEO accession
def find_sra_ids(geo_accession):
    query = f"{geo_accession} AND transcriptomic"
    handle = Entrez.esearch(db="sra", term=query, retmax=100)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

# Search GEO for RNA-Seq datasets related to autoimmune disease
query = "RNA-Seq AND Systemic Lupus Erythematosus AND human[ORGN]"
geo_ids = search_geo(query, retmax=100)
# Save GEO IDs to a file
with open("geo_ids.txt", "w") as f:
    for geo_id in geo_ids:
        f.write(geo_id + "\n")

print(f"Found {len(geo_ids)} GEO dataset IDs. Saved to geo_ids.txt")
# Fetch dataset details and find SRA IDs
geo_details_list = []
sra_mappings = {}
for geo_id in geo_ids:
    details = fetch_geo_details(geo_id)
    if details:
        accession = details["Accession"]
        sra_ids = find_sra_ids(accession)
        sra_mappings[accession] = sra_ids

        geo_details_list.append(f"Dataset ID: {geo_id}\n"
                                f"Title: {details['title']}\n"
                                f"Summary: {details['summary']}\n"
                                f"Accession: {accession}\n"
                                f"Type: {details['gdsType']}\n"
                                f"SRA IDs: {', '.join(sra_ids) if sra_ids else 'None'}\n"
                                "--------------------------------------------------\n")
    
    time.sleep(1)  # To avoid NCBI rate limits

# Save dataset details to a file
with open("geo_details.txt", "w", encoding="utf-8") as f:
    f.writelines(geo_details_list)

print(f"Saved dataset details to geo_details.txt")


"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script fetches a list of reviewed proteins associated with a given disease name 
from the UniProt database using the UniProt REST API.

Usage Instructions:
1. Make sure you have internet access and the required Python packages installed:
   - requests
   - pandas

2. Edit the 'disease_name' variable near the bottom of the script to specify the disease 
   you are interested in (e.g., "Systemic Lupus Erythematosus").

3. Run the script. It will:
   - Query UniProt for reviewed proteins related to the disease.
   - Extract relevant data including UniProt ID, protein name, gene names, and organism.
   - Save the results as a CSV file named '<disease_name>_proteins.csv'.
   - Print the resulting DataFrame to the console.

Note:
- Only manually reviewed (Swiss-Prot) entries are retrieved.
- The maximum number of returned entries is limited by the 'size' parameter (currently set to 50).
"""

import requests
import pandas as pd

# Function to fetch proteins associated with a disease
def fetch_proteins_by_disease_name(disease_name):
    """
    Fetches proteins associated with a given disease name from UniProt.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f'"{disease_name}" AND reviewed:true',  # Fetch only reviewed (Swiss-Prot) entries
        "fields": "accession,protein_name,gene_primary,organism_name",
        "format": "json",
        "size": 50  # Adjust size as needed
    }

    response = requests.get(base_url, params=params)

    if response.status_code == 200:
        data = response.json()
        results = []

        for entry in data.get("results", []):
            # Extract UniProt ID
            protein_id = entry.get("primaryAccession", "N/A")

            # Extract Protein Name correctly
            protein_name = "N/A"
            if "proteinDescription" in entry:
                recommended_name = entry["proteinDescription"].get("recommendedName", {})
                if isinstance(recommended_name, dict) and "fullName" in recommended_name:
                    full_name = recommended_name["fullName"]
                    if isinstance(full_name, dict) and "value" in full_name:
                        protein_name = full_name["value"]
                    else:
                        protein_name = full_name

            # Extract Gene Names correctly
            gene_names = "N/A"
            if "genes" in entry and entry["genes"]:
                gene_list = [gene["geneName"]["value"] for gene in entry["genes"] if "geneName" in gene]
                gene_names = ", ".join(gene_list) if gene_list else "N/A"

            # Extract Organism Name
            organism = entry.get("organism", {}).get("scientificName", "N/A")

            results.append([protein_id, protein_name, gene_names, organism])

        # Convert to Pandas DataFrame with column names using underscores
        df = pd.DataFrame(results, columns=["UniProt_ID", "Protein_Name", "Gene_Names", "Organism"])
        return df

    else:
        print("Error fetching data:", response.status_code, response.text)
        return None


# Example Usage
disease_name = "Systemic Lupus Erythematosus"  # Replace with any disease
df_proteins = fetch_proteins_by_disease_name(disease_name)

# Save results locally
if df_proteins is not None:
    df_proteins.to_csv(f"{disease_name.replace(' ', '_')}_proteins.csv", index=False)
    print(f"Data saved to {disease_name.replace(' ', '_')}_proteins.csv")
    print(df_proteins)  # Display in console


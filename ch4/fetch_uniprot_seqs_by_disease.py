"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script fetches reviewed protein sequences from the UniProt database that are associated with a specified disease 
and saves them in FASTA format.

Usage:
- Replace the value of `disease_name` with the name of the disease you want to search for (e.g., "Rheumatoid Arthritis").
- Optionally, change the `output_file` name in the function call to set a custom filename for the FASTA output.
- Run the script. If any reviewed protein sequences are found for the given disease, they will be printed and saved.

Requirements:
- Python 3.x
- The `requests` library (install via `pip install requests`)

Example:
    disease_name = "Systemic Lupus Erythematosus"
    fetch_fasta_by_disease_name(disease_name, output_file="lupus_proteins.fasta")
"""

import requests

# Function to fetch protein sequences and save them in FASTA format
def fetch_fasta_by_disease_name(disease_name, output_file="proteins.fasta"):
    """
    Fetches protein sequences for a given disease from UniProt and saves them in FASTA format.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f'"{disease_name}" AND reviewed:true',
        "fields": "accession,protein_name,gene_primary,organism_name,sequence",
        "format": "json",
        "size": 50  # Adjust as needed
    }

    response = requests.get(base_url, params=params)

    if response.status_code == 200:
        data = response.json()
        fasta_entries = []

        for entry in data.get("results", []):
            # Extract UniProt ID
            protein_id = entry.get("primaryAccession", "N/A")

            # Extract Protein Name
            protein_name = "N/A"
            if "proteinDescription" in entry:
                recommended_name = entry["proteinDescription"].get("recommendedName", {})
                if isinstance(recommended_name, dict) and "fullName" in recommended_name:
                    protein_name = recommended_name["fullName"]
                    if isinstance(protein_name, dict) and "value" in protein_name:
                        protein_name = protein_name["value"]

            # Extract Gene Name
            gene_name = "N/A"
            if "genes" in entry and entry["genes"]:
                gene_list = [gene["geneName"]["value"] for gene in entry["genes"] if "geneName" in gene]
                gene_name = ", ".join(gene_list) if gene_list else "N/A"

            # Extract Sequence and Length
            sequence = entry.get("sequence", {}).get("value", None)  # Corrected extraction
            seq_length = entry.get("sequence", {}).get("length", None)  # Extract length correctly

            # Debug: Print sequences that should be valid
            if sequence:
                print(f"FOUND SEQUENCE: {protein_id} | Length: {seq_length}")

            # Skip entries without valid sequences
            if not sequence or not seq_length:
                print(f"Skipping {protein_id} due to missing sequence.")
                continue  

            # Create FASTA defline
            defline = f">{protein_id}|{protein_name}|{gene_name}|{seq_length}"
            fasta_entries.append(f"{defline}\n{sequence}\n")

        # Save to FASTA file
        if fasta_entries:
            with open(output_file, "w") as fasta_file:
                fasta_file.writelines(fasta_entries)
            print(f"FASTA sequences saved to {output_file}")
        else:
            print("No valid protein sequences found for this disease.")

    else:
        print("Error fetching data:", response.status_code, response.text)


# Example Usage
disease_name = "Systemic Lupus Erythematosus"  # Replace with any disease
fetch_fasta_by_disease_name(disease_name, output_file="lupus_proteins.fasta")


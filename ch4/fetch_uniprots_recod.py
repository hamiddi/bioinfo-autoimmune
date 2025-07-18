"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script fetches the UniProt record for a given UniProt accession ID (e.g., HLA-A: P04439)
using the UniProt REST API and saves the result as a JSON file.

How to use:
- Ensure you have internet access and the 'requests' library installed.
- Modify the 'uniprotID' variable to the desired UniProt accession ID.
- Optionally, change the output filename in the function call.
- Run the script, and it will save the record as a JSON file in the current directory.

Example:
    uniprotID = "P04439"
    fetch_hla_a_uniprot(uniprotID, "hla_a.json")
"""
import requests
import json
def fetch_hla_a_uniprot(uniprotID, output_file="uniprot.json"):
    # UniProt API URL for HLA-A (P04439)
    url = f"https://rest.uniprot.org/uniprotkb/{uniprotID}.json"
    try:
        # Fetch data from UniProt API
        response = requests.get(url)
        response.raise_for_status()# Raise an error for HTTP failures
        # Parse JSON response
        data = response.json()
        # Save JSON data to a file
        with open(output_file, "w", encoding="utf-8") as json_file:
            json.dump(data, json_file, indent=4)
        print(f"{uniprotID} record saved successfully as '{output_file}'.")
    except requests.exceptions.RequestException as e:
        print(f"Failed to fetch data from UniProt: {e}")
# Run the function
uniprotID = "P04439"
fetch_hla_a_uniprot(uniprotID, "hla_a.json")


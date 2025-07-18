"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script retrieves Protein Data Bank (PDB) entries associated with a given UniProt ID.

Usage:
- Edit the `uniprot_id` variable at the bottom of the script with your desired UniProt protein accession.
- Run the script.
- If any PDB cross-references are found for the provided UniProt ID, they will be saved to a text file
  named <UniProt_ID>.txt in the current directory.
- If no entries are found or the ID is invalid, an appropriate message will be printed.

Example:
    uniprot_id = "P69905"  # Hemoglobin subunit alpha

Dependencies:
    - requests

Make sure the 'requests' package is installed:
    pip install requests
"""

import requests

def map_uniprot_to_pdb(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Error: Unable to retrieve data for UniProt ID {uniprot_id}")
        return None
    
    data = response.json()
    
    # Extract PDB information if available
    pdb_entries = []
    if 'uniProtKBCrossReferences' in data:
        for entry in data['uniProtKBCrossReferences']:
            if entry['database'] == 'PDB':
                pdb_entries.append(entry['id'])
    
    if pdb_entries:
        pdb_list = '\n'.join(pdb_entries)
        filename = f"{uniprot_id}.txt"
        
        with open(filename, "w") as file:
            file.write(f"PDB IDs for UniProt ID {uniprot_id}:\n{pdb_list}")
        
        print(f"PDB IDs saved to {filename}")
    else:
        print(f"No PDB entries found for UniProt ID {uniprot_id}")

uniprot_id = "P69905"  # Example: Hemoglobin subunit alpha (HBA1)
map_uniprot_to_pdb(uniprot_id)


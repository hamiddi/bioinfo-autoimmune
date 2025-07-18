"""
Book title: Bioinformatics of Autoimmune Diseases
Author: Hamid D. Ismail, Ph.D.
############################################################
This script retrieves OMIM (Online Mendelian Inheritance in Man) phenotype IDs associated with a given autoimmune disease
by querying the NCBI E-Utilities API. It saves the resulting IDs into a text file for downstream analysis.

How to use:
1. Set the `disease_name` variable near the bottom of the script to the name of the autoimmune disease you are interested in.
2. Run the script.
3. The OMIM phenotype IDs associated with the disease will be saved in a file named 'omim_phenotype_ids.txt' in the current directory.

Requirements:
- Python 3
- requests library (install via `pip install requests` if not already installed)
"""

import requests

def get_omim_phenotype_ids(disease_name):
    """
    Fetches the list of OMIM phenotype IDs associated with a specific autoimmune disease from NCBI E-Utils.
    
    Parameters:
        disease_name (str): The name of the autoimmune disease.

    Returns:
        list: A list of OMIM phenotype IDs related to the disease.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "omim",
        "term": f"{disease_name} AND (Phenotype[tiab])",
        "retmode": "json",
        "retmax": 1000  # Adjust to retrieve more results if needed
    }
    
    response = requests.get(base_url, params=params)
    
    if response.status_code == 200:
        data = response.json()
        return data.get("esearchresult", {}).get("idlist", [])
    else:
        print("Error fetching data:", response.status_code)
        return []

def save_omim_ids_to_file(omim_ids, filename="omim_phenotype_ids.txt"):
    """
    Saves the list of OMIM phenotype IDs to a file, with each ID on a new line.

    Parameters:
        omim_ids (list): List of OMIM phenotype IDs.
        filename (str): The name of the file to save the IDs.
    """
    with open(filename, "w") as file:
        for omim_id in omim_ids:
            file.write(f"{omim_id}\n")
    print(f"OMIM phenotype IDs saved to {filename}")

# Example usage
disease_name = "Celiac Disease"
omim_phenotype_ids = get_omim_phenotype_ids(disease_name)

if omim_phenotype_ids:
    save_omim_ids_to_file(omim_phenotype_ids)
else:
    print("No OMIM phenotype IDs found.")


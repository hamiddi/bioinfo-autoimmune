"""
Book title: Bioinformatics of Autoimmune Diseases
Author: Hamid D. Ismail, Ph.D.
############################################################
This script fetches structural variant data associated with a given gene symbol from the NCBI dbVar database
using NCBI E-Utilities (ESearch and ESummary). It outputs a CSV file containing detailed information about 
each variant.

Usage:
1. Set the `gene_symbol` variable at the end of the script to the gene of interest (e.g., "HLA-A").
2. Run the script. It will:
   - Search dbVar for variants linked to the gene symbol.
   - Retrieve summary details for each variant.
   - Save the data into a CSV file (default: "dbvar_details.csv").

Dependencies:
- This script requires an active internet connection.
- Make sure the 'requests' and 'csv' libraries are available in your Python environment.

Output:
- A CSV file named "dbvar_details.csv" containing variant details for the specified gene.
"""

import requests
import csv

def get_dbvar_numbers(gene_symbol):
    """
    Fetch dbVar numbers for a given gene symbol using NCBI E-Utils API.
    
    Parameters:
        gene_symbol (str): The gene symbol to search for.
    
    Returns:
        list: A list of dbVar IDs associated with the gene symbol.
    """
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "dbvar",
        "term": f"{gene_symbol}[gene]",
        "retmode": "json"
    }
    
    response = requests.get(url, params=params)
    
    if response.status_code == 200:
        data = response.json()
        return data.get("esearchresult", {}).get("idlist", [])
    else:
        return []

def clean_dict_values(value):
    """Convert dictionaries into string format, replacing internal commas with semicolons."""
    if isinstance(value, dict):
        return ";".join(f"{k}:{v}".replace(",", ";") for k, v in value.items())
    if isinstance(value, list):
        return "|".join(clean_dict_values(item) for item in value) if value else "N/A"
    return str(value).replace(",", ";") if value not in [None, ""] else "N/A"

def fetch_dbvar_details(dbvar_ids, output_csv="dbvar_details.csv"):
    """
    Fetch details of each dbVar using the dbVar ID and store in a CSV file.
    
    Parameters:
        dbvar_ids (list): List of dbVar IDs to fetch details for.
        output_csv (str): Output CSV file name.
    """
    fields = [
        "Variant_Region_ID", "Obj_Type", "Study", "Study_Type", "Variant_Count", "Publications",
        "Species", "Organism", "Submitted_Assemblies", "Remapped_Assemblies", "Placements",
        "Genes", "Methods", "Clinical_Significance", "Variant_Types", "Variant_Call_Count"
    ]
    
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(fields)
        
        for dbvar_id in dbvar_ids:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            params = {
                "db": "dbvar",
                "id": dbvar_id,
                "retmode": "json"
            }
            
            response = requests.get(url, params=params)
            
            if response.status_code == 200:
                data = response.json()
                summary = data.get("result", {}).get(dbvar_id, {})
                
                # Filter out records where Obj_Type is not VARIANT
                if summary.get("obj_type", "N/A") != "VARIANT":
                    continue
                
                row = [
                    clean_dict_values(summary.get("uid", "N/A")),
                    clean_dict_values(summary.get("obj_type", "N/A")),
                    clean_dict_values(summary.get("st", "N/A")),
                    clean_dict_values(summary.get("study_type", "N/A")),
                    clean_dict_values(summary.get("variant_count", "N/A")),
                    clean_dict_values(summary.get("dbvarpublicationlist", [])),
                    clean_dict_values(summary.get("dbvarstudyorglist", [])),
                    clean_dict_values(summary.get("organism", "N/A")),
                    clean_dict_values(summary.get("dbvarsubmittedassemblylist", [])),
                    clean_dict_values(summary.get("dbvarremappedassemblylist", [])),
                    clean_dict_values(summary.get("dbvarplacementlist", [])),
                    clean_dict_values(summary.get("dbvargenelist", [])),
                    clean_dict_values(summary.get("dbvarmethodlist", [])),
                    clean_dict_values(summary.get("dbvarclinicalsignificancelist", [])),
                    clean_dict_values(summary.get("dbvarvarianttypelist", [])),
                    clean_dict_values(summary.get("variant_call_count", "N/A"))
                ]
                writer.writerow(row)

# Example usage
gene_symbol = "HLA-A"
dbvar_ids = get_dbvar_numbers(gene_symbol)
fetch_dbvar_details(dbvar_ids)


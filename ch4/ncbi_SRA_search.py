"""
Book title: Bioinformatics of Autoimmune Diseases
Author: Hamid D. Ismail, Ph.D.
############################################################
This script allows users to search the NCBI SRA (Sequence Read Archive) database
for sequencing studies related to a specific disease and sequencing method(s),
and fetches relevant metadata for the resulting SRA entries.

Usage:
1. Make sure you have Biopython installed: pip install biopython
2. Replace 'your_email@example.com' with your actual email address as required by NCBI Entrez.
3. Modify the 'disease' and 'seq_types' variables in the `main()` function as needed.
   - Example disease: "rheumatoid arthritis"
   - Example sequencing types: ["RNA-Seq", "ChIP-Seq"]
4. Run the script using: python script_name.py
5. The program will print metadata for up to 10 matching SRA studies.

Note: This program uses NCBI E-Utilities, so ensure you are connected to the internet.
"""

from Bio import Entrez
import xml.etree.ElementTree as ET
# Set email (required by NCBI)
Entrez.email = "your_email@example.com"
def search_sra(disease, seq_types, max_results=10):
    """Search the SRA database for studies related to a disease and sequencing type."""
    #query = f'("{disease}"[All Fields]) AND ({" OR ".join(seq_types)}[All Fields])'
    query = (
         f'("{disease}"[All Fields]) AND ('
         f'{" OR ".join(seq_types)}[All Fields]'
         f')'
        )
    handle = Entrez.esearch(db="sra", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_sra_metadata(sra_id):
    """Fetch metadata for a given SRA ID and extract relevant fields."""
    handle = Entrez.efetch(db="sra", id=sra_id, rettype="full", retmode="xml")
    xml_data = handle.read()
    handle.close()
    root = ET.fromstring(xml_data)
    metadata = {
        "SRA Accession": sra_id,
        "Title": extract_text(root,
            ".//STUDY/DESCRIPTOR/STUDY_TITLE"),
        "Organism": extract_text(root,
            ".//SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME"),
        "Library Strategy": extract_text(root,
            ".//EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY"),
        "Library Source": extract_text(root,
            ".//EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE"),
        "Experiment Accession": extract_text(root,
            ".//EXPERIMENT/IDENTIFIERS/PRIMARY_ID"),
        "Run Accession": extract_text(root,
            ".//RUN_SET/RUN/IDENTIFIERS/PRIMARY_ID")
    }
    return metadata

def extract_text(root, xpath):
    """Helper function to safely extract text from an XML element."""
    element = root.find(xpath)
    return element.text if element is not None else "Not Available"

def main():
    disease = "rheumatoid arthritis"  # Modify as needed
    seq_types = ["RNA-Seq",
            "small RNA OR Non-coding RNA OR ncRNA OR MicroRNA OR miRNA",
            "ChIP-Seq"]  # Modify to include/exclude methods
    print("Searching SRA database...")
    sra_ids = search_sra(disease, seq_types)
    if not sra_ids:
        print("No results found.")
        return
    print("\nSRA Search Results:")
    for sra_id in sra_ids:
        metadata = fetch_sra_metadata(sra_id)
        print("\n=== SRA Study Metadata ===")
        for key, value in metadata.items():
            print(f"{key}: {value}")
if __name__ == "__main__":
    main()


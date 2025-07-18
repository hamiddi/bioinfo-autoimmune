"""
Book title: Bioinformatics of Autoimmune Diseases
Author: Hamid D. Ismail, Ph.D.
############################################################
This script retrieves and summarizes dbSNP variant information associated with a specific NCBI Gene ID.
It connects to the NCBI Entrez API to search for SNPs linked to the gene, then fetches detailed annotations
for each SNP, including variant type, allele composition, chromosome location, minor allele frequency (MAF),
functional consequence, gene name, and HGVS notation.

USAGE:
- Set your email address in the `Entrez.email` field to comply with NCBI's usage policy.
- Modify the `gene_id` variable in the `main()` function or uncomment the input() line to provide your own Gene ID.
- Run the script. It will generate a CSV file summarizing all retrieved variant data.

NOTE:
- Requires Biopython to be installed (`pip install biopython`).
- Ensure internet access is available as it uses NCBI Entrez web services.
- The output CSV file will be saved in the current directory with a filename like `dbSNP_report_Gene_<GeneID>.csv`.
"""
import time
import csv
import re
import xml.etree.ElementTree as ET
from Bio import Entrez

# Set your email for NCBI Entrez API requests
Entrez.email = "your_email@example.com"

# Define the namespace for XML parsing
NS = {"ns": "https://www.ncbi.nlm.nih.gov/SNP/docsum"}

def fetch_dbSNP(gene_id):
    """Fetches dbSNP variants linked to the given Gene ID."""
    try:
        print(f"Fetching dbSNP entries for Gene ID {gene_id}...")
        handle = Entrez.esearch(db="snp", term=f"{gene_id}[GeneID]", retmax=1000)
        record = Entrez.read(handle)
        handle.close()

        snp_ids = record["IdList"]
        return [f"rs{snp}" for snp in snp_ids]  

    except Exception as e:
        print(f"Error fetching dbSNP data: {e}")
        return []

def extract_text(root, path):
    """Extracts text from XML for a given path, handling namespaces and missing data."""
    element = root.find(path, NS)
    return element.text.strip() if element is not None and element.text else "N/A"

def extract_multiple_texts(root, path):
    """Extracts multiple values from XML, joins them with '|', handling namespaces."""
    elements = root.findall(path, NS)
    values = [elem.text.strip().replace(",", "|") for elem in elements if elem.text]
    return "|".join(values) if values else "N/A"

def extract_maf(root):
    """Extracts Minor Allele Frequency (MAF) from XML and uses '|' as separator."""
    maf_elements = root.findall(".//ns:GLOBAL_MAFS/ns:MAF/ns:FREQ", NS)
    maf_values = []
    for maf in maf_elements:
        match = re.search(r"=(\d+\.\d+)", maf.text)  # Extract numerical MAF
        if match:
            maf_values.append(match.group(1))
    return "|".join(maf_values) if maf_values else "N/A"

def extract_alleles(root):
    """Extracts allele information from the DOCSUM field."""
    docsum = extract_text(root, ".//ns:DOCSUM")
    match = re.search(r"SEQ=\[(.*?)\]", docsum)  # Extract allele info
    return match.group(1).replace(",", "|") if match else "N/A"

def extract_gene_name(root):
    """Extracts only the first gene name associated with the SNP."""
    genes = [gene.find("ns:NAME", NS).text for gene in root.findall(".//ns:GENE_E", NS) if gene.find("ns:NAME", NS) is not None]
    return genes[0] if genes else "N/A"  # Keep only the first gene

def extract_functional_consequence(root):
    """Extracts Functional Consequence and ensures values are separated by '|', not ','."""
    fxn_class = extract_text(root, ".//ns:FXN_CLASS")
    return fxn_class.replace(",", "|") if fxn_class != "N/A" else "N/A"

def extract_hgvs(root):
    """Extracts HGVS annotations from DOCSUM, replaces ',' with '|'."""
    docsum = extract_text(root, ".//ns:DOCSUM")
    matches = re.findall(r"HGVS=([^|]+)", docsum)  # Extract HGVS portion
    return "|".join(matches).replace(",", "|") if matches else "N/A"

def fetch_variant_details(rs_id):
    """Fetches detailed information for a given dbSNP ID (rsID)."""
    try:
        handle = Entrez.efetch(db="snp", id=rs_id, retmode="xml")
        xml_data = handle.read().decode("utf-8")
        handle.close()

        root = ET.fromstring(xml_data)

        variant_data = {
            "rsID": rs_id,
            "Variant_Type": extract_text(root, ".//ns:SNP_CLASS"),
            "Alleles": extract_alleles(root),
            "Chromosome": extract_text(root, ".//ns:CHR"),
            "Canonical_SPDI": extract_multiple_texts(root, ".//ns:SPDI"),
            "Gene": extract_gene_name(root),
            "Functional_Consequence": extract_functional_consequence(root),
            "Validated": "Yes" if root.find(".//ns:VALIDATED", NS) is not None else "No",
            "MAF": extract_maf(root),
            "HGVS": extract_hgvs(root),
        }

        return variant_data

    except Exception as e:
        print(f"Error fetching details for {rs_id}: {e}")
        return None

def generate_report(gene_id):
    """Generates a CSV report for all dbSNP IDs associated with a gene."""
    dbsnp_ids = fetch_dbSNP(gene_id)

    if not dbsnp_ids:
        print("No dbSNPs found for the given Gene ID.")
        return

    report_filename = f"dbSNP_report_Gene_{gene_id}.csv"
    
    with open(report_filename, mode="w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=[
            "rsID", "Variant_Type", "Alleles", "Chromosome", "Canonical_SPDI", "Gene",
            "Functional_Consequence", "Validated", "MAF", "HGVS"
        ])
        writer.writeheader()

        for rs_id in dbsnp_ids:
            print(f"Processing {rs_id}...")
            variant_data = fetch_variant_details(rs_id)

            if variant_data:
                writer.writerow(variant_data)
            time.sleep(1)  # Avoid API rate limits

    print(f"\n📂 Report generated: {report_filename}")

def main():
    #gene_id = input("Enter NCBI Gene ID: ").strip()
    gene_id = "3105"
    if not gene_id.isdigit():
        print("Invalid Gene ID. Please enter a valid numerical Gene ID.")
        return

    generate_report(gene_id)

if __name__ == "__main__":
    main()


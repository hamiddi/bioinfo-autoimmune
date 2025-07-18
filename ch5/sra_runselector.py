"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
This script retrieves metadata from the NCBI Sequence Read Archive (SRA) for a given BioProject ID
and saves the extracted study design information as a CSV file.

Usage:
1. Make sure you have the Biopython and pandas libraries installed.
2. Replace the Entrez.email value with your actual email address (required by NCBI).
3. Modify the last line of the script to provide your desired BioProject ID and output CSV filename.

Example:
    save_sra_study_design_csv("PRJEB52174", "PRJEB52174_study_design_metadata.csv")

Output:
    A CSV file containing SRA run accession, BioSample ID, scientific name, and sample attributes.

Note:
    This script uses the NCBI Entrez API and parses the returned XML to extract metadata.
    Be respectful of NCBI’s API usage limits by not overloading the server with rapid requests.
"""

from Bio import Entrez
import pandas as pd
import xml.etree.ElementTree as ET
import time

Entrez.email = "your_email@example.com"

def fetch_run_ids(bioproject_id):
    handle = Entrez.esearch(db="sra", term=bioproject_id, retmax=1000)
    record = Entrez.read(handle)
    return record['IdList']

def fetch_metadata_with_etree(sra_id):
    handle = Entrez.efetch(db="sra", id=sra_id, rettype="xml")
    xml_data = handle.read()
    handle.close()

    root = ET.fromstring(xml_data)
    metadata_list = []

    for experiment_package in root.findall(".//EXPERIMENT_PACKAGE"):
        row = {}
        
        # Run accession
        run = experiment_package.find(".//RUN")
        if run is not None:
            row['Run Accession'] = run.attrib.get("accession", "N/A")

        # BioSample accession
        sample = experiment_package.find(".//SAMPLE")
        if sample is not None:
            row['BioSample'] = sample.attrib.get("accession", "N/A")

            sample_name = sample.find("SAMPLE_NAME/SCIENTIFIC_NAME")
            if sample_name is not None:
                row['Scientific Name'] = sample_name.text

            # Extract study design fields from SAMPLE_ATTRIBUTE
            for attr in sample.findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
                tag = attr.find("TAG").text if attr.find("TAG") is not None else None
                value = attr.find("VALUE").text if attr.find("VALUE") is not None else None
                if tag and value:
                    row[tag] = value

        metadata_list.append(row)

    return metadata_list

def save_sra_study_design_csv(bioproject_id, output_csv):
    run_ids = fetch_run_ids(bioproject_id)
    print(f"Found {len(run_ids)} SRA records for {bioproject_id}")

    all_data = []
    for i, run_id in enumerate(run_ids):
        print(f"Fetching metadata {i+1}/{len(run_ids)}: {run_id}")
        try:
            records = fetch_metadata_with_etree(run_id)
            all_data.extend(records)
            time.sleep(0.5)
        except Exception as e:
            print(f"Error fetching run {run_id}: {e}")

    if all_data:
        df = pd.DataFrame(all_data)
        df.to_csv(output_csv, index=False)
        print(f"Saved to {output_csv}")
    else:
        print("No metadata found.")

# Run it
save_sra_study_design_csv("PRJEB52174", "PRJEB52174_study_design_metadata.csv")


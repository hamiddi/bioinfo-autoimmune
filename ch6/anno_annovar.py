"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script automates the annotation of VCF files using ANNOVAR. It performs the following steps:
1. Extracts `.vcf` files from compressed `.vcf.gz` archives.
2. Converts VCF files into ANNOVAR input format using `convert2annovar.pl`.
3. Annotates variants using multiple ANNOVAR databases, including refGene, ClinVar, dbNSFP, COSMIC, and dbSNP.
4. Saves the annotated output as a CSV file.
5. Reformats the annotated CSV file to ensure a consistent number of columns per row for downstream use.

The script is designed to handle batch processing of VCF files in a specified directory and automatically 
creates the required output structure.

Required Python Packages:
- os
- subprocess
- gzip
- shutil
- pandas
- csv

External Tools Required:
- ANNOVAR (convert2annovar.pl, table_annovar.pl)
- Perl
- ANNOVAR human database files (e.g., hg38 build)

Input Files:
- Compressed VCF files (`*.vcf.gz`) in the specified input directory (e.g., `test/`)
- ANNOVAR scripts and databases (e.g., `convert2annovar.pl`, `table_annovar.pl`, and `humandb38/`)

Output Files:
- Extracted VCF files (`*.vcf`) saved to the output directory (e.g., `test/annovar_anno/`)
- ANNOVAR input files (`*.avinput`)
- Annotated CSV files (e.g., `sample_annovar.hg38_multianno.csv`)
- Reformatted CSV files with fixed column count (e.g., `sample_annovar.hg38_multianno_fixed.csv`)
"""

import os
import subprocess
import gzip
import shutil
import pandas as pd
import csv

def extract_vcf(file_path, output_dir):
    """Extracts a .vcf file from .vcf.gz and saves it in the output directory."""
    file_name = os.path.basename(file_path).replace('.vcf.gz', '.vcf')
    extracted_file_path = os.path.join(output_dir, file_name)
    with gzip.open(file_path, 'rt') as f_in, open(extracted_file_path, 'w') as f_out:
        shutil.copyfileobj(f_in, f_out)
    return extracted_file_path

def convert_vcf_to_annovar_input(vcf_file, annovar_path, output_dir):
    """Converts a VCF file to ANNOVAR input format."""
    avinput_file = vcf_file.replace('.vcf', '.avinput')
    avinput_file_path = os.path.join(output_dir, os.path.basename(avinput_file))
    cmd = [
        "perl", os.path.join(annovar_path, "convert2annovar.pl"),
        "-format", "vcf4", vcf_file, ">", avinput_file_path
    ]
    subprocess.run(" ".join(cmd), shell=True)
    return avinput_file_path

def annotate_with_annovar(avinput_file, annovar_path, humandb_path, output_dir):
    """Annotates an ANNOVAR input file using multiple databases and saves output in CSV format."""
    output_prefix = os.path.join(output_dir, os.path.basename(avinput_file).replace('.avinput', '_annovar'))
    
    # List of ANNOVAR databases to use for annotation
    databases = ["refGene","avsnp151","clinvar_20240611","dbnsfp47c","cosmic70"]
    
    # Create the ANNOVAR annotation command
    cmd = [
        "perl", os.path.join(annovar_path, "table_annovar.pl"),
        avinput_file, humandb_path,
        "-buildver", "hg38",
        "-out", output_prefix,
        "-remove", "-protocol", ",".join(databases),
        "-operation", "g,f,f,f,f", "-nastring", ".", "-csvout"
    ]
    subprocess.run(cmd)
    #final_variants_annovar.hg38_multianno.csv
    #final_variants_annovar
    csv_file = output_prefix + ".hg38_multianno.csv"
    #print(f"Annotation complete for {avinput_file}. Output: {output_prefix}.csv")
    print(f"Annotation complete for {avinput_file}. Output: {csv_file}")


def reformat_csv(input_file, output_file, num_columns=18):
    """Reformats a CSV file to ensure each line has a specified number of columns."""
    with open(input_file, 'r', newline='', encoding='utf-8') as infile:
        reader = csv.reader(infile)
        rows = list(reader)

    with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.writer(outfile)

        for row in rows:
            # Ensure the row has the specified number of columns
            if len(row) < num_columns:
                # Add empty strings to make up the required number of columns
                row.extend([''] * (num_columns - len(row)))
            elif len(row) > num_columns:
                # Trim the row if it has more than the required number of columns
                row = row[:num_columns]

            # Write the adjusted row to the output file
            writer.writerow(row)

    print(f"Reformatted CSV file saved to {output_file}")

def main():
    # Specify paths
    input_dir = "test"  # Directory containing VCF files
    output_dir = "test/annovar_anno"  # Directory to save annotated files
    annovar_path = "RIGHT_PATH/annovar"  # Path to ANNOVAR
    humandb_path = f"{annovar_path}/dbs/humandb38"  # Path to ANNOVAR human database

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Iterate through all .vcf.gz files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".vcf.gz"):
            file_path = os.path.join(input_dir, file_name)

            print(f"Processing {file_name}...")

            # Extract the VCF file to the output directory
            vcf_file = extract_vcf(file_path, output_dir)
            # Convert VCF to ANNOVAR input format
            avinput_file = convert_vcf_to_annovar_input(vcf_file, annovar_path, output_dir)
            # Annotate the ANNOVAR input file and save output in CSV format
            annotate_with_annovar(avinput_file, annovar_path, humandb_path, output_dir)

            # Get base name for output
            output_prefix = os.path.join(output_dir, os.path.basename(avinput_file).replace('.avinput', '_annovar'))
            csv_file = output_prefix + ".hg38_multianno.csv"
            reformatted_csv = output_prefix + ".hg38_multianno_fixed.csv"

            # Check if the CSV was created before trying to reformat
            if os.path.exists(csv_file):
                reformat_csv(csv_file, reformatted_csv)
            else:
                print(f"Warning: {csv_file} not found. Skipping reformatting.")
if __name__ == "__main__":
    main()


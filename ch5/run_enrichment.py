"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################
Functional Enrichment Analysis using g:Profiler API

This script performs functional enrichment analysis on a list of gene symbols 
using the g:Profiler API. It queries multiple biological databases including 
Gene Ontology (BP, MF, CC), KEGG, and Reactome.

Usage:
- Modify the 'input_gene_file' variable or uncomment the `sys.argv[1]` line 
  to specify a CSV or TXT file containing gene symbols (one per line or column).
- The script reads the first column of the input file and submits the gene list 
  to g:Profiler.
- The results, including enriched terms and statistics, will be saved to 
  'enrichment_results.csv' by default.

Requirements:
- Python 3.x
- Packages: requests, pandas

Note:
- You can customize the organism (default is 'hsapiens') or the output filename 
  by editing the call to `functional_enrichment()`.

Example:
    python enrichment_script.py genes.csv
"""

import requests
import pandas as pd
import sys
import json

def functional_enrichment(gene_list, organism="hsapiens", output_file="enrichment_results.csv"):
    if not gene_list:
        print("[INFO] No significant genes provided. Skipping enrichment.")
        return

    print(f"[INFO] Querying g:Profiler for {len(gene_list)} genes...")
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    headers = {'Content-Type': 'application/json'}

    payload = {
        "organism": organism,
        "query": gene_list,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"]
    }

    try:
        response = requests.post(url, headers=headers, data=json.dumps(payload))
        response.raise_for_status()
        results = response.json()['result']
        if not results:
            print("[INFO] No enriched terms found.")
            return

        df = pd.DataFrame(results)
        df = df[['source', 'name', 'p_value', 'term_size', 'intersection_size']]
        df.to_csv(output_file, index=False)
        print(f"[INFO] Enrichment results saved to {output_file}")

    except requests.RequestException as e:
        print(f"[ERROR] Failed to connect to g:Profiler: {e}")
    except KeyError:
        print("[ERROR] Invalid response format from g:Profiler.")

if __name__ == "__main__":
    #input_gene_file = sys.argv[1]  # CSV or TXT with one gene symbol per line or column
    #Replace RIGHT_PATH 
    input_gene_file ="RIGHT_PATH/one_w_results/DEGs/top100_genes.csv"
    df = pd.read_csv(input_gene_file)
    gene_list = df.iloc[:, 0].dropna().tolist()
    functional_enrichment(gene_list)


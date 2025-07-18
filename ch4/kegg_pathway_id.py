"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script fetches and prints KEGG pathway information using the bioservices library.
How to use:
- Make sure you have the 'bioservices' package installed. You can install it with:
    pip install bioservices
- Set the `pathway_id` variable to the KEGG pathway ID you want to retrieve.
  For example, "hsa05323" corresponds to the rheumatoid arthritis pathway in humans.
- Run the script. It will fetch and display the full KEGG pathway entry from the KEGG database.

Example:
    pathway_id = "hsa05323"
"""
from bioservices import KEGG
def fetch_pathway_by_id(pathway_id):
  kegg = KEGG()
  pathway_info = kegg.get(pathway_id)
  print(pathway_info)
pathway_id = "hsa05323" #Example pathway for rheumatoid arthritis
fetch_pathway_by_id(pathway_id)

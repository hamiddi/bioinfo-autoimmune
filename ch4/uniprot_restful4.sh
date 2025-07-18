#!/bin/bash
# Book title: Bioinformatics of Autoimmune Diseases
# Author: Hamid D. Ismail, Ph.D.
URL="https://rest.uniprot.org/uniprotkb/search"
PARAM1="?query=autoimmune+disease+AND+reviewed:true"
PARAM2="&fields=accession,protein_name,gene_names,organism_name,keyword"
PARAM3="&format=json"
curl -X GET \
     "${URL}${PARAM1}${PARAM2}${PARAM3}" \
        -H "Accept: application/json" \
        > uniprot_prot_info_reviewed.json


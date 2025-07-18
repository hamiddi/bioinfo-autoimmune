#!/bin/bash
# Book title: Bioinformatics of Autoimmune Diseases
# Author: Hamid D. Ismail, Ph.D.
URL="https://rest.uniprot.org/uniprotkb/search"
PARAM="?query=autoimmune+disease+AND+reviewed:true&format=json"
curl -X GET \
     "${URL}${PARAM}" \
     -H "Accept: application/json" \
     > uniprot_autoimmune_reviewed.json


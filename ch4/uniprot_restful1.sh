#!/bin/bash
# Book title: Bioinformatics of Autoimmune Diseases
# Author: Hamid D. Ismail, Ph.D.
UNI_URL="https://rest.uniprot.org/uniprotkb/search"
PARAMS="?query=autoimmune+disease&format=json"
curl -X GET \
     ${UNI_URL}${PARAMS} \
     -H "Accept: application/json" \
     > uniprot_autoimmune.json


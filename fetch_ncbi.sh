#!/usr/bin/env bash
# Télécharge une séquence depuis NCBI par accession, via l'API E-utilities.
# Usage : ./fetch_ncbi.sh NC_045512.2 [db] [format]
#   db     : nuccore (défaut) | protein
#   format : fasta (défaut) | gb
# Exemple : ./fetch_ncbi.sh NC_045512.2 nuccore fasta > covid.fasta
set -euo pipefail

acc="${1:?Usage: $0 ACCESSION [db] [format]}"
db="${2:-nuccore}"
fmt="${3:-fasta}"

base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
curl -s "${base}?db=${db}&id=${acc}&rettype=${fmt}&retmode=text"

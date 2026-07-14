#!/usr/bin/env bash
# Convertit un FASTQ (.gz géré) en FASTA sur stdout.
# Usage : ./fastq_to_fasta.sh reads.fastq.gz > reads.fasta
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 fichier.fastq[.gz] > sortie.fasta" >&2
    exit 1
fi

reader="cat"
[[ "$1" == *.gz ]] && reader="zcat"

$reader "$1" | awk 'NR % 4 == 1 {printf ">%s\n", substr($0, 2)} NR % 4 == 2 {print}'

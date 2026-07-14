#!/usr/bin/env bash
# Compte le nombre de reads dans un ou plusieurs FASTQ (.gz gérés).
# Usage : ./count_reads.sh reads1.fastq.gz reads2.fastq ...
set -euo pipefail

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 fichier1.fastq[.gz] [fichier2 ...]" >&2
    exit 1
fi

for f in "$@"; do
    if [[ "$f" == *.gz ]]; then
        lines=$(zcat "$f" | wc -l)
    else
        lines=$(wc -l < "$f")
    fi
    printf "%s\t%d reads\n" "$f" "$((lines / 4))"
done

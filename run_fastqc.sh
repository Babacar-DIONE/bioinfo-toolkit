#!/usr/bin/env bash
# Lance FastQC sur un ou plusieurs FASTQ et lance MultiQC si disponible.
# Usage :
#   run_fastqc.sh *.fastq.gz
#   run_fastqc.sh -o qc/ -t 8 *.fastq.gz
set -euo pipefail

OUTDIR="fastqc_results"
THREADS=4

usage() {
    cat >&2 <<EOF
Usage: $0 [-o DIR] [-t INT] FICHIER(S).fastq.gz

  -o DIR   dossier de sortie (défaut : fastqc_results)
  -t INT   threads (défaut : 4)
  -h       aide
EOF
    exit 1
}

while getopts "o:t:h" opt; do
    case "$opt" in
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

[ "$#" -eq 0 ] && { echo "ERREUR : aucun fichier fourni." >&2; usage; }
command -v fastqc &>/dev/null || { echo "ERREUR : fastqc introuvable." >&2; exit 1; }

mkdir -p "$OUTDIR"

echo "FastQC sur $# fichier(s) → $OUTDIR" >&2
fastqc --outdir "$OUTDIR" --threads "$THREADS" "$@"

if command -v multiqc &>/dev/null; then
    echo "MultiQC..." >&2
    multiqc "$OUTDIR" --outdir "$OUTDIR" --filename multiqc_report --quiet
    echo "Rapport MultiQC : $OUTDIR/multiqc_report.html" >&2
else
    echo "MultiQC non trouvé — rapports individuels dans $OUTDIR/" >&2
fi

#!/usr/bin/env bash
# Fusionne les FASTQ de plusieurs lanes Illumina, par échantillon et par sens
# de lecture. Convention attendue : SAMPLE_S#_L00#_R#_001.fastq.gz
# Les fichiers .gz étant concaténables tels quels, la fusion se fait SANS
# recompression (cat binaire) → très rapide.
#
# Usage :
#   merge_lanes.sh -i RUN_DIR -o merged/          # fusionne
#   merge_lanes.sh -i RUN_DIR -o merged/ -n       # dry-run (aperçu seulement)
set -euo pipefail

indir="."
outdir="merged"
dryrun=0

usage() {
    cat >&2 <<EOF
Usage: $0 -i DOSSIER_ENTREE -o DOSSIER_SORTIE [-n]
  -i  dossier contenant les FASTQ par lane (défaut: .)
  -o  dossier de sortie (défaut: merged)
  -n  dry-run : affiche ce qui serait fait, sans rien écrire
EOF
    exit 1
}

while getopts "i:o:nh" opt; do
    case "$opt" in
        i) indir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        n) dryrun=1 ;;
        *) usage ;;
    esac
done

[ -d "$indir" ] || { echo "Dossier introuvable : $indir" >&2; exit 1; }
shopt -s nullglob

# Clés de regroupement uniques "<stem>__<read>" à partir des fichiers présents
keys=$(
  for f in "$indir"/*.gz; do
    base=$(basename "$f")
    if echo "$base" | grep -qE '_L[0-9]{3}_(R[12]|I[12])_[0-9]{3}\.f(ast)?q\.gz$'; then
      stem=$(echo "$base" | sed -E 's/_L[0-9]{3}_(R[12]|I[12])_[0-9]{3}\.f(ast)?q\.gz$//')
      rd=$(echo "$base" | sed -E 's/.*_L[0-9]{3}_(R[12]|I[12])_[0-9]{3}\.f(ast)?q\.gz$/\1/')
      printf '%s__%s\n' "$stem" "$rd"
    fi
  done | sort -u
)

if [ -z "$keys" ]; then
    echo "Aucun fichier au format SAMPLE_L00#_R#_001.fastq.gz dans $indir" >&2
    exit 1
fi

[ "$dryrun" -eq 0 ] && mkdir -p "$outdir"

echo "$keys" | while IFS= read -r key; do
    stem="${key%__*}"
    rd="${key##*__}"
    out="$outdir/${stem}_${rd}.fastq.gz"
    # Globs triés → L001, L002, L003... dans l'ordre
    files=( "$indir/${stem}"_L[0-9][0-9][0-9]_"${rd}"_[0-9][0-9][0-9].f*q.gz )
    printf '>> %s  (%d lane[s])\n' "$out" "${#files[@]}"
    printf '     %s\n' "${files[@]}"
    if [ "$dryrun" -eq 0 ]; then
        cat "${files[@]}" > "$out"   # concaténation binaire des .gz
    fi
done

[ "$dryrun" -eq 1 ] && echo "(dry-run : rien écrit)" >&2
echo "Terminé." >&2

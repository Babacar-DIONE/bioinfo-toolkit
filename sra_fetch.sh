#!/usr/bin/env bash
# Télécharge des données SRA (prefetch + fasterq-dump) avec reprise sur erreur.
#
# Dépendances : sra-tools >= 2.10 (prefetch, fasterq-dump), pigz ou gzip (--compress)
#
# Usage :
#   sra_fetch.sh SRR123456
#   sra_fetch.sh SRR123456 ERR456789 -o fastq/ -t 8 --compress
#   sra_fetch.sh -f accessions.txt -o fastq/ --compress
#
# Comportement :
#   - prefetch télécharge le fichier .sra dans OUTDIR/.sra/
#   - fasterq-dump extrait les FASTQ avec --split-files (R1/R2 séparés)
#   - Si le FASTQ de destination existe déjà, l'accession est ignorée
#   - Les fichiers .sra sont supprimés après extraction (sauf --keep-sra)
#   - En cas d'échec, la tentative est relancée N fois (défaut : 3)
set -euo pipefail

# ── Valeurs par défaut ──────────────────────────────────────────────────────
OUTDIR="."
THREADS=4
COMPRESS=0
KEEP_SRA=0
RETRIES=3
ACCFILE=""
ACCESSIONS=()

# ── Aide ────────────────────────────────────────────────────────────────────
usage() {
    cat >&2 <<EOF
Usage: $0 [OPTIONS] [ACCESSION...]

Télécharge des données SRA et extrait les FASTQ.

Options:
  -f FILE       fichier d'accessions (un par ligne ; lignes vides et # ignorées)
  -o DIR        dossier de sortie   (défaut : .)
  -t INT        threads pour fasterq-dump (défaut : 4)
  --compress    compresse les FASTQ en .gz (pigz si disponible, sinon gzip)
  --keep-sra    conserve les fichiers .sra après extraction
  --retries N   tentatives en cas d'erreur réseau (défaut : 3)
  -h            affiche cette aide

Exemples :
  $0 SRR18171532
  $0 SRR123456 ERR456789 -o fastq/ -t 8 --compress
  $0 -f my_accessions.txt -o fastq/ --compress --retries 5
EOF
    exit 1
}

# ── Parsing des arguments ────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        -f)        ACCFILE="$2"; shift 2 ;;
        -o)        OUTDIR="$2";  shift 2 ;;
        -t)        THREADS="$2"; shift 2 ;;
        --compress)  COMPRESS=1; shift ;;
        --keep-sra)  KEEP_SRA=1; shift ;;
        --retries)   RETRIES="$2"; shift 2 ;;
        -h|--help)   usage ;;
        -*)          echo "Option inconnue : $1" >&2; usage ;;
        *)           ACCESSIONS+=("$1"); shift ;;
    esac
done

# Ajouter les accessions depuis le fichier
if [[ -n "$ACCFILE" ]]; then
    [[ -f "$ACCFILE" ]] || { echo "ERREUR : fichier introuvable : $ACCFILE" >&2; exit 1; }
    while IFS= read -r line; do
        line="${line%%#*}"   # retirer les commentaires en fin de ligne
        line="${line// /}"   # retirer les espaces
        [[ -n "$line" ]] && ACCESSIONS+=("$line")
    done < "$ACCFILE"
fi

[[ ${#ACCESSIONS[@]} -eq 0 ]] && { echo "ERREUR : aucun accession fourni." >&2; usage; }

# ── Vérification des dépendances ─────────────────────────────────────────────
for tool in prefetch fasterq-dump; do
    command -v "$tool" &>/dev/null || {
        echo "ERREUR : '$tool' introuvable." >&2
        echo "  → Installe sra-tools : https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit" >&2
        exit 1
    }
done

if [[ $COMPRESS -eq 1 ]]; then
    COMPRESSOR=$(command -v pigz 2>/dev/null || command -v gzip 2>/dev/null || true)
    [[ -n "$COMPRESSOR" ]] || { echo "ERREUR : --compress demandé mais ni pigz ni gzip trouvé." >&2; exit 1; }
fi

mkdir -p "$OUTDIR"
SRA_DIR="$OUTDIR/.sra"
mkdir -p "$SRA_DIR"

# ── Fonction de téléchargement ───────────────────────────────────────────────
fetch_one() {
    local acc="$1"

    # Déterminer si déjà téléchargé (au moins un FASTQ présent)
    local existing
    existing=$(find "$OUTDIR" -maxdepth 1 -name "${acc}*.fastq" -o -name "${acc}*.fastq.gz" 2>/dev/null | head -1)
    if [[ -n "$existing" ]]; then
        echo ">>> $acc  déjà présent, ignoré." >&2
        return 0
    fi

    local attempt=0
    local ok=0

    # prefetch avec retry
    while [[ $attempt -lt $RETRIES ]]; do
        attempt=$((attempt + 1))
        echo "" >&2
        echo "━━━ $acc  — prefetch (tentative $attempt/$RETRIES) ━━━" >&2
        if prefetch "$acc" --output-directory "$SRA_DIR" --progress; then
            ok=1
            break
        fi
        [[ $attempt -lt $RETRIES ]] && { echo "Échec, nouvelle tentative dans 10s..." >&2; sleep 10; }
    done

    if [[ $ok -eq 0 ]]; then
        echo "ERREUR : prefetch $acc a échoué après $RETRIES tentatives." >&2
        return 1
    fi

    # fasterq-dump
    echo "━━━ $acc  — fasterq-dump ━━━" >&2
    fasterq-dump "$SRA_DIR/$acc" \
        --outdir "$OUTDIR" \
        --threads "$THREADS" \
        --split-files \
        --skip-technical \
        --progress

    # Compression
    if [[ $COMPRESS -eq 1 ]]; then
        echo "━━━ $acc  — compression ━━━" >&2
        while IFS= read -r -d '' fq; do
            "$COMPRESSOR" --force "$fq"
            echo "  compressé : $fq.gz" >&2
        done < <(find "$OUTDIR" -maxdepth 1 -name "${acc}*.fastq" -print0 2>/dev/null)
    fi

    # Nettoyage .sra
    if [[ $KEEP_SRA -eq 0 ]]; then
        rm -rf "${SRA_DIR:?}/$acc"
    fi

    echo ">>> $acc  OK" >&2
}

# ── Boucle principale ────────────────────────────────────────────────────────
FAILED=()
for acc in "${ACCESSIONS[@]}"; do
    fetch_one "$acc" || FAILED+=("$acc")
done

# Nettoyer le dossier .sra s'il est vide
rmdir "$SRA_DIR" 2>/dev/null || true

echo "" >&2
echo "══════════════════════════════════════════" >&2
echo "Terminé : $((${#ACCESSIONS[@]} - ${#FAILED[@]})) succès, ${#FAILED[@]} échec(s)." >&2
if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo "Accessions en échec : ${FAILED[*]}" >&2
    exit 1
fi

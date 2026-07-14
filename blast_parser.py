#!/usr/bin/env python3
"""Parse et filtre la sortie tabulaire BLAST (-outfmt 6 ou 7).

Colonnes standard attendues (12) :
  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

Si une 13e colonne qlen est présente (via -outfmt '6 std qlen'), --min-qcov devient
disponible. Une 14e colonne slen est également tolérée.

Filtres :
  --max-evalue   e-value maximale            (défaut : 1e-3)
  --min-pident   %% identité minimum         (défaut : 0)
  --min-length   longueur d'alignement min   (défaut : 0)
  --min-qcov     %% couverture query min     (nécessite col 13 = qlen)

Sélection par query :
  --best-hit     un seul hit par query (e-value ↑ puis bitscore ↓)
  --top N        N meilleurs hits par query

Usage:
    blast_parser.py hits.txt --max-evalue 1e-5 --min-pident 90
    blast_parser.py hits.txt.gz --best-hit
    blast_parser.py hits.txt --top 3 --no-header
    cat hits.txt | blast_parser.py - --min-qcov 80
"""
import argparse
import gzip
import sys
from collections import defaultdict

COLS = ["qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


def smart_open(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def passes_filters(col, args):
    try:
        pident = float(col[2])
        length = int(col[3])
        evalue = float(col[10])
    except (ValueError, IndexError):
        return False

    if evalue > args.max_evalue:
        return False
    if pident < args.min_pident:
        return False
    if length < args.min_length:
        return False

    if args.min_qcov > 0:
        if len(col) < 13:
            sys.exit("ERREUR : --min-qcov nécessite une 13e colonne qlen.\n"
                     "Relance BLAST avec : -outfmt '6 std qlen'")
        try:
            qlen = int(col[12])
            qcov = (abs(int(col[7]) - int(col[6])) + 1) / qlen * 100
        except (ValueError, ZeroDivisionError):
            return False
        if qcov < args.min_qcov:
            return False

    return True


def sort_key(h):
    """Tri : e-value croissante, bitscore décroissant."""
    return (float(h[10]), -float(h[11]))


def main():
    ap = argparse.ArgumentParser(
        description="Parse et filtre la sortie BLAST -outfmt 6/7",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("blast", help="Fichier BLAST outfmt 6/7 (.gz accepté) ou '-' pour stdin")
    ap.add_argument("--max-evalue", type=float, default=1e-3, metavar="FLOAT",
                    help="E-value maximale (défaut : 1e-3)")
    ap.add_argument("--min-pident", type=float, default=0.0, metavar="FLOAT",
                    help="%% identité minimum (défaut : 0)")
    ap.add_argument("--min-length", type=int, default=0, metavar="INT",
                    help="Longueur d'alignement minimum (défaut : 0)")
    ap.add_argument("--min-qcov", type=float, default=0.0, metavar="FLOAT",
                    help="%% couverture query minimum (nécessite col 13 = qlen)")
    ap.add_argument("--best-hit", action="store_true",
                    help="Garde uniquement le meilleur hit par query")
    ap.add_argument("--top", type=int, default=0, metavar="N",
                    help="Garde les N meilleurs hits par query")
    ap.add_argument("--no-header", action="store_true",
                    help="Ne pas écrire la ligne d'en-tête")
    args = ap.parse_args()

    hits = []
    n_total = 0
    with smart_open(args.blast) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            col = line.rstrip("\n").split("\t")
            if len(col) < 12:
                continue
            n_total += 1
            if passes_filters(col, args):
                hits.append(col)

    if not hits:
        print(f"[blast_parser] 0 hit conservé sur {n_total} (tous filtrés).",
              file=sys.stderr)
        sys.exit(0)

    # Sélection par query
    if args.best_hit or args.top > 0:
        by_query = defaultdict(list)
        for h in hits:
            by_query[h[0]].append(h)

        hits = []
        for q_hits in by_query.values():
            q_hits.sort(key=sort_key)
            hits.extend(q_hits[:1] if args.best_hit else q_hits[:args.top])

    # En-tête dynamique selon le nombre de colonnes présentes
    ncols = len(hits[0])
    header = COLS[:min(ncols, 12)]
    if ncols >= 13:
        header.append("qlen")
    if ncols >= 14:
        header.append("slen")

    if not args.no_header:
        print("\t".join(header))

    for h in hits:
        print("\t".join(h))

    print(f"[blast_parser] {len(hits)} hit(s) conservé(s) sur {n_total}",
          file=sys.stderr)


if __name__ == "__main__":
    main()

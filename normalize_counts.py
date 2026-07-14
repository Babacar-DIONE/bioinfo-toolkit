#!/usr/bin/env python3
"""Normalise une matrice de counts bruts en TPM, RPKM/FPKM ou CPM.

Formats d'entrée acceptés (auto-détectés) :
  - featureCounts  (Geneid Chr Start End Strand Length sample1.bam ...)
  - Matrice TSV    (gene_id  sample1  sample2 ...)
  - HTSeq-count    (gene  count, une seule colonne, lignes __* ignorées)

Méthodes :
  TPM   — Transcripts Per Million      (recommandée, comparable inter-échantillons)
  RPKM  — Reads Per Kilobase per Million (alias FPKM pour données paired-end)
  CPM   — Counts Per Million            (sans longueur, rapide)

TPM et RPKM nécessitent les longueurs de gènes. Elles sont extraites
automatiquement de featureCounts (colonne Length) ou fournies via --lengths.

Usage:
    normalize_counts.py featurecounts.txt --method tpm
    normalize_counts.py counts.tsv --method tpm rpkm --lengths gene_lengths.tsv
    normalize_counts.py counts.tsv.gz --method cpm
    normalize_counts.py counts.tsv --method tpm rpkm cpm
    cat htseq.txt | normalize_counts.py - --method tpm --lengths lengths.tsv
"""
import argparse
import gzip
import os
import sys


def smart_open(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_lengths_file(path):
    """TSV gene_id<TAB>length_bp — commentaires # tolérés."""
    lengths = {}
    with smart_open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    lengths[parts[0]] = int(float(parts[1]))
                except ValueError:
                    pass
    return lengths


def _is_numeric(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def parse_matrix(path):
    """
    Auto-détecte le format et retourne :
      genes   : list[str]
      samples : list[str]
      matrix  : list[list[float]]  — indexé [gène][échantillon]
      lengths : dict[str, int] | None
    """
    genes, samples, matrix = [], [], []
    lengths = None
    is_fc = False
    header_parsed = False

    with smart_open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue

            parts = line.split("\t")
            first = parts[0].lstrip("#").strip()

            if not header_parsed:
                # featureCounts header : Geneid  Chr  Start  End  Strand  Length  ...
                if first == "Geneid" and len(parts) >= 7 and parts[1] == "Chr":
                    is_fc = True
                    lengths = {}
                    # Strip BAM paths (e.g. /path/to/sample1.bam → sample1)
                    samples = [
                        os.path.splitext(os.path.basename(s))[0] for s in parts[6:]
                    ]
                    header_parsed = True
                    continue

                # featureCounts comment line  (# Program:featureCounts ...)
                if line.startswith("# Program:") or line.startswith("##"):
                    continue

                # Regular header : non-numeric second column
                if len(parts) >= 2 and not _is_numeric(parts[1]):
                    samples = parts[1:]
                    header_parsed = True
                    continue

                # No header found — generate sample names and process line as data
                samples = [f"Sample_{i + 1}" for i in range(len(parts) - 1)]
                header_parsed = True
                # fall through to data processing below

            # HTSeq summary lines
            if first.startswith("__"):
                continue

            gene = first

            if is_fc:
                if len(parts) < 7:
                    continue
                try:
                    gene_len = int(parts[5])
                    counts = [float(x) for x in parts[6:]]
                except ValueError:
                    continue
                lengths[gene] = gene_len
            else:
                try:
                    counts = [float(x) for x in parts[1:]]
                except ValueError:
                    continue

            if len(counts) != len(samples):
                continue

            genes.append(gene)
            matrix.append(counts)

    return genes, samples, matrix, lengths


# ── Méthodes de normalisation ────────────────────────────────────────────────

def calc_tpm(matrix, lengths, genes):
    n = len(matrix)
    ns = len(matrix[0]) if n else 0
    missing = set()

    # Étape 1 : RPK par gène × échantillon
    rpk = []
    for i, row in enumerate(matrix):
        glen = lengths.get(genes[i], 0)
        if glen == 0:
            missing.add(genes[i])
            rpk.append([0.0] * ns)
        else:
            rpk.append([c / (glen / 1000.0) for c in row])

    # Étape 2 : mise à l'échelle par échantillon → somme = 1 M
    out = [[0.0] * ns for _ in range(n)]
    for j in range(ns):
        total = sum(rpk[i][j] for i in range(n))
        if total == 0:
            continue
        sf = 1e6 / total
        for i in range(n):
            out[i][j] = rpk[i][j] * sf

    return out, missing


def calc_rpkm(matrix, lengths, genes):
    n = len(matrix)
    ns = len(matrix[0]) if n else 0
    missing = set()

    totals = [sum(matrix[i][j] for i in range(n)) for j in range(ns)]

    out = []
    for i, row in enumerate(matrix):
        glen = lengths.get(genes[i], 0)
        if glen == 0:
            missing.add(genes[i])
            out.append([0.0] * ns)
            continue
        vals = []
        for j, c in enumerate(row):
            t = totals[j]
            vals.append(c / (glen / 1000.0) / (t / 1e6) if t else 0.0)
        out.append(vals)

    return out, missing


def calc_cpm(matrix):
    n = len(matrix)
    ns = len(matrix[0]) if n else 0
    totals = [sum(matrix[i][j] for i in range(n)) for j in range(ns)]
    out = []
    for row in matrix:
        out.append([
            row[j] / totals[j] * 1e6 if totals[j] else 0.0
            for j in range(ns)
        ])
    return out


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Normalise une matrice de counts (TPM / RPKM / CPM)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("counts",
                    help="Matrice de counts (.tsv/.gz) ou '-' pour stdin")
    ap.add_argument("--method", nargs="+",
                    choices=["tpm", "rpkm", "cpm"], default=["tpm"],
                    metavar="METHOD",
                    help="Méthode(s) : tpm rpkm cpm — plusieurs possibles (défaut : tpm)")
    ap.add_argument("--lengths", metavar="FILE",
                    help="TSV gene_id<TAB>longueur_bp "
                         "(requis pour tpm/rpkm si entrée n'est pas featureCounts)")
    ap.add_argument("--round", type=int, default=3, metavar="N",
                    help="Décimales dans la sortie (défaut : 3)")
    args = ap.parse_args()

    genes, samples, matrix, auto_lengths = parse_matrix(args.counts)
    if not genes:
        sys.exit("ERREUR : aucune ligne de données trouvée.")

    lengths = dict(auto_lengths) if auto_lengths else {}
    if args.lengths:
        lengths.update(read_lengths_file(args.lengths))

    need_len = any(m in args.method for m in ("tpm", "rpkm"))
    if need_len and not lengths:
        sys.exit(
            "ERREUR : tpm/rpkm nécessitent les longueurs de gènes.\n"
            "  → Utilise un fichier featureCounts (longueurs extraites automatiquement)\n"
            "  → Ou passe --lengths gene_lengths.tsv  (format : gene_id<TAB>longueur_bp)"
        )

    # Calcul
    METHOD_LABEL = {"tpm": "TPM", "rpkm": "RPKM", "cpm": "CPM"}
    results = {}
    all_missing = set()

    for m in args.method:
        label = METHOD_LABEL[m]
        if m == "tpm":
            mat, miss = calc_tpm(matrix, lengths, genes)
        elif m == "rpkm":
            mat, miss = calc_rpkm(matrix, lengths, genes)
        else:
            mat, miss = calc_cpm(matrix), set()
        results[label] = mat
        all_missing |= miss

    if all_missing:
        ex = next(iter(all_missing))
        print(f"[normalize_counts] AVERTISSEMENT : {len(all_missing)} gène(s) sans "
              f"longueur → valeur 0 (ex. : {ex})", file=sys.stderr)

    # En-tête
    labels = list(results)
    if len(labels) == 1:
        col_headers = samples
    else:
        # Groupé par méthode : sample1_TPM sample2_TPM … sample1_RPKM …
        col_headers = [f"{s}_{lbl}" for lbl in labels for s in samples]

    print("gene_id\t" + "\t".join(col_headers))

    fmt = f"{{:.{args.round}f}}"
    ns = len(samples)
    for i, gene in enumerate(genes):
        if len(labels) == 1:
            vals = [fmt.format(results[labels[0]][i][j]) for j in range(ns)]
        else:
            vals = [
                fmt.format(results[lbl][i][j])
                for lbl in labels
                for j in range(ns)
            ]
        print(gene + "\t" + "\t".join(vals))

    print(
        f"[normalize_counts] {len(genes)} gènes × {ns} échantillons "
        f"— {', '.join(labels)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

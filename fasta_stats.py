#!/usr/bin/env python3
"""Statistiques d'un fichier FASTA : nb de séquences, longueur totale,
min / max / moyenne / médiane, N50 et %GC. Gère les fichiers .gz.

Usage:
    fasta_stats.py genome.fasta
    zcat reads.fa.gz | fasta_stats.py -
"""
import argparse
import gzip
import statistics
import sys


def smart_open(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def iter_fasta(handle):
    """Génère des tuples (header, sequence) sans charger tout le fichier."""
    name, seq = None, []
    for line in handle:
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def n50(lengths):
    if not lengths:
        return 0
    lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    cumulative = 0
    for length in lengths:
        cumulative += length
        if cumulative >= half:
            return length
    return lengths[-1]


def main():
    ap = argparse.ArgumentParser(description="Statistiques d'un fichier FASTA")
    ap.add_argument("fasta", help="Fichier FASTA (.fa/.fasta/.gz) ou '-' pour stdin")
    args = ap.parse_args()

    lengths, gc = [], 0
    with smart_open(args.fasta) as fh:
        for _, seq in iter_fasta(fh):
            lengths.append(len(seq))
            gc += seq.upper().count("G") + seq.upper().count("C")

    if not lengths:
        sys.exit("Aucune séquence trouvée.")

    total = sum(lengths)
    print(f"Séquences\t{len(lengths)}")
    print(f"Longueur_totale\t{total}")
    print(f"Min\t{min(lengths)}")
    print(f"Max\t{max(lengths)}")
    print(f"Moyenne\t{total / len(lengths):.1f}")
    print(f"Médiane\t{statistics.median(lengths):.1f}")
    print(f"N50\t{n50(lengths)}")
    print(f"GC%\t{100 * gc / total:.2f}")


if __name__ == "__main__":
    main()

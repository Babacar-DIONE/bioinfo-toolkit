#!/usr/bin/env python3
"""Statistiques d'un fichier FASTQ : nb de reads, bases totales,
longueur min/max/moyenne, qualité moyenne (Phred+33). Gère les .gz.

Usage:
    fastq_stats.py reads.fastq.gz
    zcat reads.fq.gz | fastq_stats.py -
"""
import argparse
import gzip
import sys


def smart_open(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def main():
    ap = argparse.ArgumentParser(description="Statistiques d'un fichier FASTQ")
    ap.add_argument("fastq", help="Fichier FASTQ (.fq/.fastq/.gz) ou '-' pour stdin")
    ap.add_argument("--phred", type=int, default=33, help="Offset Phred (défaut 33)")
    args = ap.parse_args()

    reads = bases = qual_sum = 0
    min_len, max_len = None, 0

    with smart_open(args.fastq) as fh:
        for i, line in enumerate(fh):
            m = i % 4
            if m == 1:  # séquence
                length = len(line.rstrip())
                reads += 1
                bases += length
                max_len = max(max_len, length)
                min_len = length if min_len is None else min(min_len, length)
            elif m == 3:  # qualité
                for c in line.rstrip():
                    qual_sum += ord(c) - args.phred

    if reads == 0:
        sys.exit("Aucun read trouvé.")

    print(f"Reads\t{reads}")
    print(f"Bases_totales\t{bases}")
    print(f"Longueur_min\t{min_len}")
    print(f"Longueur_max\t{max_len}")
    print(f"Longueur_moyenne\t{bases / reads:.1f}")
    print(f"Qualité_moyenne\t{qual_sum / bases:.2f}")


if __name__ == "__main__":
    main()

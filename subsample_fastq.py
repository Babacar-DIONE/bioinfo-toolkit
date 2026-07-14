#!/usr/bin/env python3
"""Sous-échantillonne N reads d'un FASTQ (reservoir sampling, une passe,
mémoire maîtrisée). Reproductible via --seed. Gère les .gz.

Usage:
    subsample_fastq.py reads.fastq.gz -n 100000 --seed 42 > sub.fastq
"""
import argparse
import gzip
import random
import sys


def smart_open(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def main():
    ap = argparse.ArgumentParser(description="Sous-échantillonnage de reads FASTQ")
    ap.add_argument("fastq", help="Fichier FASTQ ou '-' pour stdin")
    ap.add_argument("-n", "--number", type=int, required=True, help="Nombre de reads à garder")
    ap.add_argument("--seed", type=int, default=None, help="Graine aléatoire (reproductibilité)")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    reservoir = []
    read, count = [], 0

    with smart_open(args.fastq) as fh:
        for i, line in enumerate(fh):
            read.append(line.rstrip("\n"))
            if i % 4 == 3:  # read complet (4 lignes)
                count += 1
                if len(reservoir) < args.number:
                    reservoir.append(read)
                else:
                    j = rng.randint(0, count - 1)
                    if j < args.number:
                        reservoir[j] = read
                read = []

    for r in reservoir:
        print("\n".join(r))
    print(f"[subsample_fastq] {len(reservoir)} reads gardés sur {count}", file=sys.stderr)


if __name__ == "__main__":
    main()

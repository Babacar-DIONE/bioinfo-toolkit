#!/usr/bin/env python3
"""Filtre un FASTA par longueur de séquence. Sortie sur stdout.

Usage:
    filter_fasta.py genome.fa --min 500 > filtre.fa
    filter_fasta.py genome.fa --min 200 --max 10000 > filtre.fa
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


def iter_fasta(handle):
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


def main():
    ap = argparse.ArgumentParser(description="Filtre un FASTA par longueur")
    ap.add_argument("fasta", help="Fichier FASTA ou '-' pour stdin")
    ap.add_argument("--min", type=int, default=0, help="Longueur minimale")
    ap.add_argument("--max", type=int, default=float("inf"), help="Longueur maximale")
    ap.add_argument("-w", "--wrap", type=int, default=0,
                    help="Retour à la ligne toutes les N bases (0 = pas de wrap)")
    args = ap.parse_args()

    kept = 0
    with smart_open(args.fasta) as fh:
        for name, seq in iter_fasta(fh):
            if args.min <= len(seq) <= args.max:
                kept += 1
                print(f">{name}")
                if args.wrap:
                    for j in range(0, len(seq), args.wrap):
                        print(seq[j:j + args.wrap])
                else:
                    print(seq)
    print(f"[filter_fasta] {kept} séquence(s) conservée(s)", file=sys.stderr)


if __name__ == "__main__":
    main()

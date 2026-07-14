#!/usr/bin/env python3
"""Complément inverse de chaque séquence d'un FASTA. Gère les .gz et
les bases ambiguës (IUPAC). Sortie sur stdout.

Usage:
    reverse_complement.py seqs.fa > seqs_rc.fa
"""
import argparse
import gzip
import sys

COMP = str.maketrans("ACGTUNRYSWKMBDHVacgtunryswkmbdhv",
                     "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb")


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


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def main():
    ap = argparse.ArgumentParser(description="Complément inverse d'un FASTA")
    ap.add_argument("fasta", help="Fichier FASTA ou '-' pour stdin")
    args = ap.parse_args()

    with smart_open(args.fasta) as fh:
        for name, seq in iter_fasta(fh):
            print(f">{name} revcomp\n{revcomp(seq)}")


if __name__ == "__main__":
    main()

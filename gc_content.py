#!/usr/bin/env python3
"""%GC par séquence d'un FASTA. Sortie tabulée : id  longueur  GC%.

Usage:
    gc_content.py genome.fa
    gc_content.py genome.fa | sort -k3 -n   # trier par GC
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
    ap = argparse.ArgumentParser(description="%GC par séquence")
    ap.add_argument("fasta", help="Fichier FASTA ou '-' pour stdin")
    args = ap.parse_args()

    print("id\tlongueur\tGC%")
    with smart_open(args.fasta) as fh:
        for name, seq in iter_fasta(fh):
            seq_id = name.split()[0]
            up = seq.upper()
            gc = up.count("G") + up.count("C")
            pct = 100 * gc / len(seq) if seq else 0
            print(f"{seq_id}\t{len(seq)}\t{pct:.2f}")


if __name__ == "__main__":
    main()

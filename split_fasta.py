#!/usr/bin/env python3
"""Découpe un multi-FASTA. Deux modes :
  --per-file N : N séquences par fichier
  --parts K    : K fichiers de tailles ~égales

Usage:
    split_fasta.py genome.fa --per-file 100 -o chunks/
    split_fasta.py genome.fa --parts 8 -o parts/
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
    ap = argparse.ArgumentParser(description="Découpe un multi-FASTA")
    ap.add_argument("fasta", help="Fichier FASTA ou '-' pour stdin")
    ap.add_argument("-o", "--outdir", default=".", help="Répertoire de sortie")
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--per-file", type=int, help="N séquences par fichier")
    grp.add_argument("--parts", type=int, help="K fichiers ~égaux")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    records = list(iter_fasta(smart_open(args.fasta)))
    total = len(records)
    if total == 0:
        sys.exit("Aucune séquence.")

    if args.per_file:
        chunk = args.per_file
    else:  # --parts : taille de chunk arrondie au supérieur
        chunk = -(-total // args.parts)

    n = 0
    for start in range(0, total, chunk):
        n += 1
        path = os.path.join(args.outdir, f"part_{n:04d}.fa")
        with open(path, "w") as out:
            for name, seq in records[start:start + chunk]:
                out.write(f">{name}\n{seq}\n")
    print(f"[split_fasta] {total} séquences -> {n} fichier(s) dans {args.outdir}",
          file=sys.stderr)


if __name__ == "__main__":
    main()

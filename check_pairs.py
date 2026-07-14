#!/usr/bin/env python3
"""Vérifie l'appariement de deux FASTQ (R1 / R2) : même nombre de reads
et IDs correspondants dans le même ordre. Traitement en streaming
(compatible avec des millions de reads). Gère les .gz.

Code de sortie : 0 si tout est apparié, 1 sinon → utilisable dans un pipeline.

Usage:
    check_pairs.py SampleA_R1.fastq.gz SampleA_R2.fastq.gz
"""
import argparse
import gzip
import sys
from itertools import zip_longest


def smart_open(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def iter_headers(handle):
    """Ne renvoie que la ligne d'en-tête de chaque read (1 ligne sur 4)."""
    for i, line in enumerate(handle):
        if i % 4 == 0:
            yield line.rstrip()


def core_id(header):
    """Identifiant commun aux deux mates, sans le marqueur R1/R2.
    Gère '@ID 1:N:0:...' (Illumina récent) et '@ID/1' (ancien)."""
    h = header[1:] if header.startswith("@") else header
    h = h.split()[0]                 # coupe au 1er espace (retire ' 1:N:0:..')
    if h.endswith("/1") or h.endswith("/2"):
        h = h[:-2]
    return h


def main():
    ap = argparse.ArgumentParser(description="Vérifie l'appariement R1/R2")
    ap.add_argument("r1", help="FASTQ R1 (.gz accepté)")
    ap.add_argument("r2", help="FASTQ R2 (.gz accepté)")
    args = ap.parse_args()

    n = 0
    with smart_open(args.r1) as f1, smart_open(args.r2) as f2:
        for h1, h2 in zip_longest(iter_headers(f1), iter_headers(f2)):
            if h1 is None or h2 is None:
                extra = "R2" if h1 is None else "R1"
                print(f"ERREUR : nombre de reads différent — {extra} a plus de reads "
                      f"(divergence après {n} read[s]).", file=sys.stderr)
                sys.exit(1)
            if core_id(h1) != core_id(h2):
                print(f"ERREUR : IDs non appariés au read #{n + 1}\n"
                      f"  R1: {core_id(h1)}\n  R2: {core_id(h2)}", file=sys.stderr)
                sys.exit(1)
            n += 1

    print(f"OK : {n} paires correctement appariées.")


if __name__ == "__main__":
    main()

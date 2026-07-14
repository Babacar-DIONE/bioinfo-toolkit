#!/usr/bin/env python3
"""Extrait d'un FASTA les séquences dont l'ID figure dans une liste
(un ID par ligne). Avec -v : garde celles qui N'y sont PAS.

Usage:
    extract_seqs.py genome.fa --ids liste.txt > sous_ensemble.fa
    extract_seqs.py genome.fa --ids retirer.txt -v > restant.fa
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
    ap = argparse.ArgumentParser(description="Extrait des séquences par ID")
    ap.add_argument("fasta", help="Fichier FASTA ou '-' pour stdin")
    ap.add_argument("--ids", required=True, help="Fichier des IDs (un par ligne)")
    ap.add_argument("-v", "--invert", action="store_true",
                    help="Garde les séquences absentes de la liste")
    args = ap.parse_args()

    with open(args.ids) as fh:
        wanted = {line.strip() for line in fh if line.strip()}

    kept = 0
    with smart_open(args.fasta) as fh:
        for name, seq in iter_fasta(fh):
            seq_id = name.split()[0]  # 1er mot du header
            present = seq_id in wanted
            if present != args.invert:
                kept += 1
                print(f">{name}\n{seq}")
    print(f"[extract_seqs] {kept} séquence(s) extraite(s)", file=sys.stderr)


if __name__ == "__main__":
    main()

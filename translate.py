#!/usr/bin/env python3
"""Traduit des séquences ADN en protéines (code génétique standard).
Cadre 1 par défaut, ou les 6 cadres avec --six.

Usage:
    translate.py cds.fa > prot.fa
    translate.py contigs.fa --six > prot_6frames.fa
"""
import argparse
import gzip
import sys

CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


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


def translate(seq):
    seq = seq.upper().replace("U", "T")
    prot = []
    for i in range(0, len(seq) - 2, 3):
        prot.append(CODONS.get(seq[i:i + 3], "X"))
    return "".join(prot)


def main():
    ap = argparse.ArgumentParser(description="Traduction ADN -> protéine")
    ap.add_argument("fasta", help="Fichier FASTA ou '-' pour stdin")
    ap.add_argument("--six", action="store_true", help="Traduire les 6 cadres de lecture")
    args = ap.parse_args()

    with smart_open(args.fasta) as fh:
        for name, seq in iter_fasta(fh):
            if not args.six:
                print(f">{name}\n{translate(seq)}")
            else:
                rc = seq.upper().translate(COMP)[::-1]
                for frame in range(3):
                    print(f">{name}_frame+{frame + 1}\n{translate(seq[frame:])}")
                for frame in range(3):
                    print(f">{name}_frame-{frame + 1}\n{translate(rc[frame:])}")


if __name__ == "__main__":
    main()

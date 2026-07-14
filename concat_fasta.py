#!/usr/bin/env python3
"""Concatène plusieurs FASTA en un seul multi-FASTA (stdout). Contrairement
à un simple `cat`, garantit la séparation correcte des enregistrements
(pas de header collé à la séquence précédente si un fichier ne finit pas
par un saut de ligne). Gère les .gz.

Options utiles :
  --prefix-file : préfixe chaque header par le nom du fichier (traçabilité)
  --dedup       : ignore les IDs déjà rencontrés (garde la 1re occurrence)

Usage:
    concat_fasta.py *.fasta > all.fasta
    concat_fasta.py echantillons/*.fa --prefix-file > all.fasta
    concat_fasta.py transcrits/*.fa --dedup > all.fasta
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


def stem(path):
    b = os.path.basename(path)
    for ext in (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna", ".gz"):
        if b.endswith(ext):
            return b[: -len(ext)]
    return b


def main():
    ap = argparse.ArgumentParser(description="Concatène des FASTA en multi-FASTA")
    ap.add_argument("fasta", nargs="+", help="Fichiers FASTA (.gz acceptés)")
    ap.add_argument("--prefix-file", action="store_true",
                    help="Préfixe chaque header par le nom du fichier : >fichier|header")
    ap.add_argument("--dedup", action="store_true",
                    help="Ignore les IDs déjà vus (garde la 1re occurrence)")
    ap.add_argument("--sep", default="|", help="Séparateur pour --prefix-file (défaut: |)")
    args = ap.parse_args()

    seen = set()
    kept = dup = 0
    for path in args.fasta:
        n = 0
        with smart_open(path) as fh:
            for name, seq in iter_fasta(fh):
                seq_id = name.split()[0]
                if args.dedup and seq_id in seen:
                    dup += 1
                    continue
                if args.dedup:
                    seen.add(seq_id)
                header = f"{stem(path)}{args.sep}{name}" if args.prefix_file else name
                sys.stdout.write(f">{header}\n{seq}\n")
                kept += 1
                n += 1
        print(f"[concat_fasta] {path}: {n} séquence(s)", file=sys.stderr)

    msg = (f"[concat_fasta] TOTAL: {kept} séquence(s) depuis "
           f"{len(args.fasta)} fichier(s)")
    if args.dedup:
        msg += f" ({dup} doublon[s] ignoré[s])"
    print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()

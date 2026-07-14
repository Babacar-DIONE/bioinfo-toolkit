#!/usr/bin/env python3
"""Reconstruit les séquences de transcrits épissés à partir d'un génome
FASTA et d'une annotation GTF/GFF. Pour chaque transcrit, les exons sont
concaténés dans l'ordre puis, sur le brin '-', l'ensemble est
complémenté-inversé. Gère les .gz.

Coordonnées GTF/GFF : 1-basées et inclusives (converties correctement).
Regroupement : GTF -> transcript_id ; GFF3 -> Parent de l'exon.

Usage:
    extract_transcripts.py genome.fa annotation.gtf > transcrits.fa
    extract_transcripts.py genome.fa.gz annot.gff3.gz --format gff > tr.fa
    extract_transcripts.py genome.fa annot.gtf --feature CDS > cds.fa
"""
import argparse
import gzip
import re
import sys

COMP = str.maketrans("ACGTUNRYSWKMBDHVacgtunryswkmbdhv",
                     "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb")


def smart_open(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def load_genome(path):
    """Charge le génome en mémoire : {chrom: séquence}."""
    genome, name, seq = {}, None, []
    with smart_open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    genome[name] = "".join(seq)
                name = line[1:].split()[0]   # 1er mot du header
                seq = []
            else:
                seq.append(line)
    if name is not None:
        genome[name] = "".join(seq)
    return genome


def detect_format(attr):
    if 'transcript_id "' in attr or re.search(r'\w+ "', attr):
        return "gtf"
    if "=" in attr:
        return "gff"
    return "gtf"


def get_id(attr, fmt, id_attr):
    if fmt == "gtf":
        key = id_attr or "transcript_id"
        m = re.search(rf'{key}\s+"([^"]+)"', attr)
        return m.group(1) if m else None
    else:  # gff
        key = id_attr or "Parent"
        m = re.search(rf'(?:^|;)\s*{key}=([^;]+)', attr)
        return m.group(1).split(",")[0] if m else None


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def main():
    ap = argparse.ArgumentParser(description="Extrait les transcrits épissés")
    ap.add_argument("genome", help="Génome FASTA (.gz accepté)")
    ap.add_argument("annotation", help="Annotation GTF/GFF (.gz accepté)")
    ap.add_argument("--format", choices=["auto", "gtf", "gff"], default="auto",
                    help="Format d'annotation (défaut: auto-détecté)")
    ap.add_argument("--feature", default="exon",
                    help="Type de feature à assembler (défaut: exon ; ex: CDS)")
    ap.add_argument("--id-attr", default=None,
                    help="Attribut d'identifiant (défaut: transcript_id GTF / Parent GFF)")
    args = ap.parse_args()

    # Collecte des exons par transcrit : tid -> {chrom, strand, exons:[(start,end)]}
    tx, fmt = {}, args.format
    with smart_open(args.annotation) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            col = line.rstrip("\n").split("\t")
            if len(col) < 9 or col[2] != args.feature:
                continue
            if fmt == "auto":
                fmt = detect_format(col[8])
            tid = get_id(col[8], fmt, args.id_attr)
            if tid is None:
                continue
            chrom, start, end, strand = col[0], int(col[3]), int(col[4]), col[6]
            rec = tx.setdefault(tid, {"chrom": chrom, "strand": strand, "exons": []})
            rec["exons"].append((start, end))

    if not tx:
        sys.exit(f"Aucune feature '{args.feature}' exploitable "
                 f"(format détecté : {fmt}).")

    genome = load_genome(args.genome)

    written, missing = 0, 0
    for tid, rec in tx.items():
        chrom = rec["chrom"]
        if chrom not in genome:
            missing += 1
            print(f"[extract_transcripts] chrom absent du génome : {chrom} "
                  f"(transcrit {tid})", file=sys.stderr)
            continue
        chrom_seq = genome[chrom]
        exons = sorted(rec["exons"])          # ordre génomique croissant
        seq = "".join(chrom_seq[s - 1:e] for s, e in exons)   # 1-basé -> 0-basé
        if rec["strand"] == "-":
            seq = revcomp(seq)
        print(f">{tid}\n{seq}")
        written += 1

    msg = f"[extract_transcripts] {written} transcrit(s) écrit(s) (format {fmt})"
    if missing:
        msg += f", {missing} ignoré(s) (chrom manquant)"
    print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()

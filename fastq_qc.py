#!/usr/bin/env python3
"""Contrôle qualité d'un fichier FASTQ brut — Python stdlib, aucune dépendance.

Métriques calculées :
  - Stats globales : reads, bases, longueurs
  - Qualité par position : moyenne + Q10/Q25/médiane/Q75/Q90 (Phred)
  - Contenu par position : %A %T %G %C %N
  - Distribution du %GC par read
  - Distribution de la qualité moyenne par read
  - Détection d'adaptateurs Illumina (TruSeq, Nextera, Small RNA)
  - Séquences sur-représentées (50 pb, premiers 200 k reads)
  - Résumé PASS / WARN / FAIL par module

Gère les fichiers .gz et lit depuis stdin avec -.

Usage:
    fastq_qc.py reads.fastq.gz
    fastq_qc.py reads.fastq.gz --max-reads 500000
    fastq_qc.py R1.fq.gz > R1_qc.txt
    zcat reads.fq.gz | fastq_qc.py -
"""
import argparse
import collections
import gzip
import sys

# ── Constantes ───────────────────────────────────────────────────────────────
ADAPTERS = {
    "Illumina TruSeq":      "AGATCGGAAGAGC",
    "Nextera Transposase":  "CTGTCTCTTATACACATCT",
    "Illumina Small RNA":   "TGGAATTCTCGG",
}

MAX_Q       = 42          # Phred 0..41
OVR_LEN     = 50          # longueur du kmer pour sur-représentation
OVR_READS   = 200_000     # nb max de reads analysés pour sur-représentation
BASE_IDX    = {"A": 0, "T": 1, "G": 2, "C": 3}  # N → 4 implicitement


# ── Helpers ──────────────────────────────────────────────────────────────────

def smart_open(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def pct_from_hist(hist, p):
    """Percentile p (0-100) depuis un histogramme entier."""
    total = sum(hist)
    if total == 0:
        return 0
    target = max(1, round(total * p / 100))
    cum = 0
    for i, c in enumerate(hist):
        cum += c
        if cum >= target:
            return i
    return len(hist) - 1


def status(fail_cond, warn_cond):
    if fail_cond:
        return "FAIL"
    if warn_cond:
        return "WARN"
    return "PASS"


def print_section(title):
    print(f"\n>>  {title}  {'─' * max(0, 50 - len(title))}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Contrôle qualité FASTQ brut",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("fastq",
                    help="Fichier FASTQ (.gz accepté) ou '-' pour stdin")
    ap.add_argument("--phred", type=int, default=33, metavar="INT",
                    help="Offset Phred (défaut : 33)")
    ap.add_argument("--max-reads", type=int, default=0, metavar="N",
                    help="Analyser seulement les N premiers reads (0 = tous)")
    args = ap.parse_args()

    phred     = args.phred
    max_reads = args.max_reads

    # ── Accumulateurs ─────────────────────────────────────────────────────
    pos_q = []          # pos_q[i] = histogramme Phred [0..MAX_Q] à la position i
    pos_b = []          # pos_b[i] = [A, T, G, C, N] à la position i
    gc_hist  = [0] * 101
    mq_hist  = [0] * MAX_Q
    len_hist = collections.Counter()
    seq_cnt  = collections.Counter()
    adap_hits = {k: 0 for k in ADAPTERS}
    total_reads = total_bases = 0

    # ── Lecture ───────────────────────────────────────────────────────────
    with smart_open(args.fastq) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline().rstrip("\n")
            fh.readline()                           # ligne '+'
            qraw = fh.readline().rstrip("\n")

            if not seq or not qraw:
                break
            if len(seq) != len(qraw):
                continue                            # read corrompu

            L = len(seq)
            total_reads += 1
            total_bases += L
            len_hist[L] += 1

            # Étendre les tableaux par position si nécessaire
            while len(pos_q) < L:
                pos_q.append([0] * MAX_Q)
                pos_b.append([0, 0, 0, 0, 0])

            # Parcours base par base
            gc = q_sum = 0
            seq_up = seq.upper()
            for i in range(L):
                base = seq_up[i]
                q    = max(0, min(MAX_Q - 1, ord(qraw[i]) - phred))
                pos_q[i][q] += 1
                pos_b[i][BASE_IDX.get(base, 4)] += 1
                q_sum += q
                if base in "GC":
                    gc += 1

            gc_hist[round(gc * 100 / L)] += 1
            mq_hist[round(q_sum / L)] += 1

            # Adaptateurs — cherche dans les 40 derniers pb
            tail = seq_up[-40:]
            for name, adap in ADAPTERS.items():
                if adap in tail:
                    adap_hits[name] += 1

            # Sur-représentation — 50 premiers pb, limité aux OVR_READS premiers reads
            if total_reads <= OVR_READS and L >= OVR_LEN:
                seq_cnt[seq_up[:OVR_LEN]] += 1

            if max_reads and total_reads >= max_reads:
                break

    if total_reads == 0:
        sys.exit("ERREUR : aucun read trouvé.")

    # ── Calculs de synthèse ───────────────────────────────────────────────
    filename   = args.fastq if args.fastq != "-" else "stdin"
    len_vals   = sorted(len_hist.keys())
    min_len    = len_vals[0]
    max_len    = len_vals[-1]
    mean_gc    = sum(gc * c for gc, c in enumerate(gc_hist)) / total_reads
    median_mq  = pct_from_hist(mq_hist, 50)

    # Flags PASS/WARN/FAIL calculés lors de l'affichage des sections
    flags = {}

    # ═══════════════════════════════════════════════════════════════════════
    #  RAPPORT
    # ═══════════════════════════════════════════════════════════════════════
    print("=" * 60)
    print("  RAPPORT QC FASTQ")
    print("=" * 60)

    # ── 1. Stats globales ─────────────────────────────────────────────────
    print_section("1. Stats globales")
    print(f"Fichier            {filename}")
    print(f"Reads analysés     {total_reads:,}")
    print(f"Bases totales      {total_bases:,}")
    print(f"Longueur min       {min_len}")
    print(f"Longueur max       {max_len}")
    print(f"Longueurs uniques  {len(len_hist)}"
          + ("  (longueurs variables)" if len(len_hist) > 1 else ""))
    print(f"GC%% moyen         {mean_gc:.1f}%%")
    print(f"Qualité médiane    Q{median_mq}")
    flags["Stats globales"] = "PASS"

    # ── 2. Qualité par position ───────────────────────────────────────────
    print_section("2. Qualité par position (Phred)")
    print(f"{'Pos':>5}  {'Moy':>5}  {'Q10':>4}  {'Q25':>4}  {'Méd':>4}  {'Q75':>4}  {'Q90':>4}")
    print("-" * 42)
    fail_q = warn_q = False
    for i, h in enumerate(pos_q):
        if sum(h) == 0:
            continue
        mean_q_pos = sum(q * c for q, c in enumerate(h)) / sum(h)
        q10 = pct_from_hist(h, 10)
        q25 = pct_from_hist(h, 25)
        med = pct_from_hist(h, 50)
        q75 = pct_from_hist(h, 75)
        q90 = pct_from_hist(h, 90)
        print(f"{i + 1:>5}  {mean_q_pos:>5.1f}  {q10:>4}  {q25:>4}  {med:>4}  {q75:>4}  {q90:>4}")
        if med < 20:
            fail_q = True
        elif med < 28:
            warn_q = True
    flags["Qualité par position"] = status(fail_q, warn_q)

    # ── 3. Contenu par position ───────────────────────────────────────────
    print_section("3. Contenu par position (%%)")
    print(f"{'Pos':>5}  {'%A':>6}  {'%T':>6}  {'%G':>6}  {'%C':>6}  {'%N':>6}")
    print("-" * 44)
    fail_n = warn_n = False
    for i, b in enumerate(pos_b):
        total = sum(b)
        if total == 0:
            continue
        pct = [round(100 * b[j] / total, 1) for j in range(5)]
        print(f"{i + 1:>5}  {pct[0]:>6.1f}  {pct[1]:>6.1f}  {pct[2]:>6.1f}  {pct[3]:>6.1f}  {pct[4]:>6.1f}")
        if pct[4] >= 20:
            fail_n = True
        elif pct[4] >= 5:
            warn_n = True
    flags["%N par position"] = status(fail_n, warn_n)

    # ── 4. Distribution GC% par read ─────────────────────────────────────
    print_section("4. Distribution GC%% par read")
    print(f"{'GC%%':>5}  {'Reads':>10}  {'%%':>6}")
    print("-" * 28)
    for gc, count in enumerate(gc_hist):
        if count:
            print(f"{gc:>5}  {count:>10,}  {count * 100 / total_reads:>6.2f}")
    # Distribution anormale : trop étalée (écart-type > 15)
    variance = sum((gc - mean_gc) ** 2 * c for gc, c in enumerate(gc_hist)) / total_reads
    gc_sd = variance ** 0.5
    flags["Distribution GC"] = status(gc_sd > 20, gc_sd > 15)
    print(f"  → GC%% moyen : {mean_gc:.1f}%%   écart-type : {gc_sd:.1f}%%")

    # ── 5. Distribution qualité par read ─────────────────────────────────
    print_section("5. Distribution qualité par read (Q moyen)")
    print(f"{'Q moy':>6}  {'Reads':>10}  {'%%':>6}")
    print("-" * 28)
    for q, count in enumerate(mq_hist):
        if count:
            print(f"{q:>6}  {count:>10,}  {count * 100 / total_reads:>6.2f}")
    pct_below_q20 = sum(mq_hist[:20]) * 100 / total_reads
    flags["Qualité par read"] = status(median_mq < 20, median_mq < 27)
    print(f"  → Qualité médiane : Q{median_mq}   reads < Q20 : {pct_below_q20:.1f}%%")

    # ── 6. Adaptateurs ───────────────────────────────────────────────────
    print_section("6. Détection d'adaptateurs")
    print(f"{'Adaptateur':<30}  {'Reads':>8}  {'%%':>6}")
    print("-" * 50)
    worst_adap = 0.0
    for name, hits in adap_hits.items():
        pct = hits * 100 / total_reads
        worst_adap = max(worst_adap, pct)
        marker = "  ← DÉTECTÉ" if pct > 1 else ""
        print(f"{name:<30}  {hits:>8,}  {pct:>6.2f}{marker}")
    flags["Adaptateurs"] = status(worst_adap > 5, worst_adap > 1)

    # ── 7. Séquences sur-représentées ────────────────────────────────────
    print_section("7. Séquences sur-représentées (50 pb)")
    ref = min(total_reads, OVR_READS)
    threshold = max(1, round(ref * 0.001))       # seuil 0.1%%
    top = [(s, c) for s, c in seq_cnt.most_common(20) if c >= threshold]
    worst_ovr = 0.0
    if top:
        print(f"{'Séquence (50 pb)':<52}  {'Count':>8}  {'%%':>6}")
        print("-" * 68)
        for s, c in top:
            pct = c * 100 / ref
            worst_ovr = max(worst_ovr, pct)
            print(f"{s:<52}  {c:>8,}  {pct:>6.2f}")
    else:
        print("  Aucune séquence sur-représentée (> 0.1%% des reads)")
    flags["Sur-représentation"] = status(worst_ovr > 1, worst_ovr > 0.1)

    # ── 8. Résumé PASS / WARN / FAIL ─────────────────────────────────────
    print("\n" + "=" * 60)
    print("  RÉSUMÉ")
    print("=" * 60)
    ICONS = {"PASS": "[PASS]", "WARN": "[WARN]", "FAIL": "[FAIL]"}
    for module, s in flags.items():
        print(f"  {ICONS[s]}  {module}")
    print("=" * 60)

    n_fail = sum(1 for s in flags.values() if s == "FAIL")
    n_warn = sum(1 for s in flags.values() if s == "WARN")
    if n_fail:
        print(f"\n  {n_fail} module(s) en FAIL — action recommandée avant alignement.")
    elif n_warn:
        print(f"\n  {n_warn} module(s) en WARN — à surveiller.")
    else:
        print("\n  Tous les modules : PASS.")


if __name__ == "__main__":
    main()

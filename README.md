# bioinfo-toolkit

Boîte à outils de scripts bioinformatiques du quotidien : manipulation de FASTA/FASTQ,
statistiques, filtrage, conversion, traduction. **Aucune dépendance** : Python 3
(bibliothèque standard) + Bash. Chaque script est autonome — tu peux en copier un seul.

Tous les scripts gèrent les fichiers compressés `.gz` et lisent depuis `stdin` avec `-`,
donc ils s'enchaînent avec des pipes Unix.

## Scripts Python

| Script | Rôle | Exemple |
|--------|------|---------|
| `fasta_stats.py` | Stats FASTA (N50, GC%, longueurs) | `fasta_stats.py genome.fa` |
| `fastq_stats.py` | Stats FASTQ (reads, qualité moyenne) | `fastq_stats.py reads.fq.gz` |
| `gc_content.py` | %GC par séquence (tabulé) | `gc_content.py genome.fa` |
| `filter_fasta.py` | Filtre par longueur | `filter_fasta.py in.fa --min 500` |
| `extract_seqs.py` | Extrait/exclut par liste d'IDs | `extract_seqs.py in.fa --ids ids.txt` |
| `reverse_complement.py` | Complément inverse (IUPAC) | `reverse_complement.py in.fa` |
| `translate.py` | Traduction ADN→protéine (6 cadres) | `translate.py cds.fa --six` |
| `split_fasta.py` | Découpe un multi-FASTA | `split_fasta.py in.fa --parts 8 -o out/` |
| `subsample_fastq.py` | Sous-échantillonne N reads | `subsample_fastq.py r.fq.gz -n 1e5 --seed 42` |
| `concat_fasta.py` | Fusionne plusieurs FASTA en multi-FASTA | `concat_fasta.py *.fa --prefix-file` |
| `check_pairs.py` | Vérifie l'appariement R1/R2 | `check_pairs.py A_R1.fq.gz A_R2.fq.gz` |
| `extract_transcripts.py` | Transcrits épissés depuis génome + GTF/GFF | `extract_transcripts.py genome.fa a.gtf` |
| `blast_parser.py` | Parse/filtre BLAST -outfmt 6/7 (e-value, pident, qcov, best-hit) | `blast_parser.py hits.txt --max-evalue 1e-5 --best-hit` |
| `normalize_counts.py` | Normalise counts bruts en TPM / RPKM / CPM (featureCounts + matrice TSV) | `normalize_counts.py featurecounts.txt --method tpm` |

## Scripts Bash

| Script | Rôle |
|--------|------|
| `count_reads.sh` | Compte les reads d'un/plusieurs FASTQ |
| `fastq_to_fasta.sh` | Conversion FASTQ → FASTA |
| `fetch_ncbi.sh` | Télécharge une séquence NCBI par accession |
| `merge_lanes.sh` | Fusionne les FASTQ des lanes Illumina par échantillon |
| `sra_fetch.sh` | Télécharge des SRA via prefetch + fasterq-dump (retry, compress) |
| `run_fastqc.sh` | Lance FastQC (+ MultiQC si dispo) sur un ou plusieurs FASTQ |

## Installation

```bash
git clone https://github.com/Babacar-DIONE/bioinfo-toolkit.git
cd bioinfo-toolkit
bash install.sh        # rend les scripts exécutables + les ajoute au PATH
```

Ou manuellement : `chmod +x *.py *.sh` puis appelle-les directement.

## Exemples de pipelines

```bash
# Stats des contigs > 1 kb
filter_fasta.py assembly.fa --min 1000 | fasta_stats.py -

# GC des séquences, triées
gc_content.py genome.fa | sort -k3 -n

# Extraire des gènes puis les traduire
extract_seqs.py genes.fa --ids liste.txt | translate.py -

# Transcriptomique : fusionner les lanes d'un run, vérifier l'appariement
merge_lanes.sh -i run_brut/ -o merged/
check_pairs.py merged/SampleA_S1_R1.fastq.gz merged/SampleA_S1_R2.fastq.gz

# Regrouper des transcriptomes assemblés en une référence unique
concat_fasta.py assemblies/*.fa --prefix-file --dedup > reference.fa

# Reconstruire les transcrits épissés depuis un génome annoté
extract_transcripts.py genome.fa annotation.gtf > transcrits.fa

# Télécharger et compresser plusieurs SRA
sra_fetch.sh -f accessions.txt -o fastq/ -t 8 --compress

# Parser un BLAST : meilleur hit par query, identité ≥ 90 %, e-value ≤ 1e-5
blast_parser.py blast_out.txt --max-evalue 1e-5 --min-pident 90 --best-hit

# Top 3 hits, avec couverture query ≥ 80 % (nécessite -outfmt '6 std qlen')
blast_parser.py blast_out.txt --top 3 --min-qcov 80

# FastQC + MultiQC sur tous les FASTQ d'un run
run_fastqc.sh -o qc/ -t 8 raw_data/*.fastq.gz

# Normalisation TPM depuis featureCounts (longueurs extraites automatiquement)
normalize_counts.py featurecounts.txt --method tpm > tpm.tsv

# TPM + CPM en une passe depuis une matrice simple + fichier de longueurs
normalize_counts.py counts.tsv --method tpm cpm --lengths gene_lengths.tsv
```

## Licence

MIT — voir `LICENSE`.

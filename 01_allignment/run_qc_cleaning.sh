#!/usr/bin/env bash
# ============================================================
# BIFS619 – QC + READ CLEANING ONLY (run in 01_allignment/)
# - Downloads 3 paired-end replicates (SRR9613403/04/05)
# - FastQC (raw) + MultiQC (pre)
# - fastp trimming (q20, length>=50) using 10 threads
# - FastQC (trimmed) + MultiQC (post)
# - Tables: raw counts + duplicates; raw vs cleaned counts
# ============================================================
set -euo pipefail

THREADS=10
SAMPLES=("SRR9613403" "SRR9613404" "SRR9613405")

# ---------- 1) DOWNLOAD READS ----------
for S in "${SAMPLES[@]}"; do
  case "$S" in
    SRR9613403) BASE="003";;
    SRR9613404) BASE="004";;
    SRR9613405) BASE="005";;
  esac
  [ -f ${S}_1.fastq.gz ] || wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/${BASE}/${S}/${S}_1.fastq.gz
  [ -f ${S}_2.fastq.gz ] || wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/${BASE}/${S}/${S}_2.fastq.gz
done

# ---------- 2) QC ON RAW READS ----------
# --extract makes a folder <sample>_fastqc/ with fastqc_data.txt (easy to parse later)
fastqc --extract -t "$THREADS" \
  SRR9613403_1.fastq.gz SRR9613403_2.fastq.gz \
  SRR9613404_1.fastq.gz SRR9613404_2.fastq.gz \
  SRR9613405_1.fastq.gz SRR9613405_2.fastq.gz

multiqc . -n multiqc_raw.html

# ---------- 3) READ CLEANING (fastp) ----------
for S in "${SAMPLES[@]}"; do
  fastp \
    -i ${S}_1.fastq.gz -I ${S}_2.fastq.gz \
    -o ${S}_trimmed_1.fastq.gz -O ${S}_trimmed_2.fastq.gz \
    -q 20 -u 30 -n 5 -l 50 \
    --detect_adapter_for_pe \
    -w "$THREADS" \
    --html ${S}_fastp.html --json ${S}_fastp.json
done

# ---------- 4) QC ON TRIMMED READS ----------
fastqc --extract -t "$THREADS" \
  SRR9613403_trimmed_1.fastq.gz SRR9613403_trimmed_2.fastq.gz \
  SRR9613404_trimmed_1.fastq.gz SRR9613404_trimmed_2.fastq.gz \
  SRR9613405_trimmed_1.fastq.gz SRR9613405_trimmed_2.fastq.gz

multiqc . -n multiqc_post_trim.html

# ---------- 5) TABLE 1: RAW READ COUNTS + DUPLICATE RATE ----------
# raw read pairs from fastp 'before_filtering.total_reads'
# duplicate % = average of FastQC "Total Deduplicated Percentage" from R1 and R2 (raw)
echo -e "sample\traw_read_pairs\tduplication_rate_percent" > qc_raw_counts_duplicates.tsv
for S in "${SAMPLES[@]}"; do
  RAW=$(grep -A2 '"before_filtering"' ${S}_fastp.json | grep -m1 '"total_reads"' | sed 's/[^0-9]//g')

  D1=$(awk -F'\t' '/^Total Deduplicated Percentage/{print $2}' ${S}_1_fastqc/fastqc_data.txt | tr -d '%')
  D2=$(awk -F'\t' '/^Total Deduplicated Percentage/{print $2}' ${S}_2_fastqc/fastqc_data.txt | tr -d '%')
  # duplication rate = 100 - deduplicated%; average R1/R2
  if [[ -n "$D1" && -n "$D2" ]]; then
    DUP=$(awk -v a="$D1" -v b="$D2" 'BEGIN{printf "%.3f", 100 - ((a+b)/2)}')
  elif [[ -n "$D1" ]]; then
    DUP=$(awk -v a="$D1" 'BEGIN{printf "%.3f", 100 - a}')
  elif [[ -n "$D2" ]]; then
    DUP=$(awk -v b="$D2" 'BEGIN{printf "%.3f", 100 - b}')
  else
    DUP="NA"
  fi

  echo -e "${S}\t${RAW}\t${DUP}" >> qc_raw_counts_duplicates.tsv
done

# ---------- 6) TABLE 2: RAW vs CLEANED READ COUNTS ----------
echo -e "sample\traw_read_pairs\tcleaned_read_pairs\tkept_percent" > cleaning_raw_vs_trimmed.tsv
for S in "${SAMPLES[@]}"; do
  RAW=$(grep -A2 '"before_filtering"' ${S}_fastp.json | grep -m1 '"total_reads"' | sed 's/[^0-9]//g')
  CLEAN=$(grep -A2 '"after_filtering"'  ${S}_fastp.json | grep -m1 '"total_reads"' | sed 's/[^0-9]//g')
  KEPT=$(awk -v r="$RAW" -v c="$CLEAN" 'BEGIN{if(r>0) printf "%.2f", (c/r)*100; else print "0.00"}')
  echo -e "${S}\t${RAW}\t${CLEAN}\t${KEPT}" >> cleaning_raw_vs_trimmed.tsv
done

echo "DONE. Key outputs:"
echo "  - multiqc_raw.html, multiqc_post_trim.html"
echo "  - *_fastp.html / *_fastp.json"
echo "  - *_fastqc.html and *_fastqc/ folders (extracted) + *_fastqc.zip"
echo "  - qc_raw_counts_duplicates.tsv"
echo "  - cleaning_raw_vs_trimmed.tsv"

#!/usr/bin/env bash
# ==============================================
# RNA-seq Alignment Script using HISAT2
# ==============================================
# Author: Mike A.
# Date: $(date)
# Description:
#   Builds HISAT2 index and aligns all trimmed paired-end reads.
#   Converts SAM → sorted BAM → indexed BAM, and summarizes stats.
# ==============================================

set -euo pipefail  # safer scripting

# === Paths ===
BASE_DIR="/home/StudentFirst/rnaseq_project/mike_code/01_alignment"
REF="/home/StudentFirst/rnaseq_project/mike_code/00_raw_data/fastq_data/reference/GCF_000005845.2_ASM584v2_genomic.fna"
TRIM_DIR="${BASE_DIR}/QC/trimmed"

# Organized subfolders
INDEX_DIR="${BASE_DIR}/alignment/hisat2_index"
BAM_DIR="${BASE_DIR}/alignment/bam"
LOG_DIR="${BASE_DIR}/alignment/logs"
TABLE_DIR="${BASE_DIR}/alignment/tables"
PLOT_DIR="${BASE_DIR}/alignment/plots"

# Create directories if not existing
mkdir -p "$INDEX_DIR" "$BAM_DIR" "$LOG_DIR" "$TABLE_DIR" "$PLOT_DIR"

# HISAT2 index prefix
INDEX_PREFIX="${INDEX_DIR}/Ecoli_index"

echo "=== [1/4] Building HISAT2 index ==="
hisat2-build "$REF" "$INDEX_PREFIX" 2>&1 | tee "${LOG_DIR}/hisat2_build.log"

echo "=== [2/4] Aligning reads ==="
for R1 in ${TRIM_DIR}/*_1P.fastq; do
  SAMPLE=$(basename "$R1" _1P.fastq)
  R2=${TRIM_DIR}/${SAMPLE}_2P.fastq

  echo "Aligning sample: $SAMPLE"

  hisat2 -p 6 \
    -x "$INDEX_PREFIX" \
    -1 "$R1" -2 "$R2" \
    -S ${BAM_DIR}/${SAMPLE}.sam \
    2> ${LOG_DIR}/${SAMPLE}_hisat2.log

  echo "Converting and sorting BAM for $SAMPLE ..."
  samtools view -bS ${BAM_DIR}/${SAMPLE}.sam > ${BAM_DIR}/${SAMPLE}.bam
  samtools sort ${BAM_DIR}/${SAMPLE}.bam -o ${BAM_DIR}/${SAMPLE}_sorted.bam
  samtools index ${BAM_DIR}/${SAMPLE}_sorted.bam
  rm ${BAM_DIR}/${SAMPLE}.sam ${BAM_DIR}/${SAMPLE}.bam

  echo "Collecting flagstats for $SAMPLE ..."
  samtools flagstat ${BAM_DIR}/${SAMPLE}_sorted.bam > ${TABLE_DIR}/${SAMPLE}_flagstat.txt
done

echo "=== [3/4] Generating mapping summary table ==="
echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapping_Percentage" > ${TABLE_DIR}/mapping_summary.tsv
for F in ${TABLE_DIR}/*_flagstat.txt; do
  SAMPLE=$(basename "$F" _flagstat.txt)
  TOTAL=$(grep "in total" $F | awk '{print $1}')
  MAPPED=$(grep "mapped (" $F | awk '{print $1}')
  PCT=$(grep "mapped (" $F | awk '{print $5}' | tr -d '()%')
  echo -e "${SAMPLE}\t${TOTAL}\t${MAPPED}\t${PCT}" >> ${TABLE_DIR}/mapping_summary.tsv
done

echo "=== [4/4] Alignment complete! ==="
echo "All outputs are saved under:"
echo "  - Index files:   ${INDEX_DIR}"
echo "  - BAM files:     ${BAM_DIR}"
echo "  - Log files:     ${LOG_DIR}"
echo "  - Tables:        ${TABLE_DIR}"
echo "  - Plots:         ${PLOT_DIR}"


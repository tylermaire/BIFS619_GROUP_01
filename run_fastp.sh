#!/usr/bin/env bash
set -euo pipefail

# Ensure outputs exist
mkdir -p trimmed reports qc_post

# Loop through samples listed one-per-line in samples.txt
# (SRR9613403, SRR9613404, SRR9613405)
while read -r sample; do
  [[ -z "${sample}" ]] && continue
  echo "==== Trimming ${sample} ===="

  # RAW INPUTS ARE PLAIN .fastq (NOT .gz)
  fastp \
    -i 00_rawdata/${sample}_1.fastq \
    -I 00_rawdata/${sample}_2.fastq \
    -o trimmed/${sample}_1.trim.fastq.gz \
    -O trimmed/${sample}_2.trim.fastq.gz \
    --detect_adapter_for_pe \
    --cut_front --cut_tail --cut_mean_quality 20 \
    --length_required 50 \
    --thread 4 \
    --html reports/${sample}.fastp.html \
    --json reports/${sample}.fastp.json

done < samples.txt

echo "==== Running FastQC on trimmed reads ===="
fastqc -t 4 trimmed/*.trim.fastq.gz -o qc_post

echo "==== Summarizing with MultiQC ===="
multiqc -o qc_post qc_post

echo "All trimming + post-trim QC completed."
echo "All trimming + post-trim QC completed

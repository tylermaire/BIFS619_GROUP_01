#!/bin/bash
# Script to run featureCounts for annotation and quantification

# Paths
GTF="/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/00_rawdata/fastq_data/reference/NZ_CP076404.1.gtf"
BAM_DIR="/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/01_allignment/alignment/bam"
OUT_DIR="/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/03_annotation/counts"
OUT_FILE="$OUT_DIR/raw_counts.txt"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Run featureCounts
featureCounts \
  -a "$GTF" \
  -o "$OUT_FILE" \
  -T 8 \
  -p \
  -B \
  -t CDS \
  -g gene_id \
  "$BAM_DIR"/*.bam


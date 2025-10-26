#!/bin/bash
# ====================================================================================
# Gene Expression Quantification with featureCounts
# ====================================================================================

# Paths
REFERENCE_GFF="/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/00_rawdata/fastq_data/reference/GCF_000005845.2.gff"
BAM_DIR="/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/01_allignment/alignment/bam"
OUT_DIR="/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/02_annotation/counts"
OUT_FILE="$OUT_DIR/raw_counts.txt"
THREADS=8

echo "=== Gene Expression Quantification with featureCounts ==="
echo "Reference GFF: $REFERENCE_GFF"
echo "BAM directory: $BAM_DIR"
echo "Output directory: $OUT_DIR"
echo ""

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Check if reference GFF exists and is not empty
if [ ! -s "$REFERENCE_GFF" ]; then
    echo "ERROR: Reference GFF file is missing or empty."
    echo "Please add a valid GFF file to: $REFERENCE_GFF"
    exit 1
fi

# Check if BAM files exist
if ! ls ${BAM_DIR}/*.bam 1> /dev/null 2>&1; then
    echo "ERROR: No BAM files found in: $BAM_DIR"
    exit 1
fi

# Count BAM files
bam_count=$(ls -1 ${BAM_DIR}/*.bam 2>/dev/null | wc -l)
echo "Found $bam_count BAM files to process"
echo ""

# Run featureCounts with GFF file
echo "Starting gene expression quantification with full genome annotation..."
echo "Using annotation file: $REFERENCE_GFF"
echo ""

featureCounts \
  -a "$REFERENCE_GFF" \
  -o "$OUT_FILE" \
  -t gene \
  -g ID \
  -p \
  -T "$THREADS" \
  ${BAM_DIR}/*.bam

# Check if counts were generated
if [ -s "$OUT_FILE" ]; then
    # Count how many genes were quantified
    gene_count=$(tail -n +3 "$OUT_FILE" | wc -l)
    echo ""
    echo "✓ FeatureCounts completed successfully."
    echo "  Quantified $gene_count genes across $bam_count samples."
    echo ""
    echo "Output files:"
    echo "  - Raw counts: $OUT_FILE"
    echo "  - Summary: ${OUT_FILE}.summary"
else
    echo ""
    echo "ERROR: FeatureCounts failed to generate count data."
    exit 1
fi

echo ""
echo "=== Quantification complete ==="

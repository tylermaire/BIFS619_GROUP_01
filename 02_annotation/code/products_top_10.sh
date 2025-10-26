#!/bin/bash
# ====================================================================================
# Gene Expression Quantification with featureCounts
# ====================================================================================

# Input/Output paths
BAM_DIR="../01_alignment/alignment/bam"
REFERENCE_GFF="../00_rawdata/fastq_data/reference/GCF_000005845.2.gff"
OUTPUT_DIR="../02_annotation/counts"
THREADS=8

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=== Starting gene expression quantification with featureCounts ==="
echo "BAM directory: $BAM_DIR"
echo "Reference GFF: $REFERENCE_GFF"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Check if GFF exists
if [ ! -f "$REFERENCE_GFF" ]; then
    echo "ERROR: Reference GFF file not found at: $REFERENCE_GFF"
    exit 1
fi

# Check if BAM files exist
if ! ls ${BAM_DIR}/*.bam 1> /dev/null 2>&1; then
    echo "ERROR: No BAM files found in: $BAM_DIR"
    exit 1
fi

# Count BAM files
bam_count=$(ls -1 ${BAM_DIR}/*.bam 2>/dev/null | wc -l)
echo "Found $bam_count BAM files"
echo ""

# Run featureCounts
echo "Running featureCounts..."
featureCounts \
  -a "$REFERENCE_GFF" \
  -o "${OUTPUT_DIR}/raw_counts.txt" \
  -t gene \
  -g ID \
  -p \
  -T "$THREADS" \
  ${BAM_DIR}/*.bam

# Check if counts were generated
if [ -s "${OUTPUT_DIR}/raw_counts.txt" ]; then
    gene_count=$(tail -n +3 "${OUTPUT_DIR}/raw_counts.txt" | wc -l)
    echo ""
    echo "✓ FeatureCounts completed successfully."
    echo "  Quantified $gene_count genes"
    echo "  Output: ${OUTPUT_DIR}/raw_counts.txt"
    echo "  Summary: ${OUTPUT_DIR}/raw_counts.txt.summary"
else
    echo ""
    echo "ERROR: FeatureCounts failed to generate count data."
    exit 1
fi

echo ""
echo "=== Gene quantification complete ==="

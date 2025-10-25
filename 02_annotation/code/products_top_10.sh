#!/bin/bash

# Paths
CSV="/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/03_annotation/counts/top_10_genes.csv"
GTF="/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/00_rawdata/fastq_data/reference/NZ_CP076404.1.gtf"
OUTPUT="/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/03_annotation/counts/top_10_genes_products.csv"

# Prepare output file with header
echo "Geneid,Product" > "$OUTPUT"

# Loop through each Geneid in the CSV (skip header)
tail -n +2 "$CSV" | cut -d',' -f1 | while read gene; do
    # Search GTF for gene_id and extract product field
    product=$(grep -F "$gene" "$GTF" | grep -m1 "CDS" | sed -n 's/.*product "\([^"]*\)".*/\1/p')
    echo "$gene,$product" >> "$OUTPUT"
done

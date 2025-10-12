#!/bin/bash
# Align cleaned reads for all samples using HISAT2

ref_dir="/home/tylermaire/Group_Project/BIFS619_GROUP_01/00_rawdata/fastq_data/reference"
ref_fa="$ref_dir/NZ_CP076404.1.fasta"
index_dir="/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/alignment/hisat2_index"
bam_dir="/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/alignment/bam"
log_dir="/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/alignment/logs"

mkdir -p "$index_dir" "$bam_dir" "$log_dir"

hisat2-build "$ref_fa" "$index_dir/genome"

for sample in SRR9613403 SRR9613404 SRR9613405; do
    hisat2 \
        -x "$index_dir/genome" \
        -1 "/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/QC/cleaned_fastq/${sample}_1.clean.fastq.gz" \
        -2 "/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/QC/cleaned_fastq/${sample}_2.clean.fastq.gz" \
        -S "$bam_dir/${sample}.sam" \
        --summary-file "$log_dir/${sample}_hisat2_summary.txt" \
        2> "$log_dir/${sample}_hisat2.log"

    samtools view -@ 8 -bS "$bam_dir/${sample}.sam" | \
        samtools sort -@ 8 -o "$bam_dir/${sample}.bam"
    samtools index "$bam_dir/${sample}.bam"
    rm "$bam_dir/${sample}.sam"
done
#!/bin/bash
# Read cleaning for all paired-end samples using fastp

input_dir="/home/tylermaire/Group_Project/BIFS619_GROUP_01/00_rawdata/fastq_data/samples"
output_dir="/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/QC/cleaned_fastq"
mkdir -p "$output_dir"

for sample in SRR9613403 SRR9613404 SRR9613405; do
    fastp \
        -i "$input_dir/${sample}_1.fastq.gz" \
        -I "$input_dir/${sample}_2.fastq.gz" \
        -o "$output_dir/${sample}_1.clean.fastq.gz" \
        -O "$output_dir/${sample}_2.clean.fastq.gz" \
        --html "$output_dir/${sample}_fastp.html" \
        --json "$output_dir/${sample}_fastp.json" \
        --thread 10
done
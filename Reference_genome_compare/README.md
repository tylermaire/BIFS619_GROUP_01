# Reference Genome Validation — Mike

## Overview
This analysis compares alignments and downstream results using the **incorrect** and **correct** reference genomes for the *E. coli* WO153 strain.  
It aims to demonstrate the impact of reference genome selection on alignment accuracy and gene expression metrics.

## Background
- Original group alignment used: **NZ_CP076404.1**
- Correct reference based on the source paper: **GCF_000005845.2**
- RNA-seq data source: WO153 *E. coli* strain (triplicates, paired-end)

## Methods
1. Quality Control: FastQC  
2. Trimming: Trimmomatic  
3. Alignment: Bowtie2
4. Deduplication and Indexing: Picard  
5. Comparison Metrics: alignment rate, mapped reads, and gene coverage.


## Key Findings
- The correct reference genome (GCF_000005845.2) yielded higher alignment precision.
- Mismatches and unmapped reads were significantly reduced.
- Downstream differential expression analysis results shifted accordingly.

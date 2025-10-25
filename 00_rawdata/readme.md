# BIFS619 Group 01 - Step 00: Data Setup and Raw QC

## Overview
This script is the first part of the BIFS619 RNA-Seq analysis workflow. Its purpose is to prepare the project environment by downloading all necessary raw data and performing an initial quality control assessment.

This script corresponds to the tasks that generate files for the `00_rawdata` directory in the master pipeline.

## What this script does:
1.  **Creates Project Directories**: Sets up the base project folder and the necessary subdirectories (`00_rawdata/fastq_data/samples`, `00_rawdata/fastq_data/reference`, `00_rawdata/fastQC`).
2.  **Downloads Raw Data**: Fetches paired-end FASTQ files for the specified samples from the European Nucleotide Archive (ENA).
3.  **Downloads Reference Genome**: Downloads the *E. coli* K-12 (ASM584v2) reference genome.
4.  **Performs Quality Control**: Runs `FastQC` on the raw FASTQ files to assess their quality.
5.  **Summarizes QC**: Runs `MultiQC` to create a single, comprehensive quality report from the FastQC outputs.

## Dataset Information
This pipeline analyzes RNA-Seq data for *Escherichia coli*. The reference genome is *E. coli* K-12 (ASM584v2), and the raw sequencing data is from the following study:

-   O'Rourke, A., et al. (2020). Mechanism-of-Action Classification of Antibiotics by Global Transcriptome Profiling. *Antimicrobial agents and chemotherapy*, 64(3). https://doi.org/10.1128/AAC.01207-19

The default samples configured in this script are:
-   `SRR9613403`
-   `SRR9613404`
-   `SRR9613405`

## How to Run

### Prerequisites
-   A Linux environment (like Ubuntu on WSL).
-   The following tools must be installed: `wget`, `gunzip`, `fastqc`, `multiqc`. The script will check for these before running.

### Instructions

1.  **Navigate to the directory containing the script.**

2.  **Make the script executable:**
    ```bash
    chmod +x 00_setup_and_qc.sh
    ```

3.  **Run the script with a project directory path:**

    The script requires one argument: a path for a **new project directory** that it will create for all output files.

    ```bash
    # Example: Create a directory named 'rnaseq_project' in your home folder
    ./00_setup_and_qc.sh ~/rnaseq_project
    ```
    This command will create the `rnaseq_project` folder and place all downloaded data and QC reports inside it.

## Output Structure
After running, the script will produce the following structure within the directory you specified:

```
<your_project_directory>/
└── 00_rawdata/
    ├── fastQC/
    │   ├── SRR9613403_1_fastqc.html
    │   ├── SRR9613403_1_fastqc.zip
    │   ├── SRR9613403_2_fastqc.html
    │   ├── SRR9613403_2_fastqc.zip
    │   ├── ... (reports for other samples)
    │   └── multiqc_report.html  <-- Main QC Summary
    ├── fastq_data/
    │   ├── reference/
    │   │   └── GCF_000005845.2.fna
    │   └── samples/
    │       ├── SRR9613403_1.fastq.gz
    │       ├── SRR9613403_2.fastq.gz
    │       ├── SRR9613404_1.fastq.gz
    │       ├── SRR9613404_2.fastq.gz
    │       ├── SRR9613405_1.fastq.gz
    │       └── SRR9613405_2.fastq.gz
    └── ...
```

## Contributors
-   BIFS619 Group 01 Team Members
-   Uploaded and README by Tyler Maire

# BIFS619 Group 01 - Step 01: Read Cleaning, Alignment, and QC

## Overview
This script is the second part of the BIFS619 RNA-Seq analysis workflow. It takes the raw data prepared by the `00_setup_and_qc.sh` script and performs read cleaning, alignment to the reference genome, and generation of summary statistics for both processes.

This script corresponds to the tasks that generate files for the `01_allignment` directory in the master pipeline.

## What this script does:
1.  **Cleans Raw Reads**: Uses `fastp` to remove low-quality bases and adapter sequences from the raw FASTQ files.
2.  **Summarizes Cleaning Results**: Runs a Python script to parse the `fastp` output and generate a summary table and a bar plot showing read counts before and after cleaning.
3.  **Builds Genome Index**: Creates a `HISAT2` index from the reference genome if one does not already exist.
4.  **Aligns Reads**: Aligns the cleaned reads to the *E. coli* reference genome using `HISAT2`.
5.  **Processes Alignments**: Converts the SAM alignment files to sorted, indexed BAM files using `samtools`.
6.  **Summarizes Alignment**: Runs a Python script to parse the `HISAT2` logs and generate a summary table and a bar plot showing the mapping rate for each sample.

## How to Run

### Prerequisites
-   A Linux environment (like Ubuntu on WSL).
-   The output from the `00_setup_and_qc.sh` script must exist in the specified project directory. This includes the raw FASTQ files and the reference genome.
-   The following tools must be installed: `fastp`, `hisat2`, `samtools`, and `python3` (with `pandas` and `matplotlib` libraries).

### Instructions

1.  **Navigate to the directory containing the script.**

2.  **Make the script executable:**
    ```bash
    chmod +x 01_alignment_and_qc.sh
    ```

3.  **Run the script with the *same* project directory path you used for the previous step:**

    ```bash
    # Example: Use the 'rnaseq_project' folder created in the first step
    ./01_alignment_and_qc.sh ~/rnaseq_project
    ```
    The script will find the raw data in `~/rnaseq_project/00_rawdata` and place its output in `~/rnaseq_project/01_allignment`.

## Expected Results

The script generates summary tables and plots to assess the quality of the cleaning and alignment steps.

### Read Cleaning Summary

The following table shows the number of reads before and after cleaning with `fastp`. A small percentage of removed reads (1-3%) is typical and indicates successful removal of low-quality data.

| Sample     | Raw Reads | Cleaned Reads | Removed Reads | Removed % |
| :--------- | :-------- | :------------ | :------------ | :-------- |
| SRR9613403 | 28299586  | 27784644      | 514942        | 1.82      |
| SRR9613404 | 22861028  | 22335460      | 525568        | 2.30      |
| SRR9613405 | 20736074  | 20294582      | 441492        | 2.13      |

This data is visualized in the following chart, which shows the total number of raw reads and the proportion of cleaned reads that are carried forward for alignment.

<img width="3600" height="1800" alt="pre_post_cleaning_barplot" src="https://github.com/user-attachments/assets/5f2e838d-3508-4e87-9f45-fcbb8ad23609" />


### Alignment Summary

After alignment, the mapping rate indicates how well the cleaned reads matched the reference genome. High mapping rates (>80%) are a good indicator of a successful alignment.

<img width="3000" height="1800" alt="mapping_percent_barplot" src="https://github.com/user-attachments/assets/3b341e3e-0f7c-463a-9335-4470f3efef09" />


## Output Structure
After running, the script will populate the `01_allignment` directory:

```
<your_project_directory>/
└── 01_allignment/
    ├── QC/
    │   ├── cleaned_fastq/  # Cleaned FASTQ files from fastp
    │   ├── plots/
    │   │   └── pre_post_cleaning_barplot.png
    │   └── tables/
    │       └── pre_post_cleaning_table.csv
    ├── alignment/
    │   ├── bam/            # Final sorted and indexed BAM files
    │   ├── hisat2_index/   # Genome index files
    │   ├── logs/           # HISAT2 log files
    │   ├── plots/
    │   │   └── mapping_percent_barplot.png
    │   └── tables/
    │       └── alignment_summary.csv
    └── code/
        ├── summarize_qc.py
        └── summarize_alignment.py
```

## Contributors
-   BIFS619 Group 01 Team Members
-   Tyler Maire (tylermaire, last updated 2025-10-25)


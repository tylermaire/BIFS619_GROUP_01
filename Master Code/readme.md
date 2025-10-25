# BIFS619 Group 01 Project

## Overview
This repository contains the group project for BIFS619, primarily focused on RNA-seq analysis with directories for raw data, alignment, and annotation workflows.

## Structure
The repository is organized as follows:
- `00_rawdata`: Contains the initial datasets used for analysis
- `01_allignment`: Contains alignment scripts and results
- `02_annotation`: Contains annotation files and processes
- `Master Code`: Contains the main pipeline script for the project workflow

## How to Run the Master Code

### Prerequisites
- Git installed on your system
- Ubuntu on Windows (WSL) or a Linux system
- Basic knowledge of bioinformatics principles

### Instructions

1. For Windows users, open Ubuntu:
   - Search for "Ubuntu" in the Start menu and open it
   - Or open Command Prompt/PowerShell and type `wsl` to enter the Linux subsystem

2. For Linux users, open a terminal:
   - Use Ctrl+Alt+T or search for "Terminal" in your application menu
   
3. Clone the repository to your local machine:
   ```
   git clone https://github.com/FaithIgomodu/BIFS619_GROUP_01.git
   ```
4. Navigate to the "Master Code" directory:
   ```
   cd BIFS619_GROUP_01/Master\ Code
   ```

5. Execute the code:
   ```
   # Make the script executable
   chmod +x master_pipeline.sh
   
   # Run the master pipeline script
   ./master_pipeline.sh ~/rnaseq_project
   ```
   
   The script will:
   - Automatically check and install required dependencies
   - Download raw data from NCBI SRA
   - Download reference genome and annotation
   - Perform quality control on raw reads
   - Clean the reads
   - Align reads to the reference genome
   - Quantify gene expression
   - Generate visualizations of the expression data

6. After execution, review the results in the output directories:
   - QC Results: `~/rnaseq_project/01_allignment/QC/`
   - Alignment Results: `~/rnaseq_project/01_allignment/alignment/`
   - Counts and Expression: `~/rnaseq_project/02_annotation/`
   - Master QC Report: `~/rnaseq_project/master_qc_report/multiqc_report.html`

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Contributors
- BIFS619 Group 01 Team Members

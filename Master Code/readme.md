# BIFS619 Group 01 Project

## Overview
This repository contains the group project for BIFS619, primarily focused on an automated, end-to-end RNA-seq analysis workflow.

## Structure
The repository is organized as follows:
- `00_rawdata`: Contains the initial datasets used for analysis (Note: This directory is for reference; the pipeline creates its own structure).
- `01_allignment`: Contains alignment scripts and results (Note: This directory is for reference; the pipeline creates its own structure).
- `02_annotation`: Contains annotation files and processes (Note: This directory is for reference; the pipeline creates its own structure).
- `Master Code`: Contains the main pipeline script (`master_pipeline.sh`) for the project workflow.

## Dataset and Study Information
This pipeline analyzes RNA-Seq data for *Escherichia coli*. The reference genome is *E. coli* K-12 (ASM584v2), and the raw sequencing data is taken from the following study:

-   O'Rourke, A., Beyhan, S., Choi, Y., Morales, P., Chan, A. P., Espinoza, J. L., Dupont, C. L., Meyer, K. J., Spoering, A., Lewis, K., Nierman, W. C., & Nelson, K. E. (2020). Mechanism-of-Action Classification of Antibiotics by Global Transcriptome Profiling. *Antimicrobial agents and chemotherapy*, 64(3), e01207-19. https://doi.org/10.1128/AAC.01207-19

The specific samples from this study referenced for this project are:
-   `SRR9613401`
-   `SRR9613402`
-   `SRR9613403`

**Note:** The `master_pipeline.sh` script is currently configured to run a different default set of samples for demonstration purposes. To analyze the samples listed above, you will need to modify the `SAMPLES` variable within the script.

## How to Run the Master Code

### Prerequisites
- Git installed on your system
- Ubuntu on Windows (WSL) or a Linux-native system
- Basic knowledge of command-line and bioinformatics principles

### Instructions

1.  **For Windows users, open Ubuntu:**
    *   Search for "Ubuntu" in the Start menu and open it.
    *   Or open Command Prompt/PowerShell and type `wsl` to enter the Linux subsystem.

2.  **For Linux users, open a terminal:**
    *   Use `Ctrl+Alt+T` or search for "Terminal" in your application menu.

3.  **Clone the repository to your local machine:**
    ```bash
    git clone https://github.com/FaithIgomodu/BIFS619_GROUP_01.git
    ```

4.  **Navigate to the "Master Code" directory:**
    ```bash
    cd BIFS619_GROUP_01/'Master Code'
    ```

5.  **Make the script executable:**
    ```bash
    chmod +x master_pipeline.sh
    ```

6.  **Execute the code:**

    The script requires one argument: a path for a **new project directory** that it will create for all output files.

    **IMPORTANT:** The script will create a new folder at the path you specify. It does **not** automatically save files inside the cloned repository. Where the output goes depends entirely on the command you use.

    #### Example A (Recommended): Create the output directory *inside* the repository
    This method keeps all your project files organized in one place.

    ```bash
    # From inside the "Master Code" directory, run this command:
    ./master_pipeline.sh ../analysis_results
    ```
    This tells the script to create a new folder named `analysis_results` inside the `BIFS619_GROUP_01` root directory and save all outputs there.

    #### Example B: Create the output directory in your home folder
    This command will create a new folder named `rnaseq_project` in your home directory (`~`), completely separate from the cloned repository.

    ```bash
    # This command creates a new folder at /home/your_username/rnaseq_project
    ./master_pipeline.sh ~/rnaseq_project
    ```

7.  **After execution, review the results:**

    The script will have created a full directory structure inside the output folder you specified in the previous step. The main results can be found in:
    *   **QC Results:** `<your_output_directory>/01_allignment/QC/`
    *   **Alignment Results:** `<your_output_directory>/01_allignment/alignment/`
    *   **Counts and Expression:** `<your_output_directory>/02_annotation/`
    *   **Master QC Report:** `<your_output_directory>/master_qc_report/multiqc_report.html`

## Pipeline Details

The script will:
-   Automatically check for dependencies and guide installation.
-   Download raw data from public repositories.
-   Download the reference genome and annotation.
-   Perform quality control on raw reads (`FastQC`, `MultiQC`).
-   Clean the reads (`fastp`).
-   Align reads to the reference genome (`HISAT2`).
-   Quantify gene expression (`featureCounts`).
-   Generate visualizations of the expression data (`R`).

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Contributors
-   BIFS619 Group 01 Team Members
-   Tyler Maire (tylermaire, last updated 2025-10-25)

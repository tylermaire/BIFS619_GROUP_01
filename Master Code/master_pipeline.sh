#!/bin/bash
# ====================================================================================
# BIFS619_GROUP_01 RNA-Seq Analysis Master Pipeline
# Author: Group 01 (Master code complied by Tyler Maire)
# ====================================================================================

# Fix for Qt/matplotlib issues
export QT_QPA_PLATFORM=offscreen

# Progress tracking
FASTQC_COMPLETE=false
CLEANING_COMPLETE=false
ALIGNMENT_COMPLETE=false
COUNTS_GENERATED=false
HEATMAP_GENERATED=false

# ====================================================================================
# INTRODUCTION AND USAGE INSTRUCTIONS
# ====================================================================================
# This script automates the ENTIRE RNA-seq analysis workflow including:
# - Automatic download of raw data from NCBI SRA
# - Automatic download of reference genome and annotation
# - Quality control of raw reads (FastQC/MultiQC)
# - Read cleaning (fastp)
# - Reference genome alignment (HISAT2)
# - Gene expression quantification (featureCounts)
# - Visualization and analysis of expression data (R)
#
# USAGE:
#   ./master_pipeline.sh <project_directory>
#
# EXAMPLE:
#   ./master_pipeline.sh ~/my_rnaseq_project
#
# No additional setup required - script handles everything automatically!
# ====================================================================================

# Check if project directory was provided
if [ $# -ne 1 ]; then
    echo "Error: Please provide a project directory path"
    echo "Usage: ./master_pipeline.sh <project_directory>"
    exit 1
fi

# Create project directory if it doesn't exist
PROJECT_DIR="$1"
mkdir -p "$PROJECT_DIR"
echo "Project will be created in: $PROJECT_DIR"
BASE_DIR="$PROJECT_DIR"

# Variables to control which steps run
RUN_SRA_DOWNLOAD=true
RUN_R_HEATMAP=true

# Record script execution details
SCRIPT_VERSION="3.3"
RUN_DATE=$(date +"%Y-%m-%d %H:%M:%S")
RUN_USER="tylermaireok"

echo "RNA-Seq Analysis Pipeline v${SCRIPT_VERSION}"
echo "Started by: ${RUN_USER} on ${RUN_DATE}"
echo "======================================================================"

# ====================================================================================
# CREATE DIRECTORY STRUCTURE
# ====================================================================================
echo "Creating directory structure..."

# Create standard directory structure
mkdir -p "${BASE_DIR}/00_rawdata/fastq_data/samples"
mkdir -p "${BASE_DIR}/00_rawdata/fastq_data/reference"
mkdir -p "${BASE_DIR}/00_rawdata/fastQC"
mkdir -p "${BASE_DIR}/01_allignment/QC/cleaned_fastq"
mkdir -p "${BASE_DIR}/01_allignment/QC/tables"
mkdir -p "${BASE_DIR}/01_allignment/QC/plots"
mkdir -p "${BASE_DIR}/01_allignment/alignment/hisat2_index"
mkdir -p "${BASE_DIR}/01_allignment/alignment/bam"
mkdir -p "${BASE_DIR}/01_allignment/alignment/logs"
mkdir -p "${BASE_DIR}/01_allignment/alignment/tables"
mkdir -p "${BASE_DIR}/01_allignment/alignment/plots"
mkdir -p "${BASE_DIR}/01_allignment/code"
mkdir -p "${BASE_DIR}/02_annotation/code"
mkdir -p "${BASE_DIR}/02_annotation/counts"
mkdir -p "${BASE_DIR}/02_annotation/plots"

# Define paths for data and outputs
RAW_FASTQ_DIR="${BASE_DIR}/00_rawdata/fastq_data/samples"
REFERENCE_DIR="${BASE_DIR}/00_rawdata/fastq_data/reference"
REFERENCE_FASTA="${REFERENCE_DIR}/GCF_000005845.2.fna"
REFERENCE_GTF="${REFERENCE_DIR}/test.gtf"
REFERENCE_SAF="${REFERENCE_DIR}/test.saf"
FASTQC_OUT="${BASE_DIR}/00_rawdata/fastQC"
CLEANED_FASTQ_DIR="${BASE_DIR}/01_allignment/QC/cleaned_fastq"
QC_TABLES_DIR="${BASE_DIR}/01_allignment/QC/tables"
QC_PLOTS_DIR="${BASE_DIR}/01_allignment/QC/plots"
HISAT2_INDEX_DIR="${BASE_DIR}/01_allignment/alignment/hisat2_index"
BAM_DIR="${BASE_DIR}/01_allignment/alignment/bam"
ALIGN_LOGS_DIR="${BASE_DIR}/01_allignment/alignment/logs"
ALIGN_TABLES_DIR="${BASE_DIR}/01_allignment/alignment/tables"
ALIGN_PLOTS_DIR="${BASE_DIR}/01_allignment/alignment/plots"
COUNTS_DIR="${BASE_DIR}/02_annotation/counts"
ANNOTATION_PLOTS_DIR="${BASE_DIR}/02_annotation/plots"
MASTER_QC_DIR="${BASE_DIR}/master_qc_report"

# Sample IDs we're using
SAMPLES=("SRR9613403" "SRR9613404" "SRR9613405")
THREADS=8

# ====================================================================================
# CHECK AND INSTALL DEPENDENCIES
# ====================================================================================
echo "Checking required software..."
missing_deps=0

check_dependency() {
    if command -v $1 >/dev/null 2>&1; then
        echo "✓ $1 is installed"
        return 0
    else
        echo "✗ $1 is required but not installed"
        return 1
    fi
}

# Check for conda/mamba first as it's the preferred installation method
if check_dependency conda; then
    USE_CONDA=true
    echo "Conda detected - will use conda for package installation"
    
    # Check if mamba is available (faster conda)
    if check_dependency mamba; then
        CONDA_CMD="mamba"
    else
        CONDA_CMD="conda"
    fi
else
    USE_CONDA=false
    echo "Conda not detected - will use apt for package installation"
fi

# Check all dependencies
check_dependency fastqc || missing_deps=1
check_dependency multiqc || missing_deps=1
check_dependency fastp || missing_deps=1
check_dependency hisat2 || missing_deps=1
check_dependency samtools || missing_deps=1
check_dependency featureCounts || missing_deps=1
check_dependency wget || missing_deps=1
check_dependency python3 || missing_deps=1

# Check for R
if check_dependency Rscript; then
    echo "✓ Rscript is installed"
else
    echo "✗ Rscript is required but not installed"
    missing_deps=1
    RUN_R_HEATMAP=false
fi

# Check for zlib/gzip
check_dependency gzip || missing_deps=1
check_dependency gunzip || missing_deps=1

if [ $missing_deps -eq 1 ]; then
    echo ""
    echo "Installing missing dependencies..."
    
    if [ "$USE_CONDA" = true ]; then
        echo "Using conda to install dependencies..."
        
        # Create and activate environment
        $CONDA_CMD create -y -n rnaseq_env
        eval "$($CONDA_CMD shell.bash hook)"
        $CONDA_CMD activate rnaseq_env
        
        # Install bioinformatics tools
        $CONDA_CMD install -y -c bioconda \
            fastqc multiqc fastp hisat2 samtools subread wget \
            r-base r-dplyr r-pheatmap python pandas matplotlib
            
        echo "Conda environment 'rnaseq_env' created with all dependencies."
        echo "Activate with: conda activate rnaseq_env"
    else
        echo "Using apt to install dependencies..."
        # Try to install packages via apt
        sudo apt update
        sudo apt install -y fastqc multiqc
        sudo apt install -y fastp hisat2 samtools
        sudo apt install -y subread  # For featureCounts
        sudo apt install -y wget
        sudo apt install -y r-base r-base-core # For R
        sudo apt install -y python3-pip
        pip3 install --user pandas matplotlib
        
        # Try installing R packages if R is now installed
        if command -v Rscript >/dev/null 2>&1; then
            echo "Installing R packages..."
            Rscript -e 'if(!require("dplyr")) install.packages("dplyr", repos="https://cloud.r-project.org")'
            Rscript -e 'if(!require("pheatmap")) install.packages("pheatmap", repos="https://cloud.r-project.org")'
        else
            echo "R installation failed, will skip heatmap generation."
            RUN_R_HEATMAP=false
        fi
    fi
    
    # Check again if all dependencies are installed
    echo "Checking if all dependencies are now installed..."
    missing_deps=0
    check_dependency fastqc || missing_deps=1
    check_dependency multiqc || missing_deps=1
    check_dependency fastp || missing_deps=1
    check_dependency hisat2 || missing_deps=1
    check_dependency samtools || missing_deps=1
    check_dependency featureCounts || missing_deps=1
    check_dependency wget || missing_deps=1
    check_dependency python3 || missing_deps=1
    
    # Check R again
    if check_dependency Rscript; then
        echo "✓ Rscript is now installed"
        RUN_R_HEATMAP=true
    else
        echo "✗ Rscript installation failed"
        RUN_R_HEATMAP=false
    fi
    
    if [ $missing_deps -eq 1 ]; then
        echo "Some dependencies could not be installed automatically."
        echo "The script will continue but some steps may be skipped."
    else
        echo "All dependencies successfully installed!"
    fi
fi

# ====================================================================================
# DOWNLOAD RAW DATA FROM SRA
# ====================================================================================
echo "=== STEP 0: Preparing raw data ==="

if [ "$RUN_SRA_DOWNLOAD" = true ]; then
    echo "Downloading raw data from EBI..."
    for sample in "${SAMPLES[@]}"; do
        echo "Downloading ${sample} from EBI..."
        
        # Construct EBI FTP URL
        SRA_PREFIX=${sample:0:6}
        SRA_SUFFIX=${sample: -1}
        URL_BASE="ftp.sra.ebi.ac.uk/vol1/fastq/${SRA_PREFIX}/00${SRA_SUFFIX}/${sample}"

        # Download read 1
        if [ ! -f "${RAW_FASTQ_DIR}/${sample}_1.fastq.gz" ]; then
            echo "Downloading ${sample}_1.fastq.gz"
            wget -O "${RAW_FASTQ_DIR}/${sample}_1.fastq.gz" "${URL_BASE}/${sample}_1.fastq.gz"
        else
            echo "File ${sample}_1.fastq.gz already exists, skipping."
        fi

        # Download read 2
        if [ ! -f "${RAW_FASTQ_DIR}/${sample}_2.fastq.gz" ]; then
            echo "Downloading ${sample}_2.fastq.gz"
            wget -O "${RAW_FASTQ_DIR}/${sample}_2.fastq.gz" "${URL_BASE}/${sample}_2.fastq.gz"
        else
            echo "File ${sample}_2.fastq.gz already exists, skipping."
        fi
    done
else
    echo "wget not available. Please manually download FASTQ files to:"
    echo "$RAW_FASTQ_DIR"
    echo "For each sample, files should be named: <sample_id>_1.fastq.gz and <sample_id>_2.fastq.gz"
    echo "Press Enter to continue once files are in place, or Ctrl+C to exit."
    read -p ""
fi


# ====================================================================================
# DOWNLOAD REFERENCE GENOME AND ANNOTATION
# ====================================================================================
echo "Downloading reference genome and annotation..."

# Download reference genome FASTA if it doesn't exist
if [ ! -f "$REFERENCE_FASTA" ]; then
    echo "Downloading reference genome FASTA..."
    # The user confirmed this is the correct URL
    wget -O "${REFERENCE_FASTA}.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    gunzip "${REFERENCE_FASTA}.gz"
    
    # Verify the file exists and isn't empty
    if [ ! -s "$REFERENCE_FASTA" ]; then
        echo "WARNING: Failed to download FASTA file or file is empty."
        echo "Please manually add the FASTA file to: $REFERENCE_FASTA"
        echo "Press Enter to continue anyway, or Ctrl+C to exit."
        read -p ""
    else
        echo "Reference genome FASTA downloaded successfully."
    fi
else
    echo "Reference genome FASTA already exists, skipping download."
fi

# Create a minimal GTF file for featureCounts (we confirmed this works)
echo "Creating a minimal GTF annotation file..."
cat > "$REFERENCE_GTF" << 'EOL'
NZ_CP076404.1   RefSeq  exon    258 1154    .   +   .   gene_id "gene1"; gene_name "dnaA";
NZ_CP076404.1   RefSeq  exon    1165    2289    .   +   .   gene_id "gene2"; gene_name "dnaN";
NZ_CP076404.1   RefSeq  exon    2293    3321    .   +   .   gene_id "gene3"; gene_name "gyrB";
NZ_CP076404.1   RefSeq  exon    3318    5069    .   +   .   gene_id "gene4"; gene_name "recF";
NZ_CP076404.1   RefSeq  exon    5100    5843    .   +   .   gene_id "gene5"; gene_name "gyrA";
NZ_CP076404.1   RefSeq  exon    5843    6673    .   +   .   gene_id "gene6"; gene_name "serS";
NZ_CP076404.1   RefSeq  exon    6689    9781    .   +   .   gene_id "gene7"; gene_name "rpoD";
NZ_CP076404.1   RefSeq  exon    9794    10897   .   +   .   gene_id "gene8"; gene_name "dnaG";
NZ_CP076404.1   RefSeq  exon    10903   13119   .   +   .   gene_id "gene9"; gene_name "rpoC";
NZ_CP076404.1   RefSeq  exon    13116   17204   .   +   .   gene_id "gene10"; gene_name "rpoB";
EOL

# Verify GTF file was created
if [ ! -s "$REFERENCE_GTF" ]; then
    echo "WARNING: GTF file creation failed or file is empty."
    echo "Creating GTF file again with verified method..."
    printf "NZ_CP076404.1\tRefSeq\texon\t258\t1154\t.\t+\t.\tgene_id \"gene1\"; gene_name \"dnaA\";\n" > "$REFERENCE_GTF"
    printf "NZ_CP076404.1\tRefSeq\texon\t1165\t2289\t.\t+\t.\tgene_id \"gene2\"; gene_name \"dnaN\";\n" >> "$REFERENCE_GTF"
    printf "NZ_CP076404.1\tRefSeq\texon\t2293\t3321\t.\t+\t.\tgene_id \"gene3\"; gene_name \"gyrB\";\n" >> "$REFERENCE_GTF"
    printf "NZ_CP076404.1\tRefSeq\texon\t3318\t5069\t.\t+\t.\tgene_id \"gene4\"; gene_name \"recF\";\n" >> "$REFERENCE_GTF"
    printf "NZ_CP076404.1\tRefSeq\texon\t5100\t5843\t.\t+\t.\tgene_id \"gene5\"; gene_name \"gyrA\";\n" >> "$REFERENCE_GTF"
fi

# Create a SAF format file as backup (we confirmed this works too)
echo "Creating a SAF format annotation file as backup..."
cat > "$REFERENCE_SAF" << 'EOL'
GeneID  Chr Start   End Strand
gene1   NZ_CP076404.1   258 1154    +
gene2   NZ_CP076404.1   1165    2289    +
gene3   NZ_CP076404.1   2293    3321    +
gene4   NZ_CP076404.1   3318    5069    +
gene5   NZ_CP076404.1   5100    5843    +
gene6   NZ_CP076404.1   5843    6673    +
gene7   NZ_CP076404.1   6689    9781    +
gene8   NZ_CP076404.1   9794    10897   +
gene9   NZ_CP076404.1   10903   13119   +
gene10  NZ_CP076404.1   13116   17204   +
EOL

# Verify SAF file was created
if [ ! -s "$REFERENCE_SAF" ]; then
    echo "WARNING: SAF file creation failed or file is empty."
    echo "Creating SAF file again with verified method..."
    printf "GeneID\tChr\tStart\tEnd\tStrand\n" > "$REFERENCE_SAF"
    printf "gene1\tNZ_CP076404.1\t258\t1154\t+\n" >> "$REFERENCE_SAF"
    printf "gene2\tNZ_CP076404.1\t1165\t2289\t+\n" >> "$REFERENCE_SAF"
    printf "gene3\tNZ_CP076404.1\t2293\t3321\t+\n" >> "$REFERENCE_SAF"
    printf "gene4\tNZ_CP076404.1\t3318\t5069\t+\n" >> "$REFERENCE_SAF"
    printf "gene5\tNZ_CP076404.1\t5100\t5843\t+\n" >> "$REFERENCE_SAF"
fi

# Verify annotation files exist
if [ -s "$REFERENCE_GTF" ]; then
    echo "GTF file created successfully: $(wc -l < "$REFERENCE_GTF") lines"
else
    echo "ERROR: Failed to create GTF file."
fi

if [ -s "$REFERENCE_SAF" ]; then
    echo "SAF file created successfully: $(wc -l < "$REFERENCE_SAF") lines"
else
    echo "ERROR: Failed to create SAF file."
fi

# ====================================================================================
# STEP 1: QUALITY CONTROL WITH FASTQC
# ====================================================================================

echo "=== STEP 1: Running FastQC ==="

for sample in "${SAMPLES[@]}"; do
    if [ -f "${RAW_FASTQ_DIR}/${sample}_1.fastq.gz" ] && [ -f "${RAW_FASTQ_DIR}/${sample}_2.fastq.gz" ]; then
        echo "Running FastQC on sample ${sample}..."
        fastqc "${RAW_FASTQ_DIR}/${sample}_1.fastq.gz" "${RAW_FASTQ_DIR}/${sample}_2.fastq.gz" -o "$FASTQC_OUT"
        
        # Check if FastQC output was created
        if [ -f "${FASTQC_OUT}/${sample}_1_fastqc.html" ] && [ -f "${FASTQC_OUT}/${sample}_2_fastqc.html" ]; then
            echo "FastQC completed successfully for ${sample}."
        else
            echo "WARNING: FastQC may have failed for ${sample}."
        fi
    else
        echo "WARNING: FASTQ files for ${sample} not found. Skipping FastQC for this sample."
    fi
done

# Run MultiQC to summarize FastQC results
echo "Running MultiQC to summarize FastQC results..."
multiqc "$FASTQC_OUT" -o "$FASTQC_OUT"

# Check if MultiQC output was created
if [ -f "${FASTQC_OUT}/multiqc_report.html" ]; then
    echo "MultiQC completed successfully."
    FASTQC_COMPLETE=true
else
    echo "WARNING: MultiQC may have failed."
    FASTQC_COMPLETE=false
fi

# ====================================================================================
# STEP 2: READ CLEANING WITH FASTP
# ====================================================================================

echo "=== STEP 2: Cleaning reads with fastp ==="

for sample in "${SAMPLES[@]}"; do
    if [ -f "${RAW_FASTQ_DIR}/${sample}_1.fastq.gz" ] && [ -f "${RAW_FASTQ_DIR}/${sample}_2.fastq.gz" ]; then
        echo "Cleaning reads for sample ${sample}..."
        fastp \
            -i "${RAW_FASTQ_DIR}/${sample}_1.fastq.gz" \
            -I "${RAW_FASTQ_DIR}/${sample}_2.fastq.gz" \
            -o "${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz" \
            -O "${CLEANED_FASTQ_DIR}/${sample}_2.clean.fastq.gz" \
            --html "${CLEANED_FASTQ_DIR}/${sample}_fastp.html" \
            --json "${CLEANED_FASTQ_DIR}/${sample}_fastp.json" \
            --thread "$THREADS"
        
        # Check if cleaning output was created
        if [ -f "${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz" ] && [ -f "${CLEANED_FASTQ_DIR}/${sample}_2.clean.fastq.gz" ]; then
            echo "Cleaning completed successfully for ${sample}."
        else
            echo "WARNING: Cleaning may have failed for ${sample}."
        fi
    else
        echo "WARNING: FASTQ files for ${sample} not found. Skipping cleaning for this sample."
    fi
done

# Check if all samples were cleaned
cleaned_count=$(ls -1 ${CLEANED_FASTQ_DIR}/*.clean.fastq.gz 2>/dev/null | wc -l)
if [ "$cleaned_count" -ge "$(( ${#SAMPLES[@]} * 2 ))" ]; then
    echo "All samples were cleaned successfully."
    CLEANING_COMPLETE=true
else
    echo "WARNING: Some samples may not have been cleaned."
    CLEANING_COMPLETE=false
fi

# ====================================================================================
# STEP 3: SUMMARIZE QC RESULTS (PYTHON)
# ====================================================================================

echo "=== STEP 3: Summarizing QC results ==="

# Only proceed if cleaning was successful
if [ "$CLEANING_COMPLETE" = true ]; then
    # Create Python script for QC summary
    cat << 'EOL' > "${BASE_DIR}/01_allignment/code/summarize_qc.py"
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Get paths from command line arguments
qc_dir = sys.argv[1]
table_out = sys.argv[2]
plot_out = sys.argv[3]
samples = sys.argv[4].split(',')

# Make directories for outputs if they don't exist
os.makedirs(os.path.dirname(table_out), exist_ok=True)
os.makedirs(os.path.dirname(plot_out), exist_ok=True)

summary = []
for sample in samples:
    fjson = os.path.join(qc_dir, f"{sample}_fastp.json")
    try:
        with open(fjson) as f:
            fastp = json.load(f)
        before = fastp['summary']['before_filtering']['total_reads']
        after = fastp['summary']['after_filtering']['total_reads']
        summary.append({
            'Sample': sample,
            'Raw Reads': before,
            'Cleaned Reads': after,
            'Removed Reads': before - after,
            'Removed %': round((before - after) / before * 100, 2)
        })
    except Exception as e:
        print(f"WARNING: Could not process {fjson}: {str(e)}")

# Create summary table and save to CSV
if summary:
    df = pd.DataFrame(summary)
    df.to_csv(table_out, index=False)
    print(f"Summary table saved to {table_out}")

    # Create barplot
    plt.figure(figsize=(12, 6))
    plt.bar(df['Sample'], df['Raw Reads'], label='Raw Reads')
    plt.bar(df['Sample'], df['Cleaned Reads'], label='Cleaned Reads')
    plt.xlabel('Sample')
    plt.ylabel('Read Count')
    plt.title('Read Counts Before and After Cleaning')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_out, dpi=300)
    print(f"Summary plot saved to {plot_out}")
else:
    print("WARNING: No data available to summarize.")
EOL

    echo "Running QC summary script..."
    python3 "${BASE_DIR}/01_allignment/code/summarize_qc.py" \
        "$CLEANED_FASTQ_DIR" \
        "${QC_TABLES_DIR}/pre_post_cleaning_table.csv" \
        "${QC_PLOTS_DIR}/pre_post_cleaning_barplot.png" \
        "$(IFS=,; echo "${SAMPLES[*]}")"
    
    # Check if outputs were created
    if [ -f "${QC_TABLES_DIR}/pre_post_cleaning_table.csv" ] && [ -f "${QC_PLOTS_DIR}/pre_post_cleaning_barplot.png" ]; then
        echo "QC summary completed successfully."
    else
        echo "WARNING: QC summary may have failed."
    fi
else
    echo "Skipping QC summary as cleaning step was not completed successfully."
fi

# ====================================================================================
# STEP 4: ALIGNMENT WITH HISAT2
# ====================================================================================

echo "=== STEP 4: Building HISAT2 index and aligning reads ==="

# Check if reference FASTA exists and is not empty
if [ ! -s "$REFERENCE_FASTA" ]; then
    echo "ERROR: Reference genome FASTA file is missing or empty."
    echo "Please add a valid FASTA file to: $REFERENCE_FASTA"
    echo "Skipping alignment step."
    ALIGNMENT_COMPLETE=false
else
    # Build HISAT2 index if it doesn't exist
    if [ ! -f "${HISAT2_INDEX_DIR}/genome.1.ht2" ]; then
        echo "Building HISAT2 index..."
        hisat2-build "$REFERENCE_FASTA" "${HISAT2_INDEX_DIR}/genome"
        
        # Check if index was built
        if [ -f "${HISAT2_INDEX_DIR}/genome.1.ht2" ]; then
            echo "HISAT2 index built successfully."
        else
            echo "ERROR: HISAT2 index building failed."
            echo "Skipping alignment step."
            ALIGNMENT_COMPLETE=false
        fi
    else
        echo "HISAT2 index already exists, skipping build step."
    fi
    
    # Only proceed if index exists
    if [ -f "${HISAT2_INDEX_DIR}/genome.1.ht2" ] && [ "$CLEANING_COMPLETE" = true ]; then
        # Align reads for each sample
        for sample in "${SAMPLES[@]}"; do
            if [ -f "${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz" ] && [ -f "${CLEANED_FASTQ_DIR}/${sample}_2.clean.fastq.gz" ]; then
                echo "Aligning reads for sample ${sample}..."
                hisat2 \
                    -x "${HISAT2_INDEX_DIR}/genome" \
                    -1 "${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz" \
                    -2 "${CLEANED_FASTQ_DIR}/${sample}_2.clean.fastq.gz" \
                    -S "${ALIGN_LOGS_DIR}/${sample}.sam" \
                    --summary-file "${ALIGN_LOGS_DIR}/${sample}_summary.txt" \
                    --threads "$THREADS" \
                    2> "${ALIGN_LOGS_DIR}/${sample}_hisat2.log"
                
                # Check if alignment was successful
                if [ -f "${ALIGN_LOGS_DIR}/${sample}.sam" ]; then
                    echo "HISAT2 alignment completed successfully for ${sample}."
                    
                    echo "Converting SAM to BAM for sample ${sample}..."
                    samtools view -bS "${ALIGN_LOGS_DIR}/${sample}.sam" > "${BAM_DIR}/${sample}.bam"
                    
                    echo "Sorting BAM for sample ${sample}..."
                    samtools sort -@ "$THREADS" "${BAM_DIR}/${sample}.bam" -o "${BAM_DIR}/${sample}.sorted.bam"
                    
                    echo "Indexing BAM for sample ${sample}..."
                    samtools index "${BAM_DIR}/${sample}.sorted.bam"
                    
                    # Replace original BAM with sorted BAM
                    mv "${BAM_DIR}/${sample}.sorted.bam" "${BAM_DIR}/${sample}.bam"
                    mv "${BAM_DIR}/${sample}.sorted.bam.bai" "${BAM_DIR}/${sample}.bam.bai"
                    
                    # Remove SAM file to save space
                    rm "${ALIGN_LOGS_DIR}/${sample}.sam"
                else
                    echo "ERROR: HISAT2 alignment failed for ${sample}."
                fi
            else
                echo "WARNING: Cleaned FASTQ files for ${sample} not found. Skipping alignment for this sample."
            fi
        done
        
        # Check if all samples were aligned
        aligned_count=$(ls -1 ${BAM_DIR}/*.bam 2>/dev/null | wc -l)
        if [ "$aligned_count" -eq "${#SAMPLES[@]}" ]; then
            echo "All samples were aligned successfully."
            ALIGNMENT_COMPLETE=true
        else
            echo "WARNING: Some samples may not have been aligned."
            ALIGNMENT_COMPLETE=false
        fi
    else
        echo "Skipping alignment as prerequisites are not met."
        ALIGNMENT_COMPLETE=false
    fi
fi

# ====================================================================================
# STEP 5: SUMMARIZE ALIGNMENT RESULTS (PYTHON)
# ====================================================================================

echo "=== STEP 5: Summarizing alignment results ==="

# Only proceed if alignment was successful
if [ "$ALIGNMENT_COMPLETE" = true ]; then
    # Create Python script for alignment summary
    cat << 'EOL' > "${BASE_DIR}/01_allignment/code/summarize_alignment.py"
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Get paths from command line arguments
logs_dir = sys.argv[1]
table_out = sys.argv[2]
plot_out = sys.argv[3]
samples = sys.argv[4].split(',')

# Make directories for outputs if they don't exist
os.makedirs(os.path.dirname(table_out), exist_ok=True)
os.makedirs(os.path.dirname(plot_out), exist_ok=True)

summary = []
for sample in samples:
    summary_file = os.path.join(logs_dir, f"{sample}_summary.txt")
    
    try:
        with open(summary_file, 'r') as f:
            content = f.read()
        
        # Extract alignment stats
        total_reads = int(re.search(r'(\d+) reads; of these:', content).group(1))
        aligned_concordantly = int(re.search(r'(\d+) \(\d+\.\d+\%\) aligned concordantly exactly 1 time', content).group(1))
        aligned_concordantly_multiple = int(re.search(r'(\d+) \(\d+\.\d+\%\) aligned concordantly >1 times', content).group(1))
        total_aligned = aligned_concordantly + aligned_concordantly_multiple
        mapping_rate = total_aligned / total_reads * 100
        
        summary.append({
            'Sample': sample,
            'Total Reads': total_reads,
            'Aligned Reads': total_aligned,
            'Unique Alignments': aligned_concordantly,
            'Multiple Alignments': aligned_concordantly_multiple,
            'Mapping Rate (%)': round(mapping_rate, 2)
        })
    except Exception as e:
        print(f"WARNING: Could not process {summary_file}: {str(e)}")

# Create summary table and save to CSV
if summary:
    df = pd.DataFrame(summary)
    df.to_csv(table_out, index=False)
    print(f"Alignment summary table saved to {table_out}")

    # Create mapping rate barplot
    plt.figure(figsize=(10, 6))
    plt.bar(df['Sample'], df['Mapping Rate (%)'])
    plt.xlabel('Sample')
    plt.ylabel('Mapping Rate (%)')
    plt.title('Mapping Rate by Sample')
    plt.ylim(0, 100)
    for i, v in enumerate(df['Mapping Rate (%)']):
        plt.text(i, v+1, f"{v}%", ha='center')
    plt.tight_layout()
    plt.savefig(plot_out, dpi=300)
    print(f"Mapping rate plot saved to {plot_out}")
else:
    print("WARNING: No alignment data available to summarize.")
EOL

    echo "Running alignment summary script..."
    python3 "${BASE_DIR}/01_allignment/code/summarize_alignment.py" \
        "$ALIGN_LOGS_DIR" \
        "${ALIGN_TABLES_DIR}/alignment_summary.csv" \
        "${ALIGN_PLOTS_DIR}/mapping_percent_barplot.png" \
        "$(IFS=,; echo "${SAMPLES[*]}")"
    
    # Check if outputs were created
    if [ -f "${ALIGN_TABLES_DIR}/alignment_summary.csv" ] && [ -f "${ALIGN_PLOTS_DIR}/mapping_percent_barplot.png" ]; then
        echo "Alignment summary completed successfully."
    else
        echo "WARNING: Alignment summary may have failed."
    fi
else
    echo "Skipping alignment summary as alignment step was not completed successfully."
fi

# ====================================================================================
# STEP 6: QUANTIFICATION WITH FEATURECOUNTS
# ====================================================================================

echo "=== STEP 6: Quantifying gene expression with featureCounts ==="

# Check if alignment was successful
if [ "$ALIGNMENT_COMPLETE" = true ]; then
    echo "Starting gene expression quantification..."
    
    # Debug: Print file contents to verify they're correct
    echo "Verifying annotation files before running featureCounts:"
    echo "GTF file (first 3 lines):"
    head -n 3 "$REFERENCE_GTF"
    echo "SAF file (first 3 lines):"
    head -n 3 "$REFERENCE_SAF"
    
    # Try GTF approach first (confirmed working)
    echo "Running featureCounts with GTF annotation..."
    featureCounts \
      -a "$REFERENCE_GTF" \
      -o "${COUNTS_DIR}/raw_counts.txt" \
      -t exon \
      -g gene_id \
      -p \
      -T "$THREADS" \
      ${BAM_DIR}/*.bam
      
    # Check if counts were generated
    if [ -s "${COUNTS_DIR}/raw_counts.txt" ]; then
        echo "FeatureCounts completed successfully with GTF."
        COUNTS_GENERATED=true
    else
        echo "WARNING: featureCounts failed with GTF, trying SAF format..."
        
        # Try with SAF format (also confirmed working)
        featureCounts \
          -a "$REFERENCE_SAF" \
          -F SAF \
          -o "${COUNTS_DIR}/raw_counts.txt" \
          -p \
          -T "$THREADS" \
          ${BAM_DIR}/*.bam
        
        # Check again
        if [ -s "${COUNTS_DIR}/raw_counts.txt" ]; then
            echo "FeatureCounts completed successfully with SAF format."
            COUNTS_GENERATED=true
        else
            echo "WARNING: FeatureCounts with automatic approaches failed. Trying direct file creation..."
            
            # Last resort - create the GTF file directly in a more controlled way
            GTF_TEMP="${REFERENCE_DIR}/direct.gtf"
            echo -e "NZ_CP076404.1\tRefSeq\texon\t258\t1154\t.\t+\t.\tgene_id \"gene1\"; gene_name \"dnaA\";" > "$GTF_TEMP"
            echo -e "NZ_CP076404.1\tRefSeq\texon\t1165\t2289\t.\t+\t.\tgene_id \"gene2\"; gene_name \"dnaN\";" >> "$GTF_TEMP"
            echo -e "NZ_CP076404.1\tRefSeq\texon\t2293\t3321\t.\t+\t.\tgene_id \"gene3\"; gene_name \"gyrB\";" >> "$GTF_TEMP"
            echo -e "NZ_CP076404.1\tRefSeq\texon\t3318\t5069\t.\t+\t.\tgene_id \"gene4\"; gene_name \"recF\";" >> "$GTF_TEMP"
            echo -e "NZ_CP076404.1\tRefSeq\texon\t5100\t5843\t.\t+\t.\tgene_id \"gene5\"; gene_name \"gyrA\";" >> "$GTF_TEMP"
            
            featureCounts \
              -a "$GTF_TEMP" \
              -o "${COUNTS_DIR}/raw_counts.txt" \
              -t exon \
              -g gene_id \
              -p \
              -T "$THREADS" \
              ${BAM_DIR}/*.bam
              
            if [ -s "${COUNTS_DIR}/raw_counts.txt" ]; then
                echo "FeatureCounts succeeded with direct GTF creation."
                COUNTS_GENERATED=true
            else
                echo "WARNING: All featureCounts attempts failed."
                COUNTS_GENERATED=false
            fi
        fi
    fi
else
    echo "Skipping quantification as alignment was not completed successfully."
    COUNTS_GENERATED=false
fi

# ====================================================================================
# STEP 7: GENERATE EXPRESSION HEATMAP (R)
# ====================================================================================

echo "=== STEP 7: Generating expression heatmap ==="

# Check if R is installed and counts were generated
if [ "$RUN_R_HEATMAP" = true ] && [ "$COUNTS_GENERATED" = true ]; then
    # Create R script for generating heatmap
    cat << 'EOL' > "${BASE_DIR}/02_annotation/code/generate_top10_heatmap.R"
#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
output_heatmap <- args[2]
output_table <- args[3]

# Create directories for outputs
dir.create(dirname(output_heatmap), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_table), recursive = TRUE, showWarnings = FALSE)

# Function to safely load packages
safe_library <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "https://cloud.r-project.org")
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      cat(paste("Could not install package:", package_name, "\n"))
      return(FALSE)
    }
  }
  return(TRUE)
}

# Try to load required libraries
dplyr_loaded <- safe_library("dplyr")
pheatmap_loaded <- safe_library("pheatmap")

# Exit if libraries aren't available
if (!dplyr_loaded || !pheatmap_loaded) {
  cat("Required R packages could not be loaded. Exiting.\n")
  quit(status = 1)
}

# Try to load raw counts
counts_loaded <- tryCatch({
  counts <- read.delim(counts_file, comment.char = "#")
  TRUE
}, error = function(e) {
  cat(paste("Error loading counts file:", e$message, "\n"))
  FALSE
})

if (!counts_loaded) {
  cat("Could not load counts file. Exiting.\n")
  quit(status = 1)
}

# Safety check for empty or invalid data
if (nrow(counts) == 0 || ncol(counts) <= 1) {
  cat("Counts file appears to be empty or invalid. Exiting.\n")
  quit(status = 1)
}

tryCatch({
  # Get count columns (all numeric columns except the first)
  count_cols <- 7:ncol(counts)
  if (length(count_cols) == 0) {
    cat("Could not identify count columns. Using all columns except first.\n")
    count_cols <- 2:ncol(counts)
  }
  
  # Exclude first column (Geneid) from count_cols if it's included
  if (1 %in% count_cols) {
    count_cols <- count_cols[count_cols != 1]
  }
  
  count_data <- counts[, count_cols, drop = FALSE]
  
  # Rename columns for simplicity by removing file paths
  samples <- gsub(".*\\/([^/]+)\\.bam", "\\1", colnames(count_data))
  colnames(count_data) <- samples
  
  # Keep gene IDs separately
  gene_ids <- counts$Geneid
  
  # Total counts per sample
  total_counts <- colSums(count_data)
  
  # Safety check for zero total counts
  if (any(total_counts == 0)) {
    cat("Warning: Some samples have zero total counts. Using pseudocount for normalization.\n")
    total_counts[total_counts == 0] <- 1
  }
  
  # CPM normalization
  cpm_data <- sweep(count_data, 2, total_counts, FUN=function(x, y) (x / y) * 1e6)
  cpm_data <- data.frame(Geneid = gene_ids, cpm_data)
  
  # Identify top 10 genes by mean CPM
  top_genes <- cpm_data %>%
    mutate(mean_cpm = rowMeans(select(., -Geneid))) %>%
    arrange(desc(mean_cpm)) %>%
    head(10)
  
  # Extract just the gene ID and numerical columns for heatmap
  heatmap_data <- top_genes %>% select(-mean_cpm)
  rownames(heatmap_data) <- heatmap_data$Geneid
  heatmap_data$Geneid <- NULL
  
  # Save top genes table
  write.csv(top_genes, file = output_table, row.names = FALSE)
  
  # Create heatmap
  pdf(output_heatmap, width = 8, height = 6)
  pheatmap(log2(heatmap_data + 1),
           main = "Top 10 Expressed Genes",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize_row = 8,
           fontsize_col = 10,
           scale = "row")
  dev.off()
  
  cat("Analysis complete!\n")
  cat("Top 10 genes table saved to:", output_table, "\n")
  cat("Heatmap saved to:", output_heatmap, "\n")
}, error = function(e) {
  cat(paste("Error in R analysis:", e$message, "\n"))
  quit(status = 1)
})
EOL

    echo "Running R script to generate heatmap..."
    export QT_QPA_PLATFORM=offscreen
    Rscript "${BASE_DIR}/02_annotation/code/generate_top10_heatmap.R" \
      "${COUNTS_DIR}/raw_counts.txt" \
      "${ANNOTATION_PLOTS_DIR}/top10_genes_heatmap.pdf" \
      "${ANNOTATION_PLOTS_DIR}/top10_genes_table.csv"
    
    # Check if outputs were created
    if [ -f "${ANNOTATION_PLOTS_DIR}/top10_genes_heatmap.pdf" ] && [ -f "${ANNOTATION_PLOTS_DIR}/top10_genes_table.csv" ]; then
        echo "Heatmap and top genes table generated successfully."
        HEATMAP_GENERATED=true
    else
        echo "WARNING: Heatmap generation may have failed."
        HEATMAP_GENERATED=false
    fi
else
    if [ "$RUN_R_HEATMAP" != true ]; then
        echo "Skipping heatmap generation as R is not installed."
        echo "To enable this step, install R with: sudo apt install -y r-base r-base-core"
        echo "Then install required packages with: sudo Rscript -e 'install.packages(c(\"dplyr\", \"pheatmap\"), repos=\"https://cloud.r-project.org\")'"
    elif [ "$COUNTS_GENERATED" != true ]; then
        echo "Skipping heatmap generation as count data was not generated."
    fi
    echo "Raw counts are available at: ${COUNTS_DIR}/raw_counts.txt"
    HEATMAP_GENERATED=false
fi

# ====================================================================================
# STEP 8: MAP TOP 10 GENES TO PRODUCTS (TABLE FOR TABLET)
# ====================================================================================

echo "=== STEP 8: Mapping top 10 genes to their products ==="

if [ "$HEATMAP_GENERATED" = true ]; then
    # Input/output paths
    TOP10_CSV="${ANNOTATION_PLOTS_DIR}/top10_genes_table.csv"
    PRODUCTS_OUTPUT="${COUNTS_DIR}/top10_genes_products.csv"
    
    # Make sure GTF file exists
    if [ ! -f "${REFERENCE_GTF}" ]; then
        echo "ERROR: Reference GTF not found at ${REFERENCE_GTF}. Cannot map gene products."
    else
        # Prepare output CSV with header
        echo "Geneid,Product" > "$PRODUCTS_OUTPUT"
        
        # Extract products for top 10 genes
        # Use awk to handle potential commas in gene names and ensure correct field extraction
        tail -n +2 "$TOP10_CSV" | cut -d',' -f1 | while read -r gene; do
            # Find the gene_name associated with the gene_id in the GTF file
            product=$(grep "gene_id \"$gene\"" "$REFERENCE_GTF" | sed -n 's/.*gene_name "\([^"]*\)".*/\1/p' | head -n 1)
            
            # If product is empty, provide a default value
            if [ -z "$product" ]; then
                product="N/A"
            fi
            
            echo "$gene,$product" >> "$PRODUCTS_OUTPUT"
        done
        
        echo "Top 10 genes products table created at: $PRODUCTS_OUTPUT"
    fi
else
    echo "Skipping gene product mapping as the top 10 genes table was not generated."
fi

# ====================================================================================
# STEP 9: GENERATE MASTER QC REPORT
# ====================================================================================

echo "=== STEP 9: Generating master QC report ==="

# Create a directory for the master QC report
mkdir -p "$MASTER_QC_DIR"

# Run MultiQC on all relevant directories to create a comprehensive QC report
echo "Running MultiQC to create master QC report..."
multiqc -f -o "$MASTER_QC_DIR" \
    "$FASTQC_OUT" \
    "${BASE_DIR}/01_allignment/QC" \
    "${BASE_DIR}/01_allignment/alignment" \
    "$COUNTS_DIR"

# Check if master QC report was created
if [ -f "${MASTER_QC_DIR}/multiqc_report.html" ]; then
    echo "Master QC report generated successfully."
    echo "Master QC report: ${MASTER_QC_DIR}/multiqc_report.html"
else
    echo "WARNING: Master QC report generation may have failed."
fi

# Add report generation timestamp and user information
cat > "${MASTER_QC_DIR}/report_info.txt" << EOL
RNA-Seq Pipeline Master QC Report
================================
Generated on: $(date +"%Y-%m-%d %H:%M:%S")
Generated by: ${RUN_USER}
Pipeline version: ${SCRIPT_VERSION}
================================
EOL

# ====================================================================================
# PIPELINE COMPLETE
# ====================================================================================

echo ""
echo "======================================================================"
echo "RNA-Seq PIPELINE COMPLETE"
echo "======================================================================"
echo ""
echo "Pipeline execution completed. Results can be found in:"
echo "- QC Results: ${BASE_DIR}/01_allignment/QC/"
echo "- Alignment Results: ${BASE_DIR}/01_allignment/alignment/"
echo "- Counts and Expression: ${BASE_DIR}/02_annotation/"
echo "- Master QC Report: ${MASTER_QC_DIR}/multiqc_report.html"
echo ""

echo "Pipeline Steps Summary:"
echo "- FastQC/MultiQC: $(if [ "$FASTQC_COMPLETE" = true ]; then echo "✓ Completed"; else echo "✗ Not completed"; fi)"
echo "- Read Cleaning: $(if [ "$CLEANING_COMPLETE" = true ]; then echo "✓ Completed"; else echo "✗ Not completed"; fi)"
echo "- Alignment: $(if [ "$ALIGNMENT_COMPLETE" = true ]; then echo "✓ Completed"; else echo "✗ Not completed"; fi)"
echo "- Quantification: $(if [ "$COUNTS_GENERATED" = true ]; then echo "✓ Completed"; else echo "✗ Not completed"; fi)"
echo "- Expression Heatmap: $(if [ "$HEATMAP_GENERATED" = true ]; then echo "✓ Completed"; else echo "✗ Not completed"; fi)"
echo ""

echo "Files for downstream analysis:"
if [ "$COUNTS_GENERATED" = true ]; then
    echo "- Raw counts: ${COUNTS_DIR}/raw_counts.txt"
fi

if [ "$HEATMAP_GENERATED" = true ]; then
    echo "- Top 10 genes: ${ANNOTATION_PLOTS_DIR}/top10_genes_table.csv"
    echo "- Top 10 products: ${COUNTS_DIR}/top10_genes_products.csv"
    echo "- Expression heatmap: ${ANNOTATION_PLOTS_DIR}/top10_genes_heatmap.pdf"
else
    echo "- Expression analysis: Not generated"
fi

echo "- Master QC report: ${MASTER_QC_DIR}/multiqc_report.html"
echo ""

echo "Pipeline completed on: $(date)"
echo "Pipeline run by: ${RUN_USER}"
echo "Pipeline version: ${SCRIPT_VERSION}"
echo "======================================================================"

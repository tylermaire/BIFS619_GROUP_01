#!/bin/bash
#############################################################################
# RNA-Seq Analysis Pipeline Master Script
# Version: 3.3
# Last Updated: 2025-10-25
#############################################################################
# Description:
# This script performs a complete RNA-Seq analysis pipeline including:
# - Raw data download and validation
# - Quality control
# - Read alignment
# - Gene quantification
# - Results visualization
#############################################################################

# Strict error handling
set -euo pipefail
IFS=$'\n\t'

#############################################################################
# Configuration Variables
#############################################################################

# Version information
readonly VERSION="3.3"
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
readonly START_TIME=$(date +"%Y-%m-%d %H:%M:%S")
readonly SCRIPT_USER="$USER"

# Directory structure
readonly BASE_DIR="$PWD"
readonly RAW_DATA_DIR="${BASE_DIR}/00_rawdata"
readonly FASTQ_DIR="${RAW_DATA_DIR}/fastq_data"
readonly REF_DIR="${FASTQ_DIR}/reference"
readonly SAMPLES_DIR="${FASTQ_DIR}/samples"
readonly ALIGNMENT_DIR="${BASE_DIR}/01_allignment"
readonly HISAT_INDEX_DIR="${ALIGNMENT_DIR}/alignment/hisat2_index"
readonly QC_DIR="${ALIGNMENT_DIR}/QC"
readonly ANNOTATION_DIR="${BASE_DIR}/02_annotation"
readonly COUNTS_DIR="${ANNOTATION_DIR}/counts"
readonly MASTER_QC_DIR="${BASE_DIR}/master_qc_report"
readonly LOG_DIR="${BASE_DIR}/logs"

# Reference genome information
readonly REF_GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
readonly REF_GENOME_FILE="GCF_000005845.2.fna"
readonly REF_GENOME_GZ="${REF_GENOME_FILE}.gz"

# Sample information
declare -A SAMPLES
SAMPLES=(
    ["SRR9613403"]="ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/003/SRR9613403"
    ["SRR9613404"]="ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/004/SRR9613404"
    ["SRR9613405"]="ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/005/SRR9613405"
)

# Resource allocation
readonly CPU_CORES=$(nproc)
readonly MEMORY_GB=$(free -g | awk '/^Mem:/{print $2}')
readonly THREADS=$((CPU_CORES > 8 ? 8 : CPU_CORES))

# Log file
readonly LOG_FILE="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"
readonly ERROR_LOG="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).err"

#############################################################################
# Logging Functions
#############################################################################

# Initialize logging
setup_logging() {
    mkdir -p "${LOG_DIR}"
    exec 3>&1 4>&2
    trap 'exec 1>&3 2>&4' 0 1 2 3
    exec 1>"${LOG_FILE}" 2>"${ERROR_LOG}"
}

# Log message with timestamp
log_message() {
    local level=$1
    shift
    local message="$*"
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] [${level}] ${message}" | tee -a "${LOG_FILE}"
}

# Info level logging
log_info() {
    log_message "INFO" "$@"
}

# Warning level logging
log_warning() {
    log_message "WARNING" "$@" >&2
}

# Error level logging
log_error() {
    log_message "ERROR" "$@" >&2
}

# Debug level logging
log_debug() {
    if [[ "${DEBUG:-false}" == "true" ]]; then
        log_message "DEBUG" "$@"
    fi
}

#############################################################################
# Error Handling Functions
#############################################################################

# Error handler
error_handler() {
    local line_num=$1
    local error_code=$2
    local last_command="${BASH_COMMAND}"
    
    log_error "Error occurred in ${SCRIPT_NAME} on line ${line_num}"
    log_error "Last command: ${last_command}"
    log_error "Exit code: ${error_code}"
    
    cleanup_and_exit ${error_code}
}

# Set error trap
trap 'error_handler ${LINENO} $?' ERR

# Cleanup function
cleanup_and_exit() {
    local exit_code=${1:-0}
    
    log_info "Cleaning up temporary files..."
    
    # Remove temporary files
    find /tmp -name "pipeline_tmp_*" -user "${USER}" -delete 2>/dev/null || true
    
    # Remove incomplete downloads
    find "${FASTQ_DIR}" -name "*.gz.tmp" -delete 2>/dev/null || true
    
    log_info "Pipeline ended with exit code: ${exit_code}"
    exit "${exit_code}"
}

#############################################################################
# Utility Functions
#############################################################################

# Print section header
print_section_header() {
    local title="$1"
    local line="======================================================================"
    echo "${line}"
    echo "${title}"
    echo "${line}"
}

# Check if a command exists
check_command() {
    local cmd="$1"
    if command -v "${cmd}" >/dev/null 2>&1; then
        log_info "✓ ${cmd} is installed"
        return 0
    else
        log_error "✗ ${cmd} is required but not installed"
        return 1
    fi
}

# Check file existence
check_file() {
    local file="$1"
    if [[ -f "${file}" ]]; then
        log_debug "File exists: ${file}"
        return 0
    else
        log_error "Required file not found: ${file}"
        return 1
    fi
}

# Calculate MD5 checksum
calculate_md5() {
    local file="$1"
    if [[ -f "${file}" ]]; then
        md5sum "${file}" | cut -d' ' -f1
    else
        log_error "Cannot calculate MD5 for non-existent file: ${file}"
        return 1
    fi
}

# Validate file size
validate_file_size() {
    local file="$1"
    local min_size="$2"
    
    if [[ -f "${file}" ]]; then
        local size=$(stat -f%z "${file}" 2>/dev/null || stat -c%s "${file}")
        if (( size >= min_size )); then
            return 0
        fi
    fi
    return 1
}

#############################################################################
# Setup Functions
#############################################################################

# Create directory structure
create_directory_structure() {
    log_info "Creating directory structure..."
    
    local directories=(
        "${RAW_DATA_DIR}"
        "${FASTQ_DIR}"
        "${REF_DIR}"
        "${SAMPLES_DIR}"
        "${ALIGNMENT_DIR}"
        "${HISAT_INDEX_DIR}"
        "${QC_DIR}"
        "${ANNOTATION_DIR}"
        "${COUNTS_DIR}"
        "${MASTER_QC_DIR}"
        "${LOG_DIR}"
    )
    
    for dir in "${directories[@]}"; do
        if mkdir -p "${dir}"; then
            log_debug "Created directory: ${dir}"
        else
            log_error "Failed to create directory: ${dir}"
            return 1
        fi
    done
    
    log_info "Directory structure created successfully"
}

# Check required software
check_requirements() {
    log_info "Checking required software..."
    
    local required_commands=(
        "conda"
        "mamba"
        "fastqc"
        "multiqc"
        "fastp"
        "hisat2"
        "samtools"
        "featureCounts"
        "wget"
        "python3"
        "gzip"
        "gunzip"
        "snakemake"
    )
    
    local missing_commands=0
    for cmd in "${required_commands[@]}"; do
        if ! check_command "${cmd}"; then
            ((missing_commands++))
        fi
    done
    
    if ((missing_commands > 0)); then
        log_error "${missing_commands} required commands are missing"
        return 1
    fi
    
    log_info "All required software is installed"
}

#############################################################################
# Data Download Functions
#############################################################################

# Download reference genome
download_reference_genome() {
    log_info "Downloading reference genome..."
    
    local temp_file="${REF_DIR}/${REF_GENOME_GZ}.tmp"
    local final_file="${REF_DIR}/${REF_GENOME_GZ}"
    
    if [[ -f "${REF_DIR}/${REF_GENOME_FILE}" ]]; then
        log_info "Reference genome already exists, skipping download"
        return 0
    fi
    
    wget -O "${temp_file}" "${REF_GENOME_URL}" || {
        log_error "Failed to download reference genome"
        return 1
    }
    
    mv "${temp_file}" "${final_file}"
    
    gunzip -f "${final_file}" || {
        log_error "Failed to decompress reference genome"
        return 1
    }
    
    if ! validate_file_size "${REF_DIR}/${REF_GENOME_FILE}" 1000000; then
        log_error "Reference genome file appears to be incomplete"
        return 1
    }
    
    log_info "Reference genome downloaded and decompressed successfully"
}

# Download sample data
download_sample_data() {
    log_info "Downloading sample data..."
    
    local download_errors=0
    
    for sample in "${!SAMPLES[@]}"; do
        local base_url="${SAMPLES[${sample}]}"
        
        for read in 1 2; do
            local filename="${sample}_${read}.fastq.gz"
            local temp_file="${SAMPLES_DIR}/${filename}.tmp"
            local final_file="${SAMPLES_DIR}/${filename}"
            local url="${base_url}/${filename}"
            
            if [[ -f "${final_file}" ]]; then
                log_info "Sample file already exists: ${filename}"
                continue
            }
            
            log_info "Downloading ${filename}..."
            
            if ! wget -O "${temp_file}" "${url}"; then
                log_error "Failed to download ${filename}"
                ((download_errors++))
                continue
            fi
            
            mv "${temp_file}" "${final_file}"
            
            if ! validate_file_size "${final_file}" 1000000; then
                log_error "Sample file appears to be incomplete: ${filename}"
                ((download_errors++))
                continue
            }
            
            log_info "Successfully downloaded ${filename}"
        done
    done
    
    if ((download_errors > 0)); then
        log_error "${download_errors} sample files failed to download"
        return 1
    fi
    
    touch "${SAMPLES_DIR}/DOWNLOAD_COMPLETE.txt"
    log_info "All sample data downloaded successfully"
}

#############################################################################
# Quality Control Functions
#############################################################################

# Run FastQC on samples
run_fastqc() {
    log_info "Running FastQC analysis..."
    
    local fastqc_errors=0
    
    mkdir -p "${QC_DIR}/fastqc"
    
    for file in "${SAMPLES_DIR}"/*.fastq.gz; do
        if [[ -f "${file}" ]]; then
            log_info "Running FastQC on $(basename "${file}")..."
            
            if ! fastqc -o "${QC_DIR}/fastqc" -t "${THREADS}" "${file}"; then
                log_error "FastQC failed for $(basename "${file}")"
                ((fastqc_errors++))
            fi
        fi
    done
    
    if ((fastqc_errors > 0)); then
        log_error "FastQC failed for ${fastqc_errors} files"
        return 1
    fi
    
    log_info "FastQC analysis completed successfully"
}

# Run MultiQC
run_multiqc() {
    local input_dir="$1"
    local output_dir="$2"
    local title="$3"
    
    log_info "Running MultiQC for ${title}..."
    
    mkdir -p "${output_dir}"
    
    if ! multiqc -f -o "${output_dir}" "${input_dir}"; then
        log_error "MultiQC failed for ${title}"
        return 1
    fi
    
    log_info "MultiQC completed successfully for ${title}"
}

#############################################################################
# Alignment Functions
#############################################################################

# Build HISAT2 index
build_hisat2_index() {
    log_info "Building HISAT2 index..."
    
    if [[ -f "${HISAT_INDEX_DIR}/genome.1.ht2" ]]; then
        log_info "HISAT2 index already exists, skipping build"
        return 0
    fi
    
    mkdir -p "${HISAT_INDEX_DIR}"
    
    if ! hisat2-build "${REF_DIR}/${REF_GENOME_FILE}" "${HISAT_INDEX_DIR}/genome" -p "${THREADS}"; then
        log_error "HISAT2 index building failed"
        return 1
    fi
    
    log_info "HISAT2 index built successfully"
}

# Align reads with HISAT2
align_reads() {
    log_info "Aligning reads with HISAT2..."
    
    local alignment_errors=0
    
    for sample in "${!SAMPLES[@]}"; do
        log_info "Aligning sample ${sample}..."
        
        local read1="${SAMPLES_DIR}/${sample}_1.fastq.gz"
        local read2="${SAMPLES_DIR}/${sample}_2.fastq.gz"
        local sam_file="${ALIGNMENT_DIR}/alignment/${sample}.sam"
        local bam_file="${ALIGNMENT_DIR}/alignment/${sample}.sorted.bam"
        local log_file="${ALIGNMENT_DIR}/alignment/${sample}.log"
        
        # Check input files
        if ! check_file "${read1}" || ! check_file "${read2}"; then
            log_error "Input files missing for sample ${sample}"
            ((alignment_errors++))
            continue
        }
        
        # Align reads
        if ! hisat2 -x "${HISAT_INDEX_DIR}/genome" \
                    -1 "${read1}" \
                    -2 "${read2}" \
                    -S "${sam_file}" \
                    --threads "${THREADS}" \
                    2> "${log_file}"; then
            log_error "HISAT2 alignment failed for sample ${sample}"
            ((alignment_errors++))
            continue
        }
        
        # Convert SAM to BAM and sort
        if ! samtools view -bS "${sam_file}" | \
             samtools sort -@ "${THREADS}" -o "${bam_file}"; then
            log_error "SAM to BAM conversion failed for sample ${sample}"
            ((alignment_errors++))
            continue
        }
        
        # Index BAM file
        if ! samtools index "${bam_file}"; then
            log_error "BAM indexing failed for sample ${sample}"
            ((alignment_errors++))
            continue
        }
        
        # Remove SAM file to save space
        rm -f "${sam_file}"
        
        log_info "Successfully processed sample ${sample}"
    done
    
    if ((alignment_errors > 0)); then
        log_error "Alignment failed for ${alignment_errors} samples"
        return 1
    fi
    
    log_info "All samples aligned successfully"
}

#############################################################################
# Quantification Functions
#############################################################################

# Count reads with featureCounts
count_reads() {
    log_info "Counting reads with featureCounts..."
    
    local bam_files=("${ALIGNMENT_DIR}"/alignment/*.sorted.bam)
    
    if [[ ${#bam_files[@]} -eq 0 ]]; then
        log_error "No BAM files found for counting"
        return 1
    fi
    
    mkdir -p "${COUNTS_DIR}"
    
    if ! featureCounts -p \
                      -t gene \
                      -g gene_id \
                      -a "${REF_DIR}/${REF_GENOME_FILE%.fna}.gtf" \
                      -o "${COUNTS_DIR}/raw_counts.txt" \
                      -T "${THREADS}" \
                      "${bam_files[@]}"; then
        log_error "featureCounts failed"
        return 1
    fi
    
    log_info "Read counting completed successfully"
}

#############################################################################
# Main Pipeline Function
#############################################################################

run_pipeline() {
    local start_time=$(date +%s)
    
    # Print pipeline header
    print_section_header "RNA-Seq Analysis Pipeline v${VERSION}"
    log_info "Started by: ${SCRIPT_USER} on ${START_TIME}"
    
    # Initialize pipeline
    setup_logging
    create_directory_structure || exit 1
    check_requirements || exit 1
    
    # Download data
    print_section_header "=== STEP 0: Downloading data ==="
    download_reference_genome || exit 1
    download_sample_data || exit 1
    
    # Quality control
    print_section_header "=== STEP 1: Quality Control ==="
    run_fastqc || exit 1
    run_multiqc "${QC_DIR}/fastqc" "${QC_DIR}/multiqc_fastqc" "FastQC results" || exit 1
    
    # Alignment
    print_section_header "=== STEP 2: Read Alignment ==="
    build_hisat2_index || exit 1
    align_reads || exit 1
    run_multiqc "${ALIGNMENT_DIR}/alignment" "${QC_DIR}/multiqc_alignment" "Alignment results" || exit 1
    
    # Quantification
    print_section_header "=== STEP 3: Gene Quantification ==="
    count_reads || exit 1
    
    # Generate final report
    print_section_header "=== STEP 4: Generating Final Report ==="
    run_multiqc "${BASE_DIR}" "${MASTER_QC_DIR}" "Complete pipeline results" || exit 1
    
    # Calculate execution time
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    # Print completion message
    print_section_header "RNA-Seq PIPELINE COMPLETE"
    log_info "Pipeline execution completed in $(date -u -d @${duration} +'%H:%M:%S')"
    log_info "Results can be found in:"
    log_info "- QC Results: ${QC_DIR}"
    log_info "- Alignment Results: ${ALIGNMENT_DIR}/alignment"
    log_info "- Counts and Expression: ${COUNTS_DIR}"
    log_info "- Master QC Report: ${MASTER_QC_DIR}/multiqc_report.html"
    log_info "Pipeline run by: ${SCRIPT_USER}"
    log_info "Pipeline version: ${VERSION}"
    
    return 0
}

#############################################################################
# Script Entry Point
#############################################################################

# Run the pipeline
run_pipeline

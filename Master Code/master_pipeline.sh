#!/bin/bash
# ====================================================================================
# BIFS619_GROUP_01 RNA-Seq Analysis Master Pipeline
# Author: Group 01 (Master code compiled by Tyler Maire)
# Updated: improved reference download & validation; better FASTQ detection diagnostics
# ====================================================================================

set -euo pipefail

# Fix for Qt/matplotlib issues when running headless
export QT_QPA_PLATFORM=offscreen

# Progress tracking
FASTQC_COMPLETE=false
CLEANING_COMPLETE=false
ALIGNMENT_COMPLETE=false
COUNTS_GENERATED=false
HEATMAP_GENERATED=false

# ====================================================================================
# USAGE
#   ./master_pipeline.sh <project_directory>
# The script assumes your repository (including 00_rawdata/fastq_data/samples) is cloned
# into <project_directory>. It will detect FASTQs there and will download + validate
# the reference FASTA from NCBI when required.
# ====================================================================================

if [ $# -ne 1 ]; then
  echo "Error: Please provide a project directory path"
  echo "Usage: $0 <project_directory>"
  exit 1
fi

PROJECT_DIR="$1"
mkdir -p "$PROJECT_DIR"
echo "Project will be created / used at: $PROJECT_DIR"
BASE_DIR="$PROJECT_DIR"

# Threading + key paths
THREADS="${THREADS:-8}"
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

RUN_R_HEATMAP="${RUN_R_HEATMAP:-true}"
SCRIPT_VERSION="3.4-local-fixed-ref"
RUN_DATE=$(date +"%Y-%m-%d %H:%M:%S")
RUN_USER="${USER:-$(whoami)}"

echo "RNA-Seq Analysis Pipeline v${SCRIPT_VERSION}"
echo "Started by: ${RUN_USER} on ${RUN_DATE}"
echo "======================================================================"

# Create directories
mkdir -p \
  "${RAW_FASTQ_DIR}" "${REFERENCE_DIR}" "${FASTQC_OUT}" \
  "${CLEANED_FASTQ_DIR}" "${QC_TABLES_DIR}" "${QC_PLOTS_DIR}" \
  "${HISAT2_INDEX_DIR}" "${BAM_DIR}" "${ALIGN_LOGS_DIR}" "${ALIGN_TABLES_DIR}" "${ALIGN_PLOTS_DIR}" \
  "${BASE_DIR}/01_allignment/code" "${BASE_DIR}/02_annotation/code" "${COUNTS_DIR}" "${ANNOTATION_PLOTS_DIR}" \
  "${MASTER_QC_DIR}"

# Basic dependency checks
echo "Checking required software..."
check_dep() { command -v "$1" >/dev/null 2>&1 && echo "✓ $1 is installed" || { echo "✗ $1 missing"; return 1; }; }
missing=0
for bin in fastqc multiqc fastp hisat2 samtools featureCounts wget curl python3 gzip gunzip; do
  check_dep "$bin" || missing=1
done

if command -v Rscript >/dev/null 2>&1; then
  echo "✓ Rscript is installed"
else
  echo "✗ Rscript missing. Heatmap generation will be skipped."
  RUN_R_HEATMAP=false
fi

[ $missing -eq 1 ] && echo "Note: some tools are missing; pipeline will attempt available steps."

# ============================================================================
# STEP 0: detect FASTQ files in repo clone (and diagnostics)
# ============================================================================
echo "=== STEP 0: Preparing raw data (using local FASTQ files) ==="

SAMPLES=()
shopt -s nullglob nocaseglob

# Prefer paired read patterns
for f in "${RAW_FASTQ_DIR}"/*_1.fastq.gz "${RAW_FASTQ_DIR}"/*_R1.fastq.gz; do
  [ -e "$f" ] || continue
  base=$(basename "$f")
  sample="${base%_1.fastq.gz}"
  sample="${sample%_R1.fastq.gz}"
  sample="${sample%.fastq.gz}"
  if [[ ! " ${SAMPLES[*]} " =~ " ${sample} " ]]; then
    SAMPLES+=("$sample")
  fi
done

# If none, fall back to any *.fastq.gz (single-end)
if [ ${#SAMPLES[@]} -eq 0 ]; then
  for f in "${RAW_FASTQ_DIR}"/*.fastq.gz; do
    [ -e "$f" ] || continue
    base=$(basename "$f")
    sample="${base%.fastq.gz}"
    if [[ ! " ${SAMPLES[*]} " =~ " ${sample} " ]]; then
      SAMPLES+=("$sample")
    fi
  done
fi

shopt -u nullglob nocaseglob

# Diagnostic: list files present in RAW_FASTQ_DIR
echo "Files in ${RAW_FASTQ_DIR}:"
ls -lah "${RAW_FASTQ_DIR}" 2>/dev/null || echo "(directory empty or not accessible)"

if [ ${#SAMPLES[@]} -eq 0 ]; then
  echo "ERROR: No FASTQ files found in ${RAW_FASTQ_DIR}."
  echo "If your GitHub clone contains FASTQs, run the pipeline with the path to that clone."
  echo "Example: ./master_pipeline.sh ~/path/to/BIFS619_GROUP_01"
  exit 1
fi

echo "Detected samples:"
for s in "${SAMPLES[@]}"; do echo " - $s"; done

# We'll not attempt SRA downloads here
RUN_SRA_DOWNLOAD=false

# ============================================================================
# STEP 0.5: Download and validate reference FASTA (NCBI)
# ============================================================================
echo "=== STEP 0.5: Ensuring reference FASTA is present and valid ==="
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"

# Validate a FASTA: starts with '>' and has at least one sequence with alphabetic bases
validate_fasta() {
  local f="$1"
  if [ ! -s "$f" ]; then
    return 1
  fi
  # Check first non-empty line starts with '>'
  local first
  first=$(grep -m1 -v '^$' "$f" | head -n1 || true)
  if [[ "$first" != '>'* ]]; then
    return 2
  fi
  # Check that there is at least one sequence line with letters A/C/G/T/N (case-insensitive)
  if ! grep -m1 -E -i '^[ACGTURYKMSWBDHVN]+$' "$f" >/dev/null 2>&1; then
    # if no pure letter-only lines, still allow non-strict sequences but warn
    return 3
  fi
  return 0
}

download_reference() {
  mkdir -p "$REFERENCE_DIR"
  local gz_path="${REFERENCE_DIR}/GCF_000005845.2_genomic.fna.gz"
  echo "[INFO] Downloading reference from: $REF_URL"
  if command -v wget >/dev/null 2>&1; then
    wget -c -O "$gz_path" "$REF_URL"
  else
    curl -f -L -o "$gz_path" "$REF_URL"
  fi

  if [ ! -s "$gz_path" ]; then
    echo "[ERROR] Reference download produced no file."
    [ -f "$gz_path" ] && rm -f "$gz_path"
    return 1
  fi

  # Decompress to a temp file first and validate before replacing REFERENCE_FASTA
  local tmp_fasta="${REFERENCE_DIR}/tmp_reference.fna"
  if ! gunzip -c "$gz_path" > "$tmp_fasta" 2>/dev/null; then
    echo "[ERROR] gunzip -c failed to decompress $gz_path"
    rm -f "$tmp_fasta"
    return 1
  fi

  # Validate the decompressed file
  validate_fasta "$tmp_fasta"
  local v=$?
  if [ $v -eq 0 ]; then
    mv -f "$tmp_fasta" "$REFERENCE_FASTA"
    rm -f "$gz_path" || true
    echo "[OK] Reference downloaded and validated: $REFERENCE_FASTA"
    return 0
  elif [ $v -eq 2 ]; then
    echo "[ERROR] Decompressed file does not start with '>' — not a FASTA."
    rm -f "$tmp_fasta" "$gz_path"
    return 2
  elif [ $v -eq 3 ]; then
    echo "[WARN] Decompressed file has no strict sequence lines; moving it into place but HISAT2 may fail."
    mv -f "$tmp_fasta" "$REFERENCE_FASTA"
    rm -f "$gz_path" || true
    return 0
  else
    echo "[ERROR] Unknown validation failure."
    rm -f "$tmp_fasta" "$gz_path"
    return 1
  fi
}

# If file exists, validate and re-download only if invalid
if [ -s "$REFERENCE_FASTA" ]; then
  echo "[INFO] Found existing reference: $REFERENCE_FASTA — validating..."
  validate_fasta "$REFERENCE_FASTA"
  v=$?
  if [ $v -eq 0 ]; then
    echo "[OK] Reference FASTA appears valid."
  else
    echo "[WARN] Reference FASTA failed validation (code $v). Attempting re-download."
    if ! download_reference; then
      echo "[ERROR] Reference re-download failed. Please inspect $REFERENCE_FASTA or download manually."
    fi
  fi
else
  if ! download_reference; then
    echo "[WARN] Could not download a valid reference. Alignment will be skipped unless you provide a valid FASTA at:"
    echo "  $REFERENCE_FASTA"
  fi
fi

# ============================================================================
# STEP 1: FastQC
# ============================================================================
echo "=== STEP 1: Running FastQC ==="
mkdir -p "$FASTQC_OUT"

for sample in "${SAMPLES[@]}"; do
  f1="${RAW_FASTQ_DIR}/${sample}_1.fastq.gz"
  f2="${RAW_FASTQ_DIR}/${sample}_2.fastq.gz"
  single="${RAW_FASTQ_DIR}/${sample}.fastq.gz"

  if [ -f "$f1" ] && [ -f "$f2" ]; then
    echo "Running FastQC on sample ${sample} (paired)..."
    fastqc "$f1" "$f2" -o "$FASTQC_OUT" || echo "[WARN] FastQC failed for ${sample}"
  elif [ -f "$single" ]; then
    echo "Running FastQC on sample ${sample} (single-end)..."
    fastqc "$single" -o "$FASTQC_OUT" || echo "[WARN] FastQC failed for ${sample}"
  else
    echo "WARNING: FASTQ files for ${sample} not found at expected paths:"
    echo "  - ${f1}"
    echo "  - ${f2}"
    echo "  - ${single}"
    echo "  (Check that you passed the path to the cloned repository containing 00_rawdata/fastq_data/samples)"
  fi
done

echo "Running MultiQC to summarize FastQC results..."
multiqc "$FASTQC_OUT" -o "$FASTQC_OUT" || true
if [ -f "${FASTQC_OUT}/multiqc_report.html" ]; then
  FASTQC_COMPLETE=true
fi

# ============================================================================
# STEP 2: fastp cleaning
# ============================================================================
echo "=== STEP 2: Cleaning reads with fastp ==="
mkdir -p "$CLEANED_FASTQ_DIR"

for sample in "${SAMPLES[@]}"; do
  r1="${RAW_FASTQ_DIR}/${sample}_1.fastq.gz"
  r2="${RAW_FASTQ_DIR}/${sample}_2.fastq.gz"
  single="${RAW_FASTQ_DIR}/${sample}.fastq.gz"

  if [ -f "$r1" ] && [ -f "$r2" ]; then
    echo "Cleaning paired reads for ${sample}..."
    fastp -i "$r1" -I "$r2" \
      -o "${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz" \
      -O "${CLEANED_FASTQ_DIR}/${sample}_2.clean.fastq.gz" \
      --html "${CLEANED_FASTQ_DIR}/${sample}_fastp.html" \
      --json "${CLEANED_FASTQ_DIR}/${sample}_fastp.json" \
      --thread "$THREADS" || echo "[WARN] fastp failed for ${sample}"
  elif [ -f "$single" ]; then
    echo "Cleaning single-end reads for ${sample}..."
    fastp -i "$single" -o "${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz" \
      --html "${CLEANED_FASTQ_DIR}/${sample}_fastp.html" \
      --json "${CLEANED_FASTQ_DIR}/${sample}_fastp.json" \
      --thread "$THREADS" || echo "[WARN] fastp failed for ${sample}"
  else
    echo "WARNING: No input FASTQ found for ${sample}; skipping cleaning."
  fi
done

# Check cleaning results robustly (default 0)
cleaned_count=$(ls -1 ${CLEANED_FASTQ_DIR}/*.clean.fastq.gz 2>/dev/null | wc -l || echo 0)
if ! [[ "$cleaned_count" =~ ^[0-9]+$ ]]; then cleaned_count=0; fi

if [ "$cleaned_count" -ge 1 ]; then
  # ensure each sample has at least one cleaned file
  ok=0
  for sample in "${SAMPLES[@]}"; do
    if [ -f "${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz" ] || [ -f "${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz" ]; then
      ok=$((ok+1))
    fi
  done
  if [ "$ok" -eq "${#SAMPLES[@]}" ]; then
    CLEANING_COMPLETE=true
    echo "[OK] Cleaning completed for all detected samples."
  else
    echo "WARNING: Cleaning incomplete for some samples (${ok}/${#SAMPLES[@]})."
    CLEANING_COMPLETE=false
  fi
else
  echo "WARNING: No cleaned FASTQ files found after fastp."
  CLEANING_COMPLETE=false
fi

# ============================================================================
# STEP 3: Summarize QC (only if cleaning completed)
# ============================================================================
echo "=== STEP 3: Summarizing QC results ==="
if [ "$CLEANING_COMPLETE" = true ]; then
  cat << 'PY' > "${BASE_DIR}/01_allignment/code/summarize_qc.py"
#!/usr/bin/env python3
import os, json, sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

qc_dir = sys.argv[1]
table_out = sys.argv[2]
plot_out = sys.argv[3]
samples = sys.argv[4].split(',')

os.makedirs(os.path.dirname(table_out), exist_ok=True)
os.makedirs(os.path.dirname(plot_out), exist_ok=True)
summary=[]
for sample in samples:
    jf = os.path.join(qc_dir, f"{sample}_fastp.json")
    try:
        with open(jf) as fh:
            j=json.load(fh)
        before = j['summary']['before_filtering']['total_reads']
        after = j['summary']['after_filtering']['total_reads']
        summary.append({'Sample':sample,'Raw Reads':before,'Cleaned Reads':after})
    except Exception as e:
        print("WARNING:", jf, e)

if summary:
    df=pd.DataFrame(summary)
    df.to_csv(table_out,index=False)
    plt.figure(figsize=(8,4))
    x=range(len(df))
    plt.bar(x, df['Raw Reads'], label='Raw')
    plt.bar(x, df['Cleaned Reads'], label='Cleaned')
    plt.xticks(x, df['Sample'])
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_out,dpi=150)
else:
    print("No QC summary data.")
PY

  python3 "${BASE_DIR}/01_allignment/code/summarize_qc.py" \
    "$CLEANED_FASTQ_DIR" "${QC_TABLES_DIR}/pre_post_cleaning_table.csv" "${QC_PLOTS_DIR}/pre_post_cleaning_barplot.png" "$(IFS=,; echo "${SAMPLES[*]}")" || true
else
  echo "Skipping QC summary because cleaning didn't complete for all samples."
fi

# ============================================================================
# STEP 4: HISAT2 index build and alignment (only if reference valid & cleaning done)
# ============================================================================
echo "=== STEP 4: Building HISAT2 index and aligning reads ==="

# validate reference before building index
if [ -s "$REFERENCE_FASTA" ]; then
  validate_fasta "$REFERENCE_FASTA"
  v=$?
  if [ $v -ne 0 ]; then
    echo "[ERROR] Reference FASTA at $REFERENCE_FASTA failed validation (code $v). Alignment will be skipped."
    ALIGNMENT_COMPLETE=false
  else
    # Build index if missing and cleaning was done
    if [ -f "${HISAT2_INDEX_DIR}/genome.1.ht2" ]; then
      echo "[INFO] HISAT2 index already exists."
    else
      echo "Building HISAT2 index..."
      hisat2-build "$REFERENCE_FASTA" "${HISAT2_INDEX_DIR}/genome" 2>&1 | tee "${HISAT2_INDEX_DIR}/hisat2-build.log" || true
    fi

    if [ -f "${HISAT2_INDEX_DIR}/genome.1.ht2" ] && [ "$CLEANING_COMPLETE" = true ]; then
      for sample in "${SAMPLES[@]}"; do
        in1="${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz"
        in2="${CLEANED_FASTQ_DIR}/${sample}_2.clean.fastq.gz"
        single="${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz"

        if [ -f "$in1" ] && [ -f "$in2" ]; then
          echo "Aligning paired sample ${sample}..."
          hisat2 -x "${HISAT2_INDEX_DIR}/genome" -1 "$in1" -2 "$in2" -S "${ALIGN_LOGS_DIR}/${sample}.sam" --summary-file "${ALIGN_LOGS_DIR}/${sample}_summary.txt" --threads "$THREADS" 2> "${ALIGN_LOGS_DIR}/${sample}_hisat2.log" || true
        elif [ -f "$single" ]; then
          echo "Aligning single-end sample ${sample}..."
          hisat2 -x "${HISAT2_INDEX_DIR}/genome" -U "$single" -S "${ALIGN_LOGS_DIR}/${sample}.sam" --summary-file "${ALIGN_LOGS_DIR}/${sample}_summary.txt" --threads "$THREADS" 2> "${ALIGN_LOGS_DIR}/${sample}_hisat2.log" || true
        else
          echo "WARNING: No cleaned FASTQ for ${sample}; skipping."
          continue
        fi

        if [ -f "${ALIGN_LOGS_DIR}/${sample}.sam" ]; then
          samtools view -bS "${ALIGN_LOGS_DIR}/${sample}.sam" > "${BAM_DIR}/${sample}.bam" || true
          samtools sort -@ "$THREADS" -o "${BAM_DIR}/${sample}.sorted.bam" "${BAM_DIR}/${sample}.bam" || true
          samtools index "${BAM_DIR}/${sample}.sorted.bam" || true
          mv -f "${BAM_DIR}/${sample}.sorted.bam" "${BAM_DIR}/${sample}.bam" || true
          [ -f "${BAM_DIR}/${sample}.sorted.bam.bai" ] && mv -f "${BAM_DIR}/${sample}.sorted.bam.bai" "${BAM_DIR}/${sample}.bam.bai" || true
          rm -f "${ALIGN_LOGS_DIR}/${sample}.sam" || true
        fi
      done

      aligned_count=$(ls -1 ${BAM_DIR}/*.bam 2>/dev/null | wc -l || echo 0)
      if ! [[ "$aligned_count" =~ ^[0-9]+$ ]]; then aligned_count=0; fi
      if [ "$aligned_count" -ge 1 ]; then
        ALIGNMENT_COMPLETE=true
      else
        ALIGNMENT_COMPLETE=false
      fi
    else
      echo "Skipping alignment (index missing or cleaning incomplete)."
      ALIGNMENT_COMPLETE=false
    fi
  fi
else
  echo "Reference FASTA missing; skipping alignment."
  ALIGNMENT_COMPLETE=false
fi

# ============================================================================
# Remaining steps (summaries, featureCounts, heatmap, MultiQC) unchanged except for guards
# ============================================================================
# ... (rest of pipeline: alignment summary, featureCounts, heatmap, MultiQC)
# For brevity, the remainder of the pipeline remains the same as before and will
# run conditionally based on ALIGNMENT_COMPLETE and COUNTS_GENERATED flags.
# You already have those steps in your script; keep them unchanged below this point.
#
# (If you want, I can paste the full remainder here with the same guards to be fully self-contained.)
#
echo "Pipeline script has updated reference download/validation and FASTQ diagnostics."
echo "If HISAT2 still reports 'Reference file does not seem to be a FASTA file', run the following checks:"
echo "  1) Inspect the first 20 lines of the FASTA: head -n 20 ${REFERENCE_FASTA}"
echo "  2) Check file type: file ${REFERENCE_FASTA}"
echo "  3) Confirm it starts with '>' and contains sequence characters: grep -m1 -v '^$' ${REFERENCE_FASTA} | head -n1"
echo ""
echo "If you want, I can replace the remainder of the pipeline here as well (fully expanded)."

#!/bin/bash
# ====================================================================================
# BIFS619_GROUP_01 RNA-Seq Analysis Master Pipeline
# Author: Group 01 (Master code compiled by Tyler Maire)
# Updated: force-download correct reference FASTA from NCBI GCF_000005845.2_ASM584v2
# ====================================================================================

set -euo pipefail

export QT_QPA_PLATFORM=offscreen

FASTQC_COMPLETE=false
CLEANING_COMPLETE=false
ALIGNMENT_COMPLETE=false
COUNTS_GENERATED=false
HEATMAP_GENERATED=false

if [ $# -ne 1 ]; then
  echo "Error: Please provide a project directory path"
  echo "Usage: $0 <project_directory>"
  exit 1
fi

PROJECT_DIR="$1"
mkdir -p "$PROJECT_DIR"
BASE_DIR="$PROJECT_DIR"
echo "Project will be created / used at: $BASE_DIR"

THREADS="${THREADS:-8}"
RAW_FASTQ_DIR="${BASE_DIR}/00_rawdata/fastq_data/samples"
REFERENCE_DIR="${BASE_DIR}/00_rawdata/fastq_data/reference"
REFERENCE_FASTA="${REFERENCE_DIR}/GCF_000005845.2.fna"
FASTQC_OUT="${BASE_DIR}/00_rawdata/fastQC"
CLEANED_FASTQ_DIR="${BASE_DIR}/01_allignment/QC/cleaned_fastq"
HISAT2_INDEX_DIR="${BASE_DIR}/01_allignment/alignment/hisat2_index"
BAM_DIR="${BASE_DIR}/01_allignment/alignment/bam"
ALIGN_LOGS_DIR="${BASE_DIR}/01_allignment/alignment/logs"
COUNTS_DIR="${BASE_DIR}/02_annotation/counts"
MASTER_QC_DIR="${BASE_DIR}/master_qc_report"

SCRIPT_VERSION="3.4-local-reffix"
RUN_DATE=$(date +"%Y-%m-%d %H:%M:%S")
RUN_USER="${USER:-$(whoami)}"

echo "RNA-Seq Analysis Pipeline v${SCRIPT_VERSION}"
echo "Started by: ${RUN_USER} on ${RUN_DATE}"
echo "======================================================================"

# create dirs
mkdir -p "${RAW_FASTQ_DIR}" "${REFERENCE_DIR}" "${FASTQC_OUT}" "${CLEANED_FASTQ_DIR}" "${HISAT2_INDEX_DIR}" "${BAM_DIR}" "${ALIGN_LOGS_DIR}" "${COUNTS_DIR}" "${MASTER_QC_DIR}"

# quick dep check
echo "Checking required software..."
which wget >/dev/null 2>&1 && echo "✓ wget" || echo "✗ wget (required to download reference)"
which gunzip >/dev/null 2>&1 && echo "✓ gunzip" || echo "✗ gunzip"
which hisat2-build >/dev/null 2>&1 && echo "✓ hisat2-build" || echo "✗ hisat2-build (needed for index)"

# ========================================================================
# STEP 0: detect FASTQ files (from cloned repo)
# ========================================================================
echo "=== STEP 0: Detecting local FASTQ files ==="
shopt -s nullglob nocaseglob
SAMPLES=()
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

echo "Files in ${RAW_FASTQ_DIR}:"
ls -lah "${RAW_FASTQ_DIR}" 2>/dev/null || true

if [ ${#SAMPLES[@]} -eq 0 ]; then
  echo "ERROR: No FASTQ files detected in ${RAW_FASTQ_DIR}. Aborting."
  exit 1
fi
echo "Detected samples:"
for s in "${SAMPLES[@]}"; do echo " - $s"; done

# ========================================================================
# STEP 0.5: Force-download correct reference FASTA from NCBI (user-provided URL)
# ========================================================================
echo "=== STEP 0.5: Downloading reference FASTA from NCBI (GCF_000005845.2_ASM584v2) ==="

# Exact file to fetch from the directory:
REF_GZ_NAME="GCF_000005845.2_ASM584v2_genomic.fna.gz"
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/${REF_GZ_NAME}"
GZ_PATH="${REFERENCE_DIR}/${REF_GZ_NAME}"
TMP_FNA="${REFERENCE_DIR}/.tmp_reference_decompressed.fna"

# Remove old invalid REFERENCE_FASTA if present (back it up first)
if [ -f "$REFERENCE_FASTA" ]; then
  echo "[INFO] Backing up existing reference FASTA to ${REFERENCE_FASTA}.bak"
  mv -f "$REFERENCE_FASTA" "${REFERENCE_FASTA}.bak" || true
fi

echo "[INFO] Downloading: $REF_URL"
if wget -c -O "$GZ_PATH" "$REF_URL"; then
  echo "[INFO] Download finished: $GZ_PATH"
else
  echo "[ERROR] wget failed to download reference from NCBI. Check network/connectivity and the URL:"
  echo "  $REF_URL"
  [ -f "$GZ_PATH" ] && rm -f "$GZ_PATH"
  exit 1
fi

echo "[INFO] Decompressing downloaded file to temporary location..."
if gunzip -c "$GZ_PATH" > "$TMP_FNA"; then
  echo "[INFO] Decompression succeeded"
else
  echo "[ERROR] Decompression failed (gunzip -c). Inspect $GZ_PATH"
  rm -f "$TMP_FNA" "$GZ_PATH"
  exit 1
fi

# Validate the decompressed file looks like FASTA
first_line=$(grep -m1 -v '^$' "$TMP_FNA" | head -n1 || true)
if [[ "$first_line" != '>'* ]]; then
  echo "[ERROR] Decompressed file does not start with '>'; not a valid FASTA. First non-empty line:"
  echo "-----"
  echo "$first_line"
  echo "-----"
  rm -f "$TMP_FNA" "$GZ_PATH"
  # restore backup if it exists
  if [ -f "${REFERENCE_FASTA}.bak" ]; then
    mv -f "${REFERENCE_FASTA}.bak" "$REFERENCE_FASTA"
    echo "[INFO] Restored previous reference FASTA from backup."
  fi
  exit 1
fi

# Quick sequence content check: ensure there is at least one line with base letters
if ! grep -m1 -E -i '^[ACGTURYKMSWBDHVN]+$' "$TMP_FNA" >/dev/null 2>&1; then
  echo "[WARN] FASTA does not contain plain sequence lines matching A/C/G/T/N patterns. HISAT2 may fail."
fi

# Move validated file into final location
mv -f "$TMP_FNA" "$REFERENCE_FASTA"
rm -f "$GZ_PATH" || true
chmod 644 "$REFERENCE_FASTA"
echo "[OK] Reference FASTA downloaded and placed at: $REFERENCE_FASTA"

# ========================================================================
# Continue with pipeline (FastQC, fastp, HISAT2, etc.)
# For brevity we keep the rest minimal; previous code will run same steps.
# ========================================================================
echo "Proceeding with FastQC/fastp/hisat2 steps..."

# STEP 1: FastQC
mkdir -p "$FASTQC_OUT"
for sample in "${SAMPLES[@]}"; do
  f1="${RAW_FASTQ_DIR}/${sample}_1.fastq.gz"
  f2="${RAW_FASTQ_DIR}/${sample}_2.fastq.gz"
  single="${RAW_FASTQ_DIR}/${sample}.fastq.gz"
  if [ -f "$f1" ] && [ -f "$f2" ]; then
    fastqc "$f1" "$f2" -o "$FASTQC_OUT" || true
  elif [ -f "$single" ]; then
    fastqc "$single" -o "$FASTQC_OUT" || true
  else
    echo "WARNING: FASTQ not found for sample ${sample}"
  fi
done
multiqc "$FASTQC_OUT" -o "$FASTQC_OUT" || true
FASTQC_COMPLETE=true

# STEP 2: fastp cleaning
mkdir -p "$CLEANED_FASTQ_DIR"
for sample in "${SAMPLES[@]}"; do
  r1="${RAW_FASTQ_DIR}/${sample}_1.fastq.gz"
  r2="${RAW_FASTQ_DIR}/${sample}_2.fastq.gz"
  single="${RAW_FASTQ_DIR}/${sample}.fastq.gz"
  if [ -f "$r1" ] && [ -f "$r2" ]; then
    fastp -i "$r1" -I "$r2" -o "${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz" -O "${CLEANED_FASTQ_DIR}/${sample}_2.clean.fastq.gz" --thread "$THREADS" || true
  elif [ -f "$single" ]; then
    fastp -i "$single" -o "${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz" --thread "$THREADS" || true
  fi
done

cleaned_count=$(ls -1 ${CLEANED_FASTQ_DIR}/*.clean.fastq.gz 2>/dev/null | wc -l || echo 0)
if ! [[ "$cleaned_count" =~ ^[0-9]+$ ]]; then cleaned_count=0; fi
if [ "$cleaned_count" -ge 1 ]; then
  CLEANING_COMPLETE=true
else
  CLEANING_COMPLETE=false
fi

# STEP 4: HISAT2 build + alignment (only if reference validated)
if [ "$CLEANING_COMPLETE" = true ] && [ -s "$REFERENCE_FASTA" ]; then
  echo "Building HISAT2 index..."
  hisat2-build "$REFERENCE_FASTA" "${HISAT2_INDEX_DIR}/genome" 2>&1 | tee "${HISAT2_INDEX_DIR}/hisat2-build.log" || true

  if [ -f "${HISAT2_INDEX_DIR}/genome.1.ht2" ]; then
    echo "Index built. Aligning reads..."
    mkdir -p "$BAM_DIR" "$ALIGN_LOGS_DIR"
    for sample in "${SAMPLES[@]}"; do
      in1="${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz"
      in2="${CLEANED_FASTQ_DIR}/${sample}_2.clean.fastq.gz"
      single="${CLEANED_FASTQ_DIR}/${sample}_1.clean.fastq.gz"
      if [ -f "$in1" ] && [ -f "$in2" ]; then
        hisat2 -x "${HISAT2_INDEX_DIR}/genome" -1 "$in1" -2 "$in2" -S "${ALIGN_LOGS_DIR}/${sample}.sam" --threads "$THREADS" 2> "${ALIGN_LOGS_DIR}/${sample}_hisat2.log" || true
      elif [ -f "$single" ]; then
        hisat2 -x "${HISAT2_INDEX_DIR}/genome" -U "$single" -S "${ALIGN_LOGS_DIR}/${sample}.sam" --threads "$THREADS" 2> "${ALIGN_LOGS_DIR}/${sample}_hisat2.log" || true
      fi

      if [ -f "${ALIGN_LOGS_DIR}/${sample}.sam" ]; then
        samtools view -bS "${ALIGN_LOGS_DIR}/${sample}.sam" > "${BAM_DIR}/${sample}.bam" || true
        samtools sort -@ "$THREADS" -o "${BAM_DIR}/${sample}.sorted.bam" "${BAM_DIR}/${sample}.bam" || true
        samtools index "${BAM_DIR}/${sample}.sorted.bam" || true
        mv -f "${BAM_DIR}/${sample}.sorted.bam" "${BAM_DIR}/${sample}.bam" || true
      fi
    done

    aligned_count=$(ls -1 ${BAM_DIR}/*.bam 2>/dev/null | wc -l || echo 0)
    if ! [[ "$aligned_count" =~ ^[0-9]+$ ]]; then aligned_count=0; fi
    if [ "$aligned_count" -ge 1 ]; then
      ALIGNMENT_COMPLETE=true
    fi
  else
    echo "[ERROR] HISAT2 index build failed. See ${HISAT2_INDEX_DIR}/hisat2-build.log"
  fi
else
  echo "Skipping alignment: either cleaning incomplete or reference missing/invalid."
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

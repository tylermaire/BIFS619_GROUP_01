#!/bin/bash
# ====================================================================================
# BIFS619_GROUP_01 RNA-Seq Analysis Master Pipeline
# Author: Group 01 (Master code compiled by Tyler Maire)
# Updated: Use already-cloned FASTQ files and download the reference FASTA from NCBI
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
# into <project_directory>. It will:
#  - detect FASTQ files in 00_rawdata/fastq_data/samples and operate on them
#  - download the reference FASTA from NCBI (URL below) if not present
# ====================================================================================

# --- Required argument ---
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
# canonical name used in downstream code
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
SCRIPT_VERSION="3.4-local"
RUN_DATE=$(date +"%Y-%m-%d %H:%M:%S")
RUN_USER="${USER:-$(whoami)}"

echo "RNA-Seq Analysis Pipeline v${SCRIPT_VERSION}"
echo "Started by: ${RUN_USER} on ${RUN_DATE}"
echo "======================================================================"

# ====================================================================================
# CREATE DIRECTORY STRUCTURE
# ====================================================================================
echo "Creating directory structure..."
mkdir -p \
  "${RAW_FASTQ_DIR}" "${REFERENCE_DIR}" "${FASTQC_OUT}" \
  "${CLEANED_FASTQ_DIR}" "${QC_TABLES_DIR}" "${QC_PLOTS_DIR}" \
  "${HISAT2_INDEX_DIR}" "${BAM_DIR}" "${ALIGN_LOGS_DIR}" "${ALIGN_TABLES_DIR}" "${ALIGN_PLOTS_DIR}" \
  "${BASE_DIR}/01_allignment/code" "${BASE_DIR}/02_annotation/code" "${COUNTS_DIR}" "${ANNOTATION_PLOTS_DIR}" \
  "${MASTER_QC_DIR}"

# ====================================================================================
# DEPENDENCY CHECKS (minimal)
# ====================================================================================
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

[ $missing -eq 1 ] && echo "Note: some tools are missing; script will attempt the steps that can run."

# ====================================================================================
# STEP 0 — USE LOCAL FASTQ FILES (from cloned GitHub repo)
# ====================================================================================
echo "=== STEP 0: Preparing raw data (using local FASTQ files) ==="

# Detect samples automatically from filenames in RAW_FASTQ_DIR
SAMPLES=()
shopt -s nullglob nocaseglob
# Look for paired-end patterns first
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

# If none found, look for single-end *.fastq.gz
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

if [ ${#SAMPLES[@]} -eq 0 ]; then
  echo "ERROR: No FASTQ files found in ${RAW_FASTQ_DIR}."
  echo "Ensure the repository has files under: ${RAW_FASTQ_DIR}"
  exit 1
fi

echo "Detected samples:"
for s in "${SAMPLES[@]}"; do echo " - $s"; done

# We won't attempt SRA downloads; working from local clone.
RUN_SRA_DOWNLOAD=false

# ====================================================================================
# STEP 0.5 — DOWNLOAD REFERENCE FASTA (NCBI)
# ====================================================================================
echo "=== STEP 0.5: Ensuring reference FASTA is present ==="
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"

download_reference() {
  mkdir -p "$REFERENCE_DIR"
  local gz_path="${REFERENCE_DIR}/GCF_000005845.2_genomic.fna.gz"
  if [ -s "$REFERENCE_FASTA" ]; then
    echo "[OK] Reference FASTA already exists: $REFERENCE_FASTA"
    return 0
  fi

  echo "[INFO] Downloading reference FASTA from NCBI..."
  if command -v wget >/dev/null 2>&1; then
    wget -c -O "$gz_path" "$REF_URL"
    rc=$?
  else
    curl -f -L -o "$gz_path" "$REF_URL"
    rc=$?
  fi

  if [ "$rc" -ne 0 ] || [ ! -s "$gz_path" ]; then
    echo "ERROR: Reference download failed (wget/curl returned non-zero or empty file)."
    [ -f "$gz_path" ] && rm -f "$gz_path"
    return 1
  fi

  echo "[INFO] Decompressing reference FASTA to $REFERENCE_FASTA"
  if gunzip -c "$gz_path" > "$REFERENCE_FASTA"; then
    rm -f "$gz_path"
  else
    echo "ERROR: Decompression failed. Keeping gz file at $gz_path for inspection."
    return 1
  fi

  # Basic validation: file must start with '>'
  if head -n 1 "$REFERENCE_FASTA" | grep -q '^>'; then
    echo "[OK] Reference FASTA downloaded and validated: $REFERENCE_FASTA"
    return 0
  else
    echo "WARNING: Reference FASTA does not appear to be a FASTA (first line:)"
    head -n 3 "$REFERENCE_FASTA"
    return 1
  fi
}

# Try to download reference; if it fails, we'll continue but alignment will be skipped later
if ! download_reference; then
  echo "[WARN] Reference FASTA not available; alignment step will be skipped unless you provide $REFERENCE_FASTA"
fi

# ====================================================================================
# STEP 1: QUALITY CONTROL WITH FASTQC
# ====================================================================================
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
    echo "WARNING: FASTQ files for ${sample} not found. Skipping FastQC for this sample."
  fi
done

echo "Running MultiQC to summarize FastQC results..."
multiqc "$FASTQC_OUT" -o "$FASTQC_OUT" || true

if [ -f "${FASTQC_OUT}/multiqc_report.html" ]; then
  echo "[OK] MultiQC created at ${FASTQC_OUT}/multiqc_report.html"
  FASTQC_COMPLETE=true
else
  echo "WARNING: MultiQC did not produce a report."
  FASTQC_COMPLETE=false
fi

# ====================================================================================
# STEP 2: READ CLEANING WITH FASTP
# ====================================================================================
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

# Check cleaning results
cleaned_count=$(ls -1 ${CLEANED_FASTQ_DIR}/*.clean.fastq.gz 2>/dev/null | wc -l || echo 0)
if [ "$cleaned_count" -ge 1 ]; then
  # consider cleaning complete if every sample produced at least one cleaned file
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
    CLEANING_COMPLETE=false
    echo "WARNING: Cleaning incomplete for some samples (${ok}/${#SAMPLES[@]})."
  fi
else
  CLEANING_COMPLETE=false
  echo "WARNING: No cleaned FASTQ files found after fastp."
fi

# ====================================================================================
# STEP 3: Summarize QC results (Python)
# ====================================================================================
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

# ====================================================================================
# STEP 4: ALIGNMENT WITH HISAT2 (only if reference present and cleaning complete)
# ====================================================================================
echo "=== STEP 4: Building HISAT2 index and aligning reads ==="

if [ ! -s "$REFERENCE_FASTA" ]; then
  echo "Reference FASTA missing or empty; skipping alignment."
  ALIGNMENT_COMPLETE=false
else
  # Build HISAT2 index if needed
  if [ ! -f "${HISAT2_INDEX_DIR}/genome.1.ht2" ]; then
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
    if [ "$aligned_count" -ge 1 ]; then
      ALIGNMENT_COMPLETE=true
    else
      ALIGNMENT_COMPLETE=false
    fi
  else
    echo "Skipping alignment (HISAT2 index missing or cleaning incomplete)."
    ALIGNMENT_COMPLETE=false
  fi
fi

# ====================================================================================
# STEP 5: Summarize alignment results
# ====================================================================================
echo "=== STEP 5: Summarizing alignment results ==="
if [ "$ALIGNMENT_COMPLETE" = true ]; then
  cat << 'PY' > "${BASE_DIR}/01_allignment/code/summarize_alignment.py"
#!/usr/bin/env python3
import os,sys,re,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logs_dir = sys.argv[1]
table_out = sys.argv[2]
plot_out = sys.argv[3]
samples = sys.argv[4].split(',')

os.makedirs(os.path.dirname(table_out), exist_ok=True)
os.makedirs(os.path.dirname(plot_out), exist_ok=True)

rows=[]
for s in samples:
    sf=os.path.join(logs_dir, f"{s}_summary.txt")
    try:
        t=open(sf).read()
        total=int(re.search(r'(\d+) reads; of these:', t).group(1))
        uniq=int(re.search(r'(\d+) \(\d+\.\d+\%\) aligned concordantly exactly 1 time', t).group(1))
        multi=int(re.search(r'(\d+) \(\d+\.\d+\%\) aligned concordantly >1 times', t).group(1))
        rows.append({'Sample':s,'Total':total,'Unique':uniq,'Multiple':multi,'Mapped':uniq+multi})
    except Exception as e:
        print("WARN", s, e)

if rows:
    df=pd.DataFrame(rows)
    df.to_csv(table_out,index=False)
    plt.figure(figsize=(6,4))
    plt.bar(df['Sample'], df['Mapped']/df['Total']*100)
    plt.ylim(0,100)
    plt.ylabel('Mapping rate (%)')
    plt.savefig(plot_out,dpi=150)
PY

  python3 "${BASE_DIR}/01_allignment/code/summarize_alignment.py" \
    "$ALIGN_LOGS_DIR" "${ALIGN_TABLES_DIR}/alignment_summary.csv" "${ALIGN_PLOTS_DIR}/mapping_percent_barplot.png" "$(IFS=,; echo "${SAMPLES[*]}")" || true
else
  echo "Skipping alignment summary (no alignments)."
fi

# ====================================================================================
# STEP 6: Quantification with featureCounts
# ====================================================================================
echo "=== STEP 6: Quantifying gene expression with featureCounts ==="
if [ "$ALIGNMENT_COMPLETE" = true ]; then
  featureCounts -a "$REFERENCE_GTF" -o "${COUNTS_DIR}/raw_counts.txt" -t exon -g gene_id -p -T "$THREADS" ${BAM_DIR}/*.bam 2>&1 | tee "${COUNTS_DIR}/featureCounts.log" || true
  if [ -s "${COUNTS_DIR}/raw_counts.txt" ]; then
    COUNTS_GENERATED=true
  else
    COUNTS_GENERATED=false
  fi
else
  echo "Skipping featureCounts (no alignments)."
  COUNTS_GENERATED=false
fi

# ====================================================================================
# STEP 7: Generate heatmap with R
# ====================================================================================
echo "=== STEP 7: Generating expression heatmap ==="
if [ "$RUN_R_HEATMAP" = true ] && [ "$COUNTS_GENERATED" = true ]; then
  cat << 'R' > "${BASE_DIR}/02_annotation/code/generate_top10_heatmap.R"
#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
counts_file<-args[1]; heatmap_out<-args[2]; table_out<-args[3]
dat<-read.delim(counts_file, comment.char="#")
numcols<-sapply(dat, is.numeric)
if(sum(numcols)==0) { numcols <- 2:ncol(dat) }
mat<-as.matrix(dat[,numcols])
rownames(mat)<-dat$Geneid
means<-rowMeans(mat)
top<-order(means,decreasing=TRUE)[1:min(10,length(means))]
topmat<-mat[top,,drop=FALSE]
write.csv(as.data.frame(topmat), table_out)
pdf(heatmap_out, width=8, height=6)
library(pheatmap)
pheatmap(log2(topmat+1), scale="row", main="Top 10 genes")
dev.off()
R

  Rscript "${BASE_DIR}/02_annotation/code/generate_top10_heatmap.R" \
    "${COUNTS_DIR}/raw_counts.txt" "${ANNOTATION_PLOTS_DIR}/top10_genes_heatmap.pdf" "${ANNOTATION_PLOTS_DIR}/top10_genes_table.csv" || true
else
  echo "Skipping heatmap generation (R not available or counts missing)."
fi

# ====================================================================================
# STEP 8: Master QC report
# ====================================================================================
echo "=== STEP 8: Generating master QC report ==="
multiqc -f -o "$MASTER_QC_DIR" "$FASTQC_OUT" "${BASE_DIR}/01_allignment/QC" "${BASE_DIR}/01_allignment/alignment" "$COUNTS_DIR" || true

if [ -f "${MASTER_QC_DIR}/multiqc_report.html" ]; then
  echo "Master QC report: ${MASTER_QC_DIR}/multiqc_report.html"
else
  echo "WARNING: Master QC report not found."
fi

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
echo "Pipeline completed on: $(date)"
echo "Pipeline run by: ${RUN_USER}"
echo "Pipeline version: ${SCRIPT_VERSION}"
echo "======================================================================"

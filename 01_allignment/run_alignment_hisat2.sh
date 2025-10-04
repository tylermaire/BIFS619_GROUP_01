#!/usr/bin/env bash
set -euo pipefail

THREADS=10
REF="../00_rawdata/NZ_CP076404.1.fasta"
SAMPLES=(SRR9613403 SRR9613404 SRR9613405)

# Checks
[[ -f "$REF" ]] || { echo "Reference not found: $REF"; exit 1; }
for S in "${SAMPLES[@]}"; do
  [[ -f ${S}_trimmed_1.fastq.gz && -f ${S}_trimmed_2.fastq.gz ]] \
    || { echo "Missing trimmed reads for $S"; exit 1; }
done

# 1) Build HISAT2 index (outputs: hisat2_index.*)
hisat2-build -p "$THREADS" "$REF" hisat2_index

# 2) Align each sample -> sorted BAM + index + logs + flagstat
for S in "${SAMPLES[@]}"; do
  R1=${S}_trimmed_1.fastq.gz
  R2=${S}_trimmed_2.fastq.gz
  BAM=${S}.sorted.bam

  echo "[*] Aligning $S ..."
  hisat2 -p "$THREADS" -x hisat2_index -1 "$R1" -2 "$R2" \
    2> "${S}_hisat2.log" \
  | samtools sort -@ "$THREADS" -o "$BAM" -
  samtools index "$BAM"
  samtools flagstat "$BAM" > "${S}_flagstat.txt"
done

# 3) Collect metrics table (total reads, mapped reads, mapping % as number)
echo -e "sample\ttotal_reads\tmapped_reads\tmapping_percent" > alignment_summary.tsv
for S in "${SAMPLES[@]}"; do
  total=$(awk '/in total/ {print $1; exit}' "${S}_flagstat.txt")
  mapline=$(grep " mapped (" "${S}_flagstat.txt" | head -n1)
  mapped=$(awk '{print $1}' <<<"$mapline")
  pct=$(awk -F'[()%]' '{print $2}' <<<"$mapline" | tr -d '%')
  printf "%s\t%s\t%s\t%s\n" "$S" "$total" "$mapped" "$pct" >> alignment_summary.tsv
done

# 4) Plots from alignment_summary.tsv
python3 - <<'PY'
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("alignment_summary.tsv", sep="\t")
# ensure numeric types
for col in ["mapping_percent","total_reads","mapped_reads"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# Plot 1: mapping percent per sample
plt.figure()
plt.bar(df["sample"], df["mapping_percent"])
plt.ylabel("Mapping (%)")
plt.title("HISAT2 Mapping Percentage")
plt.tight_layout()
plt.savefig("plot_mapping_percent.png", dpi=150)

# Plot 2: total vs mapped reads
x = np.arange(len(df))
w = 0.4
plt.figure()
plt.bar(x - w/2, df["total_reads"],  width=w, label="Total")
plt.bar(x + w/2, df["mapped_reads"], width=w, label="Mapped")
plt.xticks(x, df["sample"])
plt.ylabel("Reads")
plt.title("Total vs Mapped Reads")
plt.legend()
plt.tight_layout()
plt.savefig("plot_reads_total_vs_mapped.png", dpi=150)
PY

echo "Done:"
echo "  - hisat2_index.*"
echo "  - *.sorted.bam + *.bai"
echo "  - *_hisat2.log, *_flagstat.txt"
echo "  - alignment_summary.tsv"
echo "  - plot_mapping_percent.png, plot_reads_total_vs_mapped.png"

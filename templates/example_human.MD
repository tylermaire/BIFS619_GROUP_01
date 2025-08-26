# Group Project Workflow: Human RNA-Seq Analysis (BIFS 619)

## Preparation
- **Data Retrieval**
  - NCBI SRA Accessions:  
    - SRRxxxxxxx1  
    - SRRxxxxxxx2  
    - SRRxxxxxxx3  
    - SRRxxxxxxx4  
    - SRRxxxxxxx5  
- **GitHub Setup**
  - Team members: __________________________  
  - Repository link: _________________________  
  - Notes on commits / collaboration: _________________________  

---

## Step 1: Quality Control (QC)
- **Tools Used:** (FastQC, MultiQC, etc.)  
- **Tasks Performed:**  
  - Generated per-sample QC reports  
  - Summarized read quality, GC content, adapter contamination, duplication  
- **Deliverables:**  
  - Insert QC plots (e.g., per base quality, adapter content)  
  - Insert Table: Raw read counts, duplication rates  
- **Interpretation:**  
  - Brief narrative of overall quality  
  - Any concerns identified in the QC step  

---

## Step 2: Read Cleaning
- **Tools Used:** (fastp OR Cutadapt/Trimmomatic)  
- **Tasks Performed:**  
  - Trimmed adapters and low-quality bases  
  - Removed short reads (<30–50 bp)  
- **Deliverables:**  
  - Insert comparison plot: Pre- vs Post-trimming quality  
  - Insert table: Raw vs Cleaned read counts  
- **Interpretation:**  
  - How much data was lost in trimming?  
  - Did trimming improve QC?  

---

## Step 3: Alignment
- **Reference Genome:** GRCh38 + GTF (GENCODE/Ensembl)  
- **Aligner Used:** (STAR, HISAT2, or other)  
- **Tasks Performed:**  
  - Aligned cleaned reads to GRCh38  
- **Deliverables:**  
  - Insert alignment summary table (total reads, % uniquely mapped, % multimapped)  
  - Insert plot (optional: read distribution, mapping quality)  
- **Interpretation:**  
  - Comment on mapping success and any issues  

---

## Step 4: Annotation & Quantification
- **Tools Used:** (featureCounts, RSEM, or Salmon)  
- **Tasks Performed:**  
  - Generated expression tables (TPM/CPM)  
  - Identified top 10–20 expressed genes  
- **Deliverables:**  
  - Insert table: Top 20 expressed genes  
  - Insert heatmap: Expression patterns of top genes  
- **Optional:** Functional enrichment results  
- **Interpretation:**  
  - What genes are most expressed?  
  - Are results consistent with tissue/cell type expected?  

---

## Step 5: Integration & Interpretation
- **Summary of Findings:**  
  - QC → Trimming → Alignment → Quantification overview  
- **Key Insights:**  
  - Which results stand out?  
  - Any biological interpretation?  
- **Limitations / Next Steps:**  
  - What could improve the analysis?  
  - What additional data would help?  

---

## Supplemental
- **GitHub Repository URL:** ___________________________  
- **Scripts Submitted:**  
  - QC script(s): __________________  
  - Trimming script(s): __________________  
  - Alignment script(s): __________________  
  - Quantification script(s): __________________  
- **Notes:**  
  - Include any difficulties, troubleshooting steps, or parameters modified

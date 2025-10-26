# BIFS619 Group 01 - RNA-Seq Analysis Pipeline

> [!IMPORTANT]
> ## 🚀 Recommended: Use the Master Pipeline Script
> 
> **For the easiest and most complete analysis experience, run the automated master pipeline:**
> 
> ```bash
> # Clone the repository
> git clone https://github.com/tylermaire/BIFS619_GROUP_01.git
> cd BIFS619_GROUP_01
> 
> # Make the script executable
> chmod +x "Master Code/master_pipeline.sh"
> 
> # Run the complete pipeline (replace path with your desired location)
> ./Master\ Code/master_pipeline.sh ~/my_rnaseq_project
> ```
> 
> **The master pipeline automatically handles:**
> - ✅ Dependency installation (via conda or apt)
> - ✅ Raw data download from NCBI SRA
> - ✅ Reference genome and annotation download
> - ✅ Quality control (FastQC/MultiQC)
> - ✅ Read cleaning (fastp)
> - ✅ Genome alignment (HISAT2)
> - ✅ Gene expression quantification (featureCounts)
> - ✅ Visualization and heatmap generation (R)
> - ✅ Master QC report generation
> 
> **📁 The individual folders and code files in this repository are provided to:**
> - Show expected outputs and results
> - Provide detailed breakdowns of each pipeline step
> - Allow manual execution of specific steps if needed
> - Serve as documentation and reference material
> 
> **Pipeline Version:** 4.2  
> **No manual setup required - the script does everything for you!**

---

## 📊 Project Overview

This repository contains a complete RNA-seq analysis pipeline for **BIFS619** developed by Group 01. The workflow processes *Escherichia coli* K-12 MG1655 RNA-seq data through quality control, read alignment, and gene expression quantification with visualization.

**Pipeline Version:** 4.3  
**Last Updated:** October 26, 2025  
**Authors:** Tyler Maire, Kenneth Roman

---

## 🧬 Dataset

- **Organism:** *Escherichia coli* K-12 MG1655
- **Reference Genome:** GCF_000005845.2
- **Samples:** 3 biological replicates
  - SRR9613403
  - SRR9613404
  - SRR9613405
- **Total Genes Quantified:** ~4,506

---

## 🎯 Two Ways to Run This Analysis

### Option 1: Automated Master Pipeline (Recommended)

See the highlighted section above for the complete automated workflow.

**Advantages:**
- ⚡ One command runs everything
- 🔧 Automatic dependency installation
- 📥 Automatic data download
- ✅ Comprehensive error checking
- 📊 Generates master QC report

### Option 2: Manual Step-by-Step Execution

Follow the individual README files in each directory:

1. **Step 1:** (https://github.com/tylermaire/BIFS619_GROUP_01/blob/main/01_allignment/README.md)
2. **Step 2:** [Gene Expression Quantification](./02_annotation/README.md)

**When to use manual execution:**
- You want to understand each step in detail
- You need to modify specific parameters
- You're troubleshooting a particular step
- You already have some intermediate files

---

## 📁 Repository Structure

```
BIFS619_GROUP_01/
├── Master Code/
│   └── master_pipeline.sh   # ⭐ AUTOMATED PIPELINE - START HERE
├── 00_rawdata/              # Raw FASTQ files and reference genome
│   └── fastq_data/
│       └── reference/
│           ├── GCF_000005845.2.fna
│           └── GCF_000005845.2.gff
├── 01_alignment/            # Read cleaning, alignment, and QC
│   ├── QC/                  # Quality control outputs
│   │   ├── cleaned_fastq/
│   │   ├── plots/
│   │   └── tables/
│   ├── alignment/           # Alignment results
│   │   ├── bam/            # BAM files (stored on Google Drive)
│   │   ├── hisat2_index/
│   │   ├── logs/
│   │   ├── plots/
│   │   └── tables/
│   ├── code/
│   └── README.md            # Manual step-by-step guide
├── 02_annotation/           # Gene expression quantification
│   ├── code/
│   │   ├── 01_generate_gene_counts.sh
│   │   ├── 02_generate_top10_heatmap.R
│   │   └── 03_generate_CPM_table.R
│   ├── counts/             # Raw and normalized counts
│   ├── plots/              # Heatmaps and visualizations
│   └── README.md            # Manual step-by-step guide
└── README.md               # This file
```

---

## 🚀 Quick Start (Automated Pipeline)

### Prerequisites

The master pipeline can automatically install dependencies, but you'll need:
- Linux environment (Ubuntu/WSL recommended)
- `conda` or `apt` package manager
- Internet connection for downloads

### Installation & Execution

```bash
# 1. Clone the repository
git clone https://github.com/tylermaire/BIFS619_GROUP_01.git
cd BIFS619_GROUP_01

# 2. Make the master script executable
chmod +x "Master Code/master_pipeline.sh"

# 3. Run the pipeline (creates project in specified directory)
./Master\ Code/master_pipeline.sh ~/my_rnaseq_project

# The script will:
# - Check and install all dependencies
# - Download raw FASTQ files from NCBI
# - Download reference genome and annotation
# - Run complete analysis workflow
# - Generate all plots and reports
```

### What the Master Pipeline Does

```
Step 0: Data Preparation
  ├── Download raw FASTQ files (SRR9613403-405)
  └── Download E. coli K-12 reference genome & GFF

Step 1: Quality Control
  ├── FastQC on raw reads
  └── MultiQC summary report

Step 2: Read Cleaning
  ├── fastp adapter trimming and quality filtering
  └── QC summary plots and tables

Step 3: Alignment
  ├── HISAT2 index building
  ├── Read alignment to reference
  └── BAM file generation and indexing

Step 4: Quantification
  ├── featureCounts gene-level quantification
  └── Raw count matrix generation

Step 5: Visualization
  ├── CPM normalization
  ├── Top 10 gene identification
  ├── Heatmap generation (regular + log2)
  └── Gene name mapping from GFF

Step 6: Master QC Report
  └── Comprehensive MultiQC report combining all steps
```

---

## 📈 Expected Results

### Read Cleaning Summary

| Sample     | Raw Reads | Cleaned Reads | Removed % |
|:-----------|:----------|:--------------|:----------|
| SRR9613403 | 28,299,586| 27,784,644    | 1.82%     |
| SRR9613404 | 22,861,028| 22,335,460    | 2.30%     |
| SRR9613405 | 20,736,074| 20,294,582    | 2.13%     |

### Alignment Summary

- **Average mapping rate:** >85% across all samples
- **Uniquely mapped reads:** >80% per sample

### Top 10 Expressed Genes

| Gene ID     | Gene Name | Mean CPM    | Function                    |
|:------------|:----------|:------------|:----------------------------|
| gene-b2911  | ssrS      | 148,490.93  | 6S RNA                      |
| gene-b3556  | cspA      | 26,698.07   | Cold shock protein A        |
| gene-b3340  | fusA      | 17,810.69   | Elongation factor G         |
| gene-b2621  | ssrA      | 14,752.66   | tmRNA                       |
| gene-b3123  | rnpB      | 12,995.40   | RNase P RNA component       |
| gene-b3320  | rplC      | 11,575.93   | 50S ribosomal protein L3    |
| gene-b3317  | rplB      | 11,281.52   | 50S ribosomal protein L2    |
| gene-b3985  | rplJ      | 11,024.63   | 50S ribosomal protein L10   |
| gene-b3295  | rpoA      | 10,434.99   | RNA polymerase alpha subunit|
| gene-b3984  | rplA      | 9,335.64    | 50S ribosomal protein L1    |

---

## 📊 Visualizations

All visualizations are automatically generated by the master pipeline:

### Pre/Post Cleaning Read Counts
<img width="720" height="360" alt="pre_post_cleaning_barplot" src="https://github.com/user-attachments/assets/5baf6924-9083-40b1-8624-c1e3f402245b" />


### Alignment Mapping Rates
<img width="720" height="360" alt="mapping_percent_barplot" src="https://github.com/user-attachments/assets/21b1ec40-032f-4e33-969b-35fc4952b8fb" />


### Top 10 Genes Expression Heatmap
<img width="600" height="500" alt="top10_genes_heatmap" src="https://github.com/user-attachments/assets/4dfc8522-dd45-424d-9195-5eafe2fa195e" />


### Log2-Scaled Expression Pattern
<img width="600" height="500" alt="top10_genes_heatmap_log2" src="https://github.com/user-attachments/assets/ed0b9185-0cb0-4a71-a6c3-ec8f17957d25" />


---

## 💾 Large File Storage

Due to GitHub's 100MB file size limit, BAM files are stored externally:

**📥 Download BAM files:**  
[Google Drive - BIFS619 BAM files](https://drive.google.com/drive/folders/1KW6iuHfbfBAplelDN6l0az1ePj9p1X8c?usp=sharing)

> **Note:** The master pipeline automatically generates BAM files, so you only need to download these if you're running manual steps or want to skip the alignment process.

**Contact for access:** t.maire@student.umgc.edu

---

## 🛠️ Troubleshooting

### Master Pipeline Issues

**"Permission denied" when running script**
```bash
chmod +x "Master Code/master_pipeline.sh"
```

**"Dependency installation failed"**
- Install conda first: `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`
- Or install manually with apt: `sudo apt install fastqc multiqc fastp hisat2 samtools subread r-base`

**"Download failed" errors**
- Check internet connection
- Try running the script again (it will skip already-downloaded files)

**"Disk space" warnings**
- Raw data requires ~15GB
- Complete analysis requires ~30GB total
- Ensure sufficient disk space before starting

### Manual Execution Issues

**"No BAM files found"**
- Run Step 1 (alignment) first, OR
- Download BAM files from Google Drive link above
- Place them in `01_alignment/alignment/bam/`

**"GFF file not found"**
- Ensure reference files are in `00_rawdata/fastq_data/reference/`
- Update paths in scripts if your directory structure differs

**"R packages not installed"**
```bash
Rscript -e 'install.packages(c("dplyr", "pheatmap"), repos="https://cloud.r-project.org")'
```

**"Sample names showing as long paths"**
- Scripts version 4.3+ handle this automatically
- Update to the latest version in `code/` folders

---

## 📚 References

- **featureCounts:** Liao Y, Smyth GK and Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7):923-930.
- **HISAT2:** Kim D, Langmead B and Salzberg SL (2015). HISAT: a fast spliced aligner with low memory requirements. *Nature Methods*, 12:357-360.
- **fastp:** Chen S, Zhou Y, Chen Y, Gu J (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17):i884-i890.
- **pheatmap:** Kolde R (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12.
- **MultiQC:** Ewels P, Magnusson M, Lundin S, Käller M (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19):3047-3048.
- **Study where data was gathered from** O'Rourke, A., et al. (2020). Mechanism-of-Action Classification of Antibiotics by Global Transcriptome Profiling. Antimicrobial agents and chemotherapy, 64(3). https://doi.org/10.1128/AAC.01207-19

---

## 👥 Contributors

- **Tyler Maire** - Pipeline development, automation, documentation, analysis
- **Kenneth Roman** - Heatmap optimization, testing, troubleshooting

---

## 📧 Contact

For questions or issues:
- 📬 Email: t.maire@student.umgc.edu
- 🐛 [Open a GitHub Issue](https://github.com/tylermaire/BIFS619_GROUP_01/issues)

---

## 📄 License

This project is for educational purposes as part of BIFS619 coursework.

---

**Last Updated:** October 26, 2025  
**Pipeline Version:** 4.3  
**Master Script Version:** 4.2

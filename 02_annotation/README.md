```
# 02_annotation

## Overview

This folder contains all scripts and outputs for **gene expression quantification and visualization** in the BIFS619_GROUP_01 RNA-seq analysis pipeline. The workflow quantifies gene expression using featureCounts and generates multiple heatmap visualizations showing the top 10 most highly expressed genes.

---

## Getting Started

### 1. Clone the Repository

Open a terminal and run:

```bash
git clone https://github.com/tylermaire/BIFS619_GROUP_01.git
cd BIFS619_GROUP_01/02_annotation

```

### 2\. Ensure 01_alignment has been run

⚠️ IMPORTANT: You must complete the alignment step before running annotation.

If you have not run the previous `01_alignment` step, please go to the `01_alignment` folder in the repository and follow the README file. Only proceed once you have successfully generated BAM files.

### 3\. Install Required Software

Make sure you have the following tools installed:

-   Bash (for shell scripts)
-   R (version 4.0+)
-   Subread (for featureCounts)
-   R packages: `dplyr`, `pheatmap`

#### Installation on Ubuntu/Linux:

bash

```
# Install R and Subread
sudo apt update
sudo apt install -y r-base r-base-dev subread

# Install required R packages
Rscript -e 'install.packages(c("dplyr", "pheatmap"), repos="https://cloud.r-project.org")'

```

#### Installation using Conda (recommended):

bash

```
conda create -n rnaseq_env -c bioconda subread r-base r-dplyr r-pheatmap
conda activate rnaseq_env

```

### 4\. Prepare Required Files

Before running the scripts, ensure you have:

✅ BAM files from alignment step (located in `../01_alignment/alignment/bam/`)\
✅ Reference GFF file (located in `../00_rawdata/fastq_data/reference/GCF_000005845.2.gff`)

* * * * *

### IMPORTANT: Update File Paths

Before running any scripts in the `code/` folder, update file paths inside each script to match your directory structure.

Edit these paths in `code/01_generate_gene_counts.sh`:

-   `REFERENCE_GFF` - Path to your GFF annotation file
-   `BAM_DIR` - Path to your sorted BAM files
-   `OUT_DIR` - Where you want output files saved

Edit these paths in `code/02_generate_top10_heatmap.R` (or use command line arguments):

-   Counts file path
-   Output paths for heatmaps and tables
-   GFF file path

* * * * *

Workflow
--------

### Step 1: Quantify Gene Expression with featureCounts

This script uses featureCounts to assign aligned reads to genes based on the GFF annotation.

bash

```
bash code/01_generate_gene_counts.sh

```

What it does:

-   Reads sorted BAM files from the alignment step
-   Uses the GFF annotation to count reads mapping to each gene
-   Generates a raw counts table with ~4,506 E. coli genes
-   Creates a summary of assignment statistics

Outputs:

-   `counts/raw_counts.txt` - Raw gene expression counts
-   `counts/raw_counts.txt.summary` - Assignment statistics

Expected output:

Code

```
✓ FeatureCounts completed successfully.
  Quantified 4506 genes across 3 samples.

```

* * * * *

### Step 2: Generate Heatmaps and Visualizations

This R script performs CPM normalization, maps gene IDs to gene names, identifies the top 10 most expressed genes, and generates three types of heatmaps.

bash

```
Rscript code/02_generate_top10_heatmap.R\
  counts/raw_counts.txt\
  plots/top10_genes_heatmap.png\
  plots/top10_genes_table.csv\
  ../00_rawdata/fastq_data/reference/GCF_000005845.2.gff

```

What it does:

1.  Extracts gene names from GFF file (e.g., `gene-b0008` → `talB`)
2.  Normalizes counts to CPM (Counts Per Million)
3.  Identifies top 10 most highly expressed genes
4.  Generates 3 heatmaps:
    -   Log2-scaled heatmap (relative expression patterns)
    -   Regular CPM heatmap (absolute expression levels)
    -   Sample correlation heatmap (sample similarity)

Outputs:

-   `plots/top10_genes_heatmap.png` - Regular CPM heatmap
-   `plots/top10_genes_heatmap_log2.png` - Log2-scaled heatmap
-   `plots/top10_genes_heatmap_correlation.png` - Sample correlation
-   `plots/top10_genes_table.csv` - Table with top 10 genes and expression values

* * * * *

Understanding the Outputs
-------------------------

### 1\. Raw Counts (`counts/raw_counts.txt`)

Tab-delimited file containing:

-   Geneid: Gene identifier (e.g., `gene-b0001`)
-   Chr: Chromosome/sequence name
-   Start, End: Gene coordinates
-   Strand: + or - strand
-   Length: Gene length in bp
-   Sample columns: Raw read counts per sample (SRR9613403, SRR9613404, SRR9613405)
-   Gene name 

### Top 10 Expressed Genes - Raw Counts

| Geneid      | Chr          | Start   | End    | Strand | Length | SRR9613403 | SRR9613404 | SRR9613405 | Gene Name |
|-------------|--------------|---------|--------|--------|--------|------------|------------|------------|-----------|
| gene-b2911  | NC_000913.3  | 3073680 | 3073863| +      | 184    | 17872      | 26152      | 34745      | ssrS      |
| gene-b3556  | NC_000913.3  | 3758241 | 3758426| -      | 186    | 3851       | 4756       | 5532       | cspA      |
| gene-b3340  | NC_000913.3  | 3508847 | 3511192| +      | 2346   | 29235      | 36433      | 29006      | fusA      |
| gene-b2621  | NC_000913.3  | 2757124 | 2757486| +      | 363    | 2336       | 2903       | 2594       | ssrA      |
| gene-b3123  | NC_000913.3  | 3274087 | 3274463| +      | 377    | 2246       | 2149       | 2453       | rnpB      |
| gene-b3320  | NC_000913.3  | 3486316 | 3486966| +      | 651    | 2005       | 2441       | 1763       | rplC      |
| gene-b3317  | NC_000913.3  | 3485473 | 3486372| +      | 900    | 1987       | 2327       | 1714       | rplB      |
| gene-b3985  | NC_000913.3  | 4168659 | 4169120| +      | 462    | 1884       | 2250       | 1719       | rplJ      |
| gene-b3295  | NC_000913.3  | 3457762 | 3458814| +      | 1053   | 1781       | 2198       | 1575       | rpoA      |
| gene-b3984  | NC_000913.3  | 4167287 | 4167958| +      | 672    | 1654       | 1875       | 1400       | rplA      |

### 2\. Top 10 Genes Table (`plots/top10_genes_table.csv`)

### Top 10 Expressed Genes

| gene_id      | gene_name | SRR9613403  | SRR9613404  | SRR9613405  | mean_cpm     |
|--------------|-----------|-------------|-------------|-------------|--------------|
| gene-b2911   | ssrS      | 108,279.83  | 135,446.53  | 201,746.44  | 148,490.93   |
| gene-b3556   | cspA      | 23,334.29   | 24,637.45   | 32,122.47   | 26,698.07    |
| gene-b3340   | fusA      | 17,712.29   | 18,874.86   | 16,844.92   | 17,810.69    |
| gene-b2621   | ssrA      | 14,149.51   | 15,039.24   | 15,069.22   | 14,752.66    |
| gene-b3123   | rnpB      | 13,611.71   | 11,121.83   | 14,252.66   | 12,995.40    |
| gene-b3320   | rplC      | 12,134.63   | 12,648.14   | 9,945.01    | 11,575.93    |
| gene-b3317   | rplB      | 12,034.20   | 12,074.31   | 9,736.04    | 11,281.52    |
| gene-b3985   | rplJ      | 11,412.62   | 11,683.54   | 9,977.71    | 11,024.63    |
| gene-b3295   | rpoA      | 10,787.53   | 11,398.34   | 9,119.11    | 10,434.99    |
| gene-b3984   | rplA      | 10,020.59   | 9,689.17    | 8,297.16    | 9,335.64     |

**Table 1:** Top 10 most highly expressed genes ranked by mean CPM (Counts Per Million) across three samples. Values represent normalized gene expression levels.

Example top genes:

-   `ssrS` - 6S RNA (very high expression ~200,000 CPM)
-   `cspA` - Cold shock protein A
-   `rplA`, `rplC`, `rplB` - Ribosomal proteins
-   `fusA` - Elongation factor G

### 3\. Regular CPM Heatmap

<img width="1200" height="1000" alt="top10_genes_heatmap" src="https://github.com/user-attachments/assets/357e440f-0087-4d85-a834-7c6837a603c0" />


Shows absolute expression levels across samples:

-   White → Yellow → Orange → Red → Dark Red (low to high CPM)
-   Easily see which genes are most abundant
-   `ssrS` appears darkest (highest expression)

### 4\. Log2-Scaled Heatmap

<img width="1200" height="1000" alt="top10_genes_heatmap_log2" src="https://github.com/user-attachments/assets/138fb749-d05d-49e4-b383-43cc2a583b9a" />

Shows relative expression patterns across samples:

-   Navy → White → Red (below average → average → above average)
-   Each gene row is scaled independently
-   Better for comparing fold-changes between samples
-   Useful for identifying sample-specific differences


Folder Structure
----------------

Code

```
02_annotation/
├── code/
│   ├── 01_generate_gene_counts.sh      # Runs featureCounts
│   ├── 02_generate_top10_heatmap.R     # Generates heatmaps (v4.3)
│   └── 03_generate_CPM_table.R         # (Optional) Standalone CPM normalization
├── counts/
│   ├── raw_counts.txt                   # Raw gene expression counts
│   └── raw_counts.txt.summary           # featureCounts summary
├── plots/
│   ├── top10_genes_heatmap.png         # Regular CPM heatmap
│   ├── top10_genes_heatmap_log2.png    # Log2-scaled heatmap
│   └── top10_genes_table.csv           # Top 10 genes table
└── README.md                            # This file

```

* * * * *

Troubleshooting
---------------

### Issue: "No BAM files found"

Solution: Make sure you've completed the alignment step and BAM files exist in `../01_alignment/alignment/bam/`

### Issue: "GFF file not found"

Solution: Check that the reference GFF file exists at the specified path and update the path in the script if needed

### Issue: "R packages not installed"

Solution:

bash

```
Rscript -e 'install.packages(c("dplyr", "pheatmap"), repos="https://cloud.r-project.org")'

```

### Issue: "Sample names showing as long paths"

Solution: The updated v4.3 script handles this automatically by extracting `SRR` IDs from full paths

### Issue: "Heatmap shows gene IDs instead of names"

Solution: Ensure the GFF file path is correct in the R script. The script extracts gene names from the GFF `gene=` attribute.

* * * * *

Pipeline Version
----------------

Current Version: 4.3\
Last Updated: October 26, 2025\
Pipeline Author: BIFS619_GROUP_01 (Tyler Maire)

* * * * *

Notes
-----

-   All scripts use E. coli K-12 MG1655 genome (GCF_000005845.2)
-   Expression is normalized as CPM (Counts Per Million) for cross-sample comparison
-   Gene names are automatically extracted from GFF annotations
-   The pipeline quantifies ~4,506 genes total
-   Top 10 genes are selected based on mean CPM across all samples

* * * * *

References
----------

-   featureCounts: Liao Y, Smyth GK and Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-930.
-   pheatmap: Kolde R (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12.

* * * * *

Contact
-------

For questions or issues:

-   Open a GitHub issue in the [BIFS619_GROUP_01 repository](https://github.com/tylermaire/BIFS619_GROUP_01)
-   Contact: Tyler Maire

* * * * *

*Last updated: October 26, 2025*


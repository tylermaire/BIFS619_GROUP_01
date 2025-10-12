# 01_allignment

This folder contains all files and outputs related to alignment, quality control, and associated analysis scripts for the BIFS619 group project.

## Folder Structure

### `QC/`
Contains outputs from the quality control step, including:
- `cleaned_fastq/`  
  Contains cleaned FASTQ files and FASTP reports (HTML/JSON) for each sample, generated during sequence quality filtering.
- `plots/`  
  Visualization files (e.g. barplots) showing metrics such as pre- and post-cleaning read counts.
- `tables/`  
  CSV tables summarizing QC metrics, including numbers of reads before and after cleaning for each sample.

### `alignment/`
Contains files related to the alignment step:
- `hisat2_index/`  
  HISAT2 index files (`genome.1.ht2` – `genome.8.ht2`) used for mapping reads to the reference genome.
- `logs/`  
  Log files and summary reports from HISAT2 alignment runs.
- `plots/`  
  Barplots visualizing mapping percentages and read counts per sample post-alignment.
- `tables/`  
  CSV tables reporting alignment metrics for all processed samples.
- `bam/`  
  **BAM files are not stored in this repository due to GitHub’s file size limits.  
  Instead, download them from our Google Drive:**  
  [Google Drive - BIFS619 BAM files](https://drive.google.com/drive/folders/1KW6iuHfbfBAplelDN6l0az1ePj9p1X8c?usp=sharing)

### `code/`
Contains analysis scripts used in the pipeline:
- `align_reads.sh`  
  Bash script to run HISAT2 alignments on all samples.
- `read_cleaning.sh`  
  Bash script for performing read cleaning and QC using FASTP.
- `summarize_alignment.py`  
  Python script for summarizing alignment metrics and generating tables/plots.
- `summarize_qc.py`  
  Python script for parsing QC outputs and producing summary tables/plots.

### `test/`
(If present) Contains test files or sample data used for validation.

## Notes

- All file and folder names are referenced directly in analysis scripts.
- For questions or issues, please contact the repository owner or open a GitHub issue.
- See each subfolder’s README, if present, for more details.

---

_Last updated: October 2025_

---

![image2](image2)
![image3](image3)
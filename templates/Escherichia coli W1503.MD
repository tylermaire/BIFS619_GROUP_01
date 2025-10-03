# Group Project: RNA-Seq Workflow (Microbial)

## Preparation
- **GitHub Setup**
  - Team members: __________________________  

## Data Retrieval
- **NCBI SRA accessions used:**
  - SRR9613403
  - SRR9613404
  - SRR9613405


- **Commands used to download data:**

    -  There are three replicates for: 1-NU-C-6-KIRROMYCI 

        - Replicate 1: 1-NU-C-6-KIRROMYCIN-NCATS-52B7-1_S1
      	   - wget "ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/003/SRR9613403/SRR9613403_1.fastq.gz"
      	   - wget "ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/003/SRR9613403/SRR9613403_2.fastq.gz"
	
        - Replicate 2: 1-NU-C-6-KIRROMYCIN-NCATS-52B7-1_S2
	         - wget "ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/004/SRR9613404/SRR9613404_1.fastq.gz"
           - wget "ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/004/SRR9613404/SRR9613404_2.fastq.gz"

       -  Replicate 3: 1-NU-C-6-KIRROMYCIN-NCATS-52B7-1_S3
          -	wget "ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/005/SRR9613405/SRR9613405_1.fastq.gz"
        	- wget "ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/005/SRR9613405/SRR9613405_2.fastq.gz"

       - The reference genome was downloaded from NCBI Genbank database 
         - wget -O NZ_CP076404.1.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=NZ_CP076404.1&extrafeat=976&conwithfeat=on&hide-cdd=on"


## Quality Control (QC)
- **Tools Used:** 
  - (Example) FastQC
  - (Example) MultiQC
- **Tasks Performed:**  
  - (Example) Generated per-sample QC reports  
  - (Example) Summarized read quality, GC content, adapter contamination, duplication  
- **Deliverables:**  
  - (Example) Insert QC plots (e.g., per base quality, adapter content)  
  - (Example) Insert Table: Raw read counts, duplication rates  
- **Interpretation:**  
  - Brief narrative of overall quality  
  - Any concerns identified in the QC step  

## Read Cleaning
- **Tools Used:** 
- **Tasks Performed:**  
- **Deliverables:**  
- **Interpretation:**  

## Read Alignment 
- **Tools Used:** 
- **Tasks Performed:**  
- **Deliverables:**  
- **Interpretation:**  



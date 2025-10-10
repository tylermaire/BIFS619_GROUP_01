rule all:
    input:
        # The final output of the download_fastq rule (the dummy file)
        "fastq_data/samples/DOWNLOAD_COMPLETE.txt",
        # The final output of the download_genome rule
        "fastq_data/reference/NZ_CP076404.1.fasta"

rule download_fastq:
    input:
        # No inputs - these are source files
    output:
        # Define a single, dummy file to signal completion. 
        # Snakemake will look for this file, not the directory.
        "fastq_data/samples/DOWNLOAD_COMPLETE.txt" 
    # REMOVED: params section since we're not using wildcards
    shell:
        # 1. Create the directory
        "mkdir -p fastq_data/samples/ && "
        "echo 'Downloading Samples' && "
        
        # 2. Each wget must explicitly specify the full file path.
        "wget -O fastq_data/samples/SRR9613403_1.fastq.gz \"ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/003/SRR9613403/SRR9613403_1.fastq.gz\" && "
        "wget -O fastq_data/samples/SRR9613403_2.fastq.gz \"ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/003/SRR9613403/SRR9613403_2.fastq.gz\" && "
        "wget -O fastq_data/samples/SRR9613404_1.fastq.gz \"ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/004/SRR9613404/SRR9613404_1.fastq.gz\" && "
        "wget -O fastq_data/samples/SRR9613404_2.fastq.gz \"ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/004/SRR9613404/SRR9613404_2.fastq.gz\" && "
        "wget -O fastq_data/samples/SRR9613405_1.fastq.gz \"ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/005/SRR9613405/SRR9613405_1.fastq.gz\" && "
        "wget -O fastq_data/samples/SRR9613405_2.fastq.gz \"ftp.sra.ebi.ac.uk/vol1/fastq/SRR961/005/SRR9613405/SRR9613405_2.fastq.gz\" && "
        
        # 3. Create the dummy output file to satisfy Snakemake
        "touch {output}"

rule download_genome:
    input:
        # No inputs
    output:
        "fastq_data/reference/NZ_CP076404.1.fasta"
    shell:
        "mkdir -p fastq_data/reference && "  # FIXED: Use explicit path instead of dirname
        "echo 'Downloading Reference genome' && "
        "wget -O {output} \"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=NZ_CP076404.1&extrafeat=976&conwithfeat=on&hide-cdd=on\""

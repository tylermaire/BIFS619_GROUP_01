# ==============================================================================
# 0. CONFIGURATION AND SAMPLES
# ==============================================================================
configfile: "BIFS619_GROUP_01/00_rawdata/config/QAConfig.yaml"

SAMPLES = config["samples"]  # e.g., ["SRR9613403", "SRR9613404", ...]

# Note: CONDA_ENV_DIR is not strictly needed if 'conda' directive is removed,
# but kept here as a comment for context if environments were used elsewhere.
# CONDA_ENV_DIR = "BIFS619_GROUP_01/00_rawdata/envs"

# ==============================================================================
# 1. TARGET RULE
# ==============================================================================
rule all:
    input:
        "BIFS619_GROUP_01/00_rawdata/fastQC/multiqc_report.html",
        expand("BIFS619_GROUP_01/00_rawdata/fastQC/{sample}_fastqc.html",
               sample=SAMPLES)

# ==============================================================================
# 2. FASTQC RULE
# ==============================================================================
rule fastqc_individual:
    input:
        "BIFS619_GROUP_01/00_rawdata/fastq_data/samples/{sample}.fastq.gz"
    output:
        html = "BIFS619_GROUP_01/00_rawdata/fastQC/{sample}_fastqc.html",
        zip = "BIFS619_GROUP_01/00_rawdata/fastQC/{sample}_fastqc.zip"
    shell:
        # The output directory for fastqc is typically the directory of the
        # specified output files, which is "BIFS619_GROUP_01/00_rawdata/fastQC"
        "fastqc {input} -o BIFS619_GROUP_01/00_rawdata/fastQC"

# ==============================================================================
# 3. MULTIQC RULE
# ==============================================================================
rule multiqc_report:
    input:
        expand("BIFS619_GROUP_01/00_rawdata/fastQC/{sample}_fastqc.zip",
               sample=SAMPLES)
    output:
        "BIFS619_GROUP_01/00_rawdata/fastQC/multiqc_report.html"
    shell:
        # MultiQC scans the input directory and writes the report to the output directory.
        "multiqc BIFS619_GROUP_01/00_rawdata/fastQC -o BIFS619_GROUP_01/00_rawdata/fastQC"

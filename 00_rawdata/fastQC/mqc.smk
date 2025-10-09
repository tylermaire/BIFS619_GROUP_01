# ==============================================================================
# 0. CONFIGURATION AND SAMPLES
# ==============================================================================
configfile: "BIFS619_GROUP_01/00_rawdata/config/QAConfig.yaml"

SAMPLES = config["samples"]  

PROJECT_ROOT = "/home/newuser/BIFS619/group_project/" 


# ==============================================================================
# 1. TARGET RULE
# ==============================================================================
rule all:
    input:
        "BIFS619_GROUP_01/00_rawdata/fastQC/multiqc_report.html"

# ==============================================================================
# 2. MULTIQC RULE
# ==============================================================================
rule multiqc_report:
    input:
        # MultiQC scans directories, so we need to ensure the directory exists
        # OR list some expected files that should be there
        # Let's use a dummy input or directory
        directory("BIFS619_GROUP_01/00_rawdata/fastQC")
    output:
        "BIFS619_GROUP_01/00_rawdata/fastQC/multiqc_report.html"
    shell:
        "multiqc BIFS619_GROUP_01/00_rawdata/fastQC -o BIFS619_GROUP_01/00_rawdata/fastQC"

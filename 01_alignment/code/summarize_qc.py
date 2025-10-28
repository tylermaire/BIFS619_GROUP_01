import os
import json
import pandas as pd
import matplotlib.pyplot as plt

# Define paths
samples = ["SRR9613403", "SRR9613404", "SRR9613405"]
qc_dir = "/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/QC/cleaned_fastq"
table_out = "/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/QC/tables/pre_post_cleaning_table.csv"
plot_out = "/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/QC/plots/pre_post_cleaning_barplot.png"

summary = []
for sample in samples:
    fjson = os.path.join(qc_dir, f"{sample}_fastp.json")
    with open(fjson) as f:
        fastp = json.load(f)
    before = fastp['summary']['before_filtering']['total_reads']
    after = fastp['summary']['after_filtering']['total_reads']
    summary.append({"sample": sample, "raw_reads": before, "cleaned_reads": after})

df = pd.DataFrame(summary)
os.makedirs(os.path.dirname(table_out), exist_ok=True)
os.makedirs(os.path.dirname(plot_out), exist_ok=True)
df.to_csv(table_out, index=False)
print(df)

plt.figure(figsize=(7,5))
plt.bar(df['sample'], df['raw_reads'], label="Raw", alpha=0.7)
plt.bar(df['sample'], df['cleaned_reads'], label="Cleaned", alpha=0.7)
plt.ylabel("Read Count")
plt.title("Reads Before and After Cleaning")
plt.legend()
plt.tight_layout()
plt.savefig(plot_out)
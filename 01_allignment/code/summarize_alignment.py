import os, pandas as pd, re, matplotlib.pyplot as plt

samples = ["SRR9613403", "SRR9613404", "SRR9613405"]
metrics = []

for sample in samples:
    log = f"/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/alignment/logs/{sample}_hisat2_summary.txt"
    total, mapped, pct = None, None, None
    with open(log) as f:
        for line in f:
            if "reads; of these:" in line:
                total = int(re.findall(r"\d+", line)[0])
            if "aligned exactly 1 time" in line:
                mapped = int(re.findall(r"\d+", line)[0])
            if "%" in line and "overall alignment rate" in line:
                pct = float(re.findall(r"(\d+\.\d+)%", line)[0])
    metrics.append({"sample": sample, "total_reads": total, "mapped_reads": mapped, "mapping_percent": pct})

df = pd.DataFrame(metrics)
os.makedirs("/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/alignment/tables", exist_ok=True)
os.makedirs("/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/alignment/plots", exist_ok=True)
df.to_csv("/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/alignment/tables/alignment_metrics.csv", index=False)
print(df)

# Graph 1: Mapping percentage bar plot
plt.figure(figsize=(6,4))
plt.bar(df["sample"], df["mapping_percent"], color="mediumslateblue")
plt.ylabel("Mapping %")
plt.title("Mapping Percentage per Sample (HISAT2)")
plt.tight_layout()
plt.savefig("/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/alignment/plots/mapping_percent_barplot.png")

# Graph 2: Total reads and mapped reads per sample
plt.figure(figsize=(7,5))
plt.bar(df["sample"], df["total_reads"], label="Total Reads", alpha=0.7)
plt.bar(df["sample"], df["mapped_reads"], label="Mapped Reads", alpha=0.7)
plt.ylabel("Read Count")
plt.title("Total and Mapped Reads per Sample")
plt.legend()
plt.tight_layout()
plt.savefig("/home/tylermaire/Group_Project/BIFS619_GROUP_01/01_allignment/alignment/plots/read_counts_barplot.png")
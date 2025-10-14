#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(pheatmap)

# Load raw counts using full path
counts <- read.delim("/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/03_annotation/counts/raw_counts.txt",
                     comment.char = "#")

# Keep only BAM count columns
count_data <- counts %>% select(starts_with("X"))

# Rename columns for simplicity
colnames(count_data) <- c("SRR9613403","SRR9613404","SRR9613405")

# Keep gene IDs separately
gene_ids <- counts$Geneid

# Total counts per sample
total_counts <- colSums(count_data)

# CPM normalization
cpm_data <- sweep(count_data, 2, total_counts, FUN=function(x, y) (x / y) * 1e6)
cpm_data <- cbind(Geneid = gene_ids, cpm_data)

# Identify top 10 genes by mean CPM
top_genes <- cpm_data %>%
  rowwise() %>%
  mutate(mean_cpm = mean(c_across(-Geneid))) %>%
  ungroup() %>%
  arrange(desc(mean_cpm)) %>%
  slice(1:10) %>%
  select(-mean_cpm)

# Save top 10 genes to CSV with full path
write.csv(top_genes, 
          "/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/03_annotation/counts/top_10_genes.csv",
          row.names = FALSE)

# Prepare matrix for heatmap
top_cpm_mat <- as.matrix(top_genes[,-1])
rownames(top_cpm_mat) <- top_genes$Geneid

# Generate heatmap for only these top 10 genes
pheatmap(top_cpm_mat,
         scale = "none",         # do not standardize rows
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 10,
         fontsize_col = 10,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Top 10 Expressed Genes (CPM)",
         filename = "/home/StudentFirst/BIFS619/Group_Project/BIFS619_GROUP_01/03_annotation/counts/top_10_genes_heatmap.png")

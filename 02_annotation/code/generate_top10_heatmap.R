#!/usr/bin/env Rscript

# ====================================================================================
# Top 10 Genes Heatmap Generation Script
# Author: BIFS619_GROUP_01
# Description: Generates heatmaps showing top 10 expressed genes with gene names
# Version: 4.3
# ====================================================================================

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
output_heatmap <- args[2]
output_table <- args[3]
gff_file <- args[4]

cat("=== R Heatmap Generation v4.3 ===\n")
cat("Counts file:", counts_file, "\n")
cat("GFF file:", gff_file, "\n\n")

# Create directories for outputs
dir.create(dirname(output_heatmap), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_table), recursive = TRUE, showWarnings = FALSE)

# Function to safely load packages
safe_library <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "https://cloud.r-project.org")
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      cat(paste("Could not install package:", package_name, "\n"))
      return(FALSE)
    }
  }
  return(TRUE)
}

# Load required libraries
dplyr_loaded <- safe_library("dplyr")
pheatmap_loaded <- safe_library("pheatmap")

if (!dplyr_loaded || !pheatmap_loaded) {
  cat("Required R packages could not be loaded. Exiting.\n")
  quit(status = 1)
}

# Function to extract gene names from GFF - CORRECTED FOR E. coli K-12 GFF FORMAT
extract_gene_names <- function(gff_file) {
  cat("\nExtracting gene names from GFF file...\n")
  if (!file.exists(gff_file)) {
    cat("GFF file not found. Using gene IDs as names.\n")
    return(NULL)
  }
  
  gff_lines <- readLines(gff_file)
  gene_map <- data.frame(gene_id = character(), gene_name = character(), stringsAsFactors = FALSE)
  
  for (line in gff_lines) {
    # Skip comment lines
    if (grepl("^#", line)) next
    
    # Look for gene features only
    if (grepl("\tgene\t", line)) {
      # Extract ID - in E. coli K-12 GFF it's like: ID=gene-b0001
      # featureCounts uses this full ID: gene-b0001
      id_match <- regmatches(line, regexpr('ID=[^;]+', line))
      if (length(id_match) > 0) {
        # Keep the full ID as featureCounts uses it
        gene_id <- gsub('ID=', '', id_match)
        
        # Extract gene name from gene= attribute
        # Example: gene=thrL
        gene_name <- gene_id  # default to ID if no gene name found
        
        # Try to get gene name from ";gene=" attribute
        gene_match <- regmatches(line, regexpr(';gene=[^;]+', line))
        if (length(gene_match) > 0) {
          gene_name <- gsub(';gene=', '', gene_match)
        } else {
          # If no ;gene=, try gene= at start of attributes (though unlikely in this format)
          gene_match <- regmatches(line, regexpr('gene=[^;]+', line))
          if (length(gene_match) > 0) {
            gene_name <- gsub('gene=', '', gene_match)
          }
        }
        
        # Only add if we haven't seen this gene_id
        if (!gene_id %in% gene_map$gene_id) {
          gene_map <- rbind(gene_map, data.frame(gene_id = gene_id, gene_name = gene_name, stringsAsFactors = FALSE))
        }
      }
    }
  }
  
  cat("Found", nrow(gene_map), "genes in GFF file\n")
  cat("Example gene mappings (first 15):\n")
  print(head(gene_map, 15))
  return(gene_map)
}

# Load raw counts
cat("\nLoading counts file...\n")
counts_loaded <- tryCatch({
  counts <- read.delim(counts_file, comment.char = "#", stringsAsFactors = FALSE)
  cat("✓ Counts file loaded successfully\n")
  cat("Dimensions:", nrow(counts), "genes x", ncol(counts), "columns\n")
  cat("First 10 Geneid values from counts file:\n")
  print(head(counts$Geneid, 10))
  TRUE
}, error = function(e) {
  cat(paste("Error loading counts file:", e$message, "\n"))
  FALSE
})

if (!counts_loaded) {
  quit(status = 1)
}

if (nrow(counts) == 0 || ncol(counts) <= 1) {
  cat("Counts file appears to be empty or invalid. Exiting.\n")
  quit(status = 1)
}

tryCatch({
  # Load gene name mapping from GFF
  gene_map <- extract_gene_names(gff_file)
  
  # Get count columns (skip first 6 metadata columns from featureCounts)
  count_cols <- 7:ncol(counts)
  count_data <- counts[, count_cols, drop = FALSE]
  
  # Extract sample names from BAM file paths - ROBUST VERSION
  cat("\nExtracting sample names...\n")
  cat("Original column names:\n")
  print(colnames(count_data))
  
  # More robust extraction - handles both / and . as path separators
  samples <- sapply(colnames(count_data), function(x) {
    # Extract SRR followed by digits (most reliable method)
    m <- regexpr("SRR[0-9]+", x)
    if (m > 0) {
      return(regmatches(x, m))
    } else {
      return(x)  # If no SRR found, keep original
    }
  }, USE.NAMES = FALSE)
  
  cat("Extracted sample names:", paste(samples, collapse = ", "), "\n")
  colnames(count_data) <- as.character(samples)
  cat("Final column names:", paste(colnames(count_data), collapse = ", "), "\n")
  
  # Keep gene IDs
  gene_ids <- counts$Geneid
  
  # CPM normalization
  cat("\nPerforming CPM normalization...\n")
  total_counts <- colSums(count_data)
  cat("Total counts per sample:\n")
  print(total_counts)
  
  if (any(total_counts == 0)) {
    cat("Warning: Some samples have zero total counts. Using pseudocount.\n")
    total_counts[total_counts == 0] <- 1
  }
  
  cpm_data <- sweep(count_data, 2, total_counts, FUN=function(x, y) (x / y) * 1e6)
  cpm_data <- data.frame(Geneid = gene_ids, cpm_data, stringsAsFactors = FALSE)
  
  # Map gene IDs to gene names
  cat("\nMapping gene IDs to gene names...\n")
  if (!is.null(gene_map) && nrow(gene_map) > 0) {
    cat("Merging counts with gene map...\n")
    cat("Before merge - CPM data rows:", nrow(cpm_data), "\n")
    cat("Gene map rows:", nrow(gene_map), "\n")
    
    # Merge and keep all genes from cpm_data
    cpm_data <- merge(gene_map, cpm_data, by.x = "gene_id", by.y = "Geneid", all.y = TRUE)
    
    cat("After merge - CPM data rows:", nrow(cpm_data), "\n")
    
    # Fill in missing gene names with gene IDs
    cpm_data$gene_name[is.na(cpm_data$gene_name)] <- cpm_data$gene_id[is.na(cpm_data$gene_name)]
    
    cat("✓ Gene names mapped successfully\n")
    cat("Example mappings after merge (first 15):\n")
    print(head(cpm_data[, c("gene_id", "gene_name")], 15))
  } else {
    cat("Using gene IDs as gene names (no gene map available)\n")
    cpm_data$gene_name <- cpm_data$Geneid
    cpm_data$gene_id <- cpm_data$Geneid
  }
  
  # Identify top 10 genes by mean CPM
  cat("\nIdentifying top 10 expressed genes...\n")
  cat("CPM data columns:", paste(colnames(cpm_data), collapse = ", "), "\n")
  
  # Calculate mean CPM across samples (columns starting with SRR)
  sample_cols <- grep("^SRR", colnames(cpm_data))
  cat("Found", length(sample_cols), "sample columns\n")
  
  if (length(sample_cols) > 0) {
    top_genes <- cpm_data %>%
      mutate(mean_cpm = rowMeans(select(., all_of(sample_cols)))) %>%
      arrange(desc(mean_cpm)) %>%
      head(10)
  } else {
    cat("ERROR: No sample columns found starting with SRR!\n")
    quit(status = 1)
  }
  
  cat("\n=== TOP 10 GENES ===\n")
  print(top_genes %>% select(gene_id, gene_name, mean_cpm))
  
  # Create table for output
  top_genes_output <- top_genes %>%
    select(gene_id, gene_name, everything())
  
  # Extract numerical columns for heatmap with gene names as row names
  heatmap_data <- top_genes %>% 
    select(all_of(sample_cols))
  
  # Use gene names as row names
  rownames(heatmap_data) <- top_genes$gene_name
  
  cat("\nHeatmap data structure:\n")
  cat("Row names (gene names):", paste(rownames(heatmap_data), collapse = ", "), "\n")
  cat("Column names (samples):", paste(colnames(heatmap_data), collapse = ", "), "\n")
  cat("Data dimensions:", nrow(heatmap_data), "x", ncol(heatmap_data), "\n")
  
  # Verify data is numeric
  if (!all(sapply(heatmap_data, is.numeric))) {
    cat("ERROR: Heatmap data contains non-numeric values!\n")
    quit(status = 1)
  }
  
  # Save table
  write.csv(top_genes_output, file = output_table, row.names = FALSE)
  cat("\n✓ Table saved to:", output_table, "\n")
  
  # Create HEATMAP 1: Log2-scaled with row scaling (shows relative expression patterns)
  cat("✓ Generating log2-scaled heatmap...\n")
  output_heatmap_log2 <- sub("\\.png$", "_log2.png", output_heatmap)
  png(output_heatmap_log2, width = 1200, height = 1000, res = 150)
  pheatmap(log2(heatmap_data + 1),
           main = "Top 10 Expressed Genes (CPM, log2 scale)",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize_row = 11,
           fontsize_col = 12,
           scale = "row",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           border_color = "grey60",
           cellwidth = 40,
           cellheight = 25)
  dev.off()
  cat("✓ Log2 heatmap saved to:", output_heatmap_log2, "\n")
  
  # Create HEATMAP 2: Regular CPM values (shows absolute expression levels)
  cat("✓ Generating regular CPM heatmap...\n")
  png(output_heatmap, width = 1200, height = 1000, res = 150)
  pheatmap(heatmap_data,
           main = "Top 10 Expressed Genes (CPM)",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize_row = 11,
           fontsize_col = 12,
           scale = "none",
           color = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(50),
           border_color = "grey60",
           cellwidth = 40,
           cellheight = 25)
  dev.off()
  cat("✓ Regular CPM heatmap saved to:", output_heatmap, "\n")
  
  # Create HEATMAP 3: Sample Correlation Heatmap
  cat("✓ Generating sample correlation heatmap...\n")
  
  # Use ALL genes (not just top 10) for correlation
  # Get all sample data with proper column names
  all_sample_data <- cpm_data %>% 
    select(all_of(sample_cols))
  colnames(all_sample_data) <- colnames(heatmap_data)  # Use cleaned sample names
  
  # Calculate correlation matrix
  cor_matrix <- cor(all_sample_data, method = "pearson")
  
  cat("Correlation matrix:\n")
  print(cor_matrix)
  
  # Create correlation heatmap
  output_correlation <- sub("\\.png$", "_correlation.png", output_heatmap)
  png(output_correlation, width = 800, height = 800, res = 150)
  pheatmap(cor_matrix,
           main = "Sample-to-Sample Correlation (All Genes)",
           display_numbers = TRUE,
           number_format = "%.3f",
           fontsize_number = 14,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize_row = 12,
           fontsize_col = 12,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           breaks = seq(0.8, 1.0, length.out = 51),
           border_color = "grey60",
           cellwidth = 80,
           cellheight = 80)
  dev.off()
  cat("✓ Correlation heatmap saved to:", output_correlation, "\n")
  
  cat("\n=== Analysis complete! ===\n")
  cat("Generated 3 heatmaps:\n")
  cat("  1. Log2-scaled:", output_heatmap_log2, "\n")
  cat("  2. Regular CPM:", output_heatmap, "\n")
  cat("  3. Sample correlation:", output_correlation, "\n")
}, error = function(e) {
  cat(paste("\nERROR:", e$message, "\n"))
  cat("Traceback:\n")
  print(traceback())
  quit(status = 1)
})

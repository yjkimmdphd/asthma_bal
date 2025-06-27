### ========================================================================
# compare the DEG results between eos% > 1 + batch vs eos% > 1 + batch + sex
### ========================================================================
library(tidyverse)

# == load data===
# Define the path for the output folder
output_folder_cell_batch <- file.path("./reports/local_only/deg_bal_bronch~cell2025-01-03")

# Check if the folder exists; if not, create it
if (dir.exists(output_folder_cell_batch)) {
  cell_batch<-read.csv(file.path(output_folder_cell_batch,"deg_bronch_res_sig_16_~ bal_Eos_p_more_1 + Batch_2025-01-03_.csv"),row.names = 1)
}

output_folder_cell_batch_sex <- file.path("./reports/local_only/deg_bal_bronch~cell2025-04-10")

# Check if the folder exists; if not, create it
if (dir.exists(output_folder_cell_batch_sex)) {
  cell_batch_sex<-read.csv(file.path(output_folder_cell_batch_sex,"deg_bronch_res_sig_16_~ bal_Eos_p_more_1 + Batch + Sex_2025-04-10_.csv"),row.names = 1)
}


# == Separate upregulated (log2FC >= 1) and downregulated (log2FC <= -1) genes == 
# For cell_batch
up_genes_batch <- rownames(cell_batch)[cell_batch$log2FoldChange >= 1]
down_genes_batch <- rownames(cell_batch)[cell_batch$log2FoldChange <= -1]
length(up_genes_batch) # Number of upregulated genes in batch model
length(down_genes_batch) # Number of downregulated genes in batch model

# For cell_batch_sex
up_genes_batch_sex <- rownames(cell_batch_sex)[cell_batch_sex$log2FoldChange >= 1]
down_genes_batch_sex <- rownames(cell_batch_sex)[cell_batch_sex$log2FoldChange <= -1]
length(up_genes_batch_sex) # Number of upregulated genes in batch+sex model
length(down_genes_batch_sex) # Number of downregulated genes in batch+sex model

# Find common genes for each category
common_up <- intersect(up_genes_batch, up_genes_batch_sex)
length(common_up) # Number of common upregulated genes

common_down <- intersect(down_genes_batch, down_genes_batch_sex)
length(common_down) # Number of common downregulated genes

# Create a summary data frame
summary_df <- data.frame(
  Regulation = c("Upregulated", "Downregulated"),
  Batch_Only = c(length(setdiff(up_genes_batch, up_genes_batch_sex)), 
                 length(setdiff(down_genes_batch, down_genes_batch_sex))),
  BatchSex_Only = c(length(setdiff(up_genes_batch_sex, up_genes_batch)), 
                    length(setdiff(down_genes_batch_sex, down_genes_batch))),
  Common = c(length(common_up), length(common_down)),
  Total_Batch = c(length(up_genes_batch), length(down_genes_batch)),
  Total_BatchSex = c(length(up_genes_batch_sex), length(down_genes_batch_sex))
)
print(summary_df)

# Create visualizations

# 1. Bar plot for gene counts
# Install and load required packages if not already loaded
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# Reshape data for plotting
plot_data <- data.frame(
  Model = rep(c("Batch", "Batch+Sex"), each = 2),
  Regulation = rep(c("Upregulated", "Downregulated"), 2),
  Count = c(length(up_genes_batch), length(down_genes_batch),
            length(up_genes_batch_sex), length(down_genes_batch_sex))
)

# Create bar plot
p1 <- ggplot(plot_data, aes(x = Model, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Gene Count by Regulation and Model",
       x = "Model", y = "Number of Genes") +
  scale_fill_manual(values = c("Upregulated" = "#FF9999", "Downregulated" = "#9999FF"))

print(p1)

# 2. Venn diagrams for up and down regulated genes
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Upregulated genes Venn diagram
venn.diagram(
  x = list(Batch = up_genes_batch, BatchSex = up_genes_batch_sex),
  filename = "upregulated_genes_venn.png",
  category.names = c("Batch model", "Batch+Sex model"),
  main = "Upregulated Genes (log2FC >= 1)",
  fill = c("#FF9999", "#FFCCCC"),
  output = TRUE
)

# Downregulated genes Venn diagram
venn.diagram(
  x = list(Batch = down_genes_batch, BatchSex = down_genes_batch_sex),
  filename = "downregulated_genes_venn.png",
  category.names = c("Batch model", "Batch+Sex model"),
  main = "Downregulated Genes (log2FC <= -1)",
  fill = c("#9999FF", "#CCCCFF"),
  output = TRUE
)

# 3. Scatter plot comparing log2FoldChange between models
# Get all genes that are |log2FC| >= 1 in either model
all_sig_genes <- unique(c(up_genes_batch, down_genes_batch, 
                          up_genes_batch_sex, down_genes_batch_sex))

# Create data frame for comparison
comparison_df <- data.frame(
  Gene = all_sig_genes,
  Batch_L2FC = cell_batch[all_sig_genes, "log2FoldChange"],
  BatchSex_L2FC = cell_batch_sex[all_sig_genes, "log2FoldChange"]
)

# Add regulation category
comparison_df$Regulation <- "Neither"
comparison_df$Regulation[comparison_df$Batch_L2FC >= 1 & comparison_df$BatchSex_L2FC >= 1] <- "Both Up"
comparison_df$Regulation[comparison_df$Batch_L2FC <= -1 & comparison_df$BatchSex_L2FC <= -1] <- "Both Down"
comparison_df$Regulation[comparison_df$Batch_L2FC >= 1 & comparison_df$BatchSex_L2FC < 1] <- "Up in Batch only"
comparison_df$Regulation[comparison_df$Batch_L2FC <= -1 & comparison_df$BatchSex_L2FC > -1] <- "Down in Batch only"
comparison_df$Regulation[comparison_df$Batch_L2FC < 1 & comparison_df$BatchSex_L2FC >= 1] <- "Up in BatchSex only"
comparison_df$Regulation[comparison_df$Batch_L2FC > -1 & comparison_df$BatchSex_L2FC <= -1] <- "Down in BatchSex only"

# Create scatter plot
p2 <- ggplot(comparison_df, aes(x = Batch_L2FC, y = BatchSex_L2FC, color = Regulation)) +
  geom_point(aes(size = Regulation), alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Comparison of log2FoldChange Between Models",
       x = "Batch Model log2FC", y = "Batch+Sex Model log2FC") +
  scale_color_manual(values = c("Both Up" = "#FF9999", 
                                "Both Down" = "#9999FF",
                                "Up in Batch only" = "#FFCC99", 
                                "Down in Batch only" = "#99CCFF",
                                "Up in BatchSex only" = "#FF99CC", 
                                "Down in BatchSex only" = "#99FFCC",
                                "Neither" = "#CCCCCC")) +
  coord_fixed() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = -1, linetype = "dotted") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = -1, linetype = "dotted") +
  scale_size_manual(values = c("Both Up" = 3, 
                               "Both Down" = 3,
                               "Up in Batch only" = 5, 
                               "Down in Batch only" = 5,
                               "Up in BatchSex only" = 5, 
                               "Down in BatchSex only" = 5,
                               "Neither" = 1))

print(p2)

# 4. Save the most different genes between models
# Calculate absolute difference
comparison_df$Abs_Difference <- abs(comparison_df$Batch_L2FC - comparison_df$BatchSex_L2FC)

# Sort by absolute difference
comparison_df_sorted <- comparison_df[order(comparison_df$Abs_Difference, decreasing = TRUE), ]

# Print top 20 genes with largest differences
print("Top 20 genes with largest log2FoldChange differences between models:")
print(head(comparison_df_sorted, 20))

# Save results to CSV
write.csv(summary_df, "gene_regulation_summary.csv", row.names = FALSE)
write.csv(comparison_df, "gene_comparison_detailed.csv", row.names = FALSE)




# ============ VEnn diatram ====
# Create Venn Diagrams for comparison of upregulated and downregulated genes

# Load required libraries for Venn diagrams
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Suppress grid.newpage message from VennDiagram package
futile.logger::flog.threshold(futile.logger::ERROR)

# Extract upregulated and downregulated genes
# For cell_batch
up_genes_batch <- rownames(cell_batch)[cell_batch$log2FoldChange >= 1]
down_genes_batch <- rownames(cell_batch)[cell_batch$log2FoldChange <= -1]

# For cell_batch_sex
up_genes_batch_sex <- rownames(cell_batch_sex)[cell_batch_sex$log2FoldChange >= 1]
down_genes_batch_sex <- rownames(cell_batch_sex)[cell_batch_sex$log2FoldChange <= -1]

# Count genes in each category
cat("Upregulated genes (log2FC >= 1):\n")
cat("  Batch model:", length(up_genes_batch), "\n")
cat("  Batch+Sex model:", length(up_genes_batch_sex), "\n")
cat("  Common:", length(intersect(up_genes_batch, up_genes_batch_sex)), "\n\n")

cat("Downregulated genes (log2FC <= -1):\n")
cat("  Batch model:", length(down_genes_batch), "\n")
cat("  Batch+Sex model:", length(down_genes_batch_sex), "\n")
cat("  Common:", length(intersect(down_genes_batch, down_genes_batch_sex)), "\n\n")

# Create directory for output if it doesn't exist
if (!dir.exists("venn_diagrams")) {
  dir.create("venn_diagrams")
}

# ============ BASIC VENN DIAGRAMS ============
# 1. Upregulated genes Venn diagram
venn.plot.up <- venn.diagram(
  x = list(
    "Batch" = up_genes_batch, 
    "Batch+Sex" = up_genes_batch_sex
  ),
  filename = "venn_diagrams/upregulated_genes_venn.png",
  output = TRUE,
  
  # Titles
  main = "Upregulated Genes (log2FC >= 1)",
  sub = paste0("Batch: ", length(up_genes_batch), " genes, Batch+Sex: ", length(up_genes_batch_sex), " genes"),
  
  # Colors
  fill = c("#FF9999", "#FFCCCC"),
  col = c("#FF0000", "#FF6666"),
  
  # Area labels
  cat.col = c("#FF0000", "#FF6666"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  
  # Set size
  height = 3000,
  width = 3000,
  resolution = 300,
  
  # Formatting
  euler.d = FALSE,
  scaled = TRUE
)

# 2. Downregulated genes Venn diagram
venn.plot.down <- venn.diagram(
  x = list(
    "Batch" = down_genes_batch, 
    "Batch+Sex" = down_genes_batch_sex
  ),
  filename = "venn_diagrams/downregulated_genes_venn.png",
  output = TRUE,
  
  # Titles
  main = "Downregulated Genes (log2FC <= -1)",
  sub = paste0("Batch: ", length(down_genes_batch), " genes, Batch+Sex: ", length(down_genes_batch_sex), " genes"),
  
  # Colors
  fill = c("#9999FF", "#CCCCFF"),
  col = c("#0000FF", "#6666FF"),
  
  # Area labels
  cat.col = c("#0000FF", "#6666FF"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  
  # Set size
  height = 3000,
  width = 3000,
  resolution = 300,
  
  # Formatting
  euler.d = FALSE,
  scaled = TRUE
)

# ============ ENHANCED VENN DIAGRAMS WITH LABELS ============
# If you have the ggVennDiagram package, it provides nicer visualizations
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  cat("For better Venn diagrams, consider installing: install.packages('ggVennDiagram')\n")
} else {
  library(ggVennDiagram)
  library(ggplot2)
  
  # 3. Enhanced upregulated genes Venn diagram
  venn_data_up <- list(
    "Batch" = up_genes_batch, 
    "Batch+Sex" = up_genes_batch_sex
  )
  
  p3 <- ggVennDiagram(venn_data_up) +
    scale_fill_gradient(low = "#FFCCCC", high = "#FF0000") +
    labs(title = "Upregulated Genes (log2FC >= 1)") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  
  # 4. Enhanced downregulated genes Venn diagram
  venn_data_down <- list(
    "Batch" = down_genes_batch, 
    "Batch+Sex" = down_genes_batch_sex
  )
  
  p4 <- ggVennDiagram(venn_data_down) +
    scale_fill_gradient(low = "#CCCCFF", high = "#0000FF") +
    labs(title = "Downregulated Genes (log2FC <= -1)") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
}
print(p3)
print(p4)
ggsave("venn_diagrams/upregulated_genes_ggvenn.png", p3, width = 8, height = 6, dpi = 300)
ggsave("venn_diagrams/downregulated_genes_ggvenn.png", p4, width = 8, height = 6, dpi = 300)

# ============ COMBINED VENN DIAGRAM ============
# 5. Create a combined Venn diagram showing all four sets
# Install and load the required library
if (!requireNamespace("ggVennDiagram", quietly = TRUE) && !requireNamespace("patchwork", quietly = TRUE)) {
  cat("For combined Venn diagrams, consider installing: install.packages(c('ggVennDiagram', 'patchwork'))\n")
} else if (requireNamespace("ggVennDiagram", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)) {
  library(ggVennDiagram)
  library(patchwork)
  
  # Combined diagram with all four categories
  p5 <- p3 + p4 +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Comparison of Differentially Expressed Genes Between Models",
      subtitle = "Genes with |log2FoldChange| >= 1",
      theme = theme(
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5)
      )
    )
  
}

ggsave("venn_diagrams/combined_venn_diagram.png", p5, width = 12, height = 6, dpi = 300)
# ============ EXPORT GENE LISTS ============
# 6. Export lists of genes in each category
# Create lists of genes in each category
up_only_batch <- setdiff(up_genes_batch, up_genes_batch_sex)
up_only_batch_sex <- setdiff(up_genes_batch_sex, up_genes_batch)
up_both <- intersect(up_genes_batch, up_genes_batch_sex)

down_only_batch <- setdiff(down_genes_batch, down_genes_batch_sex)
down_only_batch_sex <- setdiff(down_genes_batch_sex, down_genes_batch)
down_both <- intersect(down_genes_batch, down_genes_batch_sex)

# Export to CSV files
write.csv(data.frame(Gene = up_only_batch), "venn_diagrams/up_only_batch.csv", row.names = FALSE)
write.csv(data.frame(Gene = up_only_batch_sex), "venn_diagrams/up_only_batch_sex.csv", row.names = FALSE)
write.csv(data.frame(Gene = up_both), "venn_diagrams/up_both_models.csv", row.names = FALSE)

write.csv(data.frame(Gene = down_only_batch), "venn_diagrams/down_only_batch.csv", row.names = FALSE)
write.csv(data.frame(Gene = down_only_batch_sex), "venn_diagrams/down_only_batch_sex.csv", row.names = FALSE)
write.csv(data.frame(Gene = down_both), "venn_diagrams/down_both_models.csv", row.names = FALSE)

# Create summary table
summary_df <- data.frame(
  Category = c("Upregulated in Batch only", 
               "Upregulated in Batch+Sex only", 
               "Upregulated in both models", 
               "Downregulated in Batch only", 
               "Downregulated in Batch+Sex only", 
               "Downregulated in both models"),
  Count = c(length(up_only_batch), 
            length(up_only_batch_sex), 
            length(up_both), 
            length(down_only_batch), 
            length(down_only_batch_sex), 
            length(down_both))
)

# Add percentage columns
total_batch <- length(up_genes_batch) + length(down_genes_batch)
total_batch_sex <- length(up_genes_batch_sex) + length(down_genes_batch_sex)

summary_df$Percent_of_Model <- NA
summary_df$Percent_of_Model[1:3] <- c(
  round(length(up_only_batch) / total_batch * 100, 1),
  round(length(up_only_batch_sex) / total_batch_sex * 100, 1),
  round(length(up_both) / total_batch * 100, 1)
)
summary_df$Percent_of_Model[4:6] <- c(
  round(length(down_only_batch) / total_batch * 100, 1),
  round(length(down_only_batch_sex) / total_batch_sex * 100, 1),
  round(length(down_both) / total_batch * 100, 1)
)

# Export summary table
write.csv(summary_df, "venn_diagrams/gene_overlap_summary.csv", row.names = FALSE)

# Print completion message
cat("\nVenn diagrams created and saved to 'venn_diagrams' folder\n")
cat("Gene lists exported to CSV files\n")




# ============check overlap with GSEA WNT module
## GSEA found Eos% >=1% associated with downregulation of WNT pathway. Check if these are downregulated in both models

# Define the list of genes to check
gene_list <- c("RBPJ", "ADAM17", "HDAC2", "GNAI1", "CTNNB1", "CUL1", "SKP2", "CSNK1E", "NUMB", "KAT2A", 
               "LEF1", "AXIN1", "AXIN2", "HEY1", "DLL1", "NCOR2", "TP53", "HDAC5", "PPARD", "PTCH1", 
               "FZD8", "CCND2", "HDAC11", "NKD1", "PSEN2", "FRAT1", "MAML1", "NOTCH1", "FZD1", "DVL2", 
               "DKK1", "NCSTN", "JAG1", "TCF7", "NOTCH4", "JAG2", "HEY2", "WNT6", "MYC", "WNT5B")

# Function to check if a gene is in a dataframe and get its log2FoldChange
get_log2FC <- function(gene, df) {
  if(gene %in% rownames(df)) {
    return(df[gene, "log2FoldChange"])
  } else {
    return(NA)
  }
}

# Create a results dataframe
results <- data.frame(
  Gene = gene_list,
  Batch_log2FC = sapply(gene_list, function(gene) get_log2FC(gene, cell_batch)),
  BatchSex_log2FC = sapply(gene_list, function(gene) get_log2FC(gene, cell_batch_sex)),
  stringsAsFactors = FALSE
)

# Determine regulation group
results$Regulation_Batch <- ifelse(is.na(results$Batch_log2FC), "Not in dataset",
                                   ifelse(results$Batch_log2FC >= 1, "Upregulated",
                                          ifelse(results$Batch_log2FC <= -1, "Downregulated", "Not significant")))

results$Regulation_BatchSex <- ifelse(is.na(results$BatchSex_log2FC), "Not in dataset",
                                      ifelse(results$BatchSex_log2FC >= 1, "Upregulated",
                                             ifelse(results$BatchSex_log2FC <= -1, "Downregulated", "Not significant")))

# Define the combined regulation group
results$Combined_Group <- "Unknown"

for (i in 1:nrow(results)) {
  if (is.na(results$Batch_log2FC[i]) || is.na(results$BatchSex_log2FC[i])) {
    results$Combined_Group[i] <- "Missing in one or both datasets"
  } else if (results$Batch_log2FC[i] >= 1 && results$BatchSex_log2FC[i] >= 1) {
    results$Combined_Group[i] <- "Both Up"
  } else if (results$Batch_log2FC[i] <= -1 && results$BatchSex_log2FC[i] <= -1) {
    results$Combined_Group[i] <- "Both Down"
  } else if (results$Batch_log2FC[i] >= 1 && results$BatchSex_log2FC[i] < 1) {
    results$Combined_Group[i] <- "Up in Batch only"
  } else if (results$Batch_log2FC[i] <= -1 && results$BatchSex_log2FC[i] > -1) {
    results$Combined_Group[i] <- "Down in Batch only"
  } else if (results$Batch_log2FC[i] < 1 && results$BatchSex_log2FC[i] >= 1) {
    results$Combined_Group[i] <- "Up in BatchSex only"
  } else if (results$Batch_log2FC[i] > -1 && results$BatchSex_log2FC[i] <= -1) {
    results$Combined_Group[i] <- "Down in BatchSex only"
  } else {
    results$Combined_Group[i] <- "Not significant in either"
  }
}

# Calculate log2FC difference
results$Log2FC_Difference <- results$Batch_log2FC - results$BatchSex_log2FC

# Sort by regulation group
results <- results[order(results$Combined_Group, -abs(results$Log2FC_Difference)), ]

# Print the results
print(results, row.names = FALSE)

# Summarize the counts by group
group_counts <- table(results$Combined_Group)
print("Summary of regulation groups:")
print(group_counts)

# Check which genes are not in either dataset
missing_genes <- results$Gene[is.na(results$Batch_log2FC) | is.na(results$BatchSex_log2FC)]
if (length(missing_genes) > 0) {
  cat("\nThe following genes were not found in one or both datasets:\n")
  cat(paste(missing_genes, collapse = ", "), "\n")
}

# Export results to CSV
write.csv(results, "gene_regulation_analysis.csv", row.names = FALSE)
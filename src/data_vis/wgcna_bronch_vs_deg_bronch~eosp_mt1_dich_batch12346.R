library(dplyr)
wgcna_folder<-"./reports/local_only/wgcna/bronch/module-gene_list"

# ----------------
# Bronchial DEG for BAL Eos % > 1% vs <=1%
# ----------------
deg_folder<-file.path("./reports/local_only/deg_bal_bronch~cell2025-01-03")
deg_file<-if(file.exists(deg_folder)){list.files(deg_folder)[grep(".csv",list.files(deg_folder))]}else{print("folder doesn't exist")}
deg<-read.csv(file.path(deg_folder,"deg_bronch_res_sig_16_~ bal_Eos_p_more_1 + Batch_2025-01-03_.csv"), row.names = 1)

all_assessed_genes<-read.csv(file.path(deg_folder,"deg_bronch_res_all_16_~ bal_Eos_p_more_1 + Batch_2025-01-03_.csv"), row.names = 1)%>%rownames
deg_abs_lfc_mt1<-deg%>%filter(abs(log2FoldChange)>1)%>%rownames

module_list<-list.files(wgcna_folder)[grep(".txt",list.files(wgcna_folder))]

modules_files<-file.path(wgcna_folder,module_list)

m_gene_list<-lapply(modules_files,read.table)

names(m_gene_list)<-sub("batch12346.txt","",module_list)

overlap_proportion <- sapply(m_gene_list, function(x) {
  mean(unlist(x) %in% deg_abs_lfc_mt1)
})# no notable overlap between bronch down gene with WGCNA module genes 
print(overlap_proportion)

###########################################
# Example Over-representation Analysis
###########################################

# 1. Predefine a target set of genes (deg_abs_lfc_mt1)
deg_abs_lfc_mt1 

# 2. module gene sets (m_gene_list)
m_gene_list 

# 3. Define the "universe" of genes
N <- length(all_assessed_genes)  # total number of genes used in DEG analysis

# 4. Prepare a data frame to store results
results_df <- data.frame(
  gene_set       = character(),
  size_of_set    = numeric(),
  size_of_target = numeric(),
  overlap        = numeric(),
  p_value        = numeric(),
  stringsAsFactors = FALSE
)

# 5. Calculate the size of the target set
d <- length(deg_abs_lfc_mt1)

# 6. Loop over each gene set in m_gene_list
for (i in seq_along(m_gene_list)) {
  current_set <- m_gene_list[[i]]
  set_name    <- names(m_gene_list)[i]  # name of the gene set
  
  # Size of this gene set
  k <- nrow(current_set)
  
  # Overlap with the target set
  overlap <- length(intersect(unlist(current_set), deg_abs_lfc_mt1))
  
  # Construct a 2x2 contingency table for Fisher's test:
  #
  #            in_target   not_in_target
  # in_set A        i          (k - i)
  # not_in_set A   (d - i)     N - k - (d - i)
  #
  table_2x2 <- matrix(
    c(overlap,
      k - overlap,
      d - overlap,
      N - k - (d - overlap)),
    nrow = 2,
    byrow = TRUE
  )
  
  # Perform Fisher's exact test for over-representation
  fisher_res <- fisher.test(table_2x2, alternative = "greater")
  
  # Append results to results_df
  results_df <- rbind(
    results_df,
    data.frame(
      gene_set       = set_name,
      size_of_set    = k,
      size_of_target = d,
      overlap        = overlap,
      p_value        = fisher_res$p.value,
      stringsAsFactors = FALSE
    )
  )
}

# 7. Correct for multiple testing (Benjamini-Hochberg)
results_df$padj_fdr <- p.adjust(results_df$p_value, method = "fdr")

# 8. Sort by adjusted p-value and inspect top results
results_df <- results_df[order(results_df$padj_fdr), ]
print(results_df)

output_dir <- file.path(wgcna_folder, "output")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Compute the intersection between "brown" gene list and your target list
shared_genes <- intersect(unlist(m_gene_list["brown"]), deg_abs_lfc_mt1)

# Write the intersection to a tab-delimited .txt file
output_file <- file.path(output_dir, "brown_deg_abs_lfc_mt1_intersection.txt")
write.table(
  data.frame(shared_genes = shared_genes),
  file      = output_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Compute the intersection between "brown" gene list and your target list
shared_genes <- intersect(unlist(m_gene_list["ivory"]), deg_abs_lfc_mt1)

# Write the intersection to a tab-delimited .txt file
output_file <- file.path(output_dir, "ivory_deg_abs_lfc_mt1_intersection.txt")
write.table(
  data.frame(shared_genes = shared_genes),
  file      = output_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


# ----------------
# for BAL Eos % > 1 vs <=1 + BAL Neut % >4 vs <=4
# ----------------
deg_folder<-file.path("./reports/local_only/deg_bronch~eos-neut-mixed-pauci+batch123456")
deg_file<-if(file.exists(deg_folder)){list.files(deg_folder)[grep(".csv",list.files(deg_folder))]}else{print("folder doesn't exist")}
deg<-read.csv(file.path(deg_folder,"deg_bronch_res_sig_1_3_comp1_eos_vs_pauci_2024-11-05_.csv"), row.names = 1)

all_assessed_genes<-read.csv(file.path(deg_folder,"deg_bronch_res_all_1_3_comp1_eos_vs_pauci_2024-11-05_.csv"), row.names = 1)%>%rownames
deg_abs_lfc<-deg%>%filter(abs(log2FoldChange)>1)%>%rownames

module_list<-list.files(wgcna_folder)[grep(".txt",list.files(wgcna_folder))]

modules_files<-file.path(wgcna_folder,module_list)

m_gene_list<-lapply(modules_files,read.table)

names(m_gene_list)<-sub("batch12346.txt","",module_list)

overlap_proportion <- sapply(m_gene_list, function(x) {
  mean(unlist(x) %in% deg_abs_lfc)
})# no notable overlap between bronch down gene with WGCNA module genes 
print(overlap_proportion)

###########################################
# Over-representation Analysis
###########################################

# 1. Predefine a target set of genes (deg_abs_lfc)
deg_abs_lfc 

# 2. module gene sets (m_gene_list)
m_gene_list 

# 3. Define the "universe" of genes
N <- length(all_assessed_genes)  # total number of genes used in DEG analysis

# 4. Prepare a data frame to store results
results_df <- data.frame(
  gene_set       = character(),
  size_of_set    = numeric(),
  size_of_target = numeric(),
  overlap        = numeric(),
  p_value        = numeric(),
  stringsAsFactors = FALSE
)

# 5. Calculate the size of the target set
d <- length(deg_abs_lfc)

# 6. Loop over each gene set in m_gene_list
for (i in seq_along(m_gene_list)) {
  current_set <- m_gene_list[[i]]
  set_name    <- names(m_gene_list)[i]  # name of the gene set
  
  # Size of this gene set
  k <- nrow(current_set)
  
  # Overlap with the target set
  overlap <- length(intersect(unlist(current_set), deg_abs_lfc))
  
  # Construct a 2x2 contingency table for Fisher's test:
  #
  #            in_target   not_in_target
  # in_set A        i          (k - i)
  # not_in_set A   (d - i)     N - k - (d - i)
  #
  table_2x2 <- matrix(
    c(overlap,
      k - overlap,
      d - overlap,
      N - k - (d - overlap)),
    nrow = 2,
    byrow = TRUE
  )
  
  # Perform Fisher's exact test for over-representation
  fisher_res <- fisher.test(table_2x2, alternative = "greater")
  
  # Append results to results_df
  results_df <- rbind(
    results_df,
    data.frame(
      gene_set       = set_name,
      size_of_set    = k,
      size_of_target = d,
      overlap        = overlap,
      p_value        = fisher_res$p.value,
      stringsAsFactors = FALSE
    )
  )
}

# 7. Correct for multiple testing (Benjamini-Hochberg)
results_df$padj_fdr <- p.adjust(results_df$p_value, method = "fdr")

# 8. Sort by adjusted p-value and inspect top results
results_df <- results_df[order(results_df$padj_fdr), ]
print(results_df)

output_dir <- file.path(wgcna_folder, "output")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Compute the intersection between "brown" gene list and your target list
shared_genes <- intersect(unlist(m_gene_list["ivory"]), deg_abs_lfc)

# Write the intersection to a tab-delimited .txt file
output_file <- file.path(output_dir, "ivory_deg_abs_lfc_eos-pauci_intersection.txt")
write.table(
  data.frame(shared_genes = shared_genes),
  file      = output_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


library(tidyverse)
library(WGCNA)

#===============================================================================
# Setup Project Paths
#===============================================================================
project_base <- "."
resources_dir <- file.path(project_base, "resources", "processed_data")
reports_dir <- file.path(project_base, "reports", "local_only")
wgcna_dir <- file.path(reports_dir, "wgcna", "bronch")
output_dir <- file.path(wgcna_dir, "output_2025-01-18") # to analyze wGNCA from january use 'output'
deg_dir <- file.path(reports_dir, "deg_bal_bronch~cell2025-01-03")
module_list_dir <- file.path(output_dir, "module-gene_list")

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  print("output_dir doesn't exist")
  dir.create(output_dir, recursive = TRUE)
}else{
  print("output_dir exists")
}

#===============================================================================
# Load Expression Data
#===============================================================================
normalized_count_path <- file.path(resources_dir, "normalized_gene_count", 
                                   "normalized_gene_count_bronch_vsd_batch-corrected_2025-05-30.txt")
counts <- read.table(normalized_count_path, header = TRUE, row.names = 1, sep = "\t")
genes <- rownames(counts)

# Transform for expression analysis
expression <- as.data.frame(t(counts))
colnames(expression) <- genes

# Check data quality and remove outliers
gsg <- goodSamplesGenes(expression)
if (!gsg$allOK) {
  expression <- expression[gsg$goodSamples, gsg$goodGenes]
}

# Identify outlier samples using a dendrogram
sampleTree <- hclust(dist(expression), method = "average")
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 160, minSize = 7)
expression.data <- expression[cut.sampleTree == 1, ]

#===============================================================================
# Load Previously Computed WGCNA Results
#===============================================================================
# Load merged colors
mergedColors <- if(file.exists(file.path(output_dir, "mergedColors.txt"))) {
  read.delim(file.path(output_dir, "mergedColors.txt"), sep = "\t", header = TRUE, row.names = 1)
}
colnames(mergedColors) <- "modules"

# Load merged module eigengenes
mergedMEs <- if(file.exists(file.path(output_dir, "mergedMEs.txt"))) {
  read.delim(file.path(output_dir, "mergedMEs.txt"), sep = "\t", header = TRUE, row.names = 1)
}

#===============================================================================
# Load and Process Phenotype Data
#===============================================================================
phen_path <- file.path(resources_dir, "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2025-02-14.csv")
phen <- read.csv(phen_path)

# Load sampling date differences
sampling_date_path <- file.path(resources_dir, "sampling_dates", "swab-bal-cbc_differences_in_days.txt")
if (file.exists(sampling_date_path)) {
  sampling_date_diff <- read.table(sampling_date_path, row.names = NULL, header = TRUE) %>%
    filter(Comparison == "blood_bal") 
  
  colnames(sampling_date_diff)<-c("ID","Comparison","Difference")
  
  # Join with phenotype data
  phen <- left_join(phen, sampling_date_diff, by = "ID")
}

# Process phenotype data
phen$Batch <- factor(phen$Batch, levels = unique(phen$Batch))
phen_bronch <- phen[grepl("^B", phen$SampleID), ]
phen_input <- phen_bronch
phen_input$SampleID <- gsub("-", ".", phen_input$SampleID)

#===============================================================================
# Create Categorical Variables Based on Thresholds
#===============================================================================
# Define variable groups
source.cell.log <- c(
  "BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log",
  "BAL_wbc_log", "blood_eos_log", "blood_eos_p_log", "blood_neut_log",
  "blood_neut_p_log", "blood_wbc_log"
)

var_dichot_bal <- c(
  "bal_AEC_more_0", "bal_AEC_more_1", "bal_AEC_more_3", "bal_AEC_more_5",
  "bal_Eos_p_more_0", "bal_Eos_p_more_1", "bal_Eos_p_more_3",
  "bal_ANC_more_0", "bal_ANC_more_5", "bal_ANC_more_13",
  "bal_neut_p_more_0", "bal_neut_p_more_2", "bal_neut_p_more_5"
)

var_dichot_blood <- c(
  "bld_AEC_more_0", "bld_AEC_more_100",
  "bld_AEC_more_300", "bld_AEC_more_500"
)

# Create threshold variables more efficiently
threshold_variables <- list(
  BAL_eos_ct = list(prefix = "bal_AEC_more_", thresholds = c(0, 1, 1.2, 3, 5)),
  BAL_eos_p = list(prefix = "bal_Eos_p_more_", thresholds = c(0, 1, 3)),
  BAL_neut_ct = list(prefix = "bal_ANC_more_", thresholds = c(0, 5, 13)),
  BAL_neut_p = list(prefix = "bal_neut_p_more_", thresholds = c(0, 2, 5)),
  blood_eos = list(prefix = "bld_AEC_more_", thresholds = c(0, 100, 300, 500))
)

for (var_name in names(threshold_variables)) {
  var_info <- threshold_variables[[var_name]]
  for (threshold in var_info$thresholds) {
    new_var_name <- paste0(var_info$prefix, threshold)
    phen_input[[new_var_name]] <- as.numeric(phen_input[[var_name]] > threshold)
  }
}

# Convert to factors
phen_input <- phen_input %>%
  mutate(across(c(var_dichot_bal, var_dichot_blood), ~ factor(., levels = c(0, 1))))

#===============================================================================
# Prepare Cell Count Variables for Analysis
#===============================================================================
# Combine variables to test
var_to_test <- c(source.cell.log, var_dichot_bal, var_dichot_blood)

# Filter for blood-related variables
var_to_test_bld <- var_to_test[c(grep("blood", var_to_test), grep("bld", var_to_test))]

# Identify samples with CBC data and sampling date difference < 1 year
cbc_sampleID <- phen_input %>%
  # Check if any of the blood variables are NA
  filter(if_all(all_of(var_to_test_bld), ~!is.na(.)), 
         abs(Difference) < 365) %>%
  pull(SampleID)

#===============================================================================
# Prepare Trait Data for BAL and Blood Analysis
#===============================================================================
# For BAL phenotype
alltraits <- phen_input[phen_input$SampleID %in% rownames(mergedMEs), ]
rownames(alltraits) <- alltraits$SampleID

bal_vars <- c("FEV1_percent", "BAL_eos_ct_log", "BAL_eos_p_log", 
              "BAL_neut_ct_log", "BAL_neut_p_log", "BAL_wbc_log", var_dichot_bal)
alltraits_bal <- alltraits[, bal_vars]
good_bal <- !rowSums(is.na(alltraits_bal[, -1])) > 0
datTraits_bal <- alltraits_bal[good_bal, ]

# For blood phenotype
alltraits_bld <- phen_input %>%
  filter(SampleID %in% rownames(mergedMEs), 
         SampleID %in% cbc_sampleID) %>%
  select(SampleID, FEV1_percent, all_of(var_to_test_bld))
rownames(alltraits_bld) <- alltraits_bld$SampleID
alltraits_bld <- alltraits_bld[, -1]  # Remove SampleID column
good_bld <- !rowSums(is.na(alltraits_bld[, -1])) > 0
datTraits_bld <- alltraits_bld[good_bld, ]

#===============================================================================
# Create Gene-Module Table
#===============================================================================
gene.module.table <- data.frame(
  genes = colnames(expression.data),
  modules = mergedColors$modules
)

gene.module.table$modules<-factor(gene.module.table$modules,levels=unique(gene.module.table$modules))
gmt_list<-split.data.frame(gene.module.table,gene.module.table$modules)

gmt_list_folder<-file.path(output_dir,"gmt_list")
if(!file.exists(gmt_list_folder)){
  dir.create(gmt_list_folder)
}else{
  print("folder exists")
}

lapply(gmt_list,function(list){
  write.table(unlist(as.character(list$genes)),
              file.path(gmt_list_folder,
                        paste(unique(list$modules),"module_gene_list.txt",sep="_")
                        ),
              col.names = FALSE,
              row.names = FALSE,
              sep = "\t"
              )
})


#===============================================================================
# Module-Trait Correlations for BAL
#===============================================================================
# Helper function for module-trait correlation analysis
calculate_module_trait_correlation <- function(mergedMEs_subset, datTraits_subset, output_prefix, output_dir = getwd()) {
  # Convert factors to numeric (works for binary and ordinal factors)
  datTraits_numeric <- datTraits_subset
  factor_cols <- sapply(datTraits_numeric, is.factor)
  if(any(factor_cols)) {
    datTraits_numeric[factor_cols] <- lapply(datTraits_numeric[factor_cols], 
                                             function(x) as.numeric(as.character(x)))
  }
  
  # Calculate correlations
  nSamples <- nrow(mergedMEs_subset)
  correlation <- cor(mergedMEs_subset, datTraits_numeric, use = "p")
  pvalue <- corPvalueStudent(correlation, nSamples)
  
  # Format text matrix for heatmap
  textMatrix <- paste(signif(correlation, 2), "\n(",
                      signif(pvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(correlation)
  
  # Create heatmap
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  png(file.path(output_dir, paste0("module-trait-correlation_", output_prefix, "_", Sys.Date(), ".png")), 
      width = 1200, height = 1200)
  
  par(mar = c(12, 12, 8, 8))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  
  labeledHeatmap(
    Matrix = correlation,
    xLabels = names(datTraits_numeric),
    yLabels = names(mergedMEs_subset),
    ySymbols = names(mergedMEs_subset),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 1.1,
    zlim = c(-1, 1),
    main = "Module-trait correlation"
  )
  
  dev.off()
  
  return(list(correlation = correlation, pvalue = pvalue, textMatrix = textMatrix))
}

# Filter data for BAL analysis
expression.data_bal <- expression.data[rownames(expression.data) %in% rownames(datTraits_bal), ]
mergedMEs_bal <- mergedMEs[rownames(mergedMEs) %in% rownames(datTraits_bal), ]

# Calculate correlations for BAL
bal_corr_results <- calculate_module_trait_correlation(mergedMEs_bal, datTraits_bal, "bal")

# Create subset heatmap for main figure
trait_subset <- c("FEV1_percent", "BAL_eos_p_log", "bal_Eos_p_more_1")
fig_text_index <- which(colnames(bal_corr_results$correlation) %in% trait_subset)

png(file.path(output_dir, paste0("module-trait-correlation_bal_subset_", Sys.Date(), ".png")), 
    width = 400, height = 1200)

par(mar = c(10, 10, 4, 4))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

labeledHeatmap(
  Matrix = bal_corr_results$correlation[, trait_subset],
  xLabels = names(datTraits_bal[, trait_subset]),
  yLabels = names(mergedMEs_bal),
  ySymbols = names(mergedMEs_bal),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = bal_corr_results$textMatrix[, fig_text_index],
  setStdMargins = FALSE,
  cex.text = 1.1,
  # zlim = c(-1, 1),
  main = "Module-trait correlation"
)

dev.off()

#===============================================================================
# Module-Trait Correlations for Blood
#===============================================================================
# Filter data for blood analysis
expression.data_bld <- expression.data[rownames(expression.data) %in% rownames(datTraits_bld), ]
mergedMEs_bld <- mergedMEs[rownames(mergedMEs) %in% rownames(datTraits_bld), ]

# Calculate correlations for blood
bld_corr_results <- calculate_module_trait_correlation(mergedMEs_bld, datTraits_bld, "bld")

#===============================================================================
# Load DEG Results
#===============================================================================
# Fix: Use deg_dir instead of file_path
deg_folder <- deg_dir
file_names <- list.files(deg_folder)

# Load DEG results for BAL eos % > 1 
deg_file_bal_eos_p_mt1 <- file_names[grep("deg_bronch_res_sig_16_", file_names)]
deg_results_bal_eos_p_mt1 <- read.csv(file.path(deg_folder, deg_file_bal_eos_p_mt1), row.names = 1)

# Load all DEG results
deg_file <- file_names[grep("deg_bronch_res_sig", file_names)]
deg_results <- lapply(file.path(deg_folder, deg_file), function(d) read.csv(d, row.names = 1))

# Extract descriptive names from filenames
extracted_strings <- sapply(deg_file, function(x) {
  string_after_tilde <- trimws(strsplit(x, "~")[[1]][2])
  string_without_csv <- sub("\\+ Batch_2025-01-03_.csv$", "", string_after_tilde)
  return(string_without_csv)
})
names(deg_results) <- extracted_strings

#===============================================================================
# Helper Functions for Module-DEG Overlap Analysis
#===============================================================================
analyze_deg_module_overlap <- function(deg_pattern, gene_module_table, output_prefix) {
  # Filter DEG files by pattern
  input_deg_files <- deg_file[grep(deg_pattern, deg_file, ignore.case = TRUE)]
  input_deg_files <- input_deg_files[grep("sig", input_deg_files)]
  
  # Read DEG files
  deg_list <- lapply(input_deg_files, function(x) {
    read.csv(file.path(deg_folder, x), row.names = 1)
  })
  names(deg_list) <- sub(".*~\\s*(.*?)\\s*\\+.*", "\\1", input_deg_files)
  
  # Extract significant genes with |log2FC| > 1
  deg_abs_lfc <- lapply(deg_list, function(df) {
    filter(df, abs(log2FoldChange) > 1) %>% rownames()
  })
  
  # Read module gene lists
  module_gene_list <-  split(gene_module_table$genes, gene_module_table$modules)

  # Calculate overlap statistics
  overlap_proportion <- overlap_sum <- vector("list", length(deg_list))
  names(overlap_proportion) <- names(overlap_sum) <- names(deg_list)
  
  for (i in seq_along(deg_list)) {
    overlap_proportion[[i]] <- sapply(module_gene_list, function(genes) {
      round(mean(genes %in% deg_abs_lfc[[i]]), 2)
    })
    overlap_sum[[i]] <- sapply(module_gene_list, function(genes) {
      sum(genes %in% deg_abs_lfc[[i]])
    })
  }
  
  # Merge results
  merged_df <- do.call(rbind, lapply(names(overlap_proportion), function(name) {
    data.frame(
      Category = name,
      Module = names(overlap_proportion[[name]]),
      Proportion = overlap_proportion[[name]],
      Count = overlap_sum[[name]],
      stringsAsFactors = FALSE
    )
  }))
  
  # Save results
  output_file <- file.path(output_dir, paste0("module-gene_", output_prefix, "-deg-gene_overlap_", Sys.Date(), ".txt"))
  write.table(merged_df, output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  
  return(list(merged_df = merged_df, deg_abs_lfc = deg_abs_lfc, module_gene_list = module_gene_list))
}

#===============================================================================
# Analyze Blood DEG-Module Overlap
#===============================================================================
bld_overlap_results <- analyze_deg_module_overlap("bld|blood", gene.module.table, "bld")

#===============================================================================
# Analyze BAL DEG-Module Overlap
#===============================================================================
bal_overlap_results <- analyze_deg_module_overlap("bal|BAL", gene.module.table, "bal")

#===============================================================================
# Extract DEGs for BAL Eos % > 1
#===============================================================================
deg_bal_eos_p_mt1_up <- rownames(filter(deg_results$bal_Eos_p_more_1, log2FoldChange >= 1))
deg_bal_eos_p_mt1_down <- rownames(filter(deg_results$bal_Eos_p_more_1, log2FoldChange <= -1))
module_gene_list_interest <- bal_overlap_results$module_gene_list$darkgrey # darkgrey 

# Find overlap with darkslateblue module
overlap_darkgrey_up <- module_gene_list_interest[module_gene_list_interest %in% deg_bal_eos_p_mt1_up]
overlap_darkgrey_down <- module_gene_list_interest[module_gene_list_interest %in% deg_bal_eos_p_mt1_down]

#===============================================================================
# Calculate Module Membership and Gene Significance
#===============================================================================
# Extract module names
modNames <- substring(names(mergedMEs), 3)

# Calculate gene-module membership
nSamples <- nrow(expression.data)
geneModuleMembership <- as.data.frame(WGCNA::cor(expression.data, mergedMEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

#===============================================================================
# Gene Significance for BAL Eos % > 1
#===============================================================================
# Create trait vector
phen_of_interest <- as.data.frame(datTraits_bal$bal_Eos_p_more_1)

# Calculate gene-trait significance
geneTraitSignificance <- as.data.frame(WGCNA::cor(expression.data_bal, phen_of_interest, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(phen_of_interest), sep = "")
names(GSPvalue) <- paste("p.GS.", names(phen_of_interest), sep = "")

#===============================================================================
# Function for Module Membership vs Gene Significance Plots
#===============================================================================
create_membership_plot <- function(geneModuleMembership, mergedColors, modNames, 
                                   geneTraitSignificance, trait_name,
                                   deg_up, deg_down, output_file, 
                                   subset_modules = NULL) {
  
  # Determine which modules to plot
  modules_to_plot <- if (!is.null(subset_modules)) subset_modules else modNames
  
  # Calculate grid dimensions
  n_modules <- length(modules_to_plot)
  n_cols <- min(4, n_modules)
  n_rows <- ceiling(n_modules / n_cols)
  
  # Create plot
  png(output_file, width = 1200, height = 1200)
  par(mar = c(4, 4, 4, 4), mfrow = c(n_rows, n_cols))
  
  for (mod in modules_to_plot) {
    column <- match(mod, modNames)
    moduleGenes <- mergedColors == mod
    
    # Set point colors
    module_color <- rep(mod, sum(moduleGenes))
    if (mod == "ivory") module_color <- rep("grey", sum(moduleGenes))
    
    # Highlight DEGs
    gene_names <- rownames(geneModuleMembership)[moduleGenes]
    is_up_deg <- gene_names %in% deg_up
    is_down_deg <- gene_names %in% deg_down
    module_color[is_up_deg] <- "green"
    module_color[is_down_deg] <- "red"
    
    # Set point sizes
    point_sizes <- rep(2, sum(moduleGenes))
    point_sizes[is_up_deg | is_down_deg] <- 3
    
    # Create scatterplot
    verboseScatterplot(
      abs(geneModuleMembership[moduleGenes, column]),
      abs(geneTraitSignificance[moduleGenes, 1]),
      xlab = paste("Module Membership in", mod, "module"),
      ylab = paste("Gene significance for", trait_name),
      main = paste("Module membership vs. gene significance\n"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
      col = module_color, cex = point_sizes
    )
  }
  
  dev.off()
  cat("Plot saved to:", output_file, "\n")
}

#===============================================================================
# Create Full Module Membership Plot
#===============================================================================
create_membership_plot(
  geneModuleMembership = geneModuleMembership,
  mergedColors = mergedColors,
  modNames = modNames,
  geneTraitSignificance = geneTraitSignificance,
  trait_name = "eos mt1",
  deg_up = deg_bal_eos_p_mt1_up,
  deg_down = deg_bal_eos_p_mt1_down,
  output_file = file.path(output_dir, "Module_membership_vs._gene_significance_eos-p_mt1.png")
)

#===============================================================================
# Create Subset Module Membership Plot
#===============================================================================
# Find modules significantly correlated with trait
modNames_subset <- sub("ME", "", names(which(bal_corr_results$pvalue[, "bal_Eos_p_more_1"] < 0.05)))

create_membership_plot(
  geneModuleMembership = geneModuleMembership,
  mergedColors = mergedColors,
  modNames = modNames,
  geneTraitSignificance = geneTraitSignificance,
  trait_name = "eos mt1",
  deg_up = deg_bal_eos_p_mt1_up,
  deg_down = deg_bal_eos_p_mt1_down,
  output_file = file.path(output_dir, "Module_membership_vs._gene_significance_eos-p_mt1_subset_modules.png"),
  subset_modules = modNames_subset
)

#===============================================================================
# Fisher's Exact Test for Module-DEG Enrichment
#===============================================================================
# Load DEGs
deg <- read.csv(
  file.path(deg_folder, "deg_bronch_res_sig_16_~ bal_Eos_p_more_1 + Batch_2025-01-03_.csv"),
  row.names = 1
)

# Load all assessed genes
all_assessed_genes <- read.csv(
  file.path(deg_folder, "deg_bronch_res_all_16_~ bal_Eos_p_more_1 + Batch_2025-01-03_.csv"),
  row.names = 1
) %>% rownames()

# Filter DEGs by absolute log2FC > 1
deg_abs_lfc <- deg %>%
  filter(abs(log2FoldChange) >= 1) %>%
  rownames()

# Load module gene lists
m_gene_list <-split(gene.module.table$genes, gene.module.table$modules)


# Define universe of genes
all_genes <- union(unique(unlist(m_gene_list)), unique(deg_abs_lfc))
N <- length(all_genes)  # Total genes in universe
d <- length(deg_abs_lfc)  # Size of target set

# Initialize results dataframe
fisher_results_df <- data.frame(
  gene_set = character(),
  size_of_set = numeric(),
  size_of_target = numeric(),
  overlap = numeric(),
  fold_enrichment = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Perform Fisher's exact test for each module
for (i in seq_along(m_gene_list)) {
  current_set <- unlist(m_gene_list[[i]])
  set_name <- names(m_gene_list)[i]
  
  k <- length(current_set)  # Size of gene set
  overlap <- length(intersect(current_set, deg_abs_lfc))  # Overlap
  
  # Create contingency table
  table_2x2 <- matrix(
    c(overlap, k - overlap, d - overlap, N - k - (d - overlap)),
    nrow = 2, byrow = TRUE
  )
  
  # Fisher's exact test
  fisher_res <- fisher.test(table_2x2, alternative = "greater")
  
  # Calculate fold enrichment
  fold_enrichment <- if ((k * d) == 0) NA else overlap / ((k * d) / N)
  
  # Store results
  fisher_results_df <- rbind(
    fisher_results_df,
    data.frame(
      gene_set = set_name,
      size_of_set = k,
      size_of_target = d,
      overlap = overlap,
      fold_enrichment = round(fold_enrichment, 2),
      p_value = fisher_res$p.value,
      stringsAsFactors = FALSE
    )
  )
}

# Adjust p-values and sort results
fisher_results_df$padj_fdr <- p.adjust(fisher_results_df$p_value, method = "fdr")
fisher_results_df <- fisher_results_df[order(fisher_results_df$padj_fdr), ]

# Save results
write.table(
  fisher_results_df,
  file.path(output_dir, "wgcna_bronch_deg_overlap_fishers-exact_enrichment.txt"),
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
)

#===============================================================================
# Network Visualization of Eigengenes
#===============================================================================
# Create subset for visualization
datTraits_bal_subset <- datTraits_bal %>% filter(FEV1_percent > 0)

# Create trait dataframes
bal_eos_mt1 <- as.data.frame(datTraits_bal_subset$bal_Eos_p_more_1)
bal_fev1_perc <- as.data.frame(datTraits_bal_subset$FEV1_percent)
names(bal_eos_mt1) <- "eos_mt1"
names(bal_fev1_perc) <- "fev1_p"

# Combine eigengenes with traits
mergedMEs_subset <- mergedMEs[rownames(mergedMEs) %in% rownames(datTraits_bal_subset), ]
MET <- orderMEs(cbind(mergedMEs_subset, bal_eos_mt1, bal_fev1_perc))

# Create network plots
par(cex = 0.9)
plotEigengeneNetworks(
  MET, "", marDendro = c(0, 4, 1, 2), 
  marHeatmap = c(5, 4, 1, 2), cex.lab = 0.8, 
  xLabelsAngle = 90
)

# Plot dendrogram only
png(file.path(output_dir, paste0("eigengene_dendrogram_", Sys.Date(), ".png")), 
    width = 1000, height = 800)
par(cex = 1.0)
plotEigengeneNetworks(
  MET, "Eigengene dendrogram", 
  marDendro = c(0, 4, 2, 0),
  plotHeatmaps = FALSE
)
dev.off()

# Plot heatmap only
png(file.path(output_dir, paste0("eigengene_heatmap_", Sys.Date(), ".png")), 
    width = 1000, height = 800)
par(cex = 1.0, mar = c(1, 1, 1, 1))
plotEigengeneNetworks(
  MET, "Eigengene adjacency heatmap", 
  marHeatmap = c(10, 10, 2, 2),
  plotDendrograms = FALSE, 
  xLabelsAngle = 90
)
dev.off()

#===============================================================================
# Hub Gene Identification
#===============================================================================
# Calculate adjacency matrix
softPower <- 4  # Determined previously
adjacency <- adjacency(expression.data_bal, power = softPower)

# Identify hub genes
HubGenes <- chooseTopHubInEachModule(expression.data_bal, mergedColors, power = 4, type = "signed")
hubgene_table<-data.frame(Module=names(HubGenes),gene=HubGenes,hubgene="Hub_gene")

# Calculate connectivity
connectivity_allClusters <- intramodularConnectivity(adjacency, mergedColors$modules, scaleByMax = FALSE)

# Add module and gene information
matched <- match(rownames(connectivity_allClusters), gene.module.table$genes)
connectivity_allClusters$module <- gene.module.table$modules[matched]
connectivity_allClusters$gene <- gene.module.table$genes[matched]

# Find top connected genes in each module
top_kWithin_by_module <- connectivity_allClusters %>%
  group_by(module) %>%
  slice_max(order_by = kWithin, n = 20) %>%
  ungroup() %>%
  arrange(module, desc(kWithin))

# Find overlap with DEGs
top_kWithin_by_module_deg_overlap <- top_kWithin_by_module %>%
  filter(gene %in% deg_abs_lfc)%>%
  arrange(module, desc(kWithin))

# Save results
write.table(hubgene_table,file.path(output_dir,"hubgene_table.txt"), sep="\t",row.names = FALSE)
write.table(top_kWithin_by_module,file.path(output_dir,"top_kWithin_by_module.txt"), sep="\t",row.names = FALSE)
write.table(top_kWithin_by_module_deg_overlap,file.path(output_dir,"top_kWithin_by_module_deg_overlap.txt"), sep="\t",row.names = FALSE)

library(tidyverse)
library(WGCNA)

########################
## Load Readcount Data ##
########################

# Load cell count table
normalized_count_table_path <- "./resources/processed_data/normalized_gene_count/normalized_gene_count_bronch_vsd_batch-corrected.txt"
if (file.exists(normalized_count_table_path)) {
  counts <- read.table(normalized_count_table_path, 
                       header = TRUE, 
                       row.names = 1, 
                       sep = "\t")
}
genes <- rownames(counts)

# Define the path for the output folder
output_folder <- file.path("./reports/local_only/wgcna/bronch/output")

# Check if the folder exists; if not, create it
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE) # recursive = TRUE ensures all intermediate directories are created
}

# Now the folder is guaranteed to exist, and the path is set
output_folder

###
# 1. Create a new format for expression data
###
expression <- as.data.frame(t(counts))
colnames(expression) <- genes

# Check if genes/samples are good
gsg <- goodSamplesGenes(expression)

###
# 2. Remove outlier genes (and samples if necessary)
###

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(
      paste("Removing genes:", 
            paste(names(expression)[!gsg$goodGenes], collapse = ", "))
    )
  }
  
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(
      paste("Removing samples:", 
            paste(rownames(expression)[!gsg$goodSamples], collapse = ", "))
    )
  }
  
  expression <- expression[gsg$goodSamples, gsg$goodGenes]
}
###
# 3. Identify outlier samples using a dendrogram
###
sampleTree <- hclust(dist(expression), method = "average")

# png(file.path(output_folder,"sampleTree_nasal.png"), width = 800, height = 600)
# par(cex = 0.6)
# par(mar = c(0, 4, 2, 0))
# plot(
#   sampleTree,
#   main = "Sample clustering to detect outliers",
#   sub = "",
#   xlab = "",
#   cex.lab = 1.5,
#   cex.axis = 1.5,
#   cex.main = 2
# )
# abline(h = 160, col = "red")
# dev.off()

# Remove outliers
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 160, minSize = 7)
expression.data <- expression[cut.sampleTree == 1, ]

###
# 4. load merged colors
### 

mergedColors<-if(file.exists(file.path(output_folder,"mergedColors.txt"))){read.delim(file.path(output_folder,"mergedColors.txt"), sep="\t",header=TRUE,row.names = 1)}
colnames(mergedColors)<-"modules"

mergedMEs<-if(file.exists(file.path(output_folder,"mergedMEs.txt"))){read.delim(file.path(output_folder,"mergedMEs.txt"), sep="\t",header=TRUE, row.names=1)}


# -----------------------------------------------------------------------------
# Load phenotype data
# -----------------------------------------------------------------------------
phen_path <- file.path(
  "./resources/processed_data",
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv"
)
phen <- read.csv(phen_path)

# load data for sampling date differences 
sampling_date_diff<-"./resources/processed_data/sampling_dates/swab-bal-cbc_differences_in_days.txt"
sampling_date_diff<-if(file.exists(sampling_date_diff)){read.table(sampling_date_diff,row.names = NULL,header = TRUE)}else{print("sampling date data doesn't exist")}
sampling_date_diff<-sampling_date_diff%>%filter(Comparison=="blood_bal")
colnames(sampling_date_diff)[1:3]<-c("ID","sampling_date_comp","sampling_date_diff_days")

# left join sampling date data and phenotype data 
phen<-left_join(phen,sampling_date_diff,by="ID")

phen$Batch<-factor(phen$Batch,levels=unique(phen$Batch))

phen_bronch<-phen[grepl("^B",phen$SampleID),]
phen_input<-phen_bronch
phen_input$SampleID <- gsub("-", ".", phen_input$SampleID)

# Create continuous and categorical variables based on various thresholds
## Log-transformed/scaled/centered
source.cell.log <- c(
  "BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log",
  "BAL_wbc_log",    "blood_eos_log", "blood_eos_p_log", "blood_neut_log",
  "blood_neut_p_log", "blood_wbc_log"
)

## Categorical variables to test (BAL)
var_dichot_bal <- c(
  "bal_AEC_more_0", "bal_AEC_more_1", "bal_AEC_more_3", "bal_AEC_more_5",
  "bal_Eos_p_more_0", "bal_Eos_p_more_1", "bal_Eos_p_more_3",
  "bal_ANC_more_0", "bal_ANC_more_5", "bal_ANC_more_13",
  "bal_neut_p_more_0", "bal_neut_p_more_2", "bal_neut_p_more_5"
)

## Categorical variables to test (Blood)
var_dichot_blood <- c(
  "bld_AEC_more_0",   "bld_AEC_more_100",
  "bld_AEC_more_300", "bld_AEC_more_500"
)

phen_input <- phen_input %>%
  mutate(
    bal_AEC_more_0   = as.numeric(BAL_eos_ct > 0),
    bal_AEC_more_1   = as.numeric(BAL_eos_ct > 1),
    bal_AEC_more_1.2 = as.numeric(BAL_eos_ct > 1.2),
    bal_AEC_more_3   = as.numeric(BAL_eos_ct > 3),
    bal_AEC_more_5   = as.numeric(BAL_eos_ct > 5),
    
    bal_Eos_p_more_0 = as.numeric(BAL_eos_p > 0),
    bal_Eos_p_more_1 = as.numeric(BAL_eos_p > 1),
    bal_Eos_p_more_3 = as.numeric(BAL_eos_p > 3),
    
    bal_ANC_more_0   = as.numeric(BAL_neut_ct > 0),
    bal_ANC_more_5   = as.numeric(BAL_neut_ct > 5),
    bal_ANC_more_13  = as.numeric(BAL_neut_ct > 13),
    
    bal_neut_p_more_0 = as.numeric(BAL_neut_p > 0),
    bal_neut_p_more_2 = as.numeric(BAL_neut_p > 2),
    bal_neut_p_more_5 = as.numeric(BAL_neut_p > 5),
    
    bld_AEC_more_0   = as.numeric(blood_eos > 0),
    bld_AEC_more_100 = as.numeric(blood_eos > 100),
    bld_AEC_more_300 = as.numeric(blood_eos > 300),
    bld_AEC_more_500 = as.numeric(blood_eos > 500)
  ) %>%
  mutate(across(c(var_dichot_bal,var_dichot_blood), ~ factor(., levels = c(0, 1))))


#===============================================================================
# Identify Cell-Count Variables
#===============================================================================

# Combine continuous (log counts) and categorical variables
var_to_test <- c(source.cell.log, var_dichot_bal, var_dichot_blood)

# identify blood cell count phenotype variables
var_to_test_bld<-var_to_test[c(grep("blood",var_to_test),grep("bld",var_to_test))]

# Identify samples with CBC data and sampling date difference < 1 year
cbc_sampleID <- phen_input %>%
  filter(!is.na(vars(var_to_test_bld)), abs(sampling_date_diff_days) < 365) %>%
  pull(SampleID)

# subset alltraits for variables pertaining to BAL phenotype
alltraits<-phen_input[phen_input$SampleID%in%rownames(mergedMEs),]
rownames(alltraits)<-alltraits$SampleID

alltraits<-alltraits[,c("asthma_phen_ACT.score","BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log",
                        "BAL_wbc_log", var_dichot_bal)] # only select ACT and cell counts in the phenotype data
good<-!(sapply(alltraits[,-1],is.na)%>%rowSums()>0) # remove rows with at least one NA in cell count phenotype
datTraits_bal<-alltraits[good,]

# subset alltraits for variables pertaining to blood phenotype
alltraits<-phen_input[phen_input$SampleID%in%rownames(mergedMEs),]
alltraits<-alltraits[alltraits$SampleID%in%cbc_sampleID,]
rownames(alltraits)<-alltraits$SampleID

alltraits_bld<-alltraits[,c("asthma_phen_ACT.score",var_to_test_bld)] # only select ACT and cell counts in the phenotype data
good_bld<-!(sapply(alltraits_bld[,-1],is.na)%>%rowSums()>0) # remove rows with at least one NA in cell count phenotype
datTraits_bld<-alltraits_bld[good_bld,]

###
# 5. Module-Trait associations
###

gene.module.table<-data.frame(genes=colnames(expression.data),modules=mergedColors$modules)

# ==================
# quantify the association between the expression profile and ACT score and BAL cell count phenotypes
# ==================
# calculates the correlation of the trait with previously identified module eigengenes. 
# This pairwise correlation is known as the eigengene gene significance
expression.data_bal<-expression.data[rownames(expression.data)%in%rownames(datTraits_bal),]
mergedMEs_bal<-mergedMEs[rownames(mergedMEs)%in%rownames(datTraits_bal),]
nSamples_bal = nrow(expression.data_bal)
module.trait.correlation_bal = cor(mergedMEs_bal, datTraits_bal, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue_bal = corPvalueStudent(module.trait.correlation_bal, nSamples_bal) #calculate the p-value associated with the correlation

# Will display correlations and their p-values
textMatrix_bal = paste(signif(module.trait.correlation_bal, 2), "\n(",
                   signif(module.trait.Pvalue_bal, 1), ")", sep = "");
dim(textMatrix_bal) = dim(module.trait.correlation_bal)

# Display the correlation values within a heatmap plot
png(file.path(output_folder, "module-trait-correlation_bal.png"), width = 1200, height = 1200)

par(mar = c(12, 12, 8, 8))  # Adjust margins
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")  # Empty plot to set margins

labeledHeatmap(
  Matrix = module.trait.correlation_bal,
  xLabels = names(datTraits_bal),
  yLabels = names(mergedMEs_bal),
  ySymbols = names(mergedMEs_bal),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_bal,
  setStdMargins = FALSE,  # Prevents resetting margins
  cex.text = 1.1,
  zlim = c(-1, 1),
  main = "Module-trait correlation"
)

dev.off()

# subset specific traits for main figure
trait_subset<-c("asthma_phen_ACT.score","BAL_eos_p_log","bal_Eos_p_more_1")
fig_text_index<-which(colnames(module.trait.correlation_bal)%in%trait_subset)
# Display the correlation values within a heatmap plot
png(file.path(output_folder, "module-trait-correlation_bal_subset.png"), width = 400, height = 1200)

par(mar = c(10, 10, 4, 4))  # Adjust margins
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")  # Empty plot to set margins

labeledHeatmap(
  Matrix = module.trait.correlation_bal[,trait_subset],
  xLabels = names(datTraits_bal[,trait_subset]),
  yLabels = names(mergedMEs_bal),
  ySymbols = names(mergedMEs_bal),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_bal[,fig_text_index],
  setStdMargins = FALSE,  # Prevents resetting margins
  cex.text = 1.1,
  zlim = c(-1, 1),
  main = "Module-trait correlation"
)

dev.off()

# ==================
# quantify the association between the expression profile and ACT score and blood cell count profile phenotypes
# ==================
# calculates the correlation of the trait with previously identified module eigengenes. 
# This pairwise correlation is known as the eigengene gene significance

expression.data_bld<-expression.data[rownames(expression.data)%in%rownames(datTraits_bld),]
mergedMEs_bld<-mergedMEs[rownames(mergedMEs)%in%rownames(datTraits_bld),]
nSamples_bld = nrow(expression.data_bld)
module.trait.correlation_bld = cor(mergedMEs_bld, datTraits_bld, use="p") #p for pearson correlation coefficient 
module.trait.Pvalue_bld = corPvalueStudent(module.trait.correlation_bld, nSamples_bld) #calculate the p-value associated with the correlation

# Will display correlations and their p-values
textMatrix_bld = paste(signif(module.trait.correlation_bld, 2), "\n(",
                   signif(module.trait.Pvalue_bld, 1), ")", sep = "");
dim(textMatrix_bld) = dim(module.trait.correlation_bld)
# Display the correlation values within a heatmap plot
png(file.path(output_folder, "module-trait-correlation_bld.png"), width = 1200, height = 1200)

par(mar = c(12, 12, 8, 8))  # Adjust margins
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")  # Empty plot to set margins

labeledHeatmap(
  Matrix = module.trait.correlation_bld,
  xLabels = names(datTraits_bld),
  yLabels = names(mergedMEs_bld),
  ySymbols = names(mergedMEs_bld),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_bld,
  setStdMargins = FALSE,  # Prevents resetting margins
  cex.text = 1.1,
  zlim = c(-1, 1),
  main = "Module-trait correlation"
)

dev.off()
# ----------------
# proportion of module gene that overlaps with bronch DEG for 
# blood Eos %, AEC, Neut %, ANC
# ----------------

# Filter BAL-related files
input_deg_files <- deg_file[grep("bld|blood", deg_file, ignore.case = TRUE)]
input_deg_files <- input_deg_files[grep("sig", input_deg_files)]

# Read in DEG files
deg_list <- lapply(input_deg_files, function(x) {
  read.csv(file.path(deg_folder, x), row.names = 1)
})
names(deg_list) <- sub(".*~\\s*(.*?)\\s*\\+.*", "\\1", input_deg_files)


# Extract significant genes with absolute log2FoldChange > 1
deg_abs_lfc <- lapply(deg_list, function(df) {
  filter(df, abs(log2FoldChange) > 1) %>% rownames()
})

# Define WGCNA module folder
module_output_folder <- "./reports/local_only/wgcna/bronch/output/module-gene_list"
wgcna_files <- list.files(module_output_folder, pattern = "\\.txt$", full.names = TRUE)

# Read in module gene lists
module_gene_list <- lapply(wgcna_files, function(file) {
  read.table(file, header = TRUE, sep = "\t", row.names = 1) %>% unlist() %>% as.vector()
})
names(module_gene_list) <- sub("batch12346.txt", "", basename(wgcna_files))

# Calculate overlap proportions and counts
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

# Merge overlap results into a dataframe
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
output_file <- file.path(output_folder, "module-gene_bld-deg-gene_overlap.txt")
write.table(merged_df, output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Print first few rows
head(merged_df)
# ----------------
# bal Eos %, AEC, Neut %, ANC
# ----------------

# Filter BAL-related files
input_deg_files <- deg_file[grep("bal|BAL", deg_file, ignore.case = TRUE)]
input_deg_files <- input_deg_files[grep("sig", input_deg_files)]

# Read in DEG files
deg_list <- lapply(input_deg_files, function(x) {
  read.csv(file.path(deg_folder, x), row.names = 1)
})
names(deg_list) <- sub(".*~\\s*(.*?)\\s*\\+.*", "\\1", input_deg_files)


# Extract significant genes with absolute log2FoldChange > 1
deg_abs_lfc <- lapply(deg_list, function(df) {
  filter(df, abs(log2FoldChange) > 1) %>% rownames()
})

# Define WGCNA module folder
module_output_folder <- "./reports/local_only/wgcna/bronch/output/module-gene_list"
wgcna_files <- list.files(module_output_folder, pattern = "\\.txt$", full.names = TRUE)

# Read in module gene lists
module_gene_list <- lapply(wgcna_files, function(file) {
  read.table(file, header = TRUE, sep = "\t", row.names = 1) %>% unlist() %>% as.vector()
})
names(module_gene_list) <- sub("batch12346.txt", "", basename(wgcna_files))

# Calculate overlap proportions and counts
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

# Merge overlap results into a dataframe
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
output_file <- file.path(output_folder, "module-gene_bal-deg-gene_overlap.txt")
write.table(merged_df, output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Print first few rows
head(merged_df)

# ----------------
# subset bal Eos % > 1 vs <=1 DEG that overlap with ivory module genes
# ----------------
deg_bal_eos_p_mt1_up<- rownames(filter(deg_list$bal_Eos_p_more_1,log2FoldChange>=1))
deg_bal_eos_p_mt1_down<- rownames(filter(deg_list$bal_Eos_p_more_1,log2FoldChange<=-1))
module_gene_list_ivory<-module_gene_list$ivory

deg_abs_lfc_bal$bal_Eos_p_more_1


overlap_ivory_up<-module_gene_list_ivory[module_gene_list_ivory%in%deg_bal_eos_p_mt1_up]
overlap_ivory_down<-module_gene_list_ivory[module_gene_list_ivory%in%deg_bal_eos_p_mt1_down]

view(overlap_ivory_down)

### 
# 6. Target gene identification
### 

# Whole Network connectivity - a measure for how well the node is connected throughout the entire system
# Intramodular connectivity - a measure for how well the node is connected within its assigned module. Also an indicator for how well that node belongs to its module. This is also known as module membership.

modNames = substring(names(mergedMEs), 3) #extract module names

#Calculate the module membership and the associated p-values
if(nrow(expression.data)==nrow(mergedMEs)){nSamples<-nrow(expression.data)}else(print("# samples in expression data and mergedMEs do not match"))
geneModuleMembership = as.data.frame(WGCNA::cor(expression.data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# ----------------
# for BAL Eos % > 1 vs <=1
# ----------------

# Define variable weight containing the weight column of datTrait
phen_of_interest<-as.data.frame(datTraits_bal$bal_Eos_p_more_1)

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(WGCNA::cor(expression.data_bal, phen_of_interest, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(phen_of_interest), sep="")
names(GSPvalue) = paste("p.GS.", names(phen_of_interest), sep="")
head(GSPvalue)

# Define the output file path
output_file <- file.path(output_folder, "Module_membership_vs._gene_significance_eos-p_mt1.png")

# Open PNG device with appropriate dimensions
png(output_file, width = 1200, height = 1200)

# Set graphical parameters for a multi-panel plot
par(mar = c(4,4,4,4), mfrow = c(6,4))

# Loop through each module and plot
for(mod in modNames){
  module = mod
  column = match(module, modNames)
  moduleGenes = mergedColors == module
  
  # Default color is module color, but "ivory" should be "grey"
  module_color <- rep(mod, length(geneModuleMembership[moduleGenes, column]))
  module_color[module == "ivory"] <- "grey"
  
  # Identify DEGs
  deg_indices_up <- rownames(geneModuleMembership)[moduleGenes] %in% deg_bal_eos_p_mt1_up
  deg_indices_down <- rownames(geneModuleMembership)[moduleGenes] %in% deg_bal_eos_p_mt1_down
  module_color[deg_indices_up] <- "red"
  module_color[deg_indices_down] <- "cyan"
  
  # Set size for DEGs: larger than others
  point_sizes <- rep(2, length(moduleGenes))  # Default size
  point_sizes[deg_indices_up] <- 3  # Make DEGs larger
  point_sizes[deg_indices_down] <- 3  # Make DEGs larger
  
  # Plot the scatterplot
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for eos mt1",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
                     col = module_color, cex = point_sizes)  # Adjusted dot sizes
}

# Close the PNG device to save the file
dev.off()

# Print confirmation
cat("Plot saved successfully to:", output_file, "\n")


# ---------------------
# subsetted modules 
# ---------------------

# Define the output file path
output_file <- file.path(output_folder, "Module_membership_vs._gene_significance_eos-p_mt1_subset_modules.png")

# subset modules with significant correlation to traits
modNames_subset<-sub("ME","",names(which(module.trait.Pvalue_bal[,"bal_Eos_p_more_1"]<0.05)))
print(modNames_subset)

# suset genemodulemembership matrix that matches with the modNames_subset
geneModuleMembership_subset<-geneModuleMembership[,which(sub("MM","",colnames(geneModuleMembership))%in%modNames_subset)]

# Open PNG device with appropriate dimensions
png(output_file, width = 800, height = 800)

# Set graphical parameters for a multi-panel plot
par(mar = c(4,4,4,4), mfrow = c(3,3))

# Loop through each module and plot
for(mod in modNames_subset){
  module = mod
  column = match(module, modNames_subset)
  moduleGenes = mergedColors == module
  
  # Default color is module color, but "ivory" should be "grey"
  module_color <- rep(mod, length(geneModuleMembership_subset[moduleGenes, column]))
  module_color[module == "ivory"] <- "grey"
  
  # Identify DEGs
  deg_indices_up <- rownames(geneModuleMembership_subset)[moduleGenes] %in% deg_bal_eos_p_mt1_up
  deg_indices_down <- rownames(geneModuleMembership_subset)[moduleGenes] %in% deg_bal_eos_p_mt1_down
  module_color[deg_indices_up] <- "red"
  module_color[deg_indices_down] <- "cyan"
  
  # Set size for DEGs: larger than others
  point_sizes <- rep(2, length(moduleGenes))  # Default size
  point_sizes[deg_indices_up] <- 3  # Make DEGs larger
  point_sizes[deg_indices_down] <- 3  # Make DEGs larger
  
  # Plot the scatterplot
  verboseScatterplot(abs(geneModuleMembership_subset[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for eos mt1",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
                     col = module_color, cex = point_sizes)  # Adjusted dot sizes
}

# Close the PNG device to save the file
dev.off()

# Print confirmation
cat("Plot saved successfully to:", output_file, "\n")



## 
# 7. fisher's exact test
## 

# 1. Predefine a target set of genes
## 1.1 Define the folder and check that it exists
deg_folder <- file.path("./reports/local_only", "deg_bal_bronch~cell2025-01-03")

if (!dir.exists(deg_folder)) {
  stop("Folder doesn't exist: ", deg_folder)
} else {
  deg_file <- list.files(deg_folder, pattern = "\\.csv$")
  print(deg_file)  # Optional: inspect the CSV files found in deg_folder
}
deg <- read.csv(
  file.path(deg_folder, "deg_bronch_res_sig_16_~ bal_Eos_p_more_1 + Batch_2025-01-03_.csv"),
  row.names = 1
)
## 1.2 Create a vector of all assessed genes
all_assessed_genes <- read.csv(
  file.path(deg_folder, "deg_bronch_res_all_16_~ bal_Eos_p_more_1 + Batch_2025-01-03_.csv"),
  row.names = 1
) %>% rownames()

## 1.3 Filter DEGs by absolute log2 fold-change > 1
deg_abs_lfc <- deg %>%
  dplyr::filter(abs(log2FoldChange) >= 1) %>%
  rownames()

## 1.4 Define WGCNA folder and locate module files
module_output_folder <- file.path("./reports/local_only/wgcna/bronch/output/module-gene_list")

wgcna_folder <- module_output_folder  

module_list <- list.files(wgcna_folder, pattern = "\\.txt$")
modules_files <- file.path(wgcna_folder, module_list)

## 1.5 Read modules into a list and name each element
m_gene_list <- lapply(modules_files, function(x){read.table(x,header = TRUE, col.names = "genes",row.names = 1, sep = "\t")})
names(m_gene_list) <- sub("batch12346\\.txt", "", module_list) # Remove "batch12346.txt" from file names for cleaner naming

# 2. Define the "universe" of genes
all_genes <- union(unique(unlist(m_gene_list)), unique(deg_abs_lfc))
N <- length(all_genes)  # total number of genes in the universe

# 3. Prepare a data frame to store results
fisher_results_df <- data.frame(
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
  current_set <- unlist(m_gene_list[[i]])
  set_name    <- names(m_gene_list)[i]  # name of the gene set
  
  # Size of this gene set (A)
  k <- length(current_set)
  
  # Overlap with the target set (B)
  overlap <- length(intersect(current_set, deg_abs_lfc))
  
  # 2x2 contingency table for Fisher's exact test:
  #            in_target   not_in_target
  # in_set A        i          (k - i)
  # not_in_set A   (d - i)     N - k - (d - i)
  table_2x2 <- matrix(
    c(overlap,
      k - overlap,
      d - overlap,
      N - k - (d - overlap)),
    nrow  = 2,
    byrow = TRUE
  )
  
  # Fisher's exact test for over-representation
  fisher_res <- fisher.test(table_2x2, alternative = "greater")
  
  # Calculate fold enrichment
  # fold_enrichment = overlap / ((k * d) / N)
  if ((k * d) == 0) {
    fold_enrichment <- NA  # If one of the sets is empty
  } else {
    fold_enrichment <- overlap / ((k * d) / N)
  }
  
  # Append results to results_df
  fisher_results_df <- rbind(
    fisher_results_df,
    data.frame(
      gene_set        = set_name,
      size_of_set     = k,
      size_of_target  = d,
      overlap         = overlap,
      fold_enrichment = round(fold_enrichment, 2),
      p_value         = round(fisher_res$p.value, 3),
      stringsAsFactors = FALSE
    )
  )
}

# 7. Correct for multiple testing (fdr)
fisher_results_df$padj_fdr <- p.adjust(fisher_results_df$p_value, method = "fdr")

# 8. Sort by adjusted p-value and inspect top results
fisher_results_df <- fisher_results_df[order(fisher_results_df$padj_fdr), ]
print(fisher_results_df)

write.table(fisher_results_df, file.path(output_folder,"wgcna_bronch_deg_overlap_fishers-exact_enrichment.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

###
# 8. Network Visualization of Eigengenes
### 

# Isolate blood_eos_log from the clinical traits
bal_eos_mt1=as.data.frame(datTraits_bal$bal_Eos_p_more_1);
bal_act=as.data.frame(datTraits_bal$asthma_phen_ACT.score);

names(bal_eos_mt1) = "eos_mt1"
names(bal_act) = "act"

# Add the BAL_eos_ct_log to existing module eigengenes
mergedMEs<-mergedMEs[rownames(mergedMEs)%in%rownames(datTraits_bal),]
MET = orderMEs(cbind(mergedMEs, bal_eos_mt1, bal_act ))
# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# Plot the dendrogram
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0, mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(10,10,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

#####
# 9. choose hub genes
#####

# Adjacency matrix
softPower <- 4 # decided previously in the original analysis
adjacency <- adjacency(expression.data_bal, power = softPower)

HubGenes <- chooseTopHubInEachModule(expression.data_bal, mergedColors, power = 4, type= "signed")

connectivity_allClusters <- intramodularConnectivity(adjacency, mergedColors$modules, scaleByMax = FALSE)

colnames(expression.data_bal)%>%head

# Check if rownames(connectivity_allClusters) match gene.module.table$genes
matched <- match(rownames(connectivity_allClusters), gene.module.table$genes)

# Update connectivity_allClusters$module with the corresponding module values
connectivity_allClusters$module <- gene.module.table$modules[matched]
connectivity_allClusters$gene <- gene.module.table$genes[matched]

top_kWithin_by_module <- connectivity_allClusters %>%
  group_by(module) %>%
  slice_max(order_by = kWithin, n = 20) %>%
  ungroup() %>%
  arrange(module, desc(kWithin))

head(top_kWithin_by_module)

top_kWithin_by_module_deg_overlap <- top_kWithin_by_module %>%
  filter(gene %in% deg_abs_lfc)

head(top_kWithin_by_module_deg_overlap)

write.table(top_kWithin_by_module_deg_overlap, file=file.path(output_folder,"top_kWithin_by_module_deg_overlap.txt"),   sep = "\t",
            row.names = TRUE, 
            col.names = NA)

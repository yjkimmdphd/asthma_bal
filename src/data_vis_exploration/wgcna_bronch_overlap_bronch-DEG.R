# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyverse)
library(WGCNA)

# =============== 
# synopsis: plots Fisher's exact test results (enrichment and FDR) of DEG overlap with the modules including ivory
# explores what are the connections between ivory and orange module
# ===============

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
# the phenotype table
# "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv" had
# a serious defect in which the ACT score was row-shifted, and erroneously had
# incorrect association between Eos % and ACT score
# file changed to 2025-02-14
# -----------------------------------------------------------------------------
phen_path <- file.path(
  "./resources/processed_data",
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2025-02-14.csv"
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

alltraits<-alltraits[,c("FEV1_percent","BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log",
                        "BAL_wbc_log", var_dichot_bal)] # only select FEV1 predicted % and cell counts in the phenotype data
good<-!(sapply(alltraits[,-1],is.na)%>%rowSums()>0) # remove rows with at least one NA in cell count phenotype
datTraits_bal<-alltraits[good,]

# subset alltraits for variables pertaining to blood phenotype
alltraits<-phen_input[phen_input$SampleID%in%rownames(mergedMEs),]
alltraits<-alltraits[alltraits$SampleID%in%cbc_sampleID,]
rownames(alltraits)<-alltraits$SampleID

alltraits_bld<-alltraits[,c("FEV1_percent",var_to_test_bld)] # only select FEV1 and cell counts in the phenotype data
good_bld<-!(sapply(alltraits_bld[,-1],is.na)%>%rowSums()>0) # remove rows with at least one NA in cell count phenotype
datTraits_bld<-alltraits_bld[good_bld,]

###
# 5. Module-Trait associations
###

gene.module.table<-data.frame(genes=colnames(expression.data),modules=mergedColors$modules)

###
# 6. load DEG and explor overlaps
###
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


#####
# 6. more exploration
# this is a work in progress. I want to eventually want to look at how orange and ivory are connected
#####

# Read the table (assuming it's tab-delimited)
wgcna_folder<-file.path("./reports/local_only/wgcna/bronch/output/")
df <- read.delim("./reports/local_only/wgcna/bronch/output/bronch_wgcna_genelist_batch12346.txt", header = TRUE, check.names = FALSE, sep="\t",row.names = 1)
head(df)
# Convert from wide to long format
df_long <- df %>%
  pivot_longer(cols = -1, names_to = "module_color", values_to = "gene") %>%
  select(gene, module_color)  # Ensure correct column order

# View the result
head(df_long)
write.table(df_long, "./reports/local_only/wgcna/bronch/output/bronch_wgcna_long_format.txt", sep = "\t", row.names = FALSE, quote = FALSE)

ivory_network<-read.delim(file.path(wgcna_folder,"modules_output","ivory.txt"),header = TRUE, check.names = FALSE, sep="\t")
ivory_network%>%filter(weight>0.3)

orange_network<-read.delim(file.path(wgcna_folder,"modules_output","orange.txt"),header = TRUE, check.names = FALSE, sep="\t")

gmt_ivory<-gene.module.table%>%filter(modules=="ivory")
colnames(gmt_ivory)<-c("Target","modules")

orange_ivory_network<-left_join(orange_network, gmt_ivory, by="Target")
orange_ivory_network%>%filter(modules=="ivory")%>%arrange(desc(Connectivity))


# ----------------------
# visualize overlap between Module genes and DEG genes
# ----------------------
# Load Fisher exact test results 
fisher_data<-read.table(file.path(wgcna_folder,"wgcna_bronch_deg_overlap_fishers-exact_enrichment.txt"),sep="\t",header=TRUE,row.names = 1)
fisher_data$gene_set<-rownames(fisher_data)
# Order gene sets by fold enrichment
fisher_data$gene_set <- factor(fisher_data$gene_set, levels = fisher_data$gene_set[order(fisher_data$fold_enrichment, decreasing = TRUE)])

# Create a new column for color classification
fisher_data$significance <- ifelse(fisher_data$padj_fdr <= 0.05, "Significant (≤ 0.05)", "Not Significant (> 0.05)")

# Create the ggplot with horizontal bars, p-value labels, and x-axis limits
ggplot(fisher_data%>%filter(padj_fdr<=0.05), aes(x = gene_set, y = fold_enrichment, fill = significance)) +
  geom_bar(stat = "identity") +
  # geom_text(aes(label = sprintf("%.3f", padj_fdr)), hjust = -0.2, size = 3) +  # Add p-value labels
  scale_fill_manual(values = c("Significant (≤ 0.05)" = "blue", "Not Significant (> 0.05)" = "blue")) +
  labs(title = "Fisher's Exact Test: Gene Set Enrichment",
       x = "Gene Set",
       y = "Fold Enrichment",
       fill = "Significance") +
  ylim(0, 40) +  # Set x-axis limits
  theme_minimal() +
  coord_flip() +  # Flip the coordinates for horizontal bars
  theme(plot.title = element_text(hjust = 0.5))  # Center the title


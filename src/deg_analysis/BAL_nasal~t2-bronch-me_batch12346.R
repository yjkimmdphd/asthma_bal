################################################################################
# nasal RNA-seq Analysis with Additional Samples (Batch 6)
# Reference Genome: GRCh38
# uses eigen genes identified by WGCNA of bronchial RNA seq experiment as variables
#
# This script performs:
#   - Data loading and subsetting
#   - Creation of categorical variables
#   - DESeq2 analysis
#   - Export of results and summary tables
################################################################################

library(tidyverse)
library(DESeq2)

#===============================================================================
# 1. Load Gene Count Table
#===============================================================================
countdata <- file.path(
  "./resources/raw_data/MS_asthma",
  "MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt"
) 

counts <- if (file.exists(countdata)) {
  read.delim(countdata, check.names = FALSE)
}

genes <- counts[, "SampleID"]

# Select bronchial samples
samples_of_interest  <- grepl("^N", colnames(counts))
counts_of_interest   <- counts[, samples_of_interest]
rownames(counts_of_interest) <- genes
counts.ID       <- colnames(counts_of_interest)

head(counts_of_interest)

#===============================================================================
# 2. Load Phenotype Data
#===============================================================================
phen_path <- file.path(
  "./resources/processed_data",
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2025-02-14.csv"
)
phen <- read.csv(phen_path)

phen<-phen %>% filter(BAL_eos_p >= 0 & BAL_neut_p >= 0) %>%
  mutate(comp1 = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 4 ~ "mixed",
                                  BAL_eos_p > 1 & BAL_neut_p <= 4 ~ "eos",
                                  BAL_eos_p <= 1 & BAL_neut_p > 4 ~ "neut",
                                  BAL_eos_p <= 1 & BAL_neut_p <= 4 ~ "pauci"), levels = c("pauci", "neut", "mixed", "eos")),
         comp2 = factor(case_when(BAL_eos_p > 1 ~ "high_eos",
                                  BAL_eos_p <= 1 ~ "low_eos"), levels = c("high_eos", "low_eos")))
phen$Batch<-factor(phen$Batch,levels=unique(phen$Batch))

# load eigen genes values for T2 genes expressed in bronchial WGCNA modules
wgcna_bronch_t2_ME_dir<-file.path("./reports/local_only/wgcna/bronch/output")
if(file.exists(file.path(wgcna_bronch_t2_ME_dir,"eigengene_df_t2.txt"))){
  wgcna_bronch_t2_ME<-read.table(file.path(wgcna_bronch_t2_ME_dir,"eigengene_df_t2.txt"),header = TRUE,sep = "\t",row.names = NULL)
}
wgcna_bronch_t2_ME<-wgcna_bronch_t2_ME%>%select(ID,Eigengene_Level)
colnames(wgcna_bronch_t2_ME)<-c("ID","T2_ME")

# load eigen genes values for each of the bronchial WGCNA modules
wgcna_bronch_ME_dir<-file.path("./reports/local_only/wgcna/bronch/output")
if(file.exists(file.path(wgcna_bronch_ME_dir,"mergedMEs.txt"))){
  wgcna_bronch_ME<-read.table(file.path(wgcna_bronch_ME_dir,"mergedMEs.txt"),header = TRUE,sep = "\t",row.names = 1)
}

# add eigen genes as the phenotype
wgcna_bronch_ME<-cbind(ID=gsub("B","",gsub("\\..*","",rownames(wgcna_bronch_ME)))%>%as.numeric(),wgcna_bronch_ME)

# add eigen genes as the phenotype
phen<-left_join(phen,wgcna_bronch_t2_ME,by="ID")%>%left_join(wgcna_bronch_ME,by="ID")

phen_bronch<-phen[grepl("^B",phen$SampleID),]
phen_nasal<-phen[grepl("^N",phen$SampleID),]
phen_nasal<-phen_nasal[-grep("F", phen_nasal$SampleID),]

phen_input<-phen_nasal
head(phen_input)
#===============================================================================
# 3. Subset Phenotype for nasal RNA-seq
#===============================================================================
sampleID_exists <- phen_input$SampleID %in% counts.ID
phen_input  <- phen_input[sampleID_exists, ]
phen <- phen_input


#===============================================================================
# 4. Load Custom DEG Functions
#===============================================================================
# Functions that should be available:
#   - filter_low_expressed_genes_method2()
#   - rowgenes_counttable()
#   - run_deseq2_DEG_analysis()
#   - get_DEG_results()
#   - generate_DEG_input_summary_table()
#   - generate_DEG_summary_table()

source("./src/function/deg_custom_functions_v2.R")


#===============================================================================
# 5. Identify ME Variables
#===============================================================================

# select for columns of interest
me_var_i<-c(grep("^T2",colnames(phen_input)),grep("plum2",colnames(phen_input)))
me_var<-colnames(phen_input)[me_var_i]

#===============================================================================
# 6. Prepare ColData for DESeq2
#===============================================================================

# ME as variables 
var_to_test <- me_var

# Results variable names that will be used for DESeq2 results() function: continuous or "variableTRUE" for categorical
var_to_test_res <- c(
  var_to_test
)

# Create a list indicating non-missing values for each variable to test
valid_samples <- lapply(phen[, var_to_test], function(data) !is.na(data))

# Initialize a named list to store data frames for DESeq2 input
df_var <- setNames(vector("list", length(var_to_test)), var_to_test)

# Populate the list with data frames containing relevant samples
for (i in seq_along(var_to_test)) {
  df_var[[i]] <- phen[valid_samples[[i]], c("SampleID", var_to_test[i], "Batch")]
}

# Print the number of samples in each data frame
cat("Sample counts before filtering:\n")
print(sapply(df_var, function(x) dim(x)[1]))


#===============================================================================
# 7. Filter Counts and Prepare Count Table
#===============================================================================
# We use df (all available data) as input for DESeq2.
df.input <- df_var

# Subset the RNA-seq counts by the phenotype
id   <- phen$SampleID
cols <- colnames(counts_of_interest) %in% id
ct   <- counts_of_interest[ , cols]

# Create a list of count tables—one for each tested variable subset
count.table <- lapply(df.input, function(tbl) {
  ct[ , colnames(ct) %in% tbl$SampleID]
})


#===============================================================================
# 8. Construct DESeq2 Design Formula
#===============================================================================
deg.design <- paste("~", var_to_test, "+ Batch")
print(deg.design)

# Prepare empty lists for DESeq2 objects, results, and significant results
dds     <- vector("list", length(var_to_test))
res     <- vector("list", length(var_to_test))
res.sig <- vector("list", length(var_to_test))

names(res)     <- deg.design
names(res.sig) <- deg.design


#===============================================================================
# 9. Run DESeq2 and Obtain Results
#===============================================================================
assay_index <- seq_along(deg.design)

for (i in assay_index) {
  # Run custom DESeq2 pipeline
  dds[[i]] <- run_deseq2_DEG_analysis(
    count.table[[i]],
    df.input[[i]],
    deg.design[i],
    deg.design[i]
  )
  
  # Filter genes with <10 counts across all samples
  dds_temp <- dds[[i]]
  keep     <- rowSums(counts(dds_temp)) >= 10
  dds_temp <- dds_temp[keep, ]
  
  # Get results and significant results
  res[[i]]     <- get_DEG_results(dds_temp, var_to_test_res[i])
  # Remove rows with NA in padj and keep only padj <= 0.05
  res.sig[[i]] <- res[[i]][!is.na(res[[i]]$padj) & res[[i]]$padj <= 0.05, ]
  
  head(res.sig[[i]])
}


#===============================================================================
# 10. Write Significant and All Results to Disk
#===============================================================================
deg.folder <- paste0("deg_bal_nasal~me_", Sys.Date())
deg.dir    <- file.path("./reports", "temporary", deg.folder)

if (!dir.exists(deg.dir)) {
  dir.create(deg.dir)
}

if (dir.exists(deg.dir)) {
  for (i in assay_index) {
    a <- res.sig[[i]]
    b <- res[[i]]
    
    write.csv(
      a,
      row.names = TRUE,
      file.path(
        deg.dir,
        paste(
          "deg_nasal_res_sig", i, deg.design[[i]], Sys.Date(), ".csv",
          sep = "_"
        )
      )
    )
    write.csv(
      b,
      row.names = TRUE,
      file.path(
        deg.dir,
        paste(
          "deg_nasal_res_all", i, deg.design[[i]], Sys.Date(), ".csv",
          sep = "_"
        )
      )
    )
  }
}


#===============================================================================
# 11. Summarize Data Input
#===============================================================================
generate_DEG_input_summary_table <- function(
    original_ct, filtered_ct, dds, res, des
) {
  filter_method    <- "DESEq2 + ct>=10"
  n_filtered_genes <- paste(
    "analyzed n_genes:", nrow(filtered_ct), ",",
    "filtered n_genes:", nrow(original_ct) - nrow(filtered_ct)
  )
  
  samples <- sapply(dds, function(d) {
    colData(d)$SampleID %>% paste(collapse = ",")
  })
  
  dds_names <- paste0("dds", seq_along(dds))
  res_names <- paste0("res", seq_along(res))
  design    <- des
  
  data.frame(
    dds               = dds_names,
    results           = res_names,
    design            = design,
    samples           = samples,
    filter_method     = filter_method,
    n_filtered_genes  = n_filtered_genes
  )
}

if (dir.exists(deg.dir)) {
  a <- generate_DEG_input_summary_table(
    counts_of_interest[ , cols],
    ct, dds, res, deg.design
  )
  
  write.csv(
    a,
    row.names = FALSE,
    file.path(
      deg.dir,
      paste(
        "deg_input_file",
        Sys.Date(),
        ".csv",
        sep = "_"
      )
    )
  )
}

#===============================================================================
# 12. Generate Summary of DEG Analysis
#===============================================================================
generate_DEG_summary_table <- function(
    results_significant, deg_design, variable
) {
  res.sig   <- results_significant
  reslist   <- paste0("res.sig", seq_along(res.sig))
  n_sig_deg <- unlist(lapply(res.sig, nrow))
  design    <- deg_design
  var       <- variable
  
  data.frame(
    type      = "nasal",
    results   = reslist,
    n_sig_deg = n_sig_deg,
    design    = design,
    variable  = var,
    row.names = NULL
  )
}

if (dir.exists(deg.dir)) {
  a <- generate_DEG_summary_table(res.sig, deg.design, var_to_test)
  
  write.csv(
    a,
    row.names = FALSE,
    file.path(
      deg.dir,
      paste(
        "deg_res_summary",
        Sys.Date(),
        ".csv",
        sep = "_"
      )
    )
  )
}


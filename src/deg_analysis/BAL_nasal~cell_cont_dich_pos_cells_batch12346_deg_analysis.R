################################################################################
# nasal RNA-seq Analysis with Additional Samples (Batch 6)
# Reference Genome: GRCh38
# Updated Phenotype Table
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

# Select nasal samples
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
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv"
)
phen <- read.csv(phen_path)

phen$Batch<-factor(phen$Batch,levels=unique(phen$Batch))

# load data for sampling date differences 
sampling_date_diff<-"./resources/processed_data/sampling_dates/swab-bal-cbc_differences_in_days.txt"
sampling_date_diff<-if(file.exists(sampling_date_diff)){read.table(sampling_date_diff,row.names = NULL,header = TRUE)}else{print("sampling date data doesn't exist")}
sampling_date_diff<-sampling_date_diff%>%filter(Comparison=="blood_bal")
colnames(sampling_date_diff)[1:3]<-c("ID","sampling_date_comp","sampling_date_diff_days")

# left join sampling date data and phenotype data 
phen<-left_join(phen,sampling_date_diff,by="ID")
phen_bronch<-phen[grepl("^B",phen$SampleID),]
phen_nasal<-phen[grepl("^N",phen$SampleID),]
phen_nasal<-phen_nasal[-grep("F", phen_nasal$SampleID),]

phen_input<-phen_nasal
#===============================================================================
# 3. Subset Phenotype for nasal RNA-seq
#===============================================================================
sampleID_exists <- phen_input$SampleID %in% counts.ID
phen_input  <- phen_input[sampleID_exists, ]

# Create categorical variables based on various thresholds
phen_input <- phen_input %>%
  mutate(
    bal_AEC_more_0   = BAL_eos_ct > 0,
    bal_AEC_more_1   = BAL_eos_ct > 1,
    bal_AEC_more_1.2 = BAL_eos_ct > 1.2,
    bal_AEC_more_3   = BAL_eos_ct > 3,
    bal_AEC_more_5   = BAL_eos_ct > 5,
    
    bal_Eos_p_more_0 = BAL_eos_p > 0,
    bal_Eos_p_more_1 = BAL_eos_p > 1,
    bal_Eos_p_more_3 = BAL_eos_p > 3,
    
    bal_ANC_more_0   = BAL_neut_ct > 0,
    bal_ANC_more_5   = BAL_neut_ct > 5,
    bal_ANC_more_13  = BAL_neut_ct > 13,
    
    bal_neut_p_more_0 = BAL_neut_p > 0,
    bal_neut_p_more_2 = BAL_neut_p > 2,
    bal_neut_p_more_5 = BAL_neut_p > 5,
    
    bld_AEC_more_0   = blood_eos > 0,
    bld_AEC_more_100 = blood_eos > 100,
    bld_AEC_more_300 = blood_eos > 300,
    bld_AEC_more_500 = blood_eos > 500
  )

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
# 5. Identify Cell-Count Variables
#===============================================================================
# Log-transformed/scaled/centered
source.cell.log <- c(
  "BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log",
  "blood_eos_log", "blood_eos_p_log", "blood_neut_log",
  "blood_neut_p_log"
)

# Categorical variables to test (BAL)
var_dichot_bal <- c(
  "bal_AEC_more_1", "bal_AEC_more_3", "bal_AEC_more_5",
  "bal_Eos_p_more_1", "bal_Eos_p_more_3",
  "bal_ANC_more_5", "bal_ANC_more_13",
  "bal_neut_p_more_2", "bal_neut_p_more_5"
)

# Categorical variables to test (Blood)
var_dichot_blood <- c(
  "bld_AEC_more_100",
  "bld_AEC_more_300", "bld_AEC_more_500"
)

# define fx to select for sampleIDs meeting specific criteria
sampleID_selection <- function(phen, var) {
  sid <- phen %>%
    filter(.data[[var]] > 0) %>%
    select("SampleID") %>%
    unlist() %>%
    as.vector()
  return(sid)
}


# find sampleIDs with postive cell % values 
pos_sampleIDs<-
  list(sampleID_selection(phen_input,"BAL_eos_p"), 
       sampleID_selection(phen_input,"BAL_neut_p"),
       sampleID_selection(phen_input,"blood_eos_p"),
       sampleID_selection(phen_input,"blood_neut_p"))
names(pos_sampleIDs)<-c("pos_BAL_eos_p","pos_BAL_neut_p","pos_blood_eos_p","pos_blood_neut_p")

lapply(pos_sampleIDs,length)
#===============================================================================
# 6. Prepare ColData for DESeq2
#===============================================================================

# Combine continuous (log counts) and categorical variables
var_to_test <- c(source.cell.log, var_dichot_bal, var_dichot_blood)

# identify blood cell count phenotype variables
var_to_test_bld<-var_to_test[c(grep("blood",var_to_test),grep("bld",var_to_test))]

# Results variable names that will be used for DESeq2 results() function: continuous or "variableTRUE" for categorical
var_to_test_res <- c(
  source.cell.log,
  paste0(c(var_dichot_bal, var_dichot_blood), "TRUE")
)

# Create a list indicating non-missing values for each variable to test
valid_samples <- lapply(phen[, var_to_test], function(data) !is.na(data))

# Initialize a named list to store data frames for DESeq2 input
df_var <- setNames(vector("list", length(var_to_test)), var_to_test)

# Populate the list with data frames containing relevant samples
for (i in seq_along(var_to_test)) {
  df_var[[i]] <- phen[valid_samples[[i]], c("SampleID", var_to_test[i], "Batch")]
}

bal_samples<-names(df_var)[c(grep("BAL",names(df_var)),grep("bal",names(df_var)))]

bal_eos_samples<-bal_samples[c(grep("eos",bal_samples),grep("Eos",bal_samples),grep("AEC",bal_samples))]
bal_neut_samples<-bal_samples[c(grep("neut",bal_samples),grep("ANC",bal_samples))]

bld_samples<-names(df_var)[c(grep("bld",names(df_var)),grep("blood",names(df_var)))]

bld_eos_samples<-bld_samples[c(grep("eos",bld_samples),grep("Eos",bld_samples),grep("AEC",bld_samples))]
bld_neut_samples<-bld_samples[c(grep("neut",bld_samples),grep("ANC",bld_samples))]

for (df in bal_eos_samples) {
  d <- df_var[[df]]  # Access the dataframe
  
  # Print the first few rows of the dataframe
  print(head(d))
  
  # Extract SampleIDs and print them
  sid_pre <- d$SampleID
  print(sid_pre)
  
  # Identify indices of positive samples and filter the dataframe
  sid_pre_index <- sid_pre %in% pos_sampleIDs[["pos_BAL_eos_p"]]
  filtered_data <- d[sid_pre_index, ]
  # Print the dimensions of the filtered dataframe
  print(dim(filtered_data))
  
  df_var[[df]]<-filtered_data
}

for (df in bal_neut_samples) {
  d <- df_var[[df]]  # Access the dataframe
   # Extract SampleIDs and print them
  sid_pre <- d$SampleID
  # Identify indices of positive samples and filter the dataframe
  sid_pre_index <- sid_pre %in% pos_sampleIDs[["pos_BAL_neut_p"]]
  filtered_data <- d[sid_pre_index, ]
  
  df_var[[df]]<-filtered_data
}


for (df in bld_eos_samples) {
  d <- df_var[[df]]  # Access the dataframe
  # Extract SampleIDs and print them
  sid_pre <- d$SampleID
  # Identify indices of positive samples and filter the dataframe
  sid_pre_index <- sid_pre %in% pos_sampleIDs[["pos_blood_eos_p"]]
  filtered_data <- d[sid_pre_index, ]
  
  df_var[[df]]<-filtered_data
}


for (df in bld_neut_samples) {
  d <- df_var[[df]]  # Access the dataframe
  # Extract SampleIDs and print them
  sid_pre <- d$SampleID
  # Identify indices of positive samples and filter the dataframe
  sid_pre_index <- sid_pre %in% pos_sampleIDs[["pos_blood_neut_p"]]
  filtered_data <- d[sid_pre_index, ]
  
  df_var[[df]]<-filtered_data
}

# Print the number of samples in each data frame
cat("Sample counts before filtering:\n")
print(sapply(df_var, function(x) dim(x)[1]))

# Identify samples with CBC data and sampling date difference < 1 year
cbc_sampleID <- phen %>%
  filter(!is.na(vars(var_to_test_bld)), abs(sampling_date_diff_days) < 365) %>%
  pull(SampleID)

# Filter blood cell count data frames based on the identified samples
blood_df <- df_var[var_to_test_bld]
blood_df_filtered <- lapply(blood_df, function(df) filter(df, SampleID %in% cbc_sampleID))

# Print sample counts before and after filtering
cat("Sample counts before filtering blood data:\n")
print(sapply(blood_df, function(x) dim(x)[1]))
cat("Sample counts after filtering blood data:\n")
print(sapply(blood_df_filtered, function(x) dim(x)[1]))

# Replace original blood cell count data with filtered data
df_var[var_to_test_bld] <- blood_df_filtered

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
deg.folder <- paste0("deg_bal_nasal~cell", Sys.Date())
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

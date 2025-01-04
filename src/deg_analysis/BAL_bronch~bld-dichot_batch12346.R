########################################
# Bronchial RNA-seq Analysis with New Samples
# Additional samples: batch 6
# Alignment: GRCh38
# Updated phenotype table
# select for samples that have BAL and cbc sampling dates that are close
########################################

library(tidyverse)
library(DESeq2)

#===============================================================================
# 1. Load cell count table
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
samples_of_interest <- grepl("^B", colnames(counts))
counts_of_interest  <- counts[, samples_of_interest]
rownames(counts_of_interest) <- genes
counts.ID      <- colnames(counts_of_interest)

head(counts_of_interest)

#===============================================================================
# 2. Load phenotype data
#===============================================================================
phenofile <- file.path(
  "./resources/processed_data",
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv"
)

phenotype <- if (file.exists(phenofile)) {
  read.csv(phenofile, row.names = NULL)
}

# Load sampling date differences and merge
sampling_date_diff <- "./resources/processed_data/sampling_dates/swab-bal-cbc_differences_in_days.txt"
sampling_date_diff <- if (file.exists(sampling_date_diff)) {
  read.table(sampling_date_diff, row.names = NULL, header = TRUE)
}

sampling_date_diff <- sampling_date_diff %>%
  filter(Comparison == "blood_bal")

colnames(sampling_date_diff)[1:3] <- c("ID", "sampling_date_comp", "sampling_date_diff_days")
phenotype <- left_join(phenotype, sampling_date_diff, by = "ID")

#===============================================================================
# 3. Subset phenotype for bronchial RNA-seq
#===============================================================================
sample_ids <- phenotype$SampleID %in% counts.ID
phen_input  <- phenotype[sample_ids, ]

# Create categorical variables for blood eos thresholds
phen_input <- phen_input %>%
  mutate(
    bld_AEC_more_150 = blood_eos > 150,
    bld_AEC_more_300 = blood_eos > 300,
    bld_AEC_more_500 = blood_eos > 500
  ) %>%
  filter(BAL_eos_p >= 0)

phen <- phen_input

#===============================================================================
# 4. Custom DEG functions
#===============================================================================
# Should provide or have loaded:
#   - filter_low_expressed_genes_method2
#   - rowgenes_counttable
#   - run_deseq2_DEG_analysis
#   - get_DEG_results
#   - generate_DEG_input_summary_table
#   - generate_DEG_summary_table
source("./src/function/deg_custom_functions_v2.R")

#===============================================================================
# 5. Identify relevant phenotype columns
#===============================================================================
source.cell.log <- c(
  "BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log",
  "BAL_wbc_log", "blood_eos_log", "blood_eos_p_log", "blood_neut_log",
  "blood_neut_p_log", "blood_wbc_log"
)

source.cell <- c(
  "BAL_eos_ct", "BAL_eos_p", "BAL_neut_ct", "BAL_neut_p", "BAL_wbc",
  "blood_eos", "blood_eos_p", "blood_neut", "blood_neut_p", "blood_wbc"
)

# Categorical variables to test - BAL
var_dichot_bal <- c(
  "bal_AEC_more_0", "bal_AEC_more_1", "bal_AEC_more_3", "bal_AEC_more_5",
  "bal_Eos_p_more_0", "bal_Eos_p_more_1", "bal_Eos_p_more_3",
  "bal_ANC_more_0", "bal_ANC_more_5", "bal_ANC_more_13",
  "bal_neut_p_more_0", "bal_neut_p_more_2", "bal_neut_p_more_5"
)

# Categorical variables to test - Blood
var_dichot_blood <- c(
  "bld_AEC_more_150",
  "bld_AEC_more_300",
  "bld_AEC_more_500"
)

#===============================================================================
# 6. Prepare variables and subsets
#===============================================================================
var_to_test      <- c(var_dichot_blood)
var_to_test_bld  <- var_to_test[
  c(grep("blood", var_to_test), grep("bld", var_to_test))
]
var_to_test_res  <- paste0(var_dichot_blood, "TRUE")

# Create colData list for DESeq2
pi <- lapply(phen[, var_to_test], function(x) !is.na(x))

df <- vector("list", length(var_to_test))
names(df) <- paste0(var_to_test, "_all")

for (i in seq_along(var_to_test)) {
  df[[i]] <- phen[pi[[i]], c("SampleID", var_to_test[i], "Batch")]
}
print(sapply(df, dim)[1, ])

# Filter for nasal-bal sampling matches with date differences < 1 year
cbc_sampleID <- phen_input %>%
  filter(
    !is.na(vars(var_to_test_bld)),
    abs(sampling_date_diff_days) < 365
  ) %>%
  pull(SampleID)

blood_df          <- df[paste0(var_to_test_bld, "_all")]
blood_df_filtered <- lapply(blood_df, function(tbl) {
  filter(tbl, SampleID %in% cbc_sampleID)
})

sapply(blood_df, dim)[1, ]          # samples before filtering
sapply(blood_df_filtered, dim)[1, ] # samples after filtering

df[paste0(var_to_test_bld, "_all")] <- blood_df_filtered

#===============================================================================
# 7. Filter and prepare count table
#===============================================================================
id   <- phen$SampleID
cols <- colnames(counts_of_interest) %in% id
ct   <- counts_of_interest[, cols]

count.table <- lapply(df, function(tbl) {
  ct[, colnames(ct) %in% tbl$SampleID]
})

#===============================================================================
# 8. Run DESeq2
#===============================================================================
deg.design <- paste("~", var_to_test, "+ Batch")
print(deg.design)

dds     <- vector("list", length = length(var_to_test))
res     <- vector("list", length = length(var_to_test))
res.sig <- vector("list", length = length(var_to_test))

names(res)     <- deg.design
names(res.sig) <- deg.design

assay_index <- seq_along(deg.design)

for (i in assay_index) {
  dds[[i]] <- run_deseq2_DEG_analysis(
    count.table[[i]], 
    df[[i]], 
    deg.design[i], 
    deg.design[i]
  )
  
  dds_temp <- dds[[i]]
  keep     <- rowSums(counts(dds_temp)) >= 10
  dds_temp <- dds_temp[keep, ]
  
  res[[i]]     <- get_DEG_results(dds_temp, var_to_test_res[i])
  res.sig[[i]] <- res[[i]][which(res[[i]]$padj <= 0.05), ]
  
  head(res.sig[[i]])
}

#===============================================================================
# 9. Save results
#===============================================================================
deg.folder <- paste0("deg_dich_ge0_", Sys.Date())
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
          "deg_bronch_res_sig", i, deg.design[[i]], Sys.Date(), ".csv",
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
          "deg_bronch_res_all", i, deg.design[[i]], Sys.Date(), ".csv",
          sep = "_"
        )
      )
    )
  }
}

#===============================================================================
# 10. Summarize data input
#===============================================================================
generate_DEG_input_summary_table <- function(original_ct, filtered_ct, dds, res, des) {
  filter_method    <- "TMM normalized LCPM cutoff"
  n_filtered_genes <- paste(
    "analyzed n_genes:", nrow(filtered_ct), ",",
    "filtered n_genes:", nrow(original_ct) - nrow(filtered_ct)
  )
  
  samples <- sapply(dds, function(d) {
    colData(d)$SampleID %>% paste(collapse = ",")
  })
  
  dds_names   <- paste0("dds", seq_along(dds))
  results     <- paste0("res", seq_along(res))
  design      <- des
  
  data.frame(
    dds                = dds_names,
    results            = results,
    design             = design,
    samples            = samples,
    filter_method      = filter_method,
    n_filtered_genes   = n_filtered_genes
  )
}

if (dir.exists(deg.dir)) {
  a <- generate_DEG_input_summary_table(
    bronch.counts[, cols],
    ct, dds, res, deg.design
  )
  write.csv(
    a,
    row.names = FALSE,
    file.path(
      deg.dir,
      paste(
        "deg_bronch_analysis_input_cellcount_cont_ge0_+Batch",
        Sys.Date(),
        ".csv",
        sep = "_"
      )
    )
  )
}

#===============================================================================
# 11. Summary table of the DEG analysis
#===============================================================================
generate_DEG_summary_table <- function(results_significant, deg_design, variable) {
  res.sig   <- results_significant
  reslist   <- paste0("res.sig", seq_along(res.sig))
  n_sig_deg <- unlist(lapply(res.sig, nrow))
  design    <- deg_design
  var       <- variable
  
  data.frame(
    type       = "bronch",
    results    = reslist,
    n_sig_deg  = n_sig_deg,
    design     = design,
    variable   = var,
    row.names  = NULL
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
        "dds_bronch_res_summary_cellcount_cont_ge0+Batch",
        Sys.Date(),
        ".csv",
        sep = "_"
      )
    )
  )
}

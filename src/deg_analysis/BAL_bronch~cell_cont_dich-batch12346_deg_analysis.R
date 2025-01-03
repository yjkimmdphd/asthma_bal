################################################################################
# Bronchial RNA-seq Analysis with Additional Samples (Batch 6)
# Reference Genome: GRCh38
# Updated Phenotype Table
# 
# This script loads new bronchial RNA-seq data, subsets relevant phenotype data,
# creates categorical variables for DEG analysis, and runs DESeq2.
################################################################################

library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)

#===============================================================================
# 1. Load Cell Count Table
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
bronch.samples  <- grepl("^B", colnames(counts))
bronch.counts   <- counts[, bronch.samples]
rownames(bronch.counts) <- genes
counts.ID       <- colnames(bronch.counts)

head(bronch.counts)

#===============================================================================
# 2. Load Phenotype Data
#===============================================================================
phenofile <- file.path(
  "./resources/processed_data",
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv"
)

phenotype <- if (file.exists(phenofile)) {
  read.csv(phenofile, row.names = NULL)
}

#===============================================================================
# 3. Subset Phenotype for Bronchial RNA-seq
#===============================================================================
bexist <- phenotype$SampleID %in% counts.ID
bphen  <- phenotype[bexist, ]

# Create categorical variables for different eosinophil and neutrophil thresholds
bphen <- bphen %>%
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

phen <- bphen

#===============================================================================
# 4. Custom Functions for DEG Analysis
#===============================================================================
# The following functions should be loaded from an external script:
#   - filter_low_expressed_genes_method2()
#   - rowgenes_counttable()
#   - run_deseq2_DEG_analysis()
#   - get_DEG_results()
#   - generate_DEG_input_summary_table()
#   - generate_DEG_summary_table()
source("./src/function/deg_custom_functions_v2.R")

#===============================================================================
# 5. Define Variables of Interest
#===============================================================================
# Log-transformed (scaled/centered) cell counts
source.cell.log <- c(
  "BAL_eos_ct_log",
  "BAL_eos_p_log",
  "BAL_neut_ct_log",
  "BAL_neut_p_log",
  "BAL_wbc_log",
  "blood_eos_log",
  "blood_eos_p_log",
  "blood_neut_log",
  "blood_neut_p_log",
  "blood_wbc_log"
)

# Original cell counts
source.cell <- c(
  "BAL_eos_ct",
  "BAL_eos_p",
  "BAL_neut_ct",
  "BAL_neut_p",
  "BAL_wbc",
  "blood_eos",
  "blood_eos_p",
  "blood_neut",
  "blood_neut_p",
  "blood_wbc"
)

# Categorical thresholds (BAL)
var_dichot_bal <- c(
  "bal_AEC_more_0", "bal_AEC_more_1", "bal_AEC_more_3", "bal_AEC_more_5",
  "bal_Eos_p_more_0", "bal_Eos_p_more_1", "bal_Eos_p_more_3",
  "bal_ANC_more_0", "bal_ANC_more_5", "bal_ANC_more_13",
  "bal_neut_p_more_0", "bal_neut_p_more_2", "bal_neut_p_more_5"
)

# Categorical thresholds (Blood)
var_dichot_blood <- c(
  "bld_AEC_more_0",
  "bld_AEC_more_100",
  "bld_AEC_more_300",
  "bld_AEC_more_500"
)

# Combine continuous (log counts) and categorical
var_to_test <- c(source.cell.log, var_dichot_bal, var_dichot_blood)

# Identify blood-related variables
var_to_test_bld <- var_to_test[
  c(grep("blood", var_to_test), grep("bld", var_to_test))
]

# Variables for results
var_to_test_res <- c(
  source.cell.log,
  paste0(c(var_dichot_bal, var_dichot_blood), "TRUE")
)

#===============================================================================
# 6. Prepare ColData for DESeq2
#===============================================================================
# 6A. All non-NA values
pi <- lapply(phen[, var_to_test], function(x) !is.na(x))
df <- vector("list", length(var_to_test))
names(df) <- paste0(var_to_test, "_all")

for (i in seq_along(var_to_test)) {
  df[[i]] <- phen[pi[[i]], c("SampleID", var_to_test[i], "Batch")]
}

print(sapply(df, dim)[1, ])

# 6B. Only values > 0
pi.pos <- lapply(phen[, var_to_test], function(x) which(x > 0))
df.pos <- vector("list", length(var_to_test))
names(df.pos) <- paste0(var_to_test, "_pos")

for (i in seq_along(var_to_test)) {
  df.pos[[i]] <- phen[pi.pos[[i]], c("SampleID", var_to_test[i], "Batch")]
}

print(sapply(df.pos, dim)[1, ])

#===============================================================================
# 7. Filter and Prepare Count Table
#===============================================================================
# We use df.pos for the DESeq2 colData, focusing on samples with counts > 0.
df.input <- df

# Subset RNA-seq counts to match phenotype
id   <- phen$SampleID
cols <- colnames(bronch.counts) %in% id
ct   <- bronch.counts[, cols]

# Optional stringent filtering (example):
# c2 <- filter_low_expressed_genes_method2(ct, round(length(id) * 0.1, 0))
# ct <- rowgenes_counttable(ct, c2)

# Create list of count tables
count.table <- lapply(df.input, function(tbl) {
  ct[, colnames(ct) %in% tbl$SampleID]
})

#===============================================================================
# 8. Build Design and Run DESeq2
#===============================================================================
deg.design <- paste0("~ ", var_to_test, " + Batch")
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
    df.input[[i]],
    deg.design[i],
    deg.design[i]
  )
  
  dds_temp <- dds[[i]]
  keep     <- rowSums(counts(dds_temp)) >= 10
  dds_temp <- dds_temp[keep, ]
  
  res[[i]]     <- get_DEG_results(dds_temp, var_to_test_res[i])
  res.sig[[i]] <- res[[i]][res[[i]]$padj <= 0.05, ]
  
  head(res.sig[[i]])
}

#===============================================================================
# 9. Save Significant and All Results
#===============================================================================
deg.folder <- paste0("deg_cont_dich_", Sys.Date())
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
# 10. Summarize Data Input
#===============================================================================
generate_DEG_input_summary_table <- function(
    original_ct, filtered_ct, dds, res, des
) {
  filter_method    <- "count >=10"
  n_filtered_genes <- paste(
    "analyzed n_genes:", nrow(filtered_ct), ",",
    "filtered n_genes:", nrow(original_ct) - nrow(filtered_ct)
  )
  
  samples <- sapply(dds, function(d) {
    colData(d)$SampleID %>% paste(collapse = ",")
  })
  
  dds_names   <- paste0("dds", seq_along(dds))
  res_names   <- paste0("res", seq_along(res))
  design      <- des
  
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
    bronch.counts[, cols],
    ct, dds, res, deg.design
  )
  write.csv(
    a,
    row.names = FALSE,
    file.path(
      deg.dir,
      paste(
        "deg_bronch_analysis_input_cellcount_cont_dich_+Batch",
        Sys.Date(),
        ".csv",
        sep = "_"
      )
    )
  )
}

#===============================================================================
# 11. Summary Table of the DEG Analysis
#===============================================================================
generate_DEG_summary_table <- function(results_significant, deg_design, variable) {
  res.sig   <- results_significant
  reslist   <- paste0("res.sig", seq_along(res.sig))
  n_sig_deg <- unlist(lapply(res.sig, nrow))
  design    <- deg_design
  var       <- variable
  
  data.frame(
    type      = "bronch",
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
        "dds_bronch_res_summary_cellcount_cont_dich+Batch",
        Sys.Date(),
        ".csv",
        sep = "_"
      )
    )
  )
}

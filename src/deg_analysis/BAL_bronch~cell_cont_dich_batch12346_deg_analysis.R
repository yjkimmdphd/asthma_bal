################################################################################
# Bronchial RNA-seq Analysis with Additional Samples (Batch 6)
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

# Select bronchial samples
samples_of_interest  <- grepl("^B", colnames(counts))
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

phen<-phen %>% filter(BAL_eos_p >= 0 & BAL_neut_p >= 0) %>%
  mutate(comp1 = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 4 ~ "mixed",
                                  BAL_eos_p > 1 & BAL_neut_p <= 4 ~ "eos",
                                  BAL_eos_p <= 1 & BAL_neut_p > 4 ~ "neut",
                                  BAL_eos_p <= 1 & BAL_neut_p <= 4 ~ "pauci"), levels = c("pauci", "neut", "mixed", "eos")),
         comp2 = factor(case_when(BAL_eos_p > 1 ~ "high_eos",
                                  BAL_eos_p <= 1 ~ "low_eos"), levels = c("high_eos", "low_eos")))
phen$Batch<-factor(phen$Batch,levels=unique(phen$Batch))

phen_bronch<-phen[grepl("^B",phen$SampleID),]
phen_nasal<-phen[grepl("^N",phen$SampleID),]
phen_nasal<-phen_nasal[-grep("F", phen_nasal$SampleID),]

phen_input<-phen_bronch
#===============================================================================
# 3. Subset Phenotype for Bronchial RNA-seq
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
  "BAL_wbc_log",    "blood_eos_log", "blood_eos_p_log", "blood_neut_log",
  "blood_neut_p_log", "blood_wbc_log"
)

# Original counts
source.cell <- c(
  "BAL_eos_ct", "BAL_eos_p", "BAL_neut_ct", "BAL_neut_p", "BAL_wbc",
  "blood_eos",  "blood_eos_p", "blood_neut", "blood_neut_p", "blood_wbc"
)

# Categorical variables to test (BAL)
var_dichot_bal <- c(
  "bal_AEC_more_0", "bal_AEC_more_1", "bal_AEC_more_3", "bal_AEC_more_5",
  "bal_Eos_p_more_0", "bal_Eos_p_more_1", "bal_Eos_p_more_3",
  "bal_ANC_more_0", "bal_ANC_more_5", "bal_ANC_more_13",
  "bal_neut_p_more_0", "bal_neut_p_more_2", "bal_neut_p_more_5","comp1"
)

# Categorical variables to test (Blood)
var_dichot_blood <- c(
  "bld_AEC_more_0",   "bld_AEC_more_100",
  "bld_AEC_more_300", "bld_AEC_more_500"
)

# Combine continuous (log counts) and categorical variables
var_to_test <- c(source.cell.log, var_dichot_bal, var_dichot_blood)

# Identify blood-related variables in var_to_test
var_to_test_bld <- var_to_test[
  c(grep("blood", var_to_test), grep("bld", var_to_test))
]

# Results variable names: continuous or "variableTRUE" for categorical
var_to_test_res <- c(
  source.cell.log,
  paste0(c(var_dichot_bal, var_dichot_blood), "TRUE")
)

#===============================================================================
# 6. Prepare ColData for DESeq2
#===============================================================================
# 6.1. Filter phenotype rows that are non-NA
pi <- lapply(phen[ , var_to_test], function(x) !is.na(x))
df <- vector("list", length(var_to_test))
names(df) <- paste0(var_to_test, "_all")

for (i in seq_along(var_to_test)) {
  df[[i]] <- phen[pi[[i]], c("SampleID", var_to_test[i], "Batch")]
}

print(sapply(df, dim)[1, ])

#===============================================================================
# 7. Filter Counts and Prepare Count Table
#===============================================================================
# We use df (all available data) as input for DESeq2.
df.input <- df

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
  res.sig[[i]] <- res[[i]][res[[i]]$padj <= 0.05, ]
  
  head(res.sig[[i]])
}

#===============================================================================
# 10. Write Significant and All Results to Disk
#===============================================================================
deg.folder <- paste0("deg_bal_bronch~cell", Sys.Date())
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
        "deg_res_summary",
        Sys.Date(),
        ".csv",
        sep = "_"
      )
    )
  )
}

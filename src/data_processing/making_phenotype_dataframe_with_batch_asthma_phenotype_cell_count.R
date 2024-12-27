library(dplyr)

# -----------------------------------------------------------------------------
# 1) Define file paths
# -----------------------------------------------------------------------------
processed_data_folder <- "./resources/processed_data/"
phenotype_file <- file.path(processed_data_folder, "yjk_Nasal_Biomarkers_BAL_2024-12-26.csv")
batch_info_file <- file.path(processed_data_folder, "nb_studyID_sampleID_batch.csv")

# -----------------------------------------------------------------------------
# 2) Read phenotype data
# -----------------------------------------------------------------------------
phenotype <- read.csv(phenotype_file)

# -----------------------------------------------------------------------------
# 3) Calculate BAL absolute eos and neut counts 
#    based on BAL WBC and eos% / neut%
# -----------------------------------------------------------------------------
phenotype <- phenotype %>%
  mutate(
    BAL_eos_ct  = BAL_eos_p  / 100 * BAL_wbc,
    BAL_neut_ct = BAL_neut_p / 100 * BAL_wbc
  )

# -----------------------------------------------------------------------------
# 4) Log-transform relevant columns
# -----------------------------------------------------------------------------
phenotype <- phenotype %>%
  mutate(
    BAL_eos_ct_log   = log10(BAL_eos_ct + 0.001)   %>% round(2),
    BAL_eos_p_log    = log10(BAL_eos_p + 0.001)    %>% round(2), 
    BAL_wbc_log      = log10(BAL_wbc + 0.001)      %>% round(2), 
    BAL_neut_p_log   = log10(BAL_neut_p + 0.001)   %>% round(2), 
    BAL_neut_ct_log  = log10(BAL_neut_ct + 0.001)  %>% round(2),
    blood_eos_log    = log10(blood_eos + 0.001)    %>% round(2),
    blood_eos_p_log  = log10(blood_eos_p + 0.001)  %>% round(2),
    blood_neut_log   = log10(blood_neut + 0.001)   %>% round(2),
    blood_neut_p_log = log10(blood_neut_p + 0.001) %>% round(2),
    blood_wbc_log    = log10(blood_wbc)            %>% round(2)
  )

# -----------------------------------------------------------------------------
# 5) Reorder columns in the data frame
# -----------------------------------------------------------------------------
desired_order <- c(
  "BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log",
  "BAL_wbc_log", "blood_eos_log", "blood_eos_p_log", "blood_neut_log",
  "blood_neut_p_log", "blood_wbc_log",
  "BAL_eos_ct", "BAL_eos_p", "BAL_neut_ct", "BAL_neut_p",
  "BAL_wbc", "blood_eos", "blood_eos_p", "blood_neut", "blood_neut_p", "blood_wbc"
)

df <- phenotype %>%
  dplyr::select(everything(), all_of(desired_order))

# -----------------------------------------------------------------------------
# 6) Center and scale relevant columns for model fitting
# -----------------------------------------------------------------------------
source.cell.log <- c(
  "BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log",
  "BAL_wbc_log", "blood_eos_log", "blood_eos_p_log", "blood_neut_log",
  "blood_neut_p_log", "blood_wbc_log"
)

df <- df %>%
  mutate_at(vars(all_of(source.cell.log)), scale)

# -----------------------------------------------------------------------------
# 7) Load batch info, subset by sampleID (nasal vs. bronchial),
#    and merge with phenotype data
# -----------------------------------------------------------------------------
batch_info <- read.csv(batch_info_file)

nasal_batch_info <- batch_info %>% filter(grepl("^N", SampleID))
bronch_batch_info <- batch_info %>% filter(grepl("^B", SampleID))

n_phenotype <- inner_join(nasal_batch_info, df, by = "ID")
b_phenotype <- inner_join(bronch_batch_info, df, by = "ID")

phenotype <- rbind(n_phenotype, b_phenotype)

# -----------------------------------------------------------------------------
# 8) Write final data to disk
# -----------------------------------------------------------------------------
output_file <- file.path(
  processed_data_folder,
  paste("scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_",today(),".csv",sep="")
)

write.csv(phenotype, output_file, row.names = FALSE)

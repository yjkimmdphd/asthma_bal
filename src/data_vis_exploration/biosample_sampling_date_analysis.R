# Load necessary libraries
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(ggrepel)

# Load data
# -----------------------------------------------------------------------------
# 1) Define file paths
# -----------------------------------------------------------------------------
processed_data_folder <- "./resources/processed_data/"
phenotype_file_biomarker_proc <- file.path(processed_data_folder, "Nasal_Biomarkers_BAL_2025-02-01.csv")
phenotype_file_biomarker_org <-file.path(processed_data_folder, "BiomarkersOfAsthma_original.csv")
batch_info_file <- file.path(processed_data_folder, "nb_studyID_sampleID_batch.csv")

# -----------------------------------------------------------------------------
# 2) Read phenotype data
# -----------------------------------------------------------------------------
phenotype_biomarker_proc <- read.csv(phenotype_file_biomarker_proc)
phenotype_biomarker_org <- read.csv(phenotype_file_biomarker_org)

org_index<-which(phenotype_biomarker_org$ID%in%phenotype_biomarker_proc$ID)
phenotype_biomarker_org<-phenotype_biomarker_org[org_index,]
phenotype<-left_join(phenotype_biomarker_proc,phenotype_biomarker_org,by="ID")

date_var<-c("nasal_swab_date","BAL_date", "blood_date","BAL_Date","Blood_draw_date") # `BAL_date` and `blood_date` contain the dates verified on EMR. `BAL_Date` and `Blood_draw_date` are the original dates, with some of the incorrect. 

data<-phenotype[,c("ID",date_var)]

# Convert date columns to proper Date type
data <- data %>%
  mutate(
    nasal_swab_date = mdy(nasal_swab_date, quiet = TRUE),
    BAL_date = mdy(BAL_date, quiet = TRUE),
    blood_date = mdy(blood_date, quiet = TRUE)
  ) %>%
  select(
    ID,
    nasal_swab_date,
    BAL_date, 
    blood_date
  )  # Remove original character columns

# Calculate pairwise date differences (in days)
data <- data %>%
  mutate(
    bal_nasal = as.numeric(difftime(BAL_date, nasal_swab_date, units = "days")),
    blood_nasal = as.numeric(difftime(blood_date, nasal_swab_date, units = "days")),
    blood_bal = as.numeric(difftime(blood_date, BAL_date, units = "days"))
  )

# Reshape data for plotting
plot_data <- data %>%
  select(ID, bal_nasal, blood_nasal, blood_bal) %>%
  pivot_longer(cols = c(bal_nasal, blood_nasal, blood_bal), 
               names_to = "Comparison", 
               values_to = "Difference") %>%
  arrange(Difference) %>%
  mutate(ID = factor(ID, levels = unique(ID)))

# Plot differences
ggplot(plot_data, aes(x = reorder(ID, Difference), y = Difference / 365, color = Comparison)) +
  geom_point(size = 3) +
  geom_hline(yintercept = c(-1, 1), color = "red") +
  geom_label_repel(data = filter(plot_data, Difference / 365 > 1),
                   aes(label = ID), nudge_y = 1, show.legend = FALSE) +
  geom_label_repel(data = filter(plot_data, Difference / 365 < -1),
                   aes(label = ID), nudge_y = -1, show.legend = FALSE) +
  scale_y_continuous(breaks = seq(-3, 3, by = 0.5)) +
  labs(
    title = "Pairwise Differences Between Sampling Dates (nasal swab, BAL, CBC)\nfor Each Study_ID",
    x = "Study_ID (sorted by difference)",
    y = "Difference in Years",
    color = "Comparison"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save output table
write.table(plot_data,
            "./resources/processed_data/sampling_dates/swab-bal-cbc_differences_in_days.txt",
            sep = "\t", row.names = FALSE)

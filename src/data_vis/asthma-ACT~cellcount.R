########
# data visualization for asthma phenotype and bal cell count
########
library(dplyr)
# load biomarker phenotype file saved in  'phenotype'
phenotype <-
  file.path(
    "./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv"
  )
phenotype_big <-
  file.path(
    "./resources/processed_data/Nasal_Biomarkers_BAL_transformed_with_raceinfo.csv"
  )
phen <-
  if (file.exists(phenotype)) {
    read.csv(phenotype, row.names = NULL)
  }
phen_big <-
  if (file.exists(phenotype_big)) {
    read.csv(phenotype_big, row.names = NULL)
  } %>% dplyr::select("ID", "Age_at_visit", "sex", "ethnicity", "Race")
phen <- left_join(phen, phen_big, by = "ID")
bphen <- phen %>% filter(grepl("^B", SampleID))


# Load necessary libraries
library(dplyr)
library(ggplot2)

# ACT 
# Define function to perform t-test, plot, and summarize data
perform_analysis <- function(data, cutoff, eos_var, label) {
  # Create a new variable indicating above/below cutoff
  data <- data %>%
    mutate(above_cutoff = get(eos_var) > cutoff) %>%
    group_by(above_cutoff) %>%
    arrange(above_cutoff) %>%
    dplyr::select(SampleID, above_cutoff, asthma_phen_ACT.score)
  
  # Perform t-test
  t_result <- t.test(asthma_phen_ACT.score ~ above_cutoff, data = data)
  
  # Plot boxplot with t-test p-value
  boxplot(
    asthma_phen_ACT.score ~ above_cutoff,
    data = data,
    xlab = label,
    ylab = "ACT score",
    main = paste(
      "ACT ~", label, "; t-test p-value",
      signif(t_result$p.value, digits = 2)
    )
  )
  
  # Median summary
  data %>%
    dplyr::summarize(median_ACT = median(asthma_phen_ACT.score, na.rm = TRUE)) %>%
    print()
}

# Apply analysis for different cutoffs
par(mfrow=c(1,2))
perform_analysis(bphen, 0, "BAL_eos_p", "BAL Eos% > 0%")
perform_analysis(bphen, 1, "BAL_eos_p", "BAL Eos% > 1%")
perform_analysis(bphen, 3, "BAL_eos_p", "BAL Eos% > 3%")
perform_analysis(bphen, 1, "BAL_eos_ct", "BAL AEC > 1")
perform_analysis(bphen, 1, "BAL_neut_ct", "BAL ANC > 1")
perform_analysis(bphen, 4, "BAL_neut_p", "BAL neut % > 4%")
perform_analysis(bphen, 6, "BAL_neut_p", "BAL neut % > 6%")
perform_analysis(bphen, 100, "blood_eos", "Blood AEC > 100")
perform_analysis(bphen, 300, "blood_eos", "Blood AEC > 300")
perform_analysis(bphen, 500, "blood_eos", "Blood AEC > 500")
######################################################################
# FEV1perc
# Define function to perform t-test, plot, and summarize data
perform_analysis <- function(data, cutoff, eos_var, label, lm_vars = c("sex", "Race")) {
  # Create a new variable indicating above/below cutoff
  data <- data %>%
    mutate(above_cutoff = get(eos_var) > cutoff) %>%
    group_by(above_cutoff) %>%
    arrange(above_cutoff) %>%
    dplyr::select(SampleID, above_cutoff, asthma_phen_FEV1_perc, sex, Race)
  
  # Perform t-test
  t_result <- t.test(asthma_phen_FEV1_perc ~ above_cutoff, data = data)
  
  # Plot boxplot with t-test p-value
  boxplot(
    asthma_phen_FEV1_perc ~ above_cutoff,
    data = data,
    xlab = label,
    ylab = "FEV1perc",
    main = paste(
      "FEV1perc ~", label, "; t-test p-value",
      signif(t_result$p.value, digits = 2)
    )
  )
  
  # Median summary
  data %>%
    dplyr::summarize(median_ACT = median(asthma_phen_FEV1_perc, na.rm = TRUE)) %>%
    print()
  
  # Linear model controlling for additional variables
  lm_formula <- as.formula(paste("asthma_phen_FEV1_perc ~ above_cutoff +", paste(lm_vars, collapse = " + ")))
  lm_result <- lm(lm_formula, data = data)
  print(summary(lm_result))
}
# Apply analysis for different cutoffs
par(mfrow=c(2,3))
perform_analysis(bphen, 1, "BAL_eos_p", "BAL Eos% > 1%")
perform_analysis(bphen, 3, "BAL_eos_p", "BAL Eos% > 3%")
perform_analysis(bphen, 1, "BAL_eos_ct", "BAL AEC > 1")
perform_analysis(bphen, 1, "BAL_neut_ct", "BAL ANC > 1")
perform_analysis(bphen, 4, "BAL_neut_p", "BAL neut % > 4%")
perform_analysis(bphen, 6, "BAL_neut_p", "BAL neut % > 6%")
perform_analysis(bphen, 100, "blood_eos", "Blood AEC > 100")
perform_analysis(bphen, 300, "blood_eos", "Blood AEC > 300")
perform_analysis(bphen, 500, "blood_eos", "Blood AEC > 500")

hist(bphen$asthma_phen_FEV1_perc,breaks=seq(0,150,by=10))
######################################################################
# lm ACT ~ mixed cell profiles # of note, LM is not appropriate because

# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(dplyr)

## Initialize lists to store results and plots
thresholds <- 1:9
lm_result_list <- vector("list", length(thresholds))
plot_list <- vector("list", length(thresholds))

## Define function to handle threshold-based analysis, save boxplot, and add sample counts
analyze_threshold <- function(threshold) {
  # Filter out NA values for the ACT score to avoid issues with missing data
  mixed <- bphen %>%
    filter(BAL_eos_p >= 0 & BAL_neut_p >= 0, !is.na(asthma_phen_ACT.score)) %>%
    mutate(type1 = factor(case_when(
      BAL_eos_p > 1 & BAL_neut_p > threshold ~ "mixed",
      BAL_eos_p > 1 & BAL_neut_p <= threshold ~ "eos",
      BAL_eos_p <= 1 & BAL_neut_p > threshold ~ "neut",
      BAL_eos_p <= 1 & BAL_neut_p <= threshold ~ "pauci"
    ), levels = c("pauci", "neut", "mixed", "eos")))
  
  # Calculate sample counts per level
  sample_counts <- mixed %>%
    group_by(type1) %>%
    summarise(count = n())
  
  # Dynamically position sample count labels within the y-axis range
  max_score <- max(mixed$asthma_phen_ACT.score, na.rm = TRUE) * 0.95
  
  # Generate boxplot with sample count annotations
  plot_list[[threshold]] <<- ggplot(mixed, aes(x = type1, y = asthma_phen_ACT.score)) +
    geom_boxplot(na.rm = TRUE) +
    geom_text(data = sample_counts, aes(x = type1, y = max_score, label = paste("n =", count)),
              color = "blue", size = 3, vjust = -0.5) +  # Position near top of plot
    labs(
      x = "Cell profile types",
      y = "ACT score",
      title = paste("ACT ~ BAL profile types; \nThreshold - BAL Eos% 1%, Neut", threshold, "%")
    ) +
    theme_minimal()
  
  # Linear model and save summary
  lm_model <- lm(asthma_phen_ACT.score ~ type1, data = mixed)
  lm_result_list[[threshold]] <<- summary(lm_model)
  
}

## Apply function for thresholds 1 to 9
lapply(thresholds, analyze_threshold)


#################################################################################

# AOV with Dunnett's post test pauci vs neut, mixed, eos

## Initialize lists to store results and plots
thresholds <- 1:9
aov_result_list <- vector("list", length(thresholds))
plot_list <- vector("list", length(thresholds))

## Define function to handle threshold-based analysis, save boxplot, and add sample counts
library(multcomp)
analyze_threshold <- function(threshold) {
  # Filter out NA values for the ACT score to avoid issues with missing data
  mixed <- bphen %>%
    filter(BAL_eos_p >= 0 & BAL_neut_p >= 0, !is.na(asthma_phen_ACT.score)) %>%
    mutate(type1 = factor(case_when(
      BAL_eos_p > 1 & BAL_neut_p > threshold ~ "mixed",
      BAL_eos_p > 1 & BAL_neut_p <= threshold ~ "eos",
      BAL_eos_p <= 1 & BAL_neut_p > threshold ~ "neut",
      BAL_eos_p <= 1 & BAL_neut_p <= threshold ~ "pauci"
    ), levels = c("pauci", "neut", "mixed", "eos")))
  
  # Calculate sample counts per level
  sample_counts <- mixed %>%
    group_by(type1) %>%
    summarise(count = n())
  
  # Dynamically position sample count labels within the y-axis range
  max_score <- max(mixed$asthma_phen_ACT.score, na.rm = TRUE) * 0.95
  
  # Generate boxplot with sample count annotations
  plot_list[[threshold]] <<- ggplot(mixed, aes(x = type1, y = asthma_phen_ACT.score)) +
    geom_boxplot(na.rm = TRUE) +
    geom_text(data = sample_counts, aes(x = type1, y = max_score, label = paste("n =", count)),
              color = "blue", size = 3, vjust = -0.5) +  # Position near top of plot
    labs(
      x = "Cell profile types",
      y = "ACT score",
      title = paste("ACT ~ BAL profile types; \nThreshold - BAL Eos% 1%, Neut", threshold, "%")
    ) +
    theme_minimal()
  
  # Analysis of variance with Dunnett's post hoc test
  aov_model <- aov(asthma_phen_ACT.score ~ type1, data = mixed)
  posthoc <- multcomp::glht(aov_model, linfct = multcomp::mcp(type1 = "Dunnett"))
  aov_result_list[[threshold]] <<- summary(posthoc)
}

## Apply function for thresholds 1 to 9
lapply(thresholds, analyze_threshold)

## Create a 3x3 multi-panel plot using grid.arrange
multi_panel_plot <- grid.arrange(
  grobs = plot_list,
  ncol = 3,  # 3 columns
  nrow = 3   # 3 rows
)

multi_panel_plot  # Display the 3x3 panel plot

#################################################################################
# AOV with Tukey's post test pauci, neut, mixed, eos

## Initialize lists to store results and plots
thresholds <- 1:9
aov_result_list <- vector("list", length(thresholds))
plot_list <- vector("list", length(thresholds))

## Define function to handle threshold-based analysis, save boxplot, and add sample counts
analyze_threshold <- function(threshold) {
  # Filter out NA values for the ACT score to avoid issues with missing data
  mixed <- bphen %>%
    filter(BAL_eos_p >= 0 & BAL_neut_p >= 0, !is.na(asthma_phen_ACT.score)) %>%
    mutate(type1 = factor(case_when(
      BAL_eos_p > 1 & BAL_neut_p > threshold ~ "mixed",
      BAL_eos_p > 1 & BAL_neut_p <= threshold ~ "eos",
      BAL_eos_p <= 1 & BAL_neut_p > threshold ~ "neut",
      BAL_eos_p <= 1 & BAL_neut_p <= threshold ~ "pauci"
    ), levels = c("pauci", "neut", "mixed", "eos")))
  
  # Calculate sample counts per level
  sample_counts <- mixed %>%
    group_by(type1) %>%
    summarise(count = n())
  
  # Dynamically position sample count labels within the y-axis range
  max_score <- max(mixed$asthma_phen_ACT.score, na.rm = TRUE) * 0.95
  
  # Generate boxplot with sample count annotations
  plot_list[[threshold]] <<- ggplot(mixed, aes(x = type1, y = asthma_phen_ACT.score)) +
    geom_boxplot(na.rm = TRUE) +
    geom_text(data = sample_counts, aes(x = type1, y = max_score, label = paste("n =", count)),
              color = "blue", size = 3, vjust = -0.5) +  # Position near top of plot
    labs(
      x = "Cell profile types",
      y = "ACT score",
      title = paste("ACT ~ BAL profile types; \nThreshold - BAL Eos% 1%, Neut", threshold, "%")
    ) +
    theme_minimal()
  
  # Analysis of variance
  aov_model <- aov(asthma_phen_ACT.score ~ type1, data = mixed)
  
  # Tukey's post hoc test
  posthoc <- TukeyHSD(aov_model)
  aov_result_list[[threshold]] <<- posthoc
}

## Apply function for thresholds 1 to 9
lapply(thresholds, analyze_threshold)


## Create a 3x3 multi-panel plot using grid.arrange
multi_panel_plot <- grid.arrange(
  grobs = plot_list,
  ncol = 3,  # 3 columns
  nrow = 3   # 3 rows
)

multi_panel_plot  # Display the 3x3 panel plot


#################################################################################



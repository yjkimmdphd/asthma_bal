########
# data visualization for asthma phenotype and cell count
########
library(dplyr)
# load biomarker phenotype file saved in  'phenotype'
phenotype <-
  file.path(
    "./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2025-02-14.csv" # updated to a new phenotype file. the one from 2024-12-26 had error: some rows of ACT score shifted, causing artifactual differences
  )
phen <-
  if (file.exists(phenotype)) {
    read.csv(phenotype, row.names = NULL)
  }
bphen <- phen %>% filter(grepl("^B", SampleID))


# Load necessary libraries
library(ggplot2)

# ACT 
var_phenotype<-"ACT_score"
# Define function to perform t-test, plot, and summarize data


perform_analysis <- function(data, cutoff, var, var_phenotype, label, output_dir = NULL) {
  # Create a new variable indicating above/below cutoff
  data <- data %>%
    mutate(above_cutoff = get(var) > cutoff) %>%
    group_by(above_cutoff) %>%
    arrange(above_cutoff) %>%
    dplyr::select(SampleID, above_cutoff, all_of(var_phenotype))
  
  # Perform t-test
  t_result <- t.test(data[[var_phenotype]] ~ above_cutoff, data = data)
  
  # Plot boxplot with t-test p-value
  boxplot(
    data[[var_phenotype]] ~ data$above_cutoff,
    xlab = label,
    ylab = var_phenotype,
    main = paste(
      var_phenotype, "~", label, "; t-test p-value:",
      signif(t_result$p.value, digits = 2)
    ),
    col = c("skyblue", "lightgreen"),
    border = "black"
  )
  
  # Save plot if output directory is provided
  if (!is.null(output_dir)) {
    # Sanitize filename
    sanitized_label <- gsub("[[:space:]>%]", "_", label)
    file_name <- paste0(output_dir, "/", sanitized_label, "_boxplot.png")
    
    # Save the plot
    png(filename = file_name, width = 800, height = 600)
    boxplot(
      data[[var_phenotype]] ~ data$above_cutoff,
      xlab = label,
      ylab = var_phenotype,
      main = paste(
        var_phenotype, "~", label, "; t-test p-value:",
        signif(t_result$p.value, digits = 2)
      ),
      col = c("skyblue", "lightgreen"),
      border = "black"
    )
    dev.off()
  }
  
  # Summarize t-test results
  result_summary <- tibble(
    label = label,
    variable = var,
    phenotype = var_phenotype,
    cutoff = cutoff,
    t_statistic = t_result$statistic,
    p_value = t_result$p.value,
    mean_above = t_result$estimate[2],
    mean_below = t_result$estimate[1],
    conf_low = t_result$conf.int[1],
    conf_high = t_result$conf.int[2],
    degrees_freedom = t_result$parameter
  )
  
  return(result_summary)
}


output_dir <- "./reports/figures/phenotype~cell_counts_2025-02-14"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Define parameters for multiple analyses
analysis_params <- list(
  list(cutoff = 0, var = "BAL_eos_p", label = "BAL Eos% > 0%"),
  list(cutoff = 1, var = "BAL_eos_p", label = "BAL Eos% > 1%"),
  list(cutoff = 3, var = "BAL_eos_p", label = "BAL Eos% > 3%"),
  list(cutoff = 1, var = "BAL_eos_ct", label = "BAL AEC > 1"),
  list(cutoff = 1, var = "BAL_neut_ct", label = "BAL ANC > 1"),
  list(cutoff = 4, var = "BAL_neut_p", label = "BAL neut % > 4%"),
  list(cutoff = 6, var = "BAL_neut_p", label = "BAL neut % > 6%"),
  list(cutoff = 100, var = "blood_eos", label = "Blood AEC > 100"),
  list(cutoff = 300, var = "blood_eos", label = "Blood AEC > 300"),
  list(cutoff = 500, var = "blood_eos", label = "Blood AEC > 500")
)
##############
# Perform all analyses for ACT and save results
##############
results <- lapply(analysis_params, function(params) {
  perform_analysis(
    data = bphen,
    cutoff = params$cutoff,
    var = params$var,
    label = params$label,
    var_phenotype = "ACT_score",
    output_dir = output_dir
  )
})

# Combine all results into a single data frame
results_df <- bind_rows(results)

# Save results to a CSV file
write.csv(results_df, file = "./reports/local_only/astham-phenotype~cellcount/t_test_results_ACT_2025-02-14.csv", row.names = FALSE)

##############
# Perform all analyses for FEV1 and save results
##############
results <- lapply(analysis_params, function(params) {
  perform_analysis(
    data = bphen,
    cutoff = params$cutoff,
    var = params$var,
    label = params$label,
    var_phenotype = "FEV1_percent",
    output_dir = output_dir
  )
})


# Combine all results into a single data frame
results_df <- bind_rows(results)
write.csv(results_df, file = "./reports/local_only/astham-phenotype~cellcount/t_test_results_asthma_phen_FEV1_perc_2025-02-14.csv", row.names = FALSE)

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



# Load required libraries
library(tidyverse)
library(ggplot2)
library(corrr)
library(psych)           # for ICC calculations and Cohen's kappa
library(caret)           # for confusion matrix
library(pROC)            # for ROC analysis
library(BlandAltmanLeh)  # for Bland-Altman plots
library(gridExtra)       # for combining plots
library(dplyr)           # for data manipulation

# Load your data (assuming you already have the 'phen' dataframe from your code)
phen <- read.csv("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2025-02-14.csv")

# Load sampling date differences
sampling_date_diff <- "./resources/processed_data/sampling_dates/swab-bal-cbc_differences_in_days.txt"
sampling_date_diff <- if(file.exists(sampling_date_diff)){read.table(sampling_date_diff, row.names = NULL, header = TRUE)}
sampling_date_diff <- sampling_date_diff %>% filter(Comparison == "blood_bal")
colnames(sampling_date_diff)[1:3] <- c("ID", "sampling_date_comp", "sampling_date_diff_days")

# Check if phen has ID column, if not create from SampleID
if (!"ID" %in% colnames(phen)) {
  phen$ID <- phen$SampleID
}

# Merge sampling date differences with phenotype data
phen_with_dates <- left_join(phen, sampling_date_diff, by = "ID")

# Filter data as in your original code + sampling date filter
phen <- phen_with_dates %>%
  filter(grepl("^B", SampleID), BAL_eos_p >= 0 & BAL_neut_p >= 0) %>%
  filter(!is.na(BAL_eos_p) & !is.na(blood_eos)) %>%
  filter(abs(sampling_date_diff_days) < 365)  # Filter to within 365 days

# Report sampling date effects
cat("=== SAMPLING DATE ANALYSIS ===\n")
cat("Total patients after date filtering (<30 days difference):", nrow(phen), "\n")
cat("Sampling date difference range:", min(phen$sampling_date_diff_days, na.rm = TRUE), 
    "to", max(phen$sampling_date_diff_days, na.rm = TRUE), "days\n")
cat("Mean sampling date difference:", round(mean(phen$sampling_date_diff_days, na.rm = TRUE), 1), "days\n")
cat("Median sampling date difference:", median(phen$sampling_date_diff_days, na.rm = TRUE), "days\n")

# Sampling date distribution
sampling_summary <- phen %>%
  summarise(
    n_total = n(),
    n_same_day = sum(sampling_date_diff_days == 0, na.rm = TRUE),
    n_within_7_days = sum(abs(sampling_date_diff_days) <= 7, na.rm = TRUE),
    n_within_14_days = sum(abs(sampling_date_diff_days) <= 14, na.rm = TRUE),
    n_within_30_days = sum(abs(sampling_date_diff_days) <= 30, na.rm = TRUE),
    pct_same_day = round(n_same_day/n_total*100, 1),
    pct_within_7_days = round(n_within_7_days/n_total*100, 1),
    pct_within_14_days = round(n_within_14_days/n_total*100, 1)
  )

print("Sampling Date Distribution:")
print(sampling_summary)

# Analyze correlation between sampling date difference and eosinophil difference
phen$eos_difference_raw <- abs(scale(phen$BAL_eos_p)[,1] - scale(phen$blood_eos)[,1])
phen$eos_difference_log <- abs(phen$BAL_eos_p_log - phen$blood_eos_log)

date_cor_raw <- cor(abs(phen$sampling_date_diff_days), phen$eos_difference_raw, use = "complete.obs")
date_cor_log <- cor(abs(phen$sampling_date_diff_days), phen$eos_difference_log, use = "complete.obs")

cat("\nSampling Date Impact on Concordance:\n")
cat("Correlation between date difference and eosinophil discordance (raw):", round(date_cor_raw, 3), "\n")
cat("Correlation between date difference and eosinophil discordance (log):", round(date_cor_log, 3), "\n")

# Create sampling date visualization
p_dates <- ggplot(phen, aes(x = abs(sampling_date_diff_days), y = eos_difference_log)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    title = "Sampling Date Difference vs Eosinophil Discordance",
    x = "Absolute Days Between Blood and BAL Sampling",
    y = "Absolute Difference in Log-Transformed Eosinophils",
    subtitle = paste("Correlation =", round(date_cor_log, 3))
  ) +
  theme_minimal()

print(p_dates)

# Define cutoffs for high eosinophil levels
BAL_CUTOFF <- 1.19    # BAL eosinophil % cutoff (you mentioned >1.19%)
BLOOD_CUTOFF <- 150   # Blood eosinophil absolute count cutoff

# Create categorical variables (fixing the error in your original code)
phen <- phen %>%
  mutate(
    bal_eos_high = ifelse(BAL_eos_p > BAL_CUTOFF, "High", "Normal"),
    blood_eos_high = ifelse(blood_eos > BLOOD_CUTOFF, "High", "Normal"),
    bal_eos_high_binary = as.numeric(BAL_eos_p > BAL_CUTOFF),
    blood_eos_high_binary = as.numeric(blood_eos > BLOOD_CUTOFF)
  )

# Check data availability
cat("Dataset Summary:\n")
cat("Total patients with both BAL and blood eosinophil data:", nrow(phen), "\n")
cat("BAL eosinophil range:", round(min(phen$BAL_eos_p, na.rm = TRUE), 2), "-", 
    round(max(phen$BAL_eos_p, na.rm = TRUE), 2), "%\n")
cat("Blood eosinophil range:", round(min(phen$blood_eos, na.rm = TRUE), 0), "-", 
    round(max(phen$blood_eos, na.rm = TRUE), 0), "\n\n")

print("=== CATEGORICAL CONCORDANCE ANALYSIS ===")

# Cross-tabulation
cross_tab <- table(phen$bal_eos_high, phen$blood_eos_high)
print("Cross-tabulation (BAL vs Blood):")
print(cross_tab)
print("Cross-tabulation with margins:")
print(addmargins(cross_tab))

# Calculate agreement measures
total_n <- nrow(phen)
observed_agreement <- sum(diag(cross_tab)) / total_n

# Cohen's Kappa
kappa_result <- psych::cohen.kappa(cbind(phen$bal_eos_high_binary, phen$blood_eos_high_binary))
kappa_value <- kappa_result$kappa

cat("\nObserved Agreement:", round(observed_agreement * 100, 1), "%\n")
cat("Cohen's Kappa:", round(kappa_value, 3), "\n")

# Confusion Matrix with detailed statistics
conf_matrix <- confusionMatrix(
  factor(phen$blood_eos_high, levels = c("Normal", "High")), 
  factor(phen$bal_eos_high, levels = c("Normal", "High")),
  positive = "High"
)
print(conf_matrix)

# Calculate sensitivity, specificity, PPV, NPV manually
TP <- cross_tab["High", "High"]
TN <- cross_tab["Normal", "Normal"] 
FP <- cross_tab["Normal", "High"]
FN <- cross_tab["High", "Normal"]

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
ppv <- TP / (TP + FP)
npv <- TN / (TN + FN)

print("\n=== DIAGNOSTIC PERFORMANCE ===")
cat("Sensitivity (BAL detecting high blood eos):", round(sensitivity * 100, 1), "%\n")
cat("Specificity (BAL detecting normal blood eos):", round(specificity * 100, 1), "%\n")
cat("Positive Predictive Value:", round(ppv * 100, 1), "%\n")
cat("Negative Predictive Value:", round(npv * 100, 1), "%\n")

# 1. CORRELATION ANALYSIS
print("\n=== CORRELATION ANALYSIS ===")

# Pearson correlation 1: Raw values
cor_pearson_raw <- cor(phen$BAL_eos_p, phen$blood_eos, method = "pearson", use = "complete.obs")
cor_test_p_raw <- cor.test(phen$BAL_eos_p, phen$blood_eos, method = "pearson")

# Pearson correlation 2: Log-transformed values
cor_pearson_log <- cor(phen$BAL_eos_p_log, phen$blood_eos_log, method = "pearson", use = "complete.obs")
cor_test_p_log <- cor.test(phen$BAL_eos_p_log, phen$blood_eos_log, method = "pearson")

# Spearman correlation (non-parametric) using raw values
cor_spearman <- cor(phen$BAL_eos_p, phen$blood_eos, method = "spearman", use = "complete.obs")
cor_test_s <- suppressWarnings(cor.test(phen$BAL_eos_p, phen$blood_eos, method = "spearman"))

cat("Pearson correlation (raw values):", round(cor_pearson_raw, 3), "\n")
cat("Pearson p-value (raw):", round(cor_test_p_raw$p.value, 4), "\n")
cat("Pearson correlation (log-transformed):", round(cor_pearson_log, 3), "\n")
cat("Pearson p-value (log):", round(cor_test_p_log$p.value, 4), "\n")
cat("Spearman correlation (raw values):", round(cor_spearman, 3), "\n")
cat("Spearman p-value:", round(cor_test_s$p.value, 4), "\n")

# Compare the correlations
cat("\nCorrelation Comparison:\n")
cat("Raw vs Log-transformed Pearson difference:", round(abs(cor_pearson_log - cor_pearson_raw), 3), "\n")
if (abs(cor_pearson_log) > abs(cor_pearson_raw)) {
  cat("Log transformation STRENGTHENED the correlation\n")
} else if (abs(cor_pearson_log) < abs(cor_pearson_raw)) {
  cat("Log transformation WEAKENED the correlation\n") 
} else {
  cat("Log transformation had minimal effect on correlation\n")
}

# 2. SCATTER PLOT WITH REGRESSION LINE AND CUTOFFS
p1 <- ggplot(phen, aes(x = BAL_eos_p, y = blood_eos)) +
  geom_point(aes(color = interaction(bal_eos_high, blood_eos_high)), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_vline(xintercept = BAL_CUTOFF, color = "red", linetype = "solid", alpha = 0.7) +
  geom_hline(yintercept = BLOOD_CUTOFF, color = "blue", linetype = "solid", alpha = 0.7) +
  annotate("text", x = BAL_CUTOFF + 0.1, y = max(phen$blood_eos) * 0.9, 
           label = paste("BAL cutoff =", BAL_CUTOFF, "%"), color = "red", hjust = 0) +
  annotate("text", x = max(phen$BAL_eos_p) * 0.1, y = BLOOD_CUTOFF + 20, 
           label = paste("Blood cutoff =", BLOOD_CUTOFF), color = "blue", vjust = 0) +
  scale_color_manual(
    name = "Classification",
    values = c("Normal.Normal" = "green", "Normal.High" = "orange", 
               "High.Normal" = "purple", "High.High" = "red"),
    labels = c("Both Normal", "High Blood Only", "High BAL Only", "Both High")
  ) +
  labs(
    title = "BAL Eosinophil % vs Blood Eosinophil Count with Cutoffs",
    x = "BAL Eosinophil (%)",
    y = "Blood Eosinophil Absolute Count",
    subtitle = paste("Pearson r (raw) =", round(cor_pearson_raw, 3), 
                     "; Pearson r (log) =", round(cor_pearson_log, 3), 
                     ", Agreement =", round(observed_agreement * 100, 1), "%",
                     ", n =", total_n)
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)

# 3. STANDARDIZED VALUES COMPARISON
# Compare both raw and log-transformed correlations
cat("\nCorrelation Summary:\n")
cat("Raw Pearson correlation:", round(cor_pearson_raw, 3), "\n")
cat("Log-transformed Pearson correlation:", round(cor_pearson_log, 3), "\n")
cat("Spearman correlation (rank-based):", round(cor_spearman, 3), "\n")

# 4. BLAND-ALTMAN ANALYSIS (for standardized values)
print("\n=== BLAND-ALTMAN ANALYSIS ===")

phen$diff_std <- phen$BAL_eos_p_log - phen$blood_eos_log
phen$mean_std <- (phen$BAL_eos_p_log + phen$blood_eos_log) / 2

mean_diff <- mean(phen$diff_std, na.rm = TRUE)
sd_diff <- sd(phen$diff_std, na.rm = TRUE)
upper_loa <- mean_diff + 1.96 * sd_diff
lower_loa <- mean_diff - 1.96 * sd_diff

cat("Mean difference:", round(mean_diff, 3), "\n")
cat("95% Limits of Agreement:", round(lower_loa, 3), "to", round(upper_loa, 3), "\n")

# Bland-Altman plot
p2 <- ggplot(phen, aes(x = mean_std, y = diff_std)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = mean_diff, color = "blue", linetype = "solid") +
  geom_hline(yintercept = upper_loa, color = "red", linetype = "dashed") +
  geom_hline(yintercept = lower_loa, color = "red", linetype = "dashed") +
  labs(
    title = "Bland-Altman Plot (Standardized Values)",
    x = "Mean of Standardized Values",
    y = "Difference (BAL - Blood)",
    subtitle = paste("Mean diff =", round(mean_diff, 3), 
                     "; 95% LoA:", round(lower_loa, 3), "to", round(upper_loa, 3))
  ) +
  theme_minimal()

print(p2)

# 5. BIDIRECTIONAL ROC CURVE ANALYSIS
print("\n=== BIDIRECTIONAL ROC CURVE ANALYSIS ===")

# ROC curve 1: BAL predicting high blood eosinophils
roc_bal_to_blood <- roc(phen$blood_eos_high_binary, phen$BAL_eos_p, quiet = TRUE)
auc_bal_to_blood <- auc(roc_bal_to_blood)

cat("AUC for BAL predicting high blood eosinophils:", round(auc_bal_to_blood, 3), "\n")

# ROC curve 2: Blood predicting high BAL eosinophils
roc_blood_to_bal <- roc(phen$bal_eos_high_binary, phen$blood_eos, quiet = TRUE)
auc_blood_to_bal <- auc(roc_blood_to_bal)

cat("AUC for Blood predicting high BAL eosinophils:", round(auc_blood_to_bal, 3), "\n")

# Plot both ROC curves
library(gridExtra)

p3a <- ggroc(roc_bal_to_blood) +
  geom_abline(intercept = 1, slope = 1, color = "red", linetype = "dashed") +
  labs(
    title = "BAL >>Blood Eosinophils",
    subtitle = paste("AUC =", round(auc_bal_to_blood, 3)),
    x = "Specificity",
    y = "Sensitivity"
  ) +
  theme_minimal()

p3b <- ggroc(roc_blood_to_bal) +
  geom_abline(intercept = 1, slope = 1, color = "red", linetype = "dashed") +
  labs(
    title = "Blood >> BAL Eosinophils",
    subtitle = paste("AUC =", round(auc_blood_to_bal, 3)),
    x = "Specificity",
    y = "Sensitivity"
  ) +
  theme_minimal()

# Print individual plots first
print(p3a)
print(p3b)

# Then create combined plot
p3_combined <- grid.arrange(p3a, p3b, ncol = 2, 
                            top = "Bidirectional ROC Analysis")

print(p3_combined)

# Find optimal cutoffs using Youden's index for both directions
optimal_coords_bal <- coords(roc_bal_to_blood, "best", ret = c("threshold", "sensitivity", "specificity"))
optimal_cutoff_bal <- as.numeric(optimal_coords_bal["threshold"])
youden_sens_bal <- as.numeric(optimal_coords_bal["sensitivity"])
youden_spec_bal <- as.numeric(optimal_coords_bal["specificity"])

optimal_coords_blood <- coords(roc_blood_to_bal, "best", ret = c("threshold", "sensitivity", "specificity"))
optimal_cutoff_blood <- as.numeric(optimal_coords_blood["threshold"])
youden_sens_blood <- as.numeric(optimal_coords_blood["sensitivity"])
youden_spec_blood <- as.numeric(optimal_coords_blood["specificity"])

cat("\nOptimal Cutoffs (Youden's Index):\n")
cat("BAL cutoff for predicting blood eos:", round(optimal_cutoff_bal, 3), "%\n")
cat("  - Sensitivity:", round(youden_sens_bal, 3), ", Specificity:", round(youden_spec_bal, 3), "\n")
cat("Blood cutoff for predicting BAL eos:", round(optimal_cutoff_blood, 0), "\n")
cat("  - Sensitivity:", round(youden_sens_blood, 3), ", Specificity:", round(youden_spec_blood, 3), "\n")

# 6. CATEGORICAL VISUALIZATION
p4 <- ggplot(phen, aes(x = bal_eos_high, fill = blood_eos_high)) +
  geom_bar(position = "dodge", alpha = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(
    title = "Distribution of High vs Normal Classifications",
    x = "BAL Eosinophil Classification",
    y = "Count",
    fill = "Blood Eosinophil Classification"
  ) +
  theme_minimal()

print(p4)

# 7. BIDIRECTIONAL DIAGNOSTIC PERFORMANCE
print("\n=== BIDIRECTIONAL DIAGNOSTIC PERFORMANCE ===")

# Direction 1: BAL predicting Blood (already calculated above)
cat("Direction 1: BAL → Blood Eosinophils\n")
cat("Using BAL >", BAL_CUTOFF, "% to predict Blood >", BLOOD_CUTOFF, "\n")
cat("Sensitivity:", round(sensitivity * 100, 1), "% (BAL correctly identifies high blood eos)\n")
cat("Specificity:", round(specificity * 100, 1), "% (BAL correctly identifies normal blood eos)\n")
cat("PPV:", round(ppv * 100, 1), "% (if BAL high, probability blood is high)\n")
cat("NPV:", round(npv * 100, 1), "% (if BAL normal, probability blood is normal)\n\n")

# Direction 2: Blood predicting BAL
# Create confusion matrix for Blood → BAL direction
conf_matrix_blood_to_bal <- confusionMatrix(
  factor(phen$bal_eos_high, levels = c("Normal", "High")), 
  factor(phen$blood_eos_high, levels = c("Normal", "High")),
  positive = "High"
)

# Extract values for Blood → BAL
TP_b2bal <- cross_tab["High", "High"]    # Blood high, BAL high
TN_b2bal <- cross_tab["Normal", "Normal"] # Blood normal, BAL normal
FP_b2bal <- cross_tab["High", "Normal"]   # Blood high, BAL normal  
FN_b2bal <- cross_tab["Normal", "High"]   # Blood normal, BAL high

sensitivity_b2bal <- TP_b2bal / (TP_b2bal + FN_b2bal)
specificity_b2bal <- TN_b2bal / (TN_b2bal + FP_b2bal)
ppv_b2bal <- TP_b2bal / (TP_b2bal + FP_b2bal)
npv_b2bal <- TN_b2bal / (TN_b2bal + FN_b2bal)

cat("Direction 2: Blood → BAL Eosinophils\n")
cat("Using Blood >", BLOOD_CUTOFF, " to predict BAL >", BAL_CUTOFF, "%\n")
cat("Sensitivity:", round(sensitivity_b2bal * 100, 1), "% (Blood correctly identifies high BAL eos)\n")
cat("Specificity:", round(specificity_b2bal * 100, 1), "% (Blood correctly identifies normal BAL eos)\n")
cat("PPV:", round(ppv_b2bal * 100, 1), "% (if Blood high, probability BAL is high)\n")
cat("NPV:", round(npv_b2bal * 100, 1), "% (if Blood normal, probability BAL is normal)\n\n")

cat("Comparison Summary:\n")
cat("BAL → Blood: Sens =", round(sensitivity * 100, 1), "%, Spec =", round(specificity * 100, 1), "%\n")
cat("Blood → BAL: Sens =", round(sensitivity_b2bal * 100, 1), "%, Spec =", round(specificity_b2bal * 100, 1), "%\n")
print("\n=== SUMMARY STATISTICS BY CLASSIFICATION ===")

summary_by_group <- phen %>%
  group_by(bal_eos_high, blood_eos_high) %>%
  summarise(
    n = n(),
    bal_eos_mean = mean(BAL_eos_p, na.rm = TRUE),
    bal_eos_sd = sd(BAL_eos_p, na.rm = TRUE),
    blood_eos_mean = mean(blood_eos, na.rm = TRUE),
    blood_eos_sd = sd(blood_eos, na.rm = TRUE),
    .groups = 'drop'
  )

print(as.data.frame(summary_by_group))

# 9. ENHANCED RESULTS SUMMARY WITH BIDIRECTIONAL ANALYSIS
# Debug: Check lengths first
analysis_items <- c(
  "Sample Size",
  "Pearson Correlation", 
  "Spearman Correlation", 
  "Observed Agreement (%)", 
  "Cohen's Kappa",
  "BAL → Blood: Sensitivity (%)",
  "BAL → Blood: Specificity (%)",
  "BAL → Blood: PPV (%)",
  "BAL → Blood: NPV (%)",
  "BAL → Blood: AUC",
  "Blood → BAL: Sensitivity (%)",
  "Blood → BAL: Specificity (%)",
  "Blood → BAL: PPV (%)",
  "Blood → BAL: NPV (%)",
  "Blood → BAL: AUC",
  "Optimal BAL Cutoff (%)",
  "Optimal Blood Cutoff"
)

value_items <- c(
  total_n,
  round(cor_pearson, 3), 
  round(cor_spearman, 3), 
  round(observed_agreement * 100, 1),
  round(kappa_value, 3),
  round(sensitivity * 100, 1),
  round(specificity * 100, 1),
  round(ppv * 100, 1),
  round(npv * 100, 1),
  round(auc_bal_to_blood, 3),
  round(sensitivity_b2bal * 100, 1),
  round(specificity_b2bal * 100, 1),
  round(ppv_b2bal * 100, 1),
  round(npv_b2bal * 100, 1),
  round(auc_blood_to_bal, 3),
  round(optimal_cutoff_bal, 3),
  round(optimal_cutoff_blood, 0)
)

# Handle p-values more carefully
p_val_bal_blood <- tryCatch({
  if(roc_bal_to_blood$p.value < 0.001) "<0.001" else round(roc_bal_to_blood$p.value, 4)
}, error = function(e) "N/A")

p_val_blood_bal <- tryCatch({
  if(roc_blood_to_bal$p.value < 0.001) "<0.001" else round(roc_blood_to_bal$p.value, 4)
}, error = function(e) "N/A")

kappa_p_val <- tryCatch({
  if(is.na(kappa_result$p.value)) "N/A" else round(kappa_result$p.value, 4)
}, error = function(e) "N/A")

p_value_items <- c(
  "N/A",
  round(cor_test_p$p.value, 4), 
  round(cor_test_s$p.value, 4), 
  "N/A",
  kappa_p_val,
  "N/A", 
  "N/A", 
  "N/A", 
  "N/A",
  p_val_bal_blood,
  "N/A", 
  "N/A", 
  "N/A", 
  "N/A",
  p_val_blood_bal,
  "N/A",
  "N/A"
)

# Check lengths before creating data frame
cat("Length check:\n")
cat("Analysis items:", length(analysis_items), "\n")
cat("Value items:", length(value_items), "\n") 
cat("P-value items:", length(p_value_items), "\n")

# Create results summary only if lengths match
if(length(analysis_items) == length(value_items) && length(value_items) == length(p_value_items)) {
  results_summary <- data.frame(
    Analysis = analysis_items,
    Value = value_items,
    P_value = p_value_items,
    stringsAsFactors = FALSE
  )
} else {
  cat("ERROR: Length mismatch detected. Cannot create results summary.\n")
  # Print each vector to debug
  cat("Analysis vector:\n")
  print(analysis_items)
  cat("Value vector:\n") 
  print(value_items)
  cat("P-value vector:\n")
  print(p_value_items)
}

print("\n=== ENHANCED RESULTS SUMMARY WITH BIDIRECTIONAL ANALYSIS ===")
print(results_summary)

# 10. INTERPRETATION OF BIDIRECTIONAL RESULTS
print("\n=== INTERPRETATION OF BIDIRECTIONAL RESULTS ===")

cat("Clinical Decision Making:\n")
cat("========================\n")
if (auc_bal_to_blood > auc_blood_to_bal) {
  cat("BAL is BETTER at predicting blood eosinophils (AUC:", round(auc_bal_to_blood, 3), "vs", round(auc_blood_to_bal, 3), ")\n")
  cat("→ BAL may be more informative for systemic eosinophilia assessment\n")
} else if (auc_blood_to_bal > auc_bal_to_blood) {
  cat("Blood is BETTER at predicting BAL eosinophils (AUC:", round(auc_blood_to_bal, 3), "vs", round(auc_bal_to_blood, 3), ")\n")
  cat("→ Blood test may be sufficient for airway eosinophilia assessment\n")
} else {
  cat("Both directions show similar predictive ability (AUC:", round(auc_bal_to_blood, 3), "vs", round(auc_blood_to_bal, 3), ")\n")
  cat("→ Methods provide equivalent information\n")
}

cat("\nSensitivity Comparison:\n")
if (sensitivity > sensitivity_b2bal) {
  cat("BAL has higher sensitivity for detecting systemic eosinophilia\n")
  cat("→ BAL is better at catching patients with high blood eosinophils\n")
} else if (sensitivity_b2bal > sensitivity) {
  cat("Blood has higher sensitivity for detecting airway eosinophilia\n")
  cat("→ Blood test is better at catching patients with high BAL eosinophils\n")
} else {
  cat("Both methods have similar sensitivity\n")
}

cat("\nSpecificity Comparison:\n")
if (specificity > specificity_b2bal) {
  cat("BAL has higher specificity for ruling out systemic eosinophilia\n")
  cat("→ BAL is better at correctly identifying patients with normal blood eosinophils\n")
} else if (specificity_b2bal > specificity) {
  cat("Blood has higher specificity for ruling out airway eosinophilia\n")
  cat("→ Blood test is better at correctly identifying patients with normal BAL eosinophils\n")
} else {
  cat("Both methods have similar specificity\n")
}

# 11. CROSS-TABULATION DETAILS
print("\n=== CROSS-TABULATION SUMMARY ===")
cat("Total patients:", total_n, "\n")
cat("BAL High (>", BAL_CUTOFF, "%):", sum(phen$bal_eos_high_binary), "(", round(mean(phen$bal_eos_high_binary)*100, 1), "%)\n")
cat("Blood High (>", BLOOD_CUTOFF, "):", sum(phen$blood_eos_high_binary), "(", round(mean(phen$blood_eos_high_binary)*100, 1), "%)\n")
cat("Both High:", TP, "(", round(TP/total_n*100, 1), "%)\n")
cat("Both Normal:", TN, "(", round(TN/total_n*100, 1), "%)\n")
cat("Discordant cases:", FP + FN, "(", round((FP + FN)/total_n*100, 1), "%)\n")

# 13. MEDICATION EFFECTS ON EOSINOPHIL LEVELS
print("\n=== MEDICATION EFFECTS ON EOSINOPHIL LEVELS ===")

# Define medication variables
med_vars <- c("SABA", "ICS.LABA", "ICS", "LTR", "xolair.omalizumab", 
              "nucala.Mepolizumab", "fasenra.benralizumab", "atrovent")

# Check which medications are actually used (have any "Checked" values)
med_usage <- phen %>%
  summarise(across(all_of(med_vars), ~sum(. == "Checked", na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Medication", values_to = "Count") %>%
  mutate(Percentage = round((Count / nrow(phen)) * 100, 1)) %>%
  arrange(desc(Count))

print("Medication Usage Summary:")
print(med_usage)


# Function to perform medication analysis
analyze_medication <- function(med_name, data = phen) {
  if (!med_name %in% names(data)) {
    return(NULL)
  }
  
  # Skip if no one uses this medication
  if (sum(data[[med_name]] == "Checked", na.rm = TRUE) == 0) {
    cat("\nNo patients using", med_name, "- skipping analysis\n")
    return(NULL)
  }
  
  cat("\n", strrep("=", 50), "\n")
  cat("ANALYSIS FOR:", med_name, "\n")
  cat(strrep("=", 50), "\n")
  
  # Create medication groups
  data$med_group <- ifelse(data[[med_name]] == "Checked", "Using", "Not Using")
  
  # Summary statistics
  med_summary <- data %>%
    group_by(med_group) %>%
    summarise(
      n = n(),
      # BAL eosinophils
      bal_eos_mean = mean(BAL_eos_p, na.rm = TRUE),
      bal_eos_sd = sd(BAL_eos_p, na.rm = TRUE),
      bal_eos_median = median(BAL_eos_p, na.rm = TRUE),
      bal_high_n = sum(bal_eos_high == "High", na.rm = TRUE),
      bal_high_pct = round(mean(bal_eos_high == "High", na.rm = TRUE) * 100, 1),
      # Blood eosinophils
      blood_eos_mean = mean(blood_eos, na.rm = TRUE),
      blood_eos_sd = sd(blood_eos, na.rm = TRUE), 
      blood_eos_median = median(blood_eos, na.rm = TRUE),
      blood_high_n = sum(blood_eos_high == "High", na.rm = TRUE),
      blood_high_pct = round(mean(blood_eos_high == "High", na.rm = TRUE) * 100, 1),
      .groups = 'drop'
    )
  
  print(med_summary)
  
  # Statistical tests
  cat("\nStatistical Tests:\n")
  
  # BAL eosinophils - continuous
  bal_test <- tryCatch({
    if (length(unique(data$med_group)) == 2) {
      t.test(BAL_eos_p_log ~ med_group, data = data)
    } else {
      list(p.value = NA)
    }
  }, error = function(e) list(p.value = NA))
  
  cat("BAL eosinophils (continuous) - t-test test p-value:", 
      ifelse(is.na(bal_test$p.value), "N/A", round(bal_test$p.value, 4)), "\n")
  
  # Blood eosinophils - continuous  
  blood_test <- tryCatch({
    if (length(unique(data$med_group)) == 2) {
      t.test(blood_eos_log ~ med_group, data = data)
    } else {
      list(p.value = NA)
    }
  }, error = function(e) list(p.value = NA))
  
  cat("Blood eosinophils (continuous) - t-test p-value:", 
      ifelse(is.na(blood_test$p.value), "N/A", round(blood_test$p.value, 4)), "\n")
  
  # BAL eosinophils - categorical
  bal_cat_test <- tryCatch({
    if (length(unique(data$med_group)) == 2 && length(unique(data$bal_eos_high)) == 2) {
      chisq.test(table(data$med_group, data$bal_eos_high))
    } else {
      list(p.value = NA)
    }
  }, error = function(e) list(p.value = NA))
  
  cat("BAL eosinophils (high vs normal) - Chi-square test p-value:", 
      ifelse(is.na(bal_cat_test$p.value), "N/A", round(bal_cat_test$p.value, 4)), "\n")
  
  # Blood eosinophils - categorical
  blood_cat_test <- tryCatch({
    if (length(unique(data$med_group)) == 2 && length(unique(data$blood_eos_high)) == 2) {
      chisq.test(table(data$med_group, data$blood_eos_high))
    } else {
      list(p.value = NA)
    }
  }, error = function(e) list(p.value = NA))
  
  cat("Blood eosinophils (high vs normal) - Chi-square test p-value:", 
      ifelse(is.na(blood_cat_test$p.value), "N/A", round(blood_cat_test$p.value, 4)), "\n")
  
  # Create visualization if there are users of this medication
  if (sum(data$med_group == "Using") > 0) {
    
    # Box plots for continuous values
    p_bal_cont <- ggplot(data, aes(x = med_group, y = BAL_eos_p_log, fill = med_group)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      geom_hline(yintercept = BAL_CUTOFF, color = "red", linetype = "dashed") +
      labs(title = paste("BAL Eosinophils by", med_name, "Use"),
           x = paste(med_name, "Status"), 
           y = "log-standardized BAL Eosinophil %",
           subtitle = paste("p =", ifelse(is.na(bal_test$p.value), "N/A", round(bal_test$p.value, 4)))) +
      theme_minimal() +
      theme(legend.position = "none")
    
    p_blood_cont <- ggplot(data, aes(x = med_group, y = blood_eos_log, fill = med_group)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      geom_hline(yintercept = BLOOD_CUTOFF, color = "blue", linetype = "dashed") +
      labs(title = paste("Blood Eosinophils by", med_name, "Use"),
           x = paste(med_name, "Status"),
           y = "log-standardized Blood Eosinophil Count",
           subtitle = paste("p =", ifelse(is.na(blood_test$p.value), "N/A", round(blood_test$p.value, 4)))) +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Bar plots for categorical values
    bal_cat_data <- data %>%
      group_by(med_group, bal_eos_high) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(med_group) %>%
      mutate(percentage = round(count/sum(count)*100, 1))
    
    p_bal_cat <- ggplot(bal_cat_data, aes(x = med_group, y = percentage, fill = bal_eos_high)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
      geom_text(aes(label = paste0(count, "\n(", percentage, "%)")), 
                position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
      labs(title = paste("BAL Eosinophil Categories by", med_name, "Use"),
           x = paste(med_name, "Status"),
           y = "Percentage",
           fill = "BAL Eosinophil Level",
           subtitle = paste("p =", ifelse(is.na(bal_cat_test$p.value), "N/A", round(bal_cat_test$p.value, 4)))) +
      theme_minimal()
    
    blood_cat_data <- data %>%
      group_by(med_group, blood_eos_high) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(med_group) %>%
      mutate(percentage = round(count/sum(count)*100, 1))
    
    p_blood_cat <- ggplot(blood_cat_data, aes(x = med_group, y = percentage, fill = blood_eos_high)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
      geom_text(aes(label = paste0(count, "\n(", percentage, "%)")), 
                position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
      labs(title = paste("Blood Eosinophil Categories by", med_name, "Use"),
           x = paste(med_name, "Status"),
           y = "Percentage", 
           fill = "Blood Eosinophil Level",
           subtitle = paste("p =", ifelse(is.na(blood_cat_test$p.value), "N/A", round(blood_cat_test$p.value, 4)))) +
      theme_minimal()
    
    # Print plots
    combined_plot <- grid.arrange(p_bal_cont, p_blood_cont, p_bal_cat, p_blood_cat, 
                                  ncol = 2, nrow = 2,
                                  top = paste("Eosinophil Analysis:", med_name))
    print(combined_plot)
  }
  
  return(med_summary)
}

# Run analysis for each medication with sufficient usage
medications_to_analyze <- med_usage$Medication[med_usage$Count > 0]

medication_results <- list()
for (med in medications_to_analyze) {
  medication_results[[med]] <- analyze_medication(med)
}

# Create summary table of all medication effects
print("\n=== SUMMARY OF ALL MEDICATION EFFECTS ===")

# Collect p-values for all medications
med_pvalue_summary <- data.frame(
  Medication = character(),
  Users_n = numeric(),
  Users_pct = numeric(),
  BAL_continuous_p = character(),
  Blood_continuous_p = character(), 
  BAL_categorical_p = character(),
  Blood_categorical_p = character(),
  stringsAsFactors = FALSE
)

for (med in medications_to_analyze) {
  users_n <- sum(phen[[med]] == "Checked", na.rm = TRUE)
  users_pct <- round((users_n / nrow(phen)) * 100, 1)
  
  # Get p-values
  phen$temp_med_group <- ifelse(phen[[med]] == "Checked", "Using", "Not Using")
  
  bal_cont_p <- tryCatch({
    test <- wilcox.test(BAL_eos_p ~ temp_med_group, data = phen)
    ifelse(test$p.value < 0.001, "<0.001", round(test$p.value, 3))
  }, error = function(e) "N/A")
  
  blood_cont_p <- tryCatch({
    test <- wilcox.test(blood_eos ~ temp_med_group, data = phen)
    ifelse(test$p.value < 0.001, "<0.001", round(test$p.value, 3))
  }, error = function(e) "N/A")
  
  bal_cat_p <- tryCatch({
    test <- chisq.test(table(phen$temp_med_group, phen$bal_eos_high))
    ifelse(test$p.value < 0.001, "<0.001", round(test$p.value, 3))
  }, error = function(e) "N/A")
  
  blood_cat_p <- tryCatch({
    test <- chisq.test(table(phen$temp_med_group, phen$blood_eos_high))
    ifelse(test$p.value < 0.001, "<0.001", round(test$p.value, 3))
  }, error = function(e) "N/A")
  
  med_pvalue_summary <- rbind(med_pvalue_summary, data.frame(
    Medication = med,
    Users_n = users_n,
    Users_pct = users_pct,
    BAL_continuous_p = bal_cont_p,
    Blood_continuous_p = blood_cont_p,
    BAL_categorical_p = bal_cat_p,
    Blood_categorical_p = blood_cat_p,
    stringsAsFactors = FALSE
  ))
}

print(med_pvalue_summary)

# Clean up temporary variable
phen$temp_med_group <- NULL

# 15. SAMPLING DATE STRATIFIED ANALYSIS
print("\n=== SAMPLING DATE STRATIFIED ANALYSIS ===")

# Create sampling date groups for stratified analysis
phen <- phen %>%
  mutate(
    date_group = case_when(
      sampling_date_diff_days == 0 ~ "Same Day",
      abs(sampling_date_diff_days) <= 7 ~ "Within 7 Days",
      abs(sampling_date_diff_days) <= 14 ~ "Within 14 Days",
      abs(sampling_date_diff_days) <= 30 ~ "Within 30 Days",
      TRUE ~ "Other"
    ),
    date_group = factor(date_group, levels = c("Same Day", "Within 7 Days", "Within 14 Days", "Within 30 Days"))
  )

# Stratified correlation analysis
stratified_correlations <- phen %>%
  group_by(date_group) %>%
  summarise(
    n = n(),
    cor_raw = cor(BAL_eos_p, blood_eos, use = "complete.obs"),
    cor_log = cor(BAL_eos_p_log, blood_eos_log, use = "complete.obs"),
    cor_spearman = cor(BAL_eos_p, blood_eos, method = "spearman", use = "complete.obs"),
    .groups = 'drop'
  ) %>%
  mutate(
    cor_raw = round(cor_raw, 3),
    cor_log = round(cor_log, 3),
    cor_spearman = round(cor_spearman, 3)
  )

print("Correlations by Sampling Date Groups:")
print(stratified_correlations)

# Stratified agreement analysis
stratified_agreement <- phen %>%
  group_by(date_group) %>%
  summarise(
    n = n(),
    both_high = sum(bal_eos_high == "High" & blood_eos_high == "High"),
    both_normal = sum(bal_eos_high == "Normal" & blood_eos_high == "Normal"),
    discordant = sum((bal_eos_high == "High" & blood_eos_high == "Normal") | 
                       (bal_eos_high == "Normal" & blood_eos_high == "High")),
    agreement_pct = round((both_high + both_normal) / n * 100, 1),
    .groups = 'drop'
  )

print("Agreement by Sampling Date Groups:")
print(stratified_agreement)

# Visualization of correlations by date groups
p_strat_cor <- ggplot(stratified_correlations, aes(x = date_group)) +
  geom_col(aes(y = cor_log, fill = "Log-transformed"), alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_col(aes(y = cor_raw, fill = "Raw values"), alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_text(aes(y = cor_log + 0.02, label = paste0("n=", n)), size = 3) +
  labs(
    title = "Correlation Strength by Sampling Date Difference",
    x = "Sampling Date Group",
    y = "Pearson Correlation Coefficient",
    fill = "Data Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_strat_cor)

# 16. INTERPRETATION GUIDE
print("\n=== INTERPRETATION GUIDE ===")
cat("Cohen's Kappa interpretation:\n")
cat("  < 0.20: Poor agreement\n")
cat("  0.21-0.40: Fair agreement\n") 
cat("  0.41-0.60: Moderate agreement\n")
cat("  0.61-0.80: Good agreement\n")
cat("  0.81-1.00: Very good agreement\n\n")

cat("AUC interpretation:\n")
cat("  0.90-1.00: Excellent\n")
cat("  0.80-0.89: Good\n") 
cat("  0.70-0.79: Fair\n")
cat("  0.60-0.69: Poor\n")
cat("  0.50-0.59: Fail\n")

# Optional: Save results
write.csv(results_summary, "bal_blood_eosinophil_concordance_results.csv", row.names = FALSE)
write.csv(phen, "eosinophil_data_with_concordance_analysis.csv", row.names = FALSE)
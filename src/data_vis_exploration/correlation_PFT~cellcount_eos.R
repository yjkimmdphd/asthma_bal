# Load necessary libraries
library(dplyr)

# -----------------------------------------------------------------------------
# 1) Load phenotype data
# -----------------------------------------------------------------------------
phen_path <- file.path(
  "./resources/processed_data",
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv"
)
phen <- read.csv(phen_path)

phen$Batch<-factor(phen$Batch,levels=unique(phen$Batch))

phen_bronch<-phen[grepl("^B",phen$SampleID),]
phen_nasal<-phen[grepl("^N",phen$SampleID),]
phen_nasal<-phen_nasal[-grep("F", phen_nasal$SampleID),]

phen_input<-phen_bronch

# Define the variables to test
variables_to_test <- c("asthma_phen_ACT.score", "asthma_phen_FEV1_perc", 
                       "asthma_phen_FVC_perc", "asthma_phen_FEV1.FVC", 
                       "asthma_phen_Post.alb_FEV1_perc", 
                       "asthma_phen_Post.alb_FVC_perc", 
                       "asthma_phen_Post.alb_FEV1.FVC")

# Define the conditions and corresponding filters
conditions <- list(
  BAL_eos_ct_gt0 = list(filter = "BAL_eos_ct > 0", x_var = "BAL_eos_ct_log"),
  BAL_eos_ct_ge0 = list(filter = "BAL_eos_ct >= 0", x_var = "BAL_eos_ct_log"),
  BAL_eos_p_gt0 = list(filter = "BAL_eos_p > 0", x_var = "BAL_eos_p_log"),
  BAL_eos_p_ge0 = list(filter = "BAL_eos_p >= 0", x_var = "BAL_eos_p_log")
)

# Function to calculate, plot, and save correlations
calculate_and_plot_correlations <- function(data, condition, x_var, y_var) {
  # Filter the data
  filtered_data <- data %>% filter(eval(parse(text = condition)), eval(parse(text = y_var)) > 0)
  
  # Extract the relevant vectors
  x <- filtered_data %>% pull(!!sym(x_var))
  y <- filtered_data %>% pull(!!sym(y_var))
  
  # Plot the data
  plot(x, y, xlab = paste("log(", x_var, ")", sep = ""), ylab = y_var, main = paste(y_var, "vs", x_var))
  
  # Perform correlation test
  result <- cor.test(x, y, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Condition:", condition, "\n")
  cat("Variable:", y_var, "\n")
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n\n")
  
  # Linear model and regression line
  model <- lm(y ~ x)
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  text(x = min(x, na.rm = TRUE), 
       y = max(y, na.rm = TRUE), 
       labels = paste("r =", round(correlation_coefficient, 2), "R-squared:", round(r_squared, 2), "P-value:", round(p_value, 2)), 
       pos = 4, 
       col = "black")
}

# Loop over each variable and create a 2x2 grid of plots
for (y_var in variables_to_test) {
  # Generate plot file name
  plot_filename <- paste0("plots_2x2_", y_var, ".png")
  
  # Save the 2x2 grid of plots
  png(plot_filename, width = 800, height = 800)
  par(mfrow = c(2, 2))  # Set up a 2x2 plotting area
  
  for (condition_name in names(conditions)) {
    condition <- conditions[[condition_name]]$filter
    x_var <- conditions[[condition_name]]$x_var
    calculate_and_plot_correlations(bphen, condition, x_var, y_var)
  }
  
  dev.off()
}

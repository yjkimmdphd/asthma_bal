# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)

# Define folder path
folder_path <- "./reports/local_only/deg_bal_nasal~cell2025-01-04"

# List all .csv files in the folder with the specified pattern
file_list <- list.files(folder_path, pattern = "deg_nasal_res_sig_", full.names = TRUE)

# Initialize a list to store results
deg_list <- list()

# Process each file
for (file in file_list) {
  # Read the CSV file
  data <- read.csv(file, row.names = 1)
  
  # Filter DEGs based on criteria
  degs <- data %>%
    filter(abs(log2FoldChange) >= 1, padj < 0.05) %>%
    mutate(
      regulation = ifelse(log2FoldChange > 0, "up", "down")
    )
  
  # Get file name without path and extension, replace '.' with '_'
  file_name <- str_replace_all(basename(file), c("\\.csv$" = "", "\\." = "_","\\.."="_","\\..."="_"))
  
  # Create separate lists for up and down-regulated DEGs
  for (reg in c("up", "down")) {
    degs_reg <- degs %>% filter(regulation == reg)
    deg_list[[paste(file_name, reg, sep = "_")]] <- rownames(degs_reg)
  }
}

# Find the maximum length of any list
max_length <- max(sapply(deg_list, length))

# Normalize all lists to the same length by padding with empty strings
deg_list_padded <- lapply(deg_list, function(x) {
  c(x, rep("", max_length - length(x)))
})

# Combine results into a data frame
final_df <- as.data.frame(deg_list_padded, stringsAsFactors = FALSE)

# Remove columns where all entries are empty strings
final_df <- final_df[, colSums(final_df != "") > 0]

# Save the final result as a tab-delimited file
output_file <- file.path(folder_path, "final_deg_list.txt")
write.table(final_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Print a message indicating success
cat("Final DEG list saved to:", output_file, "\n")

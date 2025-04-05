##
# Perform MDS/PCA analysis with normalized counts
## 

# 1. Load necessary libraries
library(DESeq2)      # For RNA-seq data normalization
library(FactoMineR)  # For PCA
library(corrplot)    # For correlation visualization
library(tidyverse)   # For data manipulation
library(limma)       # For MDS
library(Glimma)      # For interactive visualization
library(edgeR)       # For RNA-seq analysis
library(grid)        # For graphics
library(gridExtra)   # For arranging plots
library(patchwork)   # For combining plots

# 2. Load normalized count data
load_count_data <- function(file_path) {
  if (file.exists(file_path)) {
    read.table(file_path, header = TRUE, row.names = 1, sep = "\t")
  } else {
    stop(paste("File not found:", file_path))
  }
}

normalized_counts_b <- load_count_data("./resources/processed_data/normalized_gene_count/normalized_gene_count_bronch_vsd_batch-corrected.txt")
normalized_counts_n <- load_count_data("./resources/processed_data/normalized_gene_count/normalized_gene_count_nasal_vsd_batch-corrected.txt")

# Display dimensions of the data
print(paste("Bronchial samples dimensions:", paste(dim(normalized_counts_b), collapse = " x ")))
print(paste("Nasal samples dimensions:", paste(dim(normalized_counts_n), collapse = " x ")))

# Find shared genes between nasal and bronchial samples
shared_rownames <- intersect(rownames(normalized_counts_b), rownames(normalized_counts_n))
print(paste("Number of shared genes:", length(shared_rownames)))

# Combine datasets, keeping only shared genes
normalized_counts <- cbind(
  normalized_counts_b[shared_rownames,],
  normalized_counts_n[shared_rownames,]
)

# 3. Load and process phenotype data
# Define variable groups for later use
source.cell.log <- c(
  "BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log",
  "BAL_wbc_log", "blood_eos_log", "blood_eos_p_log", "blood_neut_log",
  "blood_neut_p_log", "blood_wbc_log"
)

source.cell <- c(
  "BAL_eos_ct", "BAL_eos_p", "BAL_neut_ct", "BAL_neut_p", "BAL_wbc",
  "blood_eos", "blood_eos_p", "blood_neut", "blood_neut_p", "blood_wbc"
)

pft <- c("FEV1", "FEV1_percent", "FVC", "FVC_percent", "FEV1.FVC")

# Read phenotype data
phen <- read.csv("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2025-02-14.csv")

# Process phenotype data
phen <- phen %>% 
  # Filter samples with valid BAL cell counts
  filter(BAL_eos_p >= 0 & BAL_neut_p >= 0) %>%
  # Keep only samples that underwent RNA-seq
  filter(gsub("-", ".", SampleID) %in% colnames(normalized_counts)) %>%
  # Add sample type indicator
  mutate(
    n_b = ifelse(grepl("^B", SampleID), "bronch", "nasal"),
    # Handle missing values for ED visits and hospitalizations
    ED_visits = if_else(asthma_ED == "No", 0, ED_visits),
    admit_count = if_else(asthma_admit_hx == "No", 0, admit_count),
    # Create exacerbation count (sum of ED visits and admissions)
    exac = ED_visits + admit_count,
    # scale age
    Age_at_visit_scaled = as.numeric(scale(Age_at_visit)),
    # Process IgE values
    total_IgE_num = as.numeric(gsub(">", "", total_IgE_repeat)),
    total_IgE_num = if_else(grepl(">", total_IgE_repeat), 5000, total_IgE_num)
  )

print(paste("Number of samples after filtering:", nrow(phen)))

# 4. Prepare variables for MDS color coding
var <- c(source.cell.log, "Batch", "Age_at_visit", "Sex", "Race", "Ethnicity", "total_IgE_num", "n_b")
phen_var <- phen[, var]

# Convert categorical variables to factors
phen_var[, c("Batch", "Sex", "Race", "Ethnicity", "n_b")] <- lapply(
  phen_var[, c("Batch", "Sex", "Race", "Ethnicity", "n_b")],
  function(d) factor(d, levels = unique(d))
)

# Round numeric variables
for (n in var) {
  if (!is.factor(phen_var[, n])) {
    phen_var[, n] <- round(phen_var[, n], 1)
  }
}

# 5. Create color coding function for visualization
create_color_df <- function(df) {
  # Initialize result dataframe
  col_df <- data.frame(matrix(ncol = ncol(df), nrow = nrow(df)))
  names(col_df) <- names(df)
  
  # Define color palettes
  quart_colors <- c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D") # Red gradient
  na_color <- "#CCCCCC"  # Gray for NA values
  
  # Color palettes for categorical variables
  cat_palettes <- list(
    c("#1f77b4", "#ff7f0e"),                                    # 2 levels
    c("#1f77b4", "#ff7f0e", "#2ca02c"),                         # 3 levels
    c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"),              # 4 levels
    c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"),   # 5 levels
    c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")  # 6 levels
  )
  
  # Process each column
  for (col_name in names(df)) {
    col <- df[[col_name]]
    
    if (is.numeric(col) || is.integer(col)) {
      # For numeric columns, assign colors based on quartiles
      quarts <- quantile(col, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
      
      # Set default color as NA color
      col_df[[col_name]] <- na_color
      
      # Non-NA values
      non_na_idx <- !is.na(col)
      
      for (i in 1:4) {
        if (i < 4) {
          idx <- non_na_idx & col >= quarts[i] & col < quarts[i+1]
        } else {
          idx <- non_na_idx & col >= quarts[i] & col <= quarts[i+1]
        }
        col_df[idx, col_name] <- quart_colors[i]
      }
    } else {
      # For categorical columns, assign unique colors for each value
      levels <- unique(na.omit(col))  # Exclude NA from levels
      n_levels <- length(levels)
      
      # Select appropriate palette based on number of levels
      palette_idx <- min(n_levels, length(cat_palettes))
      colors <- cat_palettes[[palette_idx]]
      
      # If we have more levels than colors in our largest palette, recycle colors
      if (n_levels > length(colors)) {
        colors <- rep(colors, length.out = n_levels)
      }
      
      # Start with NA color for all
      col_df[[col_name]] <- na_color
      
      # Assign colors to non-NA values
      for (i in 1:n_levels) {
        idx <- which(col == levels[i])
        if (length(idx) > 0) {
          col_df[idx, col_name] <- colors[i]
        }
      }
    }
    
    # Convert the column to character
    col_df[[col_name]] <- as.character(col_df[[col_name]])
  }
  
  return(col_df)
}

# Create color dataframe
col <- create_color_df(phen_var)

# 6. Generate MDS plots for different variables
par(mar = c(4, 4, 4, 4))

# Function to create MDS plots in a more organized way
create_mds_plot <- function(data, color_var, label_var, title, dim1=1, dim2=2) {
  # If using default dimensions 1 and 2, use the builtin plotting
  if (dim1 == 1 && dim2 == 2) {
    return(plotMDS(data, col = color_var, labels = label_var, main = title))
  }
  
  # Otherwise, calculate MDS with enough dimensions and plot manually
  mds_result <- plotMDS(data, ndim=max(dim1, dim2), plot=FALSE)
  
  # Extract coordinates for the requested dimensions
  x_coords <- if (dim1 <= 2) {
    if (dim1 == 1) mds_result$x else mds_result$y
  } else {
    mds_result$eigen.vectors[, dim1] * sqrt(abs(mds_result$eigen.values[dim1]))
  }
  
  y_coords <- if (dim2 <= 2) {
    if (dim2 == 1) mds_result$x else mds_result$y
  } else {
    mds_result$eigen.vectors[, dim2] * sqrt(abs(mds_result$eigen.values[dim2]))
  }
  
  # Calculate variance explained percentages
  var_exp_x <- round(mds_result$var.explained[dim1] * 100, 1)
  var_exp_y <- round(mds_result$var.explained[dim2] * 100, 1)
  
  # Create the plot
  plot(x_coords, y_coords, 
       col = color_var, 
       pch = 19,
       main = title,
       xlab = paste0("MDS", dim1, " (", var_exp_x, "%)"),
       ylab = paste0("MDS", dim2, " (", var_exp_y, "%)"))
  
  # Add labels if provided and not NULL
  if (!is.null(label_var)) {
    text(x_coords, y_coords, labels = label_var, pos = 3, cex = 0.7)
  }
  
  # Return useful information
  invisible(list(
    x = x_coords,
    y = y_coords,
    dim1 = dim1,
    dim2 = dim2,
    var_explained = mds_result$var.explained
  ))
}


# Create MDS plots for cell count variables
mds_plots_cell <- list()
for (i in 1:length(source.cell)) {
  var_name <- source.cell.log[i]
  mds_plots_cell[[i]] <- create_mds_plot(
    normalized_counts, 
    col[[var_name]], 
    phen_var[[var_name]], 
    source.cell[i]
  )
}

# Create MDS plots for metadata variables
mds_batch <- create_mds_plot(normalized_counts, col$Batch, phen_var$Batch, "Batch", dim1=1, dim2=2)
mds_sex <- create_mds_plot(normalized_counts, col$Sex, phen_var$Sex, "Sex", dim1=1, dim2=2)
mds_sex <- create_mds_plot(normalized_counts, col$Sex, phen_var$Sex, "Sex", dim1=1, dim2=4)
mds_nasal_bronch <- create_mds_plot(normalized_counts, col$n_b, phen_var$n_b, "Sample type (nasal/bronch)", dim1=1, dim2=2)
mds_age <- create_mds_plot(normalized_counts, col$Age_at_visit, phen_var$Age_at_visit, "Age", dim1=1, dim2=2)
mds_race <- create_mds_plot(normalized_counts, col$Race, phen_var$Race, "Race", dim1=1, dim2=2)
mds_eth <- create_mds_plot(normalized_counts, col$Ethnicity, phen_var$Ethnicity, "Ethnicity", dim1=1, dim2=2)

# 7. Analyze MDS results
mds_result <- plotMDS(normalized_counts, plot = FALSE)

# Calculate the proportion of variance explained by each PC
var_explained_df <- data.frame(
  PC = factor(paste0("MDS", 1:length(mds_result$var.explained)), 
              levels = paste0("MDS", 1:length(mds_result$var.explained))),
  VarExplained = mds_result$var.explained,
  Percentage = round(mds_result$var.explained * 100, 2)
)

# Print the percentage of variance explained by each PC
print("Variance explained by top PCs:")
print(var_explained_df[1:10, ])

# 8. Create variance explained plot
p_var_explained <- ggplot(var_explained_df[1:10, ], aes(x = PC, y = Percentage)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5) +
  theme_minimal() +
  labs(title = "Percentage of Variance Explained by MDS Components",
       x = "MDS Component",
       y = "Percentage of Variance Explained")

print(p_var_explained)

# 9. Extract MDS scores and prepare for correlation analysis
mds_scores <- mds_result$eigen.vectors

# Convert to dataframe and add sample IDs
mds_df <- as.data.frame(mds_scores)
colnames(mds_df) <- paste0("MDS", seq_len(ncol(mds_df)))
mds_df$SampleID <- colnames(normalized_counts)

# Merge PC scores with metadata
merged_data <- merge(phen, mds_df, by = "SampleID")

# 10. Define variables for association analysis
continuous_vars <- c("Age_at_visit_scaled", "admit_count", "ED_visits", "ACT_score", 
                     "FEV1_percent", "BAL_eos_p_log", "blood_eos_p_log")
categorical_vars <- c("Sex", "Race", "Ethnicity", "bal_eos_p_mt1", "Batch", "n_b")
selected_pcs <- colnames(mds_df)[1:5]  # First 5 PCs

# 11. Function to calculate associations between variables and PCs
get_association <- function(var_name, pc_name, data) {
  if (var_name %in% categorical_vars) {
    # For categorical variables, use ANOVA
    data[[var_name]] <- as.factor(data[[var_name]])
    formula <- as.formula(paste(pc_name, "~", var_name))
    model <- aov(formula, data = data)
    anova_result <- summary(model)
    
    # Calculate eta-squared
    ss_effect <- anova_result[[1]]["Sum Sq"][[1]][1]
    ss_total <- sum(anova_result[[1]]["Sum Sq"][[1]])
    eta_squared <- ss_effect / ss_total
    
    # Get p-value
    p_value <- anova_result[[1]]["Pr(>F)"][[1]][1]
    
    return(list(
      association = eta_squared,
      p_value = p_value,
      method = "ANOVA"
    ))
  } else {
    # For continuous variables, use Pearson correlation
    test <- cor.test(data[[var_name]], data[[pc_name]], method = "pearson")
    return(list(
      association = as.numeric(test$estimate),
      p_value = test$p.value,
      method = "Pearson"
    ))
  }
}

# 12. Apply association function to all variables and PCs
all_vars <- c(continuous_vars, categorical_vars)
results_df <- data.frame()

for (var in all_vars) {
  for (pc in selected_pcs) {
    # Skip if variable doesn't exist in data
    if (!var %in% colnames(merged_data)) {
      warning(paste("Variable", var, "not found in data"))
      next
    }
    
    # Skip if too many missing values
    if (sum(is.na(merged_data[[var]])) > 0.8 * nrow(merged_data)) {
      warning(paste("Too many missing values in", var))
      next
    }
    
    result <- get_association(var, pc, merged_data)
    
    results_df <- rbind(results_df, data.frame(
      Variable = var,
      PC = pc,
      Association = result$association,
      P_value = result$p_value,
      Method = result$method,
      Type = ifelse(var %in% categorical_vars, "Categorical", "Continuous")
    ))
  }
}

# Apply multiple testing correction
results_df$Adj_P_value <- p.adjust(results_df$P_value, method = "fdr")
results_df$Significant <- results_df$Adj_P_value < 0.05

# 13. Create matrices for visualization
# For continuous variables
cont_matrix <- matrix(NA, nrow = length(continuous_vars), ncol = length(selected_pcs),
                      dimnames = list(continuous_vars, selected_pcs))
cont_p_matrix <- matrix(NA, nrow = length(continuous_vars), ncol = length(selected_pcs),
                        dimnames = list(continuous_vars, selected_pcs))

# For categorical variables
cat_matrix <- matrix(NA, nrow = length(categorical_vars), ncol = length(selected_pcs),
                     dimnames = list(categorical_vars, selected_pcs))
cat_p_matrix <- matrix(NA, nrow = length(categorical_vars), ncol = length(selected_pcs),
                       dimnames = list(categorical_vars, selected_pcs))

# Fill matrices
for (i in 1:nrow(results_df)) {
  row <- results_df[i, ]
  if (row$Type == "Continuous" && row$Variable %in% rownames(cont_matrix)) {
    cont_matrix[row$Variable, row$PC] <- row$Association
    cont_p_matrix[row$Variable, row$PC] <- row$Adj_P_value
  } else if (row$Type == "Categorical" && row$Variable %in% rownames(cat_matrix)) {
    cat_matrix[row$Variable, row$PC] <- row$Association
    cat_p_matrix[row$Variable, row$PC] <- row$Adj_P_value
  }
}

# 14. Visualize results with corrplot
# For continuous variables
if (requireNamespace("corrplot", quietly = TRUE)) {
  # Correlation plot for continuous variables
  corrplot(cont_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           p.mat = cont_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           cl.offset = -0.5,
           cl.align.text = "r",
           mar = c(6, 1, 7, 0),
           title = "Pearson Correlations: \nContinuous Variables vs PCs")
  
  corrplot(cont_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           p.mat = cont_p_matrix, sig.level = 0.05,
           insig = "label_sig", 
           cl.offset = -0.5,
           cl.align.text = "r",
           mar = c(6, 1, 7, 0),
           title = "Pearson Correlations: \nContinuous Variables vs PCs")
  
  # Effect size plot for categorical variables
  corrplot(cat_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           col = colorRampPalette(c("white", "pink", "red"))(100),
           p.mat = cat_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           mar = c(6, 1, 7, 0),
           cl.offset = -0.5,
           cl.align.text = "r",
           title = "Effect Sizes (η²): \nCategorical Variables vs PCs",
           is.corr = FALSE)
  
  corrplot(cat_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           col = colorRampPalette(c("white", "pink", "red"))(100),
           p.mat = cat_p_matrix, sig.level = 0.05,
           insig = "label_sig", 
           mar = c(6, 1, 7, 0),
           cl.offset = -0.5,
           cl.align.text = "r",
           title = "Effect Sizes (η²): \nCategorical Variables vs PCs",
           is.corr = FALSE)
} else {
  # Fallback to base R heatmap if corrplot is not available
  par(mfrow = c(1, 2))
  image(1:ncol(cont_matrix), 1:nrow(cont_matrix), t(cont_matrix),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        xlab = "Principal Component", ylab = "Variable",
        axes = FALSE,
        main = "Pearson Correlations: Continuous Variables vs PCs")
  axis(1, at = 1:ncol(cont_matrix), labels = colnames(cont_matrix))
  axis(2, at = 1:nrow(cont_matrix), labels = rownames(cont_matrix), las = 1)
  
  image(1:ncol(cat_matrix), 1:nrow(cat_matrix), t(cat_matrix),
        col = colorRampPalette(c("white", "red"))(100),
        xlab = "Principal Component", ylab = "Variable",
        axes = FALSE,
        main = "Effect Sizes (η²): Categorical Variables vs PCs")
  axis(1, at = 1:ncol(cat_matrix), labels = colnames(cat_matrix))
  axis(2, at = 1:nrow(cat_matrix), labels = rownames(cat_matrix), las = 1)
}

# 15. Save results to a CSV for external visualization
write.csv(results_df, "./reports/local_only/correlation_pca_trait/MDS_associations.csv", row.names = FALSE)

######################
# bronchial MDS results 
######################
phen_var<-filter(phen_var,n_b=="bronch")
normalized_counts<-normalized_counts_b

# Create color dataframe
col <- create_color_df(phen_var)



# Create MDS plots for cell count variables
mds_plots_cell <- list()
for (i in 1:length(source.cell)) {
  var_name <- source.cell.log[i]
  mds_plots_cell[[i]] <- create_mds_plot(
    normalized_counts, 
    col[[var_name]], 
    phen_var[[var_name]], 
    source.cell[i]
  )
}

# Create MDS plots for metadata variables
mds_batch <- create_mds_plot(normalized_counts, col$Batch, phen_var$Batch, "Batch", dim1=1, dim2=2)
mds_sex <- create_mds_plot(normalized_counts, col$Sex, phen_var$Sex, "Sex", dim1=1, dim2=2)
mds_sex <- create_mds_plot(normalized_counts, col$Sex, phen_var$Sex, "Sex", dim1=1, dim2=4)
mds_nasal_bronch <- create_mds_plot(normalized_counts, col$n_b, phen_var$n_b, "Sample type (nasal/bronch)", dim1=1, dim2=2)
mds_age <- create_mds_plot(normalized_counts, col$Age_at_visit, phen_var$Age_at_visit, "Age", dim1=1, dim2=2)
mds_race <- create_mds_plot(normalized_counts, col$Race, phen_var$Race, "Race", dim1=1, dim2=2)
mds_eth <- create_mds_plot(normalized_counts, col$Ethnicity, phen_var$Ethnicity, "Ethnicity", dim1=1, dim2=2)



# Analyze MDS results, bronchial only
mds_result <- plotMDS(normalized_counts_b, plot = FALSE)

# Calculate the proportion of variance explained by each PC
var_explained_df <- data.frame(
  PC = factor(paste0("MDS", 1:length(mds_result$var.explained)), 
              levels = paste0("MDS", 1:length(mds_result$var.explained))),
  VarExplained = mds_result$var.explained,
  Percentage = round(mds_result$var.explained * 100, 2)
)

# Print the percentage of variance explained by each PC
print("Variance explained by top PCs:")
print(var_explained_df[1:10, ])

# 8. Create variance explained plot
p_var_explained <- ggplot(var_explained_df[1:10, ], aes(x = PC, y = Percentage)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5) +
  theme_minimal() +
  labs(title = "Percentage of Variance Explained by MDS Components",
       x = "MDS Component",
       y = "Percentage of Variance Explained")

print(p_var_explained)

# 9. Extract MDS scores and prepare for correlation analysis
mds_scores <- mds_result$eigen.vectors

# Convert to dataframe and add sample IDs
mds_df <- as.data.frame(mds_scores)
colnames(mds_df) <- paste0("MDS", seq_len(ncol(mds_df)))
mds_df$SampleID <- colnames(normalized_counts)

# Merge PC scores with metadata
merged_data <- merge(phen, mds_df, by = "SampleID")

# 10. Define variables for association analysis
continuous_vars <- c("Age_at_visit_scaled", "admit_count", "ED_visits", "ACT_score", 
                     "FEV1_percent", "BAL_eos_p_log", "blood_eos_p_log")
categorical_vars <- c("Sex", "Race", "Ethnicity", "bal_eos_p_mt1", "Batch")
selected_pcs <- colnames(mds_df)[1:5]  # First 5 PCs

# 12. Apply association function to all variables and PCs
all_vars <- c(continuous_vars, categorical_vars)
results_df <- data.frame()

for (var in all_vars) {
  for (pc in selected_pcs) {
    # Skip if variable doesn't exist in data
    if (!var %in% colnames(merged_data)) {
      warning(paste("Variable", var, "not found in data"))
      next
    }
    
    # Skip if too many missing values
    if (sum(is.na(merged_data[[var]])) > 0.8 * nrow(merged_data)) {
      warning(paste("Too many missing values in", var))
      next
    }
    
    result <- get_association(var, pc, merged_data)
    
    results_df <- rbind(results_df, data.frame(
      Variable = var,
      PC = pc,
      Association = result$association,
      P_value = result$p_value,
      Method = result$method,
      Type = ifelse(var %in% categorical_vars, "Categorical", "Continuous")
    ))
  }
}

# Apply multiple testing correction
results_df$Adj_P_value <- p.adjust(results_df$P_value, method = "fdr")
results_df$Significant <- results_df$Adj_P_value < 0.05

# 13. Create matrices for visualization
# For continuous variables
cont_matrix <- matrix(NA, nrow = length(continuous_vars), ncol = length(selected_pcs),
                      dimnames = list(continuous_vars, selected_pcs))
cont_p_matrix <- matrix(NA, nrow = length(continuous_vars), ncol = length(selected_pcs),
                        dimnames = list(continuous_vars, selected_pcs))

# For categorical variables
cat_matrix <- matrix(NA, nrow = length(categorical_vars), ncol = length(selected_pcs),
                     dimnames = list(categorical_vars, selected_pcs))
cat_p_matrix <- matrix(NA, nrow = length(categorical_vars), ncol = length(selected_pcs),
                       dimnames = list(categorical_vars, selected_pcs))

# Fill matrices
for (i in 1:nrow(results_df)) {
  row <- results_df[i, ]
  if (row$Type == "Continuous" && row$Variable %in% rownames(cont_matrix)) {
    cont_matrix[row$Variable, row$PC] <- row$Association
    cont_p_matrix[row$Variable, row$PC] <- row$Adj_P_value
  } else if (row$Type == "Categorical" && row$Variable %in% rownames(cat_matrix)) {
    cat_matrix[row$Variable, row$PC] <- row$Association
    cat_p_matrix[row$Variable, row$PC] <- row$Adj_P_value
  }
}

# 14. Visualize results with corrplot
# For continuous variables
if (requireNamespace("corrplot", quietly = TRUE)) {
  # Correlation plot for continuous variables
  corrplot(cont_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           p.mat = cont_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           cl.offset = -0.5,
           cl.align.text = "r",
           mar = c(6, 1, 7, 0),
           title = "Pearson Correlations: \nContinuous Variables vs PCs")
  
  corrplot(cont_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           p.mat = cont_p_matrix, sig.level = 0.05,
           insig = "label_sig", 
           cl.offset = -0.5,
           cl.align.text = "r",
           mar = c(6, 1, 7, 0),
           title = "Pearson Correlations: \nContinuous Variables vs PCs")
  
  # Effect size plot for categorical variables
  corrplot(cat_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           col = colorRampPalette(c("white", "pink", "red"))(100),
           p.mat = cat_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           mar = c(6, 1, 7, 0),
           cl.offset = -0.5,
           cl.align.text = "r",
           title = "Effect Sizes (η²): \nCategorical Variables vs PCs",
           is.corr = FALSE)
  
  corrplot(cat_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           col = colorRampPalette(c("white", "pink", "red"))(100),
           p.mat = cat_p_matrix, sig.level = 0.05,
           insig = "label_sig", 
           mar = c(6, 1, 7, 0),
           cl.offset = -0.5,
           cl.align.text = "r",
           title = "Effect Sizes (η²): \nCategorical Variables vs PCs",
           is.corr = FALSE)
} else {
  # Fallback to base R heatmap if corrplot is not available
  par(mfrow = c(1, 2))
  image(1:ncol(cont_matrix), 1:nrow(cont_matrix), t(cont_matrix),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        xlab = "Principal Component", ylab = "Variable",
        axes = FALSE,
        main = "Pearson Correlations: Continuous Variables vs PCs")
  axis(1, at = 1:ncol(cont_matrix), labels = colnames(cont_matrix))
  axis(2, at = 1:nrow(cont_matrix), labels = rownames(cont_matrix), las = 1)
  
  image(1:ncol(cat_matrix), 1:nrow(cat_matrix), t(cat_matrix),
        col = colorRampPalette(c("white", "red"))(100),
        xlab = "Principal Component", ylab = "Variable",
        axes = FALSE,
        main = "Effect Sizes (η²): Categorical Variables vs PCs")
  axis(1, at = 1:ncol(cat_matrix), labels = colnames(cat_matrix))
  axis(2, at = 1:nrow(cat_matrix), labels = rownames(cat_matrix), las = 1)
}

# 15. Save results to a CSV for external visualization
write.csv(results_df, "./reports/local_only/correlation_pca_trait/MDS_associations_bronchial.csv", row.names = FALSE)


###################
# nasal MDS results
###################
var <- c(source.cell.log, "Batch", "Age_at_visit", "Sex", "Race", "Ethnicity", "total_IgE_num", "n_b")
phen_var <- phen[, var]
phen$SampleID<-gsub("-",".",phen$SampleID)

# Convert categorical variables to factors
phen_var[, c("Batch", "Sex", "Race", "Ethnicity", "n_b")] <- lapply(
  phen_var[, c("Batch", "Sex", "Race", "Ethnicity", "n_b")],
  function(d) factor(d, levels = unique(d))
)

# Round numeric variables
for (n in var) {
  if (!is.factor(phen_var[, n])) {
    phen_var[, n] <- round(phen_var[, n], 1)
  }
}

phen_var<-filter(phen_var,n_b=="nasal")
normalized_counts<-normalized_counts_n

# Create color dataframe
col <- create_color_df(phen_var)

# Create MDS plots for cell count variables
mds_plots_cell <- list()
for (i in 1:length(source.cell)) {
  var_name <- source.cell.log[i]
  mds_plots_cell[[i]] <- create_mds_plot(
    normalized_counts, 
    col[[var_name]], 
    phen_var[[var_name]], 
    source.cell[i]
  )
}

# Create MDS plots for metadata variables
mds_batch <- create_mds_plot(normalized_counts, col$Batch, phen_var$Batch, "Batch", dim1=1, dim2=2)
mds_sex <- create_mds_plot(normalized_counts, col$Sex, phen_var$Sex, "Sex", dim1=1, dim2=2)
mds_sex <- create_mds_plot(normalized_counts, col$Sex, phen_var$Sex, "Sex", dim1=1, dim2=4)
mds_nasal_bronch <- create_mds_plot(normalized_counts, col$n_b, phen_var$n_b, "Sample type (nasal/bronch)", dim1=1, dim2=2)
mds_age <- create_mds_plot(normalized_counts, col$Age_at_visit, phen_var$Age_at_visit, "Age", dim1=1, dim2=2)
mds_race <- create_mds_plot(normalized_counts, col$Race, phen_var$Race, "Race", dim1=1, dim2=2)
mds_eth <- create_mds_plot(normalized_counts, col$Ethnicity, phen_var$Ethnicity, "Ethnicity", dim1=1, dim2=2)



# Analyze MDS results, nasal only
mds_result <- plotMDS(normalized_counts_n, plot = FALSE)

# Calculate the proportion of variance explained by each PC
var_explained_df <- data.frame(
  PC = factor(paste0("MDS", 1:length(mds_result$var.explained)), 
              levels = paste0("MDS", 1:length(mds_result$var.explained))),
  VarExplained = mds_result$var.explained,
  Percentage = round(mds_result$var.explained * 100, 2)
)

# Print the percentage of variance explained by each PC
print("Variance explained by top PCs:")
print(var_explained_df[1:10, ])

# 8. Create variance explained plot
p_var_explained <- ggplot(var_explained_df[1:10, ], aes(x = PC, y = Percentage)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5) +
  theme_minimal() +
  labs(title = "Percentage of Variance Explained by MDS Components",
       x = "MDS Component",
       y = "Percentage of Variance Explained")

print(p_var_explained)

# 9. Extract MDS scores and prepare for correlation analysis
mds_scores <- mds_result$eigen.vectors

# Convert to dataframe and add sample IDs
mds_df <- as.data.frame(mds_scores)
colnames(mds_df) <- paste0("MDS", seq_len(ncol(mds_df)))
mds_df$SampleID <- colnames(normalized_counts)

# Merge PC scores with metadata
merged_data <- merge(phen, mds_df, by = "SampleID")

# 10. Define variables for association analysis
continuous_vars <- c("Age_at_visit_scaled", "admit_count", "ED_visits", "ACT_score", 
                     "FEV1_percent", "BAL_eos_p_log", "blood_eos_p_log")
categorical_vars <- c("Sex", "Race", "Ethnicity", "bal_eos_p_mt1", "Batch")
selected_pcs <- colnames(mds_df)[1:5]  # First 5 PCs

# 12. Apply association function to all variables and PCs
all_vars <- c(continuous_vars, categorical_vars)
results_df <- data.frame()

for (var in all_vars) {
  for (pc in selected_pcs) {
    # Skip if variable doesn't exist in data
    if (!var %in% colnames(merged_data)) {
      warning(paste("Variable", var, "not found in data"))
      next
    }
    
    # Skip if too many missing values
    if (sum(is.na(merged_data[[var]])) > 0.8 * nrow(merged_data)) {
      warning(paste("Too many missing values in", var))
      next
    }
    
    result <- get_association(var, pc, merged_data)
    
    results_df <- rbind(results_df, data.frame(
      Variable = var,
      PC = pc,
      Association = result$association,
      P_value = result$p_value,
      Method = result$method,
      Type = ifelse(var %in% categorical_vars, "Categorical", "Continuous")
    ))
  }
}

# Apply multiple testing correction
results_df$Adj_P_value <- p.adjust(results_df$P_value, method = "fdr")
results_df$Significant <- results_df$Adj_P_value < 0.05

# 13. Create matrices for visualization
# For continuous variables
cont_matrix <- matrix(NA, nrow = length(continuous_vars), ncol = length(selected_pcs),
                      dimnames = list(continuous_vars, selected_pcs))
cont_p_matrix <- matrix(NA, nrow = length(continuous_vars), ncol = length(selected_pcs),
                        dimnames = list(continuous_vars, selected_pcs))

# For categorical variables
cat_matrix <- matrix(NA, nrow = length(categorical_vars), ncol = length(selected_pcs),
                     dimnames = list(categorical_vars, selected_pcs))
cat_p_matrix <- matrix(NA, nrow = length(categorical_vars), ncol = length(selected_pcs),
                       dimnames = list(categorical_vars, selected_pcs))

# Fill matrices
for (i in 1:nrow(results_df)) {
  row <- results_df[i, ]
  if (row$Type == "Continuous" && row$Variable %in% rownames(cont_matrix)) {
    cont_matrix[row$Variable, row$PC] <- row$Association
    cont_p_matrix[row$Variable, row$PC] <- row$Adj_P_value
  } else if (row$Type == "Categorical" && row$Variable %in% rownames(cat_matrix)) {
    cat_matrix[row$Variable, row$PC] <- row$Association
    cat_p_matrix[row$Variable, row$PC] <- row$Adj_P_value
  }
}

# 14. Visualize results with corrplot
# For continuous variables
if (requireNamespace("corrplot", quietly = TRUE)) {
  # Correlation plot for continuous variables
  corrplot(cont_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           p.mat = cont_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           cl.offset = -0.5,
           cl.align.text = "r",
           mar = c(6, 1, 7, 0),
           title = "Pearson Correlations: \nContinuous Variables vs PCs")
 
  # Effect size plot for categorical variables
  corrplot(cat_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           col = colorRampPalette(c("white", "pink", "red"))(100),
           p.mat = cat_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           mar = c(6, 1, 7, 0),
           cl.offset = -0.5,
           cl.align.text = "r",
           title = "Effect Sizes (η²): \nCategorical Variables vs PCs",
           is.corr = FALSE)
  
  corrplot(cat_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           col = colorRampPalette(c("white", "pink", "red"))(100),
           p.mat = cat_p_matrix, sig.level = 0.05,
           insig = "label_sig", 
           mar = c(6, 1, 7, 0),
           cl.offset = -0.5,
           cl.align.text = "r",
           title = "Effect Sizes (η²): \nCategorical Variables vs PCs",
           is.corr = FALSE)
} else {
  # Fallback to base R heatmap if corrplot is not available
  par(mfrow = c(1, 2))
  image(1:ncol(cont_matrix), 1:nrow(cont_matrix), t(cont_matrix),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        xlab = "Principal Component", ylab = "Variable",
        axes = FALSE,
        main = "Pearson Correlations: Continuous Variables vs PCs")
  axis(1, at = 1:ncol(cont_matrix), labels = colnames(cont_matrix))
  axis(2, at = 1:nrow(cont_matrix), labels = rownames(cont_matrix), las = 1)
  
  image(1:ncol(cat_matrix), 1:nrow(cat_matrix), t(cat_matrix),
        col = colorRampPalette(c("white", "red"))(100),
        xlab = "Principal Component", ylab = "Variable",
        axes = FALSE,
        main = "Effect Sizes (η²): Categorical Variables vs PCs")
  axis(1, at = 1:ncol(cat_matrix), labels = colnames(cat_matrix))
  axis(2, at = 1:nrow(cat_matrix), labels = rownames(cat_matrix), las = 1)
}

# 15. Save results to a CSV for external visualization
write.csv(results_df, "./reports/local_only/correlation_pca_trait/MDS_associations_nasal.csv", row.names = FALSE)


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
library(pheatmap)    # For heatmaps
library(RColorBrewer) # For color palettes
library(vsn)         # For variance stabilizing transformation

# 2. Define variable groups for later use
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

# Helper function to fix p-value matrices for corrplot
# REPLACE THIS FUNCTION:
prepare_p_matrix <- function(p_matrix) {
  # Check if the p-matrix has any significant values
  has_sig <- any(p_matrix < 0.05, na.rm = TRUE)
  
  if (!has_sig) {
    # If no significant values, return a matrix with all values > 0.05
    return(matrix(0.06, nrow = nrow(p_matrix), ncol = ncol(p_matrix), 
                  dimnames = dimnames(p_matrix)))
  }
  
  # Replace NA values with 1 (non-significant)
  p_matrix[is.na(p_matrix)] <- 1
  return(p_matrix)
}

# 3. define resources folder
resources_dir<-file.path("./resources","processed_data")
resources_dir<-if(file.exists(resources_dir)){resources_dir}

# 4. Create color mappings for MDS plots
create_color_mappings <- function(phen_var) {
  col <- list()
  
  # Colors for categorical variables
  col$n_b <- ifelse(phen_var$n_b == "nasal", "blue", "red")
  col$Batch <- rainbow(length(unique(phen_var$Batch)))[match(phen_var$Batch, unique(phen_var$Batch))]
  col$Sex <- ifelse(phen_var$Sex == "Male", "cyan", "magenta")
  col$bal_eos_p_mt1 <- ifelse(phen_var$bal_eos_p_mt1 == TRUE, "green", "purple")
  
  # Colors for numeric variables (using gradients)
  col$Age_at_visit <- colorRampPalette(c("#FFF5F0", "#CB181D"))(100)[
    cut(phen_var$Age_at_visit, breaks = 100, labels = FALSE, include.lowest = TRUE)]
  
  col$BAL_eos_p_log <- colorRampPalette(c("#FFF5F0", "#CB181D"))(100)[
    cut(phen_var$BAL_eos_p_log, breaks = 100, labels = FALSE, include.lowest = TRUE)]
  
  col$blood_eos_log <- colorRampPalette(c("#FFF5F0", "#CB181D"))(100)[
    cut(phen_var$blood_eos_log, breaks = 100, labels = FALSE, include.lowest = TRUE)]
  
  return(col)
}

# 5. Function to perform MDS and create plots with legends

output_folder<-file.path("./reports/local_only/correlation_pca_trait")

output_folder_batch_corr<-file.path(output_folder,"batch_corrected")
if(!file.exists(output_folder_batch_corr)){dir.create(output_folder_batch_corr, recursive = TRUE, showWarnings = FALSE)}
output_folder_no_batch_corr<-file.path(output_folder,"no_batch_corrected")
if(!file.exists(output_folder_no_batch_corr)){dir.create(output_folder_no_batch_corr, recursive = TRUE, showWarnings = FALSE)}


create_multipanel_mds <- function(normalized_counts, phen_var, col, title_prefix = "", output_folder) {
  # Setup output PDF file
  pdf_filename <- file.path(output_folder,paste0(title_prefix, "_mds_analysis.pdf"))
  pdf(pdf_filename, width = 14, height = 14)
  
  # Set layout with wider right margin for legends
  par(mfrow = c(4, 2), mar = c(4, 4, 3, 8), oma = c(0, 0, 2, 0))
  
  # Calculate MDS
  mds <- limma::plotMDS(normalized_counts, ndim = 5, plot = FALSE)
  
  # Extract coordinates for dimensions
  x1 <- mds$eigen.vectors[, 1] * sqrt(mds$eigen.values[1])
  x2 <- mds$eigen.vectors[, 2] * sqrt(mds$eigen.values[2])
  x3 <- mds$eigen.vectors[, 3] * sqrt(mds$eigen.values[3])
  x4 <- mds$eigen.vectors[, 4] * sqrt(mds$eigen.values[4])
  
  # Calculate variance explained
  var1 <- round(mds$var.explained[1] * 100, 1)
  var2 <- round(mds$var.explained[2] * 100, 1)
  var3 <- round(mds$var.explained[3] * 100, 1)
  var4 <- round(mds$var.explained[4] * 100, 1)
  
  # NA color (grey)
  na_color <- "#CCCCCC"
  
  # Helper function to add legends for categorical variables
  add_categorical_legend <- function(var_name) {
    # Get unique values and corresponding colors
    unique_values <- unique(phen_var[[var_name]])
    unique_colors <- c()
    
    for (val in unique_values) {
      # Find first index with this value
      idx <- which(phen_var[[var_name]] == val)[1]
      if (length(idx) > 0) {
        unique_colors <- c(unique_colors, col[[var_name]][idx])
      }
    }
    
    # Include NA in the legend if there are any NA values
    if (any(is.na(phen_var[[var_name]]))) {
      legend_values <- c(as.character(unique_values), "NA")
      legend_colors <- c(unique_colors, na_color)
    } else {
      legend_values <- as.character(unique_values)
      legend_colors <- unique_colors
    }
    
    # Add legend outside the plot
    par(xpd = TRUE)  # Allow drawing outside plot region
    legend(par("usr")[2] * 1.05, par("usr")[4], 
           legend = legend_values, 
           fill = legend_colors, 
           title = var_name,
           cex = 0.7, bty = "n")
    par(xpd = FALSE)  # Reset to default
  }
  
  # Helper function to add legends for numeric variables
  add_numeric_legend <- function(var_name) {
    # Red gradient colors for numeric variables
    quart_colors <- c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D")
    
    # Calculate quartiles
    var_quants <- quantile(phen_var[[var_name]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    
    # Create legend labels
    quart_labels <- c(
      paste0("Q1: ", round(var_quants[1], 1), "-", round(var_quants[2], 1)),
      paste0("Q2: ", round(var_quants[2], 1), "-", round(var_quants[3], 1)),
      paste0("Q3: ", round(var_quants[3], 1), "-", round(var_quants[4], 1)),
      paste0("Q4: ", round(var_quants[4], 1), "-", round(var_quants[5], 1))
    )
    
    # Handle NA values in legend
    if (any(is.na(phen_var[[var_name]]))) {
      legend_labels <- c(quart_labels, "NA")
      legend_colors <- c(quart_colors, na_color)
    } else {
      legend_labels <- quart_labels
      legend_colors <- quart_colors
    }
    
    # Add legend outside plot
    par(xpd = TRUE)
    legend(par("usr")[2] * 1.05, par("usr")[4], 
           legend = legend_labels, 
           fill = legend_colors, 
           title = var_name,
           cex = 0.7, bty = "n")
    par(xpd = FALSE)
  }
  
  # Plot 1: Sample Type (if available) or Batch
  if ("n_b" %in% names(col) && length(unique(phen_var$n_b)) > 1) {
    plot(x1, x2, 
         col = col$n_b, pch = 19,
         xlab = paste0("MDS1 (", var1, "%)"),
         ylab = paste0("MDS2 (", var2, "%)"),
         main = "Sample Type (Nasal/Bronchial)")
    add_categorical_legend("n_b")
  } else {
    plot(x1, x2, 
         col = col$Batch, pch = 19,
         xlab = paste0("MDS1 (", var1, "%)"),
         ylab = paste0("MDS2 (", var2, "%)"),
         main = "Batch Effect")
    add_categorical_legend("Batch")
  }
  
  # Plot 2: Batch
  plot(x1, x2, 
       col = col$Batch, pch = 19,
       xlab = paste0("MDS1 (", var1, "%)"),
       ylab = paste0("MDS2 (", var2, "%)"),
       main = "Batch Effect")
  add_categorical_legend("Batch")
  
  # Plot 3: Sex
  plot(x1, x2, 
       col = col$Sex, pch = 19,
       xlab = paste0("MDS1 (", var1, "%)"),
       ylab = paste0("MDS2 (", var2, "%)"),
       main = "Sex")
  add_categorical_legend("Sex")
  
  # Plot 4: Age
  plot(x1, x2, 
       col = col$Age_at_visit, pch = 19,
       xlab = paste0("MDS1 (", var1, "%)"),
       ylab = paste0("MDS2 (", var2, "%)"),
       main = "Age")
  add_numeric_legend("Age_at_visit")
  
  # Plot 5: BAL Eos > 1%
  plot(x1, x2, 
       col = col$bal_eos_p_mt1, pch = 19,
       xlab = paste0("MDS1 (", var1, "%)"),
       ylab = paste0("MDS2 (", var2, "%)"),
       main = "BAL Eos% >1%")
  add_categorical_legend("bal_eos_p_mt1")
  
  # Plot 6: BAL Eos% (log)
  plot(x1, x2, 
       col = col$BAL_eos_p_log, pch = 19,
       xlab = paste0("MDS1 (", var1, "%)"),
       ylab = paste0("MDS2 (", var2, "%)"),
       main = "BAL Eosinophils % (log)")
  add_numeric_legend("BAL_eos_p_log")
  
  # Plot 7: Blood Eos (log)
  plot(x1, x2, 
       col = col$blood_eos_log, pch = 19,
       xlab = paste0("MDS1 (", var1, "%)"),
       ylab = paste0("MDS2 (", var2, "%)"),
       main = "Blood Eosinophils (log)")
  add_numeric_legend("blood_eos_log")
  
  # Plot 8: Sex (Dims 1 & 4)
  plot(x1, x4, 
       col = col$Sex, pch = 19,
       xlab = paste0("MDS1 (", var1, "%)"),
       ylab = paste0("MDS4 (", var4, "%)"),
       main = "Sex (Dims 1 & 4)")
  add_categorical_legend("Sex")
  
  # Add overall title
  mtext(paste(title_prefix, "Multi-dimensional Scaling Analysis of RNA-seq Data"), 
        outer = TRUE, cex = 1.5)
  
  # Close PDF device
  dev.off()
  
  message(paste("MDS plots saved to", pdf_filename))
  
  # Return MDS coordinates
  return(list(
    x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    var1 = var1, var2 = var2, var3 = var3, var4 = var4
  ))
}

# 6. Function to analyze MDS results and associations with variables
analyze_mds_results <- function(normalized_counts, phen, title_prefix = "") {
  # Calculate MDS
  mds_result <- plotMDS(normalized_counts, plot = FALSE)
  
  # Calculate variance explained
  var_explained_df <- data.frame(
    PC = factor(paste0("MDS", 1:length(mds_result$var.explained)), 
                levels = paste0("MDS", 1:length(mds_result$var.explained))),
    VarExplained = mds_result$var.explained,
    Percentage = round(mds_result$var.explained * 100, 2)
  )
  
  # Print variance explained
  print(paste(title_prefix, "Variance explained by top PCs:"))
  print(var_explained_df[1:10, ])
  
  # Create variance explained plot
  p_var_explained <- ggplot(var_explained_df[1:10, ], aes(x = PC, y = Percentage)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5) +
    theme_minimal() +
    labs(title = paste(title_prefix, "Percentage of Variance Explained by MDS Components"),
         x = "MDS Component",
         y = "Percentage of Variance Explained")
  
  print(p_var_explained)
  
  # Save plot
  ggsave(file.path(output_folder,paste0(title_prefix, "_variance_explained.pdf")), p_var_explained, width = 10, height = 6)
  
  # Extract MDS scores
  mds_scores <- mds_result$eigen.vectors
  
  # Convert to dataframe and add sample IDs
  mds_df <- as.data.frame(mds_scores)
  colnames(mds_df) <- paste0("MDS", seq_len(ncol(mds_df)))
  mds_df$SampleID <- colnames(normalized_counts)
  
  # Merge with metadata
  merged_data <- merge(phen, mds_df, by = "SampleID")
  
  # Define variables for association analysis
  continuous_vars <- c("Age_at_visit_scaled", "admit_count", "ED_visits", "ACT_score", 
                       "FEV1_percent", "BAL_eos_p_log", "blood_eos_p_log")
  categorical_vars <- c("Sex", "Race", "Ethnicity", "bal_eos_p_mt1", "Batch")
  
  # Add n_b to categorical_vars only if analyzing combined dataset
  if (length(unique(merged_data$n_b)) > 1) {
    categorical_vars <- c(categorical_vars, "n_b")
  }
  
  selected_pcs <- colnames(mds_df)[1:5]  # First 5 PCs
  
  # Function for calculating associations
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
  
  # Apply association function to all variables and PCs
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
  
  # Create matrices for visualization
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
  
  # Prepare p-value matrices
  cont_p_matrix_fixed <- prepare_p_matrix(cont_p_matrix)
  cat_p_matrix_fixed <- prepare_p_matrix(cat_p_matrix)
  
  # NEW CODE - Creating heatmaps with p-values and asterisks
  # Function to create annotation text with p-values and asterisks
  create_annotations <- function(association_matrix, pvalue_matrix) {
    # Format p-values and add asterisks for significance
    annotations <- matrix("", nrow = nrow(pvalue_matrix), ncol = ncol(pvalue_matrix),
                          dimnames = dimnames(pvalue_matrix))
    
    for (i in 1:nrow(pvalue_matrix)) {
      for (j in 1:ncol(pvalue_matrix)) {
        if (!is.na(pvalue_matrix[i, j])) {
          # Format p-value based on its magnitude
          p_value <- pvalue_matrix[i, j]
          
          # Use scientific notation for very small p-values
          if (p_value < 0.001) {
            p_formatted <- sprintf("%.2e", p_value)
          } else {
            p_formatted <- sprintf("%.3f", p_value)
          }
          
          # Add asterisk if significant
          if (p_value < 0.05) {
            annotations[i, j] <- paste0(p_formatted, "*")
          } else {
            annotations[i, j] <- p_formatted
          }
        }
      }
    }
    
    return(annotations)
  }
  
  # Create annotation matrices
  cont_annotations <- create_annotations(cont_matrix, cont_p_matrix)
  cat_annotations <- create_annotations(cat_matrix, cat_p_matrix)
  
  # Prepare data frames for ggplot2
  prepare_for_ggplot <- function(matrix_data, p_matrix, annotations, type) {
    df <- data.frame()
    
    for (i in 1:nrow(matrix_data)) {
      for (j in 1:ncol(matrix_data)) {
        if (!is.na(matrix_data[i, j])) {
          df <- rbind(df, data.frame(
            Variable = rownames(matrix_data)[i],
            PC = colnames(matrix_data)[j],
            Value = matrix_data[i, j],
            P_value = p_matrix[i, j],
            Annotation = annotations[i, j],
            Significant = p_matrix[i, j] < 0.05
          ))
        }
      }
    }
    
    df$Variable <- factor(df$Variable, levels = rev(rownames(matrix_data)))
    df$PC <- factor(df$PC, levels = colnames(matrix_data))
    df$Type <- type
    
    return(df)
  }
  
  cont_df <- prepare_for_ggplot(cont_matrix, cont_p_matrix, cont_annotations, "Continuous")
  cat_df <- prepare_for_ggplot(cat_matrix, cat_p_matrix, cat_annotations, "Categorical")
  
  # Create heatmaps with ggplot2
  create_heatmap <- function(df, title, legend_title) {
    ggplot(df, aes(x = PC, y = Variable, fill = Value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Annotation, fontface = ifelse(Significant, "bold", "plain"))) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           name = legend_title) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      labs(
        title = title,
        x = "MDS Component",
        y = "Variable"
      )
  }
  
  # Create and save continuous variables heatmap
  if (nrow(cont_df) > 0) {
    cont_heatmap <- create_heatmap(
      cont_df,
      paste(title_prefix, "Pearson Correlations: Continuous Variables vs PCs"),
      "Correlation"
    )
    
    ggsave(
      file.path(output_folder, paste0(title_prefix, "_continuous_heatmap.pdf")),
      cont_heatmap,
      width = 10,
      height = 8
    )
    
    print(cont_heatmap)
  } else {
    message("No continuous variable associations to plot")
  }
  
  # Create and save categorical variables heatmap
  if (nrow(cat_df) > 0) {
    cat_heatmap <- create_heatmap(
      cat_df,
      paste(title_prefix, "Effect Sizes (η²): Categorical Variables vs PCs"),
      "Effect Size (η²)"
    )
    
    ggsave(
      file.path(output_folder, paste0(title_prefix, "_categorical_heatmap.pdf")),
      cat_heatmap,
      width = 10,
      height = 8
    )
    
    print(cat_heatmap)
  } else {
    message("No categorical variable associations to plot")
  }
  
  dev.off()
  
  # Save results to CSV
  output_csv <- paste0(output_folder,"/MDS_associations", 
                       ifelse(title_prefix == "", "", paste0("_", tolower(title_prefix))), 
                       ".csv")
  write.csv(results_df, output_csv, row.names = FALSE)
  
  message(paste("Association results saved to", output_csv))
  
  # Return both results dataframe and heatmaps
  return(list(
    results = results_df,
    continuous_heatmap = if(nrow(cont_df) > 0) cont_heatmap else NULL,
    categorical_heatmap = if(nrow(cat_df) > 0) cat_heatmap else NULL
  ))
}

# 7. Read phenotype data
phen <- read.csv("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2025-02-14.csv")

# 8. Process phenotype data
phen <- phen %>% 
  # Filter samples with valid BAL cell counts
  filter(BAL_eos_p >= 0 & BAL_neut_p >= 0) %>%
  # Add inflammatory phenotype classifications
  mutate(
    comp1 = factor(case_when(
      BAL_eos_p > 1 & BAL_neut_p > 4 ~ "mixed",
      BAL_eos_p > 1 & BAL_neut_p <= 4 ~ "eos",
      BAL_eos_p <= 1 & BAL_neut_p > 4 ~ "neut",
      BAL_eos_p <= 1 & BAL_neut_p <= 4 ~ "pauci"), 
      levels = c("pauci", "neut", "mixed", "eos")
    ),
    comp2 = factor(case_when(
      BAL_eos_p > 1 ~ "high_eos",
      BAL_eos_p <= 1 ~ "low_eos"), 
      levels = c("high_eos", "low_eos")
    ),
    # Add sample type indicator
    n_b = ifelse(grepl("^B", SampleID), "bronch", "nasal"),
    # Handle missing values for ED visits and hospitalizations
    ED_visits = if_else(asthma_ED == "No", 0, ED_visits),
    admit_count = if_else(asthma_admit_hx == "No", 0, admit_count),
    # Create exacerbation count (sum of ED visits and admissions)
    exac = ED_visits + admit_count,
    # Scale age
    Age_at_visit_scaled = as.numeric(scale(Age_at_visit)),
    # Process IgE values
    total_IgE_num = as.numeric(gsub(">", "", total_IgE_repeat)),
    total_IgE_num = if_else(grepl(">", total_IgE_repeat), 5000, total_IgE_num),
    # Flag for BAL eosinophils > 1%
    bal_eos_p_mt1 = factor(BAL_eos_p > 1)
  )

# Factor the batch variable
phen$Batch <- factor(phen$Batch, levels = c("batch1", "batch2", "batch3", "batch4", "batch6"))

# Create separate dataframes for bronchial and nasal samples
phen_bronch <- phen[grepl("^B", phen$SampleID), ]
phen_nasal <- phen[grepl("^N", phen$SampleID), ]
phen_nasal <- phen_nasal[!grepl("F", phen_nasal$SampleID), ]

print(paste("Number of samples after filtering:", nrow(phen)))


# 9. Prepare variables for MDS color coding
var <- c(source.cell.log, "bal_eos_p_mt1", "Batch", "Age_at_visit", "Sex", "Race", "Ethnicity", "total_IgE_num", "n_b")
phen_var <- phen[, var]

# Convert categorical variables to factors
phen_var[, c("bal_eos_p_mt1", "Batch", "Sex", "Race", "Ethnicity", "n_b")] <- lapply(
  phen_var[, c("bal_eos_p_mt1", "Batch", "Sex", "Race", "Ethnicity", "n_b")],
  function(d) factor(d, levels = unique(d))
)

# Round numeric variables
for (n in var) {
  if (!is.factor(phen_var[, n])) {
    phen_var[, n] <- round(phen_var[, n], 1)
  }
}

#===============================
# analysis for batch-corrected normalized count data
#===============================
# 1. load count files
scaled_count_folder<-"normalized_gene_count"
scaled_count_file_n<-"normalized_gene_count_nasal_vsd_batch-corrected_2025-06-09.txt"
scaled_count_file_b<-"normalized_gene_count_bronch_vsd_batch-corrected_2025-06-09.txt"

# check if the files exist
print(file.exists(file.path(resources_dir, 
                            scaled_count_folder, 
                            scaled_count_file_n)))
print(file.exists(file.path(resources_dir, 
                            scaled_count_folder, 
                            scaled_count_file_b)))

# set the norm count file locations
normalized_count_path <- file.path(resources_dir, 
                                   scaled_count_folder, 
                                   c(scaled_count_file_n,scaled_count_file_b))


normalized_counts_n <- read.table(normalized_count_path[1], header = TRUE, row.names = 1, sep = "\t")
normalized_counts_b <- read.table(normalized_count_path[2], header = TRUE, row.names = 1, sep = "\t")

# 2. load normalized counts parameter information
normalization_parameter_file<-c("normalization_parameters_2025-06-09.txt")

norm_parameter<-read.delim(file.path(resources_dir, 
                                     scaled_count_folder, 
                                     normalization_parameter_file),header=FALSE,sep="\r")

print(norm_parameter)

# 3. set output directory
output_folder<-output_folder_batch_corr

# 13. Run MDS analysis and visualization for all samples
# message("Running MDS analysis for all samples...")
# col_all <- create_color_mappings(phen_var)
# mds_results_all <- create_multipanel_mds(normalized_counts, phen_var, col_all, "Combined")
# association_results_all <- analyze_mds_results(normalized_counts, phen, "Combined")

# 4. Run MDS analysis for bronchial samples only
message("Running MDS analysis for bronchial samples...")
phen_var_bronch <- filter(phen_var, n_b == "bronch")
col_bronch <- create_color_mappings(phen_var_bronch)
mds_results_bronch <- create_multipanel_mds(normalized_counts_b, phen_var_bronch, col_bronch, "Bronchial", output_folder)
association_results_bronch <- analyze_mds_results(normalized_counts_b, phen_bronch, "Bronchial")$results

# 5. Run MDS analysis for nasal samples only
message("Running MDS analysis for nasal samples...")
phen_var_nasal <- filter(phen_var, n_b == "nasal")
col_nasal <- create_color_mappings(phen_var_nasal)
mds_results_nasal <- create_multipanel_mds(normalized_counts_n, phen_var_nasal, col_nasal, "Nasal",output_folder)
association_results_nasal <- analyze_mds_results(normalized_counts_n, phen_nasal, "Nasal")$results

# 6. Print summary of key findings
message("Analysis completed!")
message("Summary of key findings:")
message("1. Variance explained by top MDS components:")
# message("   - Combined: MDS1: ", mds_results_all$var1, "%, MDS2: ", mds_results_all$var2, "%")
message("   - Bronchial: MDS1: ", mds_results_bronch$var1, "%, MDS2: ", mds_results_bronch$var2, "%")
message("   - Nasal: MDS1: ", mds_results_nasal$var1, "%, MDS2: ", mds_results_nasal$var2, "%")

# message("2. Top significant associations with MDS1:")
# top_assoc_all <- association_results_all[association_results_all$PC == "MDS1" & 
#                                            association_results_all$Significant == TRUE, ]
# if(nrow(top_assoc_all) > 0) {
#   top_assoc_all <- top_assoc_all[order(-abs(top_assoc_all$Association)), ]
#   message("   - Combined: ", paste(head(top_assoc_all$Variable, 3), collapse = ", "))
# } else {
#   message("   - Combined: No significant associations found")
# }

top_assoc_bronch <- association_results_bronch[association_results_bronch$Significant == TRUE, ]
if(nrow(top_assoc_bronch) > 0) {
  top_assoc_bronch <- top_assoc_bronch[order(-abs(top_assoc_bronch$Association)), ]
  message("   - Bronchial: ", paste(head(top_assoc_bronch$Variable, 3), collapse = ", "))
} else {
  message("   - Bronchial: No significant associations found")
}

top_assoc_nasal <- association_results_nasal[association_results_nasal$results$Significant == TRUE, ]
if(nrow(top_assoc_nasal) > 0) {
  top_assoc_nasal <- top_assoc_nasal[order(-abs(top_assoc_nasal$Association)), ]
  message("   - Nasal: ", paste(head(top_assoc_nasal$Variable, 3), collapse = ", "))
} else {
  message("   - Nasal: No significant associations found")
}


#===============================
# analysis for non-batch-corrected normalized count data
#===============================
# 1. load count files
scaled_count_folder<-"normalized_gene_count"
scaled_count_file_n<-"normalized_gene_count_nasal_vsd_no-batch-corrected_2025-06-09.txt"
scaled_count_file_b<-"normalized_gene_count_bronch_vsd_no-batch-corrected_2025-06-09.txt"

# check if the files exist
print(file.exists(file.path(resources_dir, 
                            scaled_count_folder, 
                            scaled_count_file_n)))
print(file.exists(file.path(resources_dir, 
                            scaled_count_folder, 
                            scaled_count_file_b)))

# set the norm count file locations
normalized_count_path <- file.path(resources_dir, 
                                   scaled_count_folder, 
                                   c(scaled_count_file_n,scaled_count_file_b))


normalized_counts_n <- read.table(normalized_count_path[1], header = TRUE, row.names = 1, sep = "\t")
normalized_counts_b <- read.table(normalized_count_path[2], header = TRUE, row.names = 1, sep = "\t")

# 2. load normalized counts parameter information
normalization_parameter_file<-c("normalization_parameters_2025-06-09.txt")

norm_parameter<-read.delim(file.path(resources_dir, 
                                     scaled_count_folder, 
                                     normalization_parameter_file),header=FALSE,sep="\r")

print(norm_parameter)

# 3. set output directory
output_folder<-output_folder_no_batch_corr

# 13. Run MDS analysis and visualization for all samples
# message("Running MDS analysis for all samples...")
# col_all <- create_color_mappings(phen_var)
# mds_results_all <- create_multipanel_mds(normalized_counts, phen_var, col_all, "Combined")
# association_results_all <- analyze_mds_results(normalized_counts, phen, "Combined")

# 4. Run MDS analysis for bronchial samples only
message("Running MDS analysis for bronchial samples...")
phen_var_bronch <- filter(phen_var, n_b == "bronch")
col_bronch <- create_color_mappings(phen_var_bronch)
mds_results_bronch <- create_multipanel_mds(normalized_counts_b, phen_var_bronch, col_bronch, "Bronchial",output_folder)
association_results_bronch <- analyze_mds_results(normalized_counts_b, phen_bronch, "Bronchial")$results

# 5. Run MDS analysis for nasal samples only
message("Running MDS analysis for nasal samples...")
phen_var_nasal <- filter(phen_var, n_b == "nasal")
col_nasal <- create_color_mappings(phen_var_nasal)
mds_results_nasal <- create_multipanel_mds(normalized_counts_n, phen_var_nasal, col_nasal, "Nasal",output_folder)
association_results_nasal <- analyze_mds_results(normalized_counts_n, phen_nasal, "Nasal")$results

# 6. Print summary of key findings
message("Analysis completed!")
message("Summary of key findings:")
message("1. Variance explained by top MDS components:")
# message("   - Combined: MDS1: ", mds_results_all$var1, "%, MDS2: ", mds_results_all$var2, "%")
message("   - Bronchial: MDS1: ", mds_results_bronch$var1, "%, MDS2: ", mds_results_bronch$var2, "%")
message("   - Nasal: MDS1: ", mds_results_nasal$var1, "%, MDS2: ", mds_results_nasal$var2, "%")

# message("2. Top significant associations with MDS1:")
# top_assoc_all <- association_results_all[association_results_all$PC == "MDS1" & 
#                                            association_results_all$Significant == TRUE, ]
# if(nrow(top_assoc_all) > 0) {
#   top_assoc_all <- top_assoc_all[order(-abs(top_assoc_all$Association)), ]
#   message("   - Combined: ", paste(head(top_assoc_all$Variable, 3), collapse = ", "))
# } else {
#   message("   - Combined: No significant associations found")
# }

top_assoc_bronch <- association_results_bronch[association_results_bronch$Significant == TRUE, ]
if(nrow(top_assoc_bronch) > 0) {
  top_assoc_bronch <- top_assoc_bronch[order(-abs(top_assoc_bronch$Association)), ]
  message("   - Bronchial: ", paste(head(top_assoc_bronch$Variable, 3), collapse = ", "))
} else {
  message("   - Bronchial: No significant associations found")
}

top_assoc_nasal <- association_results_nasal[association_results_nasal$results$Significant == TRUE, ]
if(nrow(top_assoc_nasal) > 0) {
  top_assoc_nasal <- top_assoc_nasal[order(-abs(top_assoc_nasal$Association)), ]
  message("   - Nasal: ", paste(head(top_assoc_nasal$Variable, 3), collapse = ", "))
} else {
  message("   - Nasal: No significant associations found")
}


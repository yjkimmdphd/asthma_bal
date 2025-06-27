##
# Investigate relationships among the cell counts in BAL and blood
## 

# Load required libraries for data manipulation, visualization, and statistical analysis
library(tidyverse)    # Data manipulation and visualization
library(DESeq2)       # Differential expression analysis
library(limma)        # Linear models for microarray data
library(edgeR)        # Empirical analysis of digital gene expression
library(corrplot)     # Correlation plot visualization
library(Hmisc)        # Statistical functions including correlation analysis

# Additional library for data manipulation
library(dplyr)

################################
## Data Loading and Preprocessing
################################

# Load gene expression count data from MS asthma study
countdata <- file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts <- if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}

# Set sample IDs as row names and extract column names (sample IDs)
rownames(counts) <- counts[,"SampleID"]
counts.ID <- colnames(counts)

# Load phenotype data containing asthma biomarkers and cell counts
phenotype <- file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2025-02-14.csv")
phenotype <- if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}

# Load sampling date differences between diffe
# rent sample types (swab, BAL, CBC)
sampling_date_diff <- "./resources/processed_data/sampling_dates/swab-bal-cbc_differences_in_days.txt"
sampling_date_diff <- if(file.exists(sampling_date_diff)){read.table(sampling_date_diff, row.names = NULL, header = TRUE)}

# Filter for blood-BAL sampling date comparisons only
sampling_date_diff <- sampling_date_diff %>% filter(Comparison == "blood_bal")
colnames(sampling_date_diff)[1:3] <- c("ID", "sampling_date_comp", "sampling_date_diff_days")

# Merge phenotype data with sampling date differences
phenotype <- left_join(phenotype, sampling_date_diff, by = "ID")

# Filter to include only samples where blood and BAL were collected within 30 days of each other
# This ensures temporal proximity for meaningful comparisons
phenotype_filtered <- phenotype %>% filter(abs(sampling_date_diff_days) < 30)
phenotype <- phenotype_filtered

################################
## Variable Definitions
################################

# Define vectors of cell count variables (log-transformed) for analysis
# BAL = Bronchoalveolar Lavage, measured from lung fluid
# blood = peripheral blood measurements
source.cell.log <- c(
  "BAL_eos_ct_log",     # BAL eosinophil count (log)
  "BAL_eos_p_log",      # BAL eosinophil percentage (log)
  "BAL_neut_ct_log",    # BAL neutrophil count (log)
  "BAL_neut_p_log",     # BAL neutrophil percentage (log)
  "BAL_wbc_log",        # BAL white blood cell count (log)
  "blood_eos_log",      # Blood eosinophil count (log)
  "blood_eos_p_log",    # Blood eosinophil percentage (log)
  "blood_neut_log",     # Blood neutrophil count (log)
  "blood_neut_p_log",   # Blood neutrophil percentage (log)
  "blood_wbc_log")      # Blood white blood cell count (log)

# Define vectors of cell count variables (original scale)
source.cell <- c(
  "BAL_eos_ct",         # BAL eosinophil count
  "BAL_eos_p",          # BAL eosinophil percentage
  "BAL_neut_ct",        # BAL neutrophil count
  "BAL_neut_p",         # BAL neutrophil percentage
  "BAL_wbc",            # BAL white blood cell count
  "blood_eos",          # Blood eosinophil count
  "blood_eos_p",        # Blood eosinophil percentage
  "blood_neut",         # Blood neutrophil count
  "blood_neut_p",       # Blood neutrophil percentage
  "blood_wbc")          # Blood white blood cell count

###########################################################################################
## Data Subset for Bronchial RNAseq Analysis
###########################################################################################

# Identify subjects who had both BAL sampling and bronchial RNAseq completed
bexist <- phenotype$SampleID %in% counts.ID

# Filter phenotype data for existing bronchial samples (IDs starting with "B")
bphen <- phenotype[bexist,] %>% filter(grepl("^B", SampleID))

# Standardize (z-score) the log-transformed cell count variables for comparison
bphen <- mutate_at(bphen, vars(all_of(source.cell.log)), scale)

################################
## Correlation Analysis
################################

# Extract cell count data matrix for correlation analysis
data <- bphen[, source.cell.log]
rownames(data) <- bphen$SampleID

# Calculate correlation matrix and p-values using Pearson correlation
# rcorr() provides both correlation coefficients and significance testing
correlation_results <- Hmisc::rcorr(as.matrix(data))
cor_matrix <- correlation_results$r  # Correlation coefficients matrix
p_matrix <- correlation_results$P    # P-values matrix

################################
## Significance Level Functions
################################

# Function to create significance labels based on p-values
# Converts p-values to asterisk notation for visualization
get_sig_labels <- function(p_matrix, sig.level = c(0.001, 0.01, 0.05)) {
  # Initialize empty matrix for significance labels
  sig_labels <- matrix("", nrow = nrow(p_matrix), ncol = ncol(p_matrix))
  
  # Assign significance levels: *** (p≤0.001), ** (p≤0.01), * (p≤0.05)
  sig_labels[p_matrix <= sig.level[1]] <- "***"
  sig_labels[p_matrix <= sig.level[2] & p_matrix > sig.level[1]] <- "**"
  sig_labels[p_matrix <= sig.level[3] & p_matrix > sig.level[2]] <- "*"
  sig_labels
}

# Generate significance labels for the correlation matrix
sig_labels <- get_sig_labels(p_matrix)

################################
## Correlation Visualization
################################

# Create correlation heatmap showing all correlations (regardless of significance)
cell_count_cor_heatmap <- corrplot::corrplot(cor_matrix, 
                                             type = "upper",           # Show upper triangle only
                                             order = "hclust",         # Order by hierarchical clustering
                                             tl.col = "black",         # Text label color
                                             tl.srt = 45,             # Text label rotation
                                             method = "ellipse",       # Use ellipses to show correlation
                                             addCoef.col = "black",    # Add correlation coefficients
                                             diag = FALSE,             # Don't show diagonal
                                             title = "Correlation Matrix, blood-bal sampling date difference < 30d, no significance",
                                             mar = c(5,5,5,5))         # Plot margins

# Create correlation heatmap showing only significant correlations (p < 0.05)
corrplot::corrplot(cor_matrix, 
                   type = "upper", 
                   order = "hclust",
                   tl.col = "black", 
                   tl.srt = 45,
                   method = "circle",            # Use circles instead of ellipses
                   addCoef.col = "black",        # Add correlation coefficients
                   p.mat = p_matrix,             # Provide p-value matrix
                   sig.level = 0.05,             # Significance threshold
                   insig = "n",                  # Hide non-significant correlations
                   addCoefasPercent = FALSE,     # Show coefficients as decimals
                   cl.cex = 0.8,                # Color legend text size
                   diag = FALSE,
                   title = "Correlation Matrix with Significance Levels",
                   mar = c(5,5,5,5))

################################
## Individual Scatter Plot Generation
################################

# Load additional libraries for plot arrangement
library(gridExtra)
library(ggpubr)

# Initialize list to store individual scatter plots
plot_list <- list()
num_columns <- ncol(data)

# Generate scatter plots for all pairwise combinations of variables
for (i in 1:(num_columns - 1)) {
  for (j in (i + 1):num_columns) {
    # Extract correlation coefficient and p-value for this pair
    cor_value <- cor_matrix[i, j]
    p_value <- p_matrix[i, j]
    
    # Create scatter plot with regression line
    p <- ggplot(data, aes_string(x = names(data)[i], y = names(data)[j])) +
      geom_point(alpha = 0.4) +                    # Add semi-transparent points
      geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add regression line
      theme_minimal() +                            # Clean theme
      labs(title = paste("Scatter plot of", names(data)[i], "vs", names(data)[j]),
           x = names(data)[i], 
           y = names(data)[j]) +
      # Annotate with correlation coefficient and p-value
      annotate("text", x = Inf, y = Inf, 
               label = sprintf("r = %.2f, p = %.3f", cor_value, p_value),
               hjust = 1.1, vjust = 1.1, size = 3, color = "red")
    
    # Add plot to list with descriptive name
    plot_list[[paste(names(data)[i], names(data)[j], sep = "_vs_")]] <- p
  }
}

################################
## Plot Display and Arrangement
################################

# Display first plot as example
print(plot_list[[1]])

# Load patchwork for plot arrangement
library(patchwork)

# Arrange all plots in a grid layout
plot_layout <- wrap_plots(plot_list, ncol = 3)
plot_layout

# Display plots in groups of 15 for better visualization
# Group 1: First 15 scatter plots
wrap_plots(plot_list[1:15], ncol = 3)

# Group 2: Next 15 scatter plots  
wrap_plots(plot_list[16:30], ncol = 3)

# Group 3: Final 15 scatter plots
wrap_plots(plot_list[31:45], ncol = 3)

# Display the correlation heatmap
cell_count_cor_heatmap
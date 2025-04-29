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

# 2. make normalized count data

# -----------------------------------------------------------------------------
# Load phenotype data
# -----------------------------------------------------------------------------
phen_path <- file.path(
  "./resources/processed_data",
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv"
)
phen <- read.csv(phen_path)%>% filter(BAL_eos_p >= 0 & BAL_neut_p >= 0) 
phen$Batch<-factor(phen$Batch,levels=unique(phen$Batch))

phen_bronch<-phen[grepl("^B",phen$SampleID),]
phen_nasal<-phen[grepl("^N",phen$SampleID),]
phen_nasal<-phen_nasal[-grep("F", phen_nasal$SampleID),]

# load count table
countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
genes<-counts[,"SampleID"]

# 
phen_input<-phen_bronch

# select bronchial samples 
sample_id<-phen_input$SampleID
counts_selected<-counts[,sample_id]
rownames(counts_selected)<-genes
head(counts_selected)

# Assuming your DESeq2 object is called 'dds'
library(DESeq2)
library(vsn)
countdata<-counts_selected
coldata_cols<-c("comp1","Batch")
coldata<-phen_input[,coldata_cols]
rownames(coldata)<-sample_id

dds<-DESeqDataSetFromMatrix(countData = countdata,colData=coldata, design= ~ comp1+Batch)

# prefilter low count genes
smallestGroupSize <- min(apply(table(coldata),1,sum)) 
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

vsd <- vst(dds, blind=FALSE) 

meanSdPlot(assay(vsd))

## vsd with vs without batch effect removal 

### vsd without batch effect removal 
normalized_counts <- assay(vsd)  # This is now your transformed expression matrix
sampleDists <- dist(t(assay(vsd)))

library(pheatmap)
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Batch
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(vsd, intgroup=c("Batch"))
plotPCA(vsd, intgroup=c("comp1"))

### vsd with batch effect removed 
library(limma)
mat <- assay(vsd)
mm <- model.matrix(~comp1, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
assay(vsd) <- mat
meanSdPlot(assay(vsd))

normalized_counts <- assay(vsd)  # This is now your transformed expression matrix
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Batch
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(vsd, intgroup=c("Batch"))
plotPCA(vsd, intgroup=c("comp1"))

### will only export vsd with batch effect removal





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
var <- c(source.cell.log,"bal_eos_p_mt1", "Batch", "Age_at_visit", "Sex", "Race", "Ethnicity", "total_IgE_num", "n_b")
phen_var <- phen[, var]

# Convert categorical variables to factors
phen_var[, c("bal_eos_p_mt1","Batch", "Sex", "Race", "Ethnicity", "n_b")] <- lapply(
  phen_var[, c("bal_eos_p_mt1","Batch", "Sex", "Race", "Ethnicity", "n_b")],
  function(d) factor(d, levels = unique(d))
)

# Round numeric variables
for (n in var) {
  if (!is.factor(phen_var[, n])) {
    phen_var[, n] <- round(phen_var[, n], 1)
  }
}

# 5. Create multipanel MDS plot for RNA-seq data analysis
# This script should be run after your initial MDS analysis
# Create a MDS multipanel plot with outside legends using Base R
# This works directly with your existing color dataframe 'col'

# Function to perform MDS and create plots with outside legends
create_multipanel_mds <- function(normalized_counts, phen_var, col) {
  # Setup 4x2 panel layout with extra space for legends
  pdf("mds_analysis_with_legends.pdf", width=14, height=14)
  # Set layout with wider right margin for legends
  par(mfrow=c(4, 2), mar=c(4, 4, 3, 8), oma=c(0, 0, 2, 0))
  
  # Calculate MDS
  mds <- limma::plotMDS(normalized_counts, ndim=5, plot=FALSE)
  
  # Extract coordinates for all dimensions we'll use
  x1 <- mds$eigen.vectors[,1] * sqrt(mds$eigen.values[1])
  x2 <- mds$eigen.vectors[,2] * sqrt(mds$eigen.values[2])
  x3 <- mds$eigen.vectors[,3] * sqrt(mds$eigen.values[3])
  x4 <- mds$eigen.vectors[,4] * sqrt(mds$eigen.values[4])
  
  # Get variance explained
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
    
    for(val in unique_values) {
      # Find first index with this value
      idx <- which(phen_var[[var_name]] == val)[1]
      unique_colors <- c(unique_colors, col[[var_name]][idx])
    }
    
    # Include NA in the legend if there are any NA values
    if(any(is.na(phen_var[[var_name]]))) {
      legend_values <- c(as.character(unique_values), "NA")
      legend_colors <- c(unique_colors, na_color)
    } else {
      legend_values <- as.character(unique_values)
      legend_colors <- unique_colors
    }
    
    # Add legend outside the plot
    par(xpd=TRUE)  # Allow drawing outside the plot region
    legend(par("usr")[2] * 1.05, par("usr")[4], 
           legend=legend_values, 
           fill=legend_colors, 
           title=var_name,
           cex=0.7, bty="n")
    par(xpd=FALSE)  # Reset to default
  }
  
  # Helper function to add legends for numeric variables
  add_numeric_legend <- function(var_name) {
    # Use red gradient colors for numeric variables
    quart_colors <- c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D")
    
    # Calculate quartiles for the variable
    var_quants <- quantile(phen_var[[var_name]], probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE)
    
    # Create legend labels (rounded to 1 decimal place)
    quart_labels <- c(
      paste0("Q1: ", round(var_quants[1], 1), "-", round(var_quants[2], 1)),
      paste0("Q2: ", round(var_quants[2], 1), "-", round(var_quants[3], 1)),
      paste0("Q3: ", round(var_quants[3], 1), "-", round(var_quants[4], 1)),
      paste0("Q4: ", round(var_quants[4], 1), "-", round(var_quants[5], 1))
    )
    
    # Include NA in the legend if there are any NA values
    if(any(is.na(phen_var[[var_name]]))) {
      legend_labels <- c(quart_labels, "NA")
      legend_colors <- c(quart_colors, na_color)
    } else {
      legend_labels <- quart_labels
      legend_colors <- quart_colors
    }
    
    # Add legend outside the plot
    par(xpd=TRUE)  # Allow drawing outside the plot region
    legend(par("usr")[2] * 1.05, par("usr")[4], 
           legend=legend_labels, 
           fill=legend_colors, 
           title=var_name,
           cex=0.7, bty="n")
    par(xpd=FALSE)  # Reset to default
  }
  
  # Plot 1: Sample Type
  plot(x1, x2, 
       col=col$n_b, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Sample Type (Nasal/Bronchial)")
  add_categorical_legend("n_b")
  
  # Plot 2: Batch
  plot(x1, x2, 
       col=col$Batch, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Batch Effect")
  add_categorical_legend("Batch")
  
  # Plot 3: Sex
  plot(x1, x2, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Sex")
  add_categorical_legend("Sex")
  
  # Plot 4: Age
  plot(x1, x2, 
       col=col$Age_at_visit, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Age")
  add_numeric_legend("Age_at_visit")
  
  # Plot 5: BAL Eos > 1%
  plot(x1, x2, 
       col=col$bal_eos_p_mt1, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eos% >1%")
  add_categorical_legend("bal_eos_p_mt1")
  
  # Plot 6: BAL Eos% (log)
  plot(x1, x2, 
       col=col$BAL_eos_p_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eosinophils % (log)")
  add_numeric_legend("BAL_eos_p_log")
  
  # Plot 7: Blood Eos (log)
  plot(x1, x2, 
       col=col$blood_eos_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Blood Eosinophils (log)")
  add_numeric_legend("blood_eos_log")
  
  # Plot 8: Sex (Dims 1 & 4)
  plot(x1, x4, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS4 (", var4, "%)"),
       main="Sex (Dims 1 & 4)")
  add_categorical_legend("Sex")
  
  # Add overall title
  mtext("Multi-dimensional Scaling Analysis of RNA-seq Data", 
        outer=TRUE, cex=1.5)
  
  # Close PDF device
  dev.off()
  
  # Also display on screen
  par(mfrow=c(4, 2), mar=c(4, 4, 3, 8), oma=c(0, 0, 2, 0))
  
  # Plot 1: Sample Type
  plot(x1, x2, 
       col=col$n_b, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Sample Type (Nasal/Bronchial)")
  add_categorical_legend("n_b")
  
  # Plot 2: Batch
  plot(x1, x2, 
       col=col$Batch, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Batch Effect")
  add_categorical_legend("Batch")
  
  # Plot 3: Sex
  plot(x1, x2, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Sex")
  add_categorical_legend("Sex")
  
  # Plot 4: Age
  plot(x1, x2, 
       col=col$Age_at_visit, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Age")
  add_numeric_legend("Age_at_visit")
  
  # Plot 5: BAL Eos > 1%
  plot(x1, x2, 
       col=col$bal_eos_p_mt1, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eos% >1%")
  add_categorical_legend("bal_eos_p_mt1")
  
  # Plot 6: BAL Eos% (log)
  plot(x1, x2, 
       col=col$BAL_eos_p_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eosinophils % (log)")
  add_numeric_legend("BAL_eos_p_log")
  
  # Plot 7: Blood Eos (log)
  plot(x1, x2, 
       col=col$blood_eos_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Blood Eosinophils (log)")
  add_numeric_legend("blood_eos_log")
  
  # Plot 8: Sex (Dims 1 & 4)
  plot(x1, x4, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS4 (", var4, "%)"),
       main="Sex (Dims 1 & 4)")
  add_categorical_legend("Sex")
  
  # Add overall title
  mtext("Multi-dimensional Scaling Analysis of RNA-seq Data", 
        outer=TRUE, cex=1.5)
  
  # Return MDS coordinates in case they're needed elsewhere
  return(list(
    x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    var1 = var1, var2 = var2, var3 = var3, var4 = var4
  ))
}

# Run the function
mds_results <- create_multipanel_mds(normalized_counts, phen_var, col)
# To save the plot (uncomment to use)
# ggsave("multipanel_mds_analysis.pdf", grid.arrange(grobs=multipanel_plots), width=12, height=10)

# Alternative using patchwork package for more control:
# If you prefer patchwork, uncomment and install if needed:
# install.packages("patchwork")
# library(patchwork)
# combined_plot <- (multipanel_plots$p_sample_type + multipanel_plots$p_batch) / 
#                  (multipanel_plots$p_sex + multipanel_plots$p_age) /
#                  (multipanel_plots$p_blood_eos + multipanel_plots$p_blood_neut) /
#                  (multipanel_plots$p_sex_dim14 + multipanel_plots$p_sample_dim23) +
#                  plot_annotation(title = "Multi-dimensional Scaling Analysis of RNA-seq Data")
# 
# ggsave("multipanel_mds_analysis_patchwork.pdf", combined_plot, width=12, height=16)
# If you want to create just the PDF without displaying on screen:
# create_multipanel_pdf <- function() {
#   pdf("mds_analysis_multipanel.pdf", width=12, height=14)
#   par(mfrow=c(4,2), mar=c(4,4,3,1), oma=c(0,0,2,0))
#   
#   # [Copy all the plot commands from create_multipanel_mds_base here]
#   
#   dev.off()
# }
# create_multipanel_pdf()
### ===============================

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
library(corrplot)
par(mfrow=c(2,2))
if (requireNamespace("corrplot", quietly = TRUE)) {
  # Correlation plot for continuous variables
  corrplot(cont_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           p.mat = cont_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           cl.offset = -0.5,
           cl.align.text = "r",
           mar = c(4, 1, 4, 0),
           title = "Pearson Correlations: \nContinuous Variables vs PCs")
  
  corrplot(cont_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           p.mat = cont_p_matrix, sig.level = 0.05,
           insig = "label_sig", 
           cl.offset = -0.5,
           cl.align.text = "r",
           mar = c(4, 1, 4, 0),
           title = "Pearson Correlations: \nContinuous Variables vs PCs")
  
  # Effect size plot for categorical variables
  corrplot(cat_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           col = colorRampPalette(c("white", "pink", "red"))(100),
           p.mat = cat_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           mar = c(4, 1, 4, 0),
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
           mar = c(4, 1, 4, 0),
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

# Create multipanel MDS plot for RNA-seq data analysis

# Create a MDS multipanel plot with outside legends using Base R
# This works directly with your existing color dataframe 'col'

# Function to perform MDS and create plots with outside legends
create_multipanel_mds <- function(normalized_counts, phen_var, col) {
  # Setup 4x2 panel layout with extra space for legends
  pdf("mds_analysis_with_legends.pdf", width=14, height=14)
  # Set layout with wider right margin for legends
  par(mfrow=c(4, 2), mar=c(4, 4, 3, 8), oma=c(0, 0, 2, 0))
  
  # Calculate MDS
  mds <- limma::plotMDS(normalized_counts, ndim=5, plot=FALSE)
  
  # Extract coordinates for all dimensions we'll use
  x1 <- mds$eigen.vectors[,1] * sqrt(mds$eigen.values[1])
  x2 <- mds$eigen.vectors[,2] * sqrt(mds$eigen.values[2])
  x3 <- mds$eigen.vectors[,3] * sqrt(mds$eigen.values[3])
  x4 <- mds$eigen.vectors[,4] * sqrt(mds$eigen.values[4])
  
  # Get variance explained
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
    
    for(val in unique_values) {
      # Find first index with this value
      idx <- which(phen_var[[var_name]] == val)[1]
      unique_colors <- c(unique_colors, col[[var_name]][idx])
    }
    
    # Include NA in the legend if there are any NA values
    if(any(is.na(phen_var[[var_name]]))) {
      legend_values <- c(as.character(unique_values), "NA")
      legend_colors <- c(unique_colors, na_color)
    } else {
      legend_values <- as.character(unique_values)
      legend_colors <- unique_colors
    }
    
    # Add legend outside the plot
    par(xpd=TRUE)  # Allow drawing outside the plot region
    legend(par("usr")[2] * 1.05, par("usr")[4], 
           legend=legend_values, 
           fill=legend_colors, 
           title=var_name,
           cex=0.7, bty="n")
    par(xpd=FALSE)  # Reset to default
  }
  
  # Helper function to add legends for numeric variables
  add_numeric_legend <- function(var_name) {
    # Use red gradient colors for numeric variables
    quart_colors <- c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D")
    
    # Calculate quartiles for the variable
    var_quants <- quantile(phen_var[[var_name]], probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE)
    
    # Create legend labels (rounded to 1 decimal place)
    quart_labels <- c(
      paste0("Q1: ", round(var_quants[1], 1), "-", round(var_quants[2], 1)),
      paste0("Q2: ", round(var_quants[2], 1), "-", round(var_quants[3], 1)),
      paste0("Q3: ", round(var_quants[3], 1), "-", round(var_quants[4], 1)),
      paste0("Q4: ", round(var_quants[4], 1), "-", round(var_quants[5], 1))
    )
    
    # Include NA in the legend if there are any NA values
    if(any(is.na(phen_var[[var_name]]))) {
      legend_labels <- c(quart_labels, "NA")
      legend_colors <- c(quart_colors, na_color)
    } else {
      legend_labels <- quart_labels
      legend_colors <- quart_colors
    }
    
    # Add legend outside the plot
    par(xpd=TRUE)  # Allow drawing outside the plot region
    legend(par("usr")[2] * 1.05, par("usr")[4], 
           legend=legend_labels, 
           fill=legend_colors, 
           title=var_name,
           cex=0.7, bty="n")
    par(xpd=FALSE)  # Reset to default
  }
  
  # Plot 1: Batch
  plot(x1, x2, 
       col=col$Batch, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Batch Effect")
  add_categorical_legend("Batch")
  
  # Plot 2: Sex
  plot(x1, x2, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Sex")
  add_categorical_legend("Sex")
  
  # Plot 3: Age
  plot(x1, x2, 
       col=col$Age_at_visit, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Age")
  add_numeric_legend("Age_at_visit")
  
  # Plot 4: BAL Eos > 1%
  plot(x1, x2, 
       col=col$bal_eos_p_mt1, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eos% >1%")
  add_categorical_legend("bal_eos_p_mt1")
  
  # Plot 5: BAL Eos% (log)
  plot(x1, x2, 
       col=col$BAL_eos_p_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eosinophils % (log)")
  add_numeric_legend("BAL_eos_p_log")
  
  # Plot 6: Blood Eos (log)
  plot(x1, x2, 
       col=col$blood_eos_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Blood Eosinophils (log)")
  add_numeric_legend("blood_eos_log")
  
  # Plot 7: Sex (Dims 1 & 4)
  plot(x1, x4, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS4 (", var4, "%)"),
       main="Sex (Dims 1 & 4)")
  add_categorical_legend("Sex")
  
  # Add overall title
  mtext("Multi-dimensional Scaling Analysis of RNA-seq Data", 
        outer=TRUE, cex=1.5)
  
  # Close PDF device
  dev.off()
  
  # Also display on screen
  par(mfrow=c(4, 2), mar=c(4, 4, 3, 8), oma=c(0, 0, 2, 0))
  
  # Plot 1: Batch
  plot(x1, x2, 
       col=col$Batch, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Batch Effect")
  add_categorical_legend("Batch")
  
  # Plot 2: Sex
  plot(x1, x2, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Sex")
  add_categorical_legend("Sex")
  
  # Plot 3: Age
  plot(x1, x2, 
       col=col$Age_at_visit, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Age")
  add_numeric_legend("Age_at_visit")
  
  # Plot 4: BAL Eos > 1%
  plot(x1, x2, 
       col=col$bal_eos_p_mt1, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eos% >1%")
  add_categorical_legend("bal_eos_p_mt1")
  
  # Plot 5: BAL Eos% (log)
  plot(x1, x2, 
       col=col$BAL_eos_p_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eosinophils % (log)")
  add_numeric_legend("BAL_eos_p_log")
  
  # Plot 6: Blood Eos (log)
  plot(x1, x2, 
       col=col$blood_eos_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Blood Eosinophils (log)")
  add_numeric_legend("blood_eos_log")
  
  # Plot 7: Sex (Dims 1 & 4)
  plot(x1, x4, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS4 (", var4, "%)"),
       main="Sex (Dims 1 & 4)")
  add_categorical_legend("Sex")
  
  # Add overall title
  mtext("Multi-dimensional Scaling Analysis of RNA-seq Data", 
        outer=TRUE, cex=1.5)
  
  # Return MDS coordinates in case they're needed elsewhere
  return(list(
    x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    var1 = var1, var2 = var2, var3 = var3, var4 = var4
  ))
}

# Run the function
mds_results <- create_multipanel_mds(normalized_counts, phen_var, col)
# To save the plot (uncomment to use)
# ggsave("multipanel_mds_analysis.pdf", grid.arrange(grobs=multipanel_plots), width=12, height=10)

# Alternative using patchwork package for more control:
# If you prefer patchwork, uncomment and install if needed:
# install.packages("patchwork")
# library(patchwork)
# combined_plot <- (multipanel_plots$p_sample_type + multipanel_plots$p_batch) / 
#                  (multipanel_plots$p_sex + multipanel_plots$p_age) /
#                  (multipanel_plots$p_blood_eos + multipanel_plots$p_blood_neut) /
#                  (multipanel_plots$p_sex_dim14 + multipanel_plots$p_sample_dim23) +
#                  plot_annotation(title = "Multi-dimensional Scaling Analysis of RNA-seq Data")
# 
# ggsave("multipanel_mds_analysis_patchwork.pdf", combined_plot, width=12, height=16)
# If you want to create just the PDF without displaying on screen:
# create_multipanel_pdf <- function() {
#   pdf("mds_analysis_multipanel.pdf", width=12, height=14)
#   par(mfrow=c(4,2), mar=c(4,4,3,1), oma=c(0,0,2,0))
#   
#   # [Copy all the plot commands from create_multipanel_mds_base here]
#   
#   dev.off()
# }
# create_multipanel_pdf()

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
par(mfrow=c(2,2))
if (requireNamespace("corrplot", quietly = TRUE)) {
  # Correlation plot for continuous variables
  corrplot(cont_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           p.mat = cont_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           cl.offset = -0.5,
           cl.align.text = "r",
           mar = c(4, 1, 4, 0),
           title = "Pearson Correlations: \nContinuous Variables vs PCs")
  
  corrplot(cont_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           p.mat = cont_p_matrix, sig.level = 0.05,
           insig = "label_sig", 
           cl.offset = -0.5,
           cl.align.text = "r",
           mar = c(4, 1, 4, 0),
           title = "Pearson Correlations: \nContinuous Variables vs PCs")
  
  # Effect size plot for categorical variables
  corrplot(cat_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           col = colorRampPalette(c("white", "pink", "red"))(100),
           p.mat = cat_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           mar = c(4, 1, 4, 0),
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
           mar = c(4, 1, 4, 0),
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
var <- c(source.cell.log,"bal_eos_p_mt1", "Batch", "Age_at_visit", "Sex", "Race", "Ethnicity", "total_IgE_num", "n_b")

phen_var <- phen[, var]

# Convert categorical variables to factors
phen_var[, c("bal_eos_p_mt1","Batch", "Sex", "Race", "Ethnicity", "n_b")] <- lapply(
  phen_var[, c("bal_eos_p_mt1","Batch", "Sex", "Race", "Ethnicity", "n_b")],
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


# Function to perform MDS and create plots with outside legends
create_multipanel_mds <- function(normalized_counts, phen_var, col) {
  # Setup 4x2 panel layout with extra space for legends
  pdf("mds_analysis_with_legends.pdf", width=14, height=14)
  # Set layout with wider right margin for legends
  par(mfrow=c(4, 2), mar=c(4, 4, 3, 8), oma=c(0, 0, 2, 0))
  
  # Calculate MDS
  mds <- limma::plotMDS(normalized_counts, ndim=5, plot=FALSE)
  
  # Extract coordinates for all dimensions we'll use
  x1 <- mds$eigen.vectors[,1] * sqrt(mds$eigen.values[1])
  x2 <- mds$eigen.vectors[,2] * sqrt(mds$eigen.values[2])
  x3 <- mds$eigen.vectors[,3] * sqrt(mds$eigen.values[3])
  x4 <- mds$eigen.vectors[,4] * sqrt(mds$eigen.values[4])
  
  # Get variance explained
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
    
    for(val in unique_values) {
      # Find first index with this value
      idx <- which(phen_var[[var_name]] == val)[1]
      unique_colors <- c(unique_colors, col[[var_name]][idx])
    }
    
    # Include NA in the legend if there are any NA values
    if(any(is.na(phen_var[[var_name]]))) {
      legend_values <- c(as.character(unique_values), "NA")
      legend_colors <- c(unique_colors, na_color)
    } else {
      legend_values <- as.character(unique_values)
      legend_colors <- unique_colors
    }
    
    # Add legend outside the plot
    par(xpd=TRUE)  # Allow drawing outside the plot region
    legend(par("usr")[2] * 1.05, par("usr")[4], 
           legend=legend_values, 
           fill=legend_colors, 
           title=var_name,
           cex=0.7, bty="n")
    par(xpd=FALSE)  # Reset to default
  }
  
  # Helper function to add legends for numeric variables
  add_numeric_legend <- function(var_name) {
    # Use red gradient colors for numeric variables
    quart_colors <- c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D")
    
    # Calculate quartiles for the variable
    var_quants <- quantile(phen_var[[var_name]], probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE)
    
    # Create legend labels (rounded to 1 decimal place)
    quart_labels <- c(
      paste0("Q1: ", round(var_quants[1], 1), "-", round(var_quants[2], 1)),
      paste0("Q2: ", round(var_quants[2], 1), "-", round(var_quants[3], 1)),
      paste0("Q3: ", round(var_quants[3], 1), "-", round(var_quants[4], 1)),
      paste0("Q4: ", round(var_quants[4], 1), "-", round(var_quants[5], 1))
    )
    
    # Include NA in the legend if there are any NA values
    if(any(is.na(phen_var[[var_name]]))) {
      legend_labels <- c(quart_labels, "NA")
      legend_colors <- c(quart_colors, na_color)
    } else {
      legend_labels <- quart_labels
      legend_colors <- quart_colors
    }
    
    # Add legend outside the plot
    par(xpd=TRUE)  # Allow drawing outside the plot region
    legend(par("usr")[2] * 1.05, par("usr")[4], 
           legend=legend_labels, 
           fill=legend_colors, 
           title=var_name,
           cex=0.7, bty="n")
    par(xpd=FALSE)  # Reset to default
  }
  
  # Plot 1: Batch
  plot(x1, x2, 
       col=col$Batch, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Batch Effect")
  add_categorical_legend("Batch")
  
  # Plot 2: Sex
  plot(x1, x2, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Sex")
  add_categorical_legend("Sex")
  
  # Plot 3: Age
  plot(x1, x2, 
       col=col$Age_at_visit, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Age")
  add_numeric_legend("Age_at_visit")
  
  # Plot 4: BAL Eos > 1%
  plot(x1, x2, 
       col=col$bal_eos_p_mt1, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eos% >1%")
  add_categorical_legend("bal_eos_p_mt1")
  
  # Plot 5: BAL Eos% (log)
  plot(x1, x2, 
       col=col$BAL_eos_p_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eosinophils % (log)")
  add_numeric_legend("BAL_eos_p_log")
  
  # Plot 6: Blood Eos (log)
  plot(x1, x2, 
       col=col$blood_eos_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Blood Eosinophils (log)")
  add_numeric_legend("blood_eos_log")
  
  # Plot 7: Sex (Dims 1 & 4)
  plot(x1, x4, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS4 (", var4, "%)"),
       main="Sex (Dims 1 & 4)")
  add_categorical_legend("Sex")
  
  # Add overall title
  mtext("Multi-dimensional Scaling Analysis of RNA-seq Data", 
        outer=TRUE, cex=1.5)
  
  # Close PDF device
  dev.off()
  
  # Also display on screen
  par(mfrow=c(4, 2), mar=c(4, 4, 3, 8), oma=c(0, 0, 2, 0))
  
  # Plot 1: Batch
  plot(x1, x2, 
       col=col$Batch, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Batch Effect")
  add_categorical_legend("Batch")
  
  # Plot 2: Sex
  plot(x1, x2, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Sex")
  add_categorical_legend("Sex")
  
  # Plot 3: Age
  plot(x1, x2, 
       col=col$Age_at_visit, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Age")
  add_numeric_legend("Age_at_visit")
  
  # Plot 4: BAL Eos > 1%
  plot(x1, x2, 
       col=col$bal_eos_p_mt1, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eos% >1%")
  add_categorical_legend("bal_eos_p_mt1")
  
  # Plot 5: BAL Eos% (log)
  plot(x1, x2, 
       col=col$BAL_eos_p_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="BAL Eosinophils % (log)")
  add_numeric_legend("BAL_eos_p_log")
  
  # Plot 6: Blood Eos (log)
  plot(x1, x2, 
       col=col$blood_eos_log, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS2 (", var2, "%)"),
       main="Blood Eosinophils (log)")
  add_numeric_legend("blood_eos_log")
  
  # Plot 7: Sex (Dims 1 & 4)
  plot(x1, x4, 
       col=col$Sex, pch=19,
       xlab=paste0("MDS1 (", var1, "%)"),
       ylab=paste0("MDS4 (", var4, "%)"),
       main="Sex (Dims 1 & 4)")
  add_categorical_legend("Sex")
  
  # Add overall title
  mtext("Multi-dimensional Scaling Analysis of RNA-seq Data", 
        outer=TRUE, cex=1.5)
  
  # Return MDS coordinates in case they're needed elsewhere
  return(list(
    x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    var1 = var1, var2 = var2, var3 = var3, var4 = var4
  ))
}

# Run the function
mds_results <- create_multipanel_mds(normalized_counts, phen_var, col)
# To save the plot (uncomment to use)
# ggsave("multipanel_mds_analysis.pdf", grid.arrange(grobs=multipanel_plots), width=12, height=10)

# Alternative using patchwork package for more control:
# If you prefer patchwork, uncomment and install if needed:
# install.packages("patchwork")
# library(patchwork)
# combined_plot <- (multipanel_plots$p_sample_type + multipanel_plots$p_batch) / 
#                  (multipanel_plots$p_sex + multipanel_plots$p_age) /
#                  (multipanel_plots$p_blood_eos + multipanel_plots$p_blood_neut) /
#                  (multipanel_plots$p_sex_dim14 + multipanel_plots$p_sample_dim23) +
#                  plot_annotation(title = "Multi-dimensional Scaling Analysis of RNA-seq Data")
# 
# ggsave("multipanel_mds_analysis_patchwork.pdf", combined_plot, width=12, height=16)
# If you want to create just the PDF without displaying on screen:
# create_multipanel_pdf <- function() {
#   pdf("mds_analysis_multipanel.pdf", width=12, height=14)
#   par(mfrow=c(4,2), mar=c(4,4,3,1), oma=c(0,0,2,0))
#   
#   # [Copy all the plot commands from create_multipanel_mds_base here]
#   
#   dev.off()
# }
# create_multipanel_pdf()

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
par(mfrow=c(2,2))
if (requireNamespace("corrplot", quietly = TRUE)) {
  # Correlation plot for continuous variables
  corrplot(cont_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           p.mat = cont_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           cl.offset = -0.5,
           cl.align.text = "r",
           mar = c(4, 1, 4, 0),
           title = "Pearson Correlations: \nContinuous Variables vs PCs")
  
  corrplot(cont_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           p.mat = cont_p_matrix, sig.level = 0.05,
           insig = "label_sig", 
           cl.offset = -0.5,
           cl.align.text = "r",
           mar = c(4, 1, 4, 0),
           title = "Pearson Correlations: \nContinuous Variables vs PCs")
  
  # Effect size plot for categorical variables
  corrplot(cat_matrix, method = "color", 
           type = "full", 
           tl.col = "black", tl.srt = 45,
           col = colorRampPalette(c("white", "pink", "red"))(100),
           p.mat = cat_p_matrix, sig.level = 0.05,
           insig = "p-value", 
           mar = c(4, 1, 4, 0),
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
           mar = c(4, 1, 4, 0),
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


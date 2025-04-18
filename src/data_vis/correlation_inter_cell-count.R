##
# Investigate relationships among the cell counts in BAL and blood
## 

library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)
library(corrplot)
library(Hmisc)

# find how many people had at least one ED visit or admission in the cohort
library(dplyr)
countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
rownames(counts)<-counts[,"SampleID"]
counts.ID<-colnames(counts)

# asthma biomarker phenotype file saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2025-02-14.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}

sampling_date_diff<-"./resources/processed_data/sampling_dates/swab-bal-cbc_differences_in_days.txt"
sampling_date_diff<-if(file.exists(sampling_date_diff)){read.table(sampling_date_diff,row.names = NULL,header = TRUE)}
sampling_date_diff<-sampling_date_diff%>%filter(Comparison=="blood_bal")
colnames(sampling_date_diff)[1:3]<-c("ID","sampling_date_comp","sampling_date_diff_days")

phenotype<-left_join(phenotype,sampling_date_diff,by="ID")
phenotype_filtered<-phenotype%>%filter(abs(sampling_date_diff_days)<30)
phenotype<-phenotype_filtered
################################
## load phenotype and batch data
################################

# make vectors of variables for later use as an input for function 'run_deseq2_DEG_analysis'

source.cell.log<-c(
  "BAL_eos_ct_log",
  "BAL_eos_p_log",
  "BAL_neut_ct_log",
  "BAL_neut_p_log",
  "BAL_wbc_log",
  "blood_eos_log",
  "blood_eos_p_log",
  "blood_neut_log",
  "blood_neut_p_log",
  "blood_wbc_log")
source.cell<-c(
  "BAL_eos_ct",
  "BAL_eos_p",
  "BAL_neut_ct",
  "BAL_neut_p",
  "BAL_wbc",
  "blood_eos",
  "blood_eos_p",
  "blood_neut",
  "blood_neut_p",
  "blood_wbc")

###########################################################################################
## subset phenotype data for which the samples exist for bronchial RNAseq experiments   
###########################################################################################
bexist<-phenotype$SampleID%in%counts.ID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
bphen<-phenotype[bexist,]%>%filter(grepl("^B",SampleID))
bphen<-mutate_at(bphen,vars(all_of(source.cell.log)),scale)

# calculating 'r' and p-val of cell count correlation 
data<-bphen[,source.cell.log]
rownames(data)<-bphen$SampleID
correlation_results <- Hmisc::rcorr(as.matrix(data))
cor_matrix <- correlation_results$r  # Correlation matrix
p_matrix <- correlation_results$P   # Matrix of p-values

# Function to create significance labels
get_sig_labels <- function(p_matrix, sig.level = c(0.001, 0.01, 0.05)) {
  # Assign significance levels as asterisks
  sig_labels <- matrix("", nrow = nrow(p_matrix), ncol = ncol(p_matrix))
  sig_labels[p_matrix <= sig.level[1]] <- "***"
  sig_labels[p_matrix <= sig.level[2] & p_matrix > sig.level[1]] <- "**"
  sig_labels[p_matrix <= sig.level[3] & p_matrix > sig.level[2]] <- "*"
  sig_labels
}

# Generate significance labels
sig_labels <- get_sig_labels(p_matrix)

# plot correlation coef and pval 
cell_count_cor_heatmap<-corrplot::corrplot(cor_matrix, type = "upper", order = "hclust", 
                   tl.col = "black", tl.srt = 45, 
                   method = "ellipse", addCoef.col = "black",
                   diag= FALSE,
                   title="Correlation Matrix, blood-bal sampling date difference < 30d, no significance",
                   mar=c(5,5,5,5))


# delete non-significant 
corrplot::corrplot(cor_matrix, type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45,
         method = "circle", addCoef.col = "black", # Add correlation coefficients
         p.mat = p_matrix, sig.level = 0.05, # Add significance levels
         insig = "n", # What to do with insignificant correlations
         addCoefasPercent = FALSE, # Show the coefficients
         cl.cex = 0.8,
         diag = FALSE,
         title="Correlation Matrix with Significance Levels",
         mar=c(5,5,5,5)
         ) # Adjust text size



# Determine the number of rows and columns dynamically
library(gridExtra)  # Ensure this is loaded
library(ggpubr)

plot_list <- list()
num_columns <- ncol(data)
for (i in 1:(num_columns - 1)) {
  for (j in (i + 1):num_columns) {
    # Fetch the correlation coefficient and p-value
    cor_value <- cor_matrix[i, j]
    p_value <- p_matrix[i, j]
    
    # Create the plot
    p <- ggplot(data, aes_string(x = names(data)[i], y = names(data)[j])) +
      geom_point(alpha = 0.4) +  # Add points
      geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add regression line
      theme_minimal() +  # Use a minimal theme
      labs(title = paste("Scatter plot of", names(data)[i], "vs", names(data)[j]),
           x = names(data)[i], y = names(data)[j]) +
      # Annotate the plot with correlation and p-value
      annotate("text", x = Inf, y = Inf, label = sprintf("r = %.2f, p = %.3f", cor_value, p_value),
               hjust = 1.1, vjust = 1.1, size = 3, color = "red")
    
    # Add the plot to the list
    plot_list[[paste(names(data)[i], names(data)[j], sep = "_vs_")]] <- p
  }
}

# Optional: Display one of the plots to check
print(plot_list[[1]])
# Combine all plots using patchwork
library(patchwork)
plot_layout <- wrap_plots(plot_list, ncol = 3)  # Adjust ncol to fit your requirement
plot_layout

wrap_plots(plot_list[1:15], ncol = 3)
wrap_plots(plot_list[16:30], ncol = 3)
wrap_plots(plot_list[31:45], ncol = 3)
cell_count_cor_heatmap
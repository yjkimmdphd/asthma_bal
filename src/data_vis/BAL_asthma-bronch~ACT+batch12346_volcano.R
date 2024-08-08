library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(gridExtra)
library(grid)
data.folder<-file.path("./reports/local_only/deg_bronch~ACTscore+batch_deg_2024-08-03")
####################################################################
# data exploration of DEG analysis using nasal/bronchial rnaseq data
# model: "Nasal DEG ~ blood AEC (>0) + Batch
# model: "Bronch DEG ~ BAL ANC(>0) + Batch
####################################################################

read.csv(file.path(data.folder,"deg_bronch_ACT-score+batch12346_res_all_1_~ asthma_phen_ACT.score + Batch_2024-08-03_.csv"),row.names = 1)

# Read the CSV file without setting row names
data <- read.csv(file.path(data.folder, "deg_bronch_ACT-score+batch12346_res_all_1_~ asthma_phen_ACT.score + Batch_2024-08-03_.csv"))

# Check for duplicates in the first column
if(anyDuplicated(data[, 1]) > 0) {
  # Make the first column unique
  data[, 1] <- make.names(data[, 1], unique = TRUE)
}

# Set the first column as row names
rownames(data) <- data[, 1]
data <- data[, -1] # Remove the first column as it's now the row names

par(mfrow=c(1,1))

# padj are BH adjusted pvalues 
p1<-EnhancedVolcano(data,
                    lab = rownames(data),
                    title='bronch ~ ACT score',
                    subtitle = "all FDR > 0.05, showing nominal p-values",
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    xlab = bquote(~Log[2]~ 'fold change'),
                    ylab = "nominal -Log(p values)",
                    xlim=c(-2,2),
                    ylim=c(0,8),
                    FCcutoff = 1,
                    cutoffLineType = 'twodash',
                    cutoffLineWidth = 0.8,
                    pointSize = 4.0,
                    labSize = 3,
                    colAlpha = 0.4,
                    legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                   'p-adj (<0.05) & Log (base 2) FC'),
                    legendPosition = 'right',
                    legendLabSize = 10,
                    legendIconSize = 5.0,    
                    drawConnectors = TRUE,
                    widthConnectors = 0.75)
p1

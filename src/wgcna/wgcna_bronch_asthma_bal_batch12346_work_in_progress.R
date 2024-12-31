##
# WGCNA of bronchial rna expression 
# reference tutorial: https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html
## 
library(limma)
library(edgeR)
library(tidyverse)
library(DESeq2)
library(WGCNA)

########################
## Load Readcount Data ##
########################

# Load cell count table
normalized_count_table_path <- "./resources/processed_data/normalized_gene_count/normalized_gene_count_bronch_vsd_batch-corrected.txt"
if (file.exists(normalized_count_table_path)) {
  ct <- read.table(normalized_count_table_path, 
                       header = TRUE, 
                       row.names = 1, 
                       sep = "\t")
}else(print("path doesn't exist"))
genes <- rownames(ct)

###############
##   WGCNA   ##
###############

options(allowParallel = TRUE)
enableWGCNAThreads()

# 1. Create a new format for expression data
expression <- as.data.frame(t(ct))
colnames(expression) <- genes

# Check if genes/samples are good
gsg <- goodSamplesGenes(expression)

# 2. Remove outlier genes (and samples if necessary)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(
      paste("Removing genes:", 
            paste(names(expression)[!gsg$goodGenes], collapse = ", "))
    )
  }
  
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(
      paste("Removing samples:", 
            paste(rownames(expression)[!gsg$goodSamples], collapse = ", "))
    )
  }
  
  expression <- expression[gsg$goodSamples, gsg$goodGenes]
}

# 3. Identify outlier samples using a dendrogram
sampleTree <- hclust(dist(expression), method = "average")

png("./reports/local_only/wgcna/bronch/sampleTree.png", width = 800, height = 600)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)
abline(h = 160, col = "red")
dev.off()

# Remove outliers
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 160, minSize = 7)
expression.data <- expression[cut.sampleTree == 1, ]

dim(expression)
dim(expression.data)

# 4. Network formulation

# Pick soft threshold
spt <- pickSoftThreshold(expression.data)

# Plot Scale Independence
png("./reports/local_only/wgcna/bronch/bronch_spt.png", width = 800, height = 600)
par(mar = c(1, 1, 1, 1))
plot(
  spt$fitIndices[, 1],
  spt$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale independence"
)
text(spt$fitIndices[, 1], spt$fitIndices[, 2], col = "red")
abline(h = 0.80, col = "red")
dev.off()

# Plot Mean Connectivity
png("./reports/local_only/wgcna/bronch/bronch_spt_vs_connectivity.png", width = 800, height = 600)
par(mar = c(1, 1, 1, 1))
plot(
  spt$fitIndices[, 1],
  spt$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)
text(spt$fitIndices[, 1], spt$fitIndices[, 5], labels = spt$fitIndices[, 1], col = "red")
dev.off()

# Adjacency matrix
softPower <- 4
adjacency <- adjacency(expression.data, power = softPower)

# 5. Module Construction

print("calculating TOM")
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1 - TOM

if (!file.exists("./resources/processed_data/Rdata/wgcna_bronch_TOM_diss_batch12346.Rdata")) {
  save(TOM.dissimilarity, 
       file = "./resources/processed_data/Rdata/wgcna_bronch_TOM_diss_batch12346.Rdata")
}

if (file.exists("./resources/processed_data/Rdata/wgcna_bronch_TOM_dissim_batch12346.Rdata")) {
  load("./resources/processed_data/Rdata/wgcna_bronch_TOM_dissim_batch12346.Rdata")
}

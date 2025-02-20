##
# WGCNA of bronchial rna expression 
# reference tutorial: https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html
## 

library(tidyverse)
library(DESeq2)
library(WGCNA)

########################
## Load Readcount Data ##
########################

# Load cell count table
normalized_count_table_path <- "./resources/processed_data/normalized_gene_count/normalized_gene_count_bronch_vsd_batch-corrected.txt"
if (file.exists(normalized_count_table_path)) {
  counts <- read.table(normalized_count_table_path, 
                       header = TRUE, 
                       row.names = 1, 
                       sep = "\t")
}
genes <- rownames(counts)

###############
##   WGCNA   ##
###############

options(allowParallel = TRUE)
enableWGCNAThreads()

# Define the path for the output folder
output_folder <- file.path("./reports/local_only/wgcna/bronch/output")

# Check if the folder exists; if not, create it
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE) # recursive = TRUE ensures all intermediate directories are created
}

# Now the folder is guaranteed to exist, and the path is set
output_folder


###
# 1. Create a new format for expression data
###
expression <- as.data.frame(t(counts))
colnames(expression) <- genes

# Check if genes/samples are good
gsg <- goodSamplesGenes(expression)

###
# 2. Remove outlier genes (and samples if necessary)
###

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
###
# 3. Identify outlier samples using a dendrogram
###
sampleTree <- hclust(dist(expression), method = "average")

png(file.path(output_folder,"sampleTree_bronch.png"), width = 800, height = 600)
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

###
# 4. Network formulation
###
# Pick soft threshold
spt <- pickSoftThreshold(expression.data)

# Plot Scale Independence
png(file.path(output_folder,"bronch_spt.png"), width = 800, height = 600)
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
png(file.path(output_folder,"bronch_spt_vs_connectivity.png"), width = 800, height = 600)
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

###
# 5. Module Construction
###

print("calculating TOM")
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1 - TOM

# ----------------------
# saving TOM.dissimilarity matrix
# ----------------------
# 
# if (!file.exists("./resources/processed_data/Rdata/wgcna_bronch_TOM_diss_batch12346.Rdata")) {
#   save(TOM.dissimilarity, 
#        file = "./resources/processed_data/Rdata/wgcna_bronch_TOM_diss_batch12346.Rdata")
# }

# --------------------------
# load TOM.dissimilarity 
# should execute block 1-3 to proceed after loading TOM.dissimilarity 
# --------------------------
# if (file.exists("./resources/processed_data/Rdata/wgcna_bronch_TOM_diss_batch12346.Rdata")) {
#   load("./resources/processed_data/Rdata/wgcna_bronch_TOM_diss_batch12346.Rdata")
# }else{print("Rdata doesn't exist")}
# --------------------------


# Hierarchical Clustering Analysis
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

# Dynamic Tree Cutting
Modules <- cutreeDynamic(
  dendro = geneTree,
  distM = TOM.dissimilarity,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 30
)

table(Modules)

ModuleColors <- labels2colors(Modules)
table(ModuleColors)

png(file.path(output_folder,"Gene_dendrogram_and_module_colors.png"), width = 800, height = 600)

plotDendroAndColors(
  geneTree,
  ModuleColors,
  "Module",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

###
# 6. module eigengene identification
# A ME (Module Eigengene) is the standardized gene expression profile for a given module.
###


MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)


write.table(MElist,file.path(output_folder,"pre-merge_MElist.txt"), sep="\t",quote=FALSE, row.names=TRUE, col.names=NA)
write.table(MEs,file.path(output_folder,"pre-merge_MEs.txt"), sep="\t",quote=FALSE, row.names=TRUE, col.names=NA)

###
# 7. Module merging
### 

# goal: condense the branches

# remove the NA values
ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic

png(file.path(output_folder,"METree.png"), width = 800, height = 600)

plot(METree)
abline(h=0.25, col = "red")  #a height of h corresponds to correlation of 1.00 - h (i.e.,  all of the modules which are more than 85% similar if h=0.15.)
dev.off()

merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = 0.25) # merge the modules which are below the threshold

# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

# plot dendrogram showing both the orginal and merged module colors

png(file.path(output_folder,"Gene_dendrogram_and_module_colors_for_original_and_merged_modules.png"), width = 800, height = 600)
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

write.table(mergedColors,file.path(output_folder,"mergedColors.txt"), sep="\t",quote=FALSE, row.names=TRUE, col.names=NA)
write.table(mergedMEs,file.path(output_folder,"mergedMEs.txt"), sep="\t",quote=FALSE, row.names=TRUE, col.names=NA)


###
# 8. construct cytoscape info
## 

# Set a threshold for adjacency to filter significant edges
threshold <- 0.1

# Get the gene names
genes <- rownames(adjacency)

# Generate the edge list
edge_list <- as.data.frame(which(adjacency > threshold, arr.ind = TRUE))
edge_list <- edge_list[edge_list$row != edge_list$col, ] # Remove self-loops
edge_list <- edge_list[!duplicated(t(apply(edge_list, 1, sort))), ] # Remove duplicates
edge_list$weight <- adjacency[cbind(edge_list$row, edge_list$col)]
edge_list <- edge_list %>% 
  mutate(Source = genes[row], Target = genes[col]) %>%
  select(Source, Target, weight)

# Save edge list to a file
write.table(edge_list, file.path(output_folder,"wgcna_bronch_edge_list.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Create the node attribute file
node_attributes <- data.frame(
  Gene = rownames(connectivity_allClusters),
  Module = connectivity_allClusters$module,
  Connectivity = connectivity_allClusters$kTotal
)

# Save node attributes to a file
write.table(node_attributes, file.path(output_folder,"wgcna_bronch_node_attributes.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# ===================
# DOWNSTREAM ANALYSES 
# ===================
# list of downstream analyses that should be performed on a different script: 
## choosing hub genes 
## Evaluate DEG and module gene overlap by fisher's exact test 
## Network Visualization of Eigengenes by checking correlation between the Eigen genes of different modules 
## Module Membership vs gene-trait significance 
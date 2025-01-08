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


# -----------------------------------------------------------------------------
# Load phenotype data
# -----------------------------------------------------------------------------
phen_path <- file.path(
  "./resources/processed_data",
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv"
)
phen <- read.csv(phen_path)

phen<-phen %>% filter(BAL_eos_p >= 0 & BAL_neut_p >= 0) %>%
  mutate(comp1_pnme = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 4 ~ 2, # "mixed"
                                  BAL_eos_p > 1 & BAL_neut_p <= 4 ~ 3, # "eos"
                                  BAL_eos_p <= 1 & BAL_neut_p > 4 ~ 1, # "neut"
                                  BAL_eos_p <= 1 & BAL_neut_p <= 4 ~ 0), # "pauci"
                        levels = c(0, 1, 2, 3)),
         comp2_eos_p_mt1 = factor(case_when(BAL_eos_p > 1 ~ 1 , #"high_eos"
                                  BAL_eos_p <= 1 ~ 0 ),# "low_eos"
                        levels = c(0,1)))
phen$Batch<-factor(phen$Batch,levels=unique(phen$Batch))

phen_bronch<-phen[grepl("^B",phen$SampleID),]
phen_nasal<-phen[grepl("^N",phen$SampleID),]
phen_nasal<-phen_nasal[-grep("F", phen_nasal$SampleID),]
phen_input<-phen_bronch
phen_input$SampleID <- gsub("-", ".", phen_input$SampleID)

alltraits<-phen_input[phen_input$SampleID%in%rownames(mergedMEs),]
rownames(alltraits)<-alltraits$SampleID

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

alltraits<-alltraits[,c("asthma_phen_ACT.score",source.cell.log, "comp1_pnme","comp2_eos_p_mt1")] # only select cell counts in the phenotype data
good<-!(sapply(alltraits,is.na)%>%rowSums()>0) # remove rows with at least one NA
datTraits<-alltraits[good,]


###
# 8. Module-Trait associations
###

gene.module.table<-data.frame(genes=colnames(expression.data),modules=mergedColors)

# Assuming list_of_dfs is your list of dataframes
list_of_dfs <- split(gene.module.table, gene.module.table$module)
genelist<-vector("list",length=length(list_of_dfs))
names(genelist)<-names(list_of_dfs)

for(i in 1:length(list_of_dfs)){
  genelist[[i]]<-list_of_dfs[[i]]$genes
}

# Finding the maximum length of vectors in the list
max_length <- max(sapply(genelist, length))

# Padding shorter vectors with NA
genelist_padded <- lapply(genelist, function(x) {
  length(x) <- max_length  # Setting the length of each vector to the maximum
  x
})

# Convert the list to a data frame
df_genelist <- as.data.frame(genelist_padded)


# Writing the data frame to a tab-delimited text file without column names
# Check if the folder exists; if not, create it
module_output_folder <- file.path("./reports/local_only/wgcna/bronch/output/module-gene_list")

if (!dir.exists(module_output_folder)) {
  dir.create(module_output_folder, recursive = TRUE)
}

write.table(df_genelist, file = file.path(output_folder,"bronch_wgcna_genelist_batch12346.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)

# Loop through each element in the list and write to a CSV file
for (module_name in names(list_of_dfs)) {
  # Extract the 'genes' column from each dataframe
  gene_data <- list_of_dfs[[module_name]]$genes
  
  # Create a file path
  file_path <- file.path(module_output_folder,paste0(module_name,"batch12346", ".txt"))
  
  # Write to a text file without row names or column names
  write.table(gene_data, file=file_path, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# ----------------
# proportion of module gene that overlaps with bronchial DEG for BAL Eos % > 1% vs <=1%
# ----------------
deg_folder<-file.path("./reports/local_only/deg_bal_bronch~cell2025-01-03")
deg_file<-if(file.exists(deg_folder)){list.files(deg_folder)[grep(".csv",list.files(deg_folder))]}else{print("folder doesn't exist")}
deg<-read.csv(file.path(deg_folder,"deg_bronch_res_sig_16_~ bal_Eos_p_more_1 + Batch_2025-01-03_.csv"), row.names = 1)
all_assessed_genes<-read.csv(file.path(deg_folder,"deg_bronch_res_all_16_~ bal_Eos_p_more_1 + Batch_2025-01-03_.csv"), row.names = 1)%>%rownames
deg_abs_lfc<-deg%>%filter(abs(log2FoldChange)>1)%>%rownames

wgcna_folder<-module_output_folder

module_list<-list.files(wgcna_folder)[grep(".txt",list.files(wgcna_folder))]

modules_files<-file.path(wgcna_folder,module_list)

m_gene_list<-lapply(modules_files,read.table)

names(m_gene_list)<-sub("batch12346.txt","",module_list)

overlap_proportion <- sapply(m_gene_list, function(x) {
  round(mean(unlist(x) %in% deg_abs_lfc),2)
})# no notable overlap between bronch down gene with WGCNA module genes 

overlap_sum <- sapply(m_gene_list, function(x) {
  sum(unlist(x) %in% deg_abs_lfc)
})

overlap_df<-data.frame(overlap_proportion,overlap_sum)
overlap_df<-overlap_df[order(-(overlap_df$overlap_sum)),]
print(overlap_df)

write.table(overlap_df,file.path(output_folder,"module-gene_deg-gene_overlap.txt"), sep="\t",quote=FALSE, row.names=TRUE, col.names=NA)

# quantify the association between the expression profile and a particular trait of interest 
# by calculating the correlation of the trait with previously identified module eigengenes. 
# This pairwise correlation is known as the eigengene gene significance
expression.data<-expression.data[rownames(expression.data)%in%rownames(datTraits),]
mergedMEs<-mergedMEs[rownames(mergedMEs)%in%rownames(datTraits),]
nSamples = nrow(expression.data)
module.trait.correlation = cor(mergedMEs, datTraits, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation

# Will display correlations and their p-values
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
                   signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
par(mar = c(6, 6, 3, 1))
# Display the correlation values within a heatmap plot
png(file.path(output_folder,"module-trait-correlation.png"), width = 1200, height = 1200)
labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait correlation"))
dev.off()



### 
# 9. Target gene identification
### 

# Whole Network connectivity - a measure for how well the node is connected throughout the entire system
# Intramodular connectivity - a measure for how well the node is connected within its assigned module. Also an indicator for how well that node belongs to its module. This is also known as module membership.

modNames = substring(names(mergedMEs), 3) #extract module names

#Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(WGCNA::cor(expression.data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# ----------------
# for log(BAL Eos %)
# ----------------

# Define variable weight containing the weight column of datTrait
BAL_eos_p_log = as.data.frame(datTraits$BAL_eos_p_log)
names(BAL_eos_p_log) = "BAL_eos_p_log"

phen_of_interest<-BAL_eos_p_log

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(WGCNA::cor(expression.data, phen_of_interest, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(phen_of_interest), sep="")
names(GSPvalue) = paste("p.GS.", names(phen_of_interest), sep="")
head(GSPvalue)

# scatter plot of gene significance vs. module membership in all the module

par(mfrow=c(6,4))

png(file.path(output_folder,"Module_membership_vs._gene_significance_log(BAL_Eos_%).png"), width = 1200, height = 1200)
for(mod in modNames){
  module = mod
  column = match(module, modNames)
  moduleGenes = mergedColors==module
  verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                     abs(geneTraitSignificance[moduleGenes,1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for eos",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}

dev.off()

# ----------------
# for BAL Eos % > 1 vs <=1 + BAL Neut % >4 vs <=4
# ----------------

# Define variable weight containing the weight column of datTrait
phen_of_interest<-as.data.frame(datTraits$comp1_pnme)

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(WGCNA::cor(expression.data, phen_of_interest, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(phen_of_interest), sep="")
names(GSPvalue) = paste("p.GS.", names(phen_of_interest), sep="")
head(GSPvalue)

# scatter plot of gene significance vs. module membership in all the module
par(mar = c(4,4,4,4),mfrow=c(6,4))

png(file.path(output_folder,"Module_membership_vs._gene_significance_mixed_cell.png"), width = 1200, height = 1200)
for(mod in modNames){
  module = mod
  column = match(module, modNames)
  moduleGenes = mergedColors==module
  verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                     abs(geneTraitSignificance[moduleGenes,1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for eos mt1",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
dev.off()
# ----------------
# for BAL Eos % > 1 vs <=1
# ----------------

# Define variable weight containing the weight column of datTrait
phen_of_interest<-as.data.frame(datTraits$comp2_eos_p_mt1)

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(WGCNA::cor(expression.data, phen_of_interest, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(phen_of_interest), sep="")
names(GSPvalue) = paste("p.GS.", names(phen_of_interest), sep="")
head(GSPvalue)

# scatter plot of gene significance vs. module membership in all the module
# DEG are colored in red
par(mar = c(4,4,4,4),mfrow=c(6,4))

png(file.path(output_folder,"Module_membership_vs._gene_significance_eos-p_mt1.png"), width = 1200, height = 1200)
for(mod in modNames){
  module = mod
  column = match(module, modNames)
  moduleGenes = mergedColors==module
  
  module_color<-rep(mod,length(geneModuleMembership[moduleGenes,column]))
  module_color[which(rownames(geneModuleMembership)[moduleGenes]%in%deg_abs_lfc)]<-"red"

  verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                     abs(geneTraitSignificance[moduleGenes,1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for eos mt1",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module_color)
}
dev.off()
###
# 10. Network Visualization of Eigengenes
### 

# Isolate blood_eos_log from the clinical traits
bal_mix=as.data.frame(datTraits$comp1_pnme);
bal_eos_mt1=as.data.frame(datTraits$comp2_eos_p_mt1);
BAL_eos_p_log = as.data.frame(datTraits$BAL_eos_p_log);
BAL_neut_p_log = as.data.frame(datTraits$BAL_neut_p_log)

names(bal_mix) = "mixed_cell"
names(bal_eos_mt1) = "eos_mt1"
names(BAL_eos_p_log) = "BAL_eos_p_log"
names(BAL_neut_p_log) = "BAL_neut_p_log"

# Add the BAL_eos_ct_log to existing module eigengenes
mergedMEs<-mergedMEs[rownames(mergedMEs)%in%rownames(datTraits),]
MET = orderMEs(cbind(mergedMEs, bal_mix, bal_eos_mt1,BAL_eos_p_log, BAL_neut_p_log ))
# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# Plot the dendrogram
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0, mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(10,10,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

###
# 11. find genes in each module that has high correlation and low p-val in gene significance vs. module membership
### 

#####
# 12. choose hub genes
#####

HubGenes <- chooseTopHubInEachModule(expression.data, mergedColors, power = 4, type= "signed")

connectivity_allClusters <- intramodularConnectivity(adjacency, mergedColors, scaleByMax = FALSE)

colnames(expression.data)%>%head

# Check if rownames(connectivity_allClusters) match gene.module.table$genes
matched <- match(rownames(connectivity_allClusters), gene.module.table$genes)

# Update connectivity_allClusters$module with the corresponding module values
connectivity_allClusters$module <- gene.module.table$modules[matched]
connectivity_allClusters$gene <- gene.module.table$genes[matched]

top_kWithin_by_module <- connectivity_allClusters %>%
  group_by(module) %>%
  slice_max(order_by = kWithin, n = 20) %>%
  ungroup() %>%
  arrange(module, desc(kWithin))

head(top_kWithin_by_module)

df_assessed_only <- top_kWithin_by_module %>%
  filter(gene %in% deg_abs_lfc)

head(df_assessed_only)


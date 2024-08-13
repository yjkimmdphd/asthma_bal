##
# WGCNA of nasal rna expression 
# reference tutorial: https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html
## 
library(limma)
library(edgeR)
library(tidyverse)
library(DESeq2)
library(WGCNA)
######################
## load readcount data
######################

countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
bronch.samples<-grepl("^B",colnames(counts))
bronch.samples<-which(bronch.samples==TRUE)
bronch.counts<-counts[,c(1,bronch.samples)]
colnames(bronch.counts)[bronch.samples]<-substr(colnames(bronch.counts)[bronch.samples],1,4)
head(bronch.counts)
counts.ID<-colnames(bronch.counts)
counts<-bronch.counts
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


# asthma biomarker phenotype file saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
numeric_id<-phenotype$ID
phenotype$ID<-sprintf("%03d",numeric_id) # adds padded zeros in front of the subject ID numbers



###########################################################################################
## subset phenotype data for which the samples exist for nasal/bronchial RNAseq experiments   
###########################################################################################
bID<-phenotype$SampleID
bexist<-bID%in%counts.ID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
bsample<-bID[bexist] # bronchial sample ID in the readcount matrix (batch 1-4,6) that has BAL phenotype data
bphen<-phenotype[phenotype$ID%in%substring(bsample,2),] # phenotype table with bsample
bphen<-mutate(bphen, SampleID=bsample)%>%relocate(SampleID, .before=1) # include sample ID for bronchial RNAseq samples


###################################
# custom functions for DEG analysis
###################################
# should load the following fx:
## filter_low_expressed_genes_method2: Filters low counts genes using TMM normalized lcpm as a cutoff point. Requires 'limma'
## rowgenes_counttable: changes the row names of the count table with gene names
## run_deseq2_DEG_analysis: takes countdata,coldata,design,des as variables, and runs DEG on DESeq2
## get_DEG_results: saves result of DESeq2 output, ordered in padj 
## generate_DEG_input_summary_table: makes a table of input information
## generate_DEG_summary_table: makes results summary (i.e., # of DEG for each analysis)

source("./src/function/deg_custom_functions.R")


##########################################################
#set colData (phenotype data) for bronchial RNAseq experiments
##########################################################

phen<-bphen

# select RNAseq counts
id<-phen$SampleID
cols<-colnames(counts)%in%id
ct<-counts[,cols] # First column is actually gene name 
genes<-counts$SampleID
rownames(ct)<-geneso

###############
#
#### WGCNA ####
#
###############

# Set the number of cores to a lower number than available
options(allowParallel = TRUE)
# If your machine has 8 cores, you might try using just 4 or even 2
enableWGCNAThreads(nThreads = 2)

##
# 1. Create a new format expression data - remove gene name column
##
expression<-ct
rownames(expression)<-NULL
expression = as.data.frame(t(expression))

#set col as gene names
colnames(expression) = genes

# gene filter
gsg <- goodSamplesGenes(expression) # gsg$allOK is false, needs some filtering

##
# 2. remove outlier genes
##
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(expression)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(expression)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
  expression <- expression[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}
##
# 3. Identifying outlier samples with dendrogram
##
sampleTree <- hclust(dist(expression), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#draw on line to show cutoff height
abline(h = 3e6, col = "red")

# height of 3e6 removes sample ID N277
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 3e6, minSize = 10) #returns numeric vector
#Remove outlier
expression.data <- expression[cut.sampleTree==1, ]
dim(expression)
dim(expression.data)

##
# 4. network formulation
##

# Determining the Soft Power Threshold
spt <- pickSoftThreshold(expression.data) 

save(spt,file="./resources/processed_data/Rdata/wgcna_bronch_spt_batch12346.Rdata") # save the spt as Rdata

# Plot the R2 values as a function of the soft thresholds
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")

#Plot mean connectivity as a function of soft thresholds
par(mar=c(1,1,1,1))
dev.off()
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")


# spt should be 4, since after it, the power is above 0.8, for a short period lol

#calling the adjacency fx
softPower <- 4
adjacency <- adjacency(expression.data, power = softPower)

###
# 5. module Construction
### 
# convert the adjacency matrix into a TOM similarity matrix we can call the WGCNA function 
# this takes the longest and needs HPC 
print("calculating TOM")
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1-TOM

save(TOM,file="./resources/processed_data/Rdata/wgcna_bronch_TOM_batch12346.Rdata") # save the TOM as Rdata

# Hierarchical Clustering Analysis

#creating the dendrogram 
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 

#plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)

table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module. 

# plot the module assignment under the gene dendrogram for visualization.
ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)

#plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


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
plot(METree)
abline(h=.25, col = "red") #a height of 0.25 corresponds to correlation of 0.75 (i.e.,  all of the modules which are more than 75% similar.)

merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25) # merge the modules which are below the threshold

# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

# plot dendrogram showing both the orgiinal and merged module colors
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
alltraits<-phen[phen$SampleID%in%rownames(mergedMEs),]
rownames(alltraits)<-alltraits$SampleID
alltraits<-alltraits[,3:12]
good<-!(sapply(alltraits,is.na)%>%rowSums()>0) # remove rows with at least one NA
datTraits<-alltraits[good,]


###
# 8. Module-Trait associations
###

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
par(mar = c(6, 8.5, 3, 1))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


### 
# 9. Target gene identification
### 

# Define variable weight containing the weight column of datTrait
blood_eos_log = as.data.frame(datTraits$blood_eos_log)
names(blood_eos_log) = "blood_eos_log"

modNames = substring(names(mergedMEs), 3) #extract module names

#Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(WGCNA::cor(expression.data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(WGCNA::cor(expression.data, blood_eos_log, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(blood_eos_log), sep="")
names(GSPvalue) = paste("p.GS.", names(blood_eos_log), sep="")
head(GSPvalue)

# scatter plot of gene significance vs. module membership in all the module

par(mfrow=c(4,4))
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


# workspace saving
save.image("wgcna_workspace_dendrogram.Rdata")

# load Rdata
# load("./resources/processed_data/wgcna_nasal/wgcna_workspace_dendrogram.Rdata")

c("lightcyan","greenyellow","cyan","darkorange","salmon","skyblue3","salmon4","pink","darkgreen","paleturquoise","plum2","violet","blue","grey")

###
# 10. Network Visualization of Eigengenes
### 

# Isolate blood_eos_log from the clinical traits
blood_eos_log = as.data.frame(datTraits$blood_eos_log);
names(blood_eos_log) = "blood_eos_log"
# Add the blood_eos_log to existing module eigengenes
MEs<-MEs[rownames(MEs)%in%rownames(datTraits),]
MET = orderMEs(cbind(MEs, blood_eos_log))
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
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

###
# 11. find genes in each module that has high correlation and low p-val in gene significance vs. module membership
### 

gene.module.table<-data.frame(genes=colnames(expression.data),modules=mergedColors)

mod=c("lightcyan","greenyellow","cyan","darkorange","salmon","skyblue3","salmon4","pink","darkgreen","paleturquoise","plum2","violet","blue","grey")

gene.module.table<-gene.module.table%>%filter(modules%in%mod)

# Assuming list_of_dfs is your list of dataframes
list_of_dfs <- split(gene.module.table, gene.module.table$module)
genelist<-vector("list",length=length(mod))
names(genelist)<-names(list_of_dfs)

for(i in 1:length(mod)){
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
write.table(df_genelist, file = "./reports/wgcna/nasal/genelist_new_spt.txt", sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)

# Loop through each element in the list and write to a CSV file
for (module_name in names(list_of_dfs)) {
  # Extract the 'genes' column from each dataframe
  gene_data <- list_of_dfs[[module_name]]$genes
  
  # Create a file path
  file_path <- paste0("./reports/wgcna/nasal/", module_name, ".txt")
  
  # Write to a text file without row names or column names
  write.table(gene_data, file=file_path, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

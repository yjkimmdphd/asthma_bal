###
# source TOM.dissimilarity from bronch_wgcna_ANC_pos_TOM
###
library(limma)
library(edgeR)
library(tidyverse)
library(DESeq2)
library(WGCNA)

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
abline(h=.10, col = "red") #a height of h corresponds to correlation of 1.00 - h (i.e.,  all of the modules which are more than 85% similar if h=0.15.)

######################################### check here 5/18/2024 progress ####################
merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .10) # merge the modules which are below the threshold

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
alltraits<-alltraits[,3:12] # only select cell counts in the phenotype data
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


###
# 10. Network Visualization of Eigengenes
### 

# Isolate blood_eos_log from the clinical traits
BAL_eos_ct_log = as.data.frame(datTraits$BAL_eos_ct_log);
BAL_neut_ct_log = as.data.frame(datTraits$BAL_neut_ct_log)
BAL_eos_p_log = as.data.frame(datTraits$BAL_eos_p_log);
BAL_neut_p_log = as.data.frame(datTraits$BAL_neut_p_log)

names(BAL_eos_ct_log) = "BAL_eos_ct_log"
names(BAL_neut_ct_log) = "BAL_neut_ct_log"
names(BAL_eos_p_log) = "BAL_eos_p_log"
names(BAL_neut_p_log) = "BAL_neut_p_log"



# Add the BAL_eos_ct_log to existing module eigengenes
MEs<-MEs[rownames(MEs)%in%rownames(datTraits),]
MET = orderMEs(cbind(MEs, BAL_eos_ct_log, BAL_neut_ct_log, BAL_eos_p_log, BAL_neut_p_log ))
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
output_folder<-file.path("./reports/wgcna/bronch")
write.table(df_genelist, file = file.path(output_folder,"bronch_wgcna_genelist_batch12346.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)

# Loop through each element in the list and write to a CSV file
for (module_name in names(list_of_dfs)) {
  # Extract the 'genes' column from each dataframe
  gene_data <- list_of_dfs[[module_name]]$genes
  
  # Create a file path
  file_path <- file.path(output_folder,paste0(module_name,"batch12346", ".txt"))
  
  # Write to a text file without row names or column names
  write.table(gene_data, file=file_path, quote=FALSE, row.names=FALSE, col.names=FALSE)
}
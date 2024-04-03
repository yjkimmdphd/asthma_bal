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

# load count data from seq batch 1-4
counts<-file.path("./resources/working_data/copy_of_batch1234_readcount_matrix_allsamples.afterQC.txt") # original RNAseqs count data table is in the MS_asthma folder 
counts<-if(file.exists(counts)){read.delim(counts)}
counts.ID<-colnames(counts)

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
phenotype<-file.path("./resources/processed_data/Nasal_Biomarkers_BAL_transformed.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
numeric_id<-phenotype$ID
phenotype$ID<-sprintf("%03d",numeric_id) # adds padded zeros in front of the subject ID numbers
phenotype<-mutate(phenotype,pos_cellcount=phenotype[,source.cell]>0)%>%arrange(ID) # check which cell counts are positive. 

###########################################################################################
## subset phenotype data for which the samples exist for nasal/bronchial RNAseq experiments   
###########################################################################################

nID<-paste0("N",phenotype$ID) # N*** indicates nasal sample ID, sequence data is not available for all as of 2023-10-04
bID<-paste0("B",phenotype$ID) # B*** indicates bronchial sample ID, sequence data is not available for all as of 2023-10-04

#find which sampleIDs have count data

nexist<-nID%in%counts.ID # find which subjects s/p nasal and had bronchial sample RNAseq completed. Nasal samples in batch 1-4 only sequenced to ID 337
bexist<-bID%in%counts.ID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
nsample<-nID[nexist] # nasal sample ID in the readcount matrix (batch 1-4) that has BAL phenotype data
bsample<-bID[bexist] # bronchial sample ID in the readcount matrix (batch 1-4) that has BAL phenotype data

# subset phenotype dataframe for which there is count data

nphen<-phenotype[phenotype$ID%in%substring(nsample,2),] # phenotype table with nsample
bphen<-phenotype[phenotype$ID%in%substring(bsample,2),] # phenotype table with bsample

nphen<-mutate(nphen, SampleID=nsample)%>%relocate(SampleID, .before=1) # include sample ID for nasal RNAseq samples
bphen<-mutate(bphen, SampleID=bsample)%>%relocate(SampleID, .before=1) # include sample ID for bronchial RNAseq samples

# left join batch info table with nasal/bronchial phenotype table  
## get batch information
batch<-file.path("./resources/processed_data/asthma_nasal_bronchial_batch1234_info.txt")
batch.info<-if(file.exists(batch)){read.delim(batch)}

## define function join_phenotype_batch_info. p is phenotype table. b is batch info table. Factorize the batch info. 
join_phenotype_batch_info<-function(p,b){
  table<-left_join(p,b, by="SampleID")
  table$Batch<-factor(table$Batch, levels=unique(table$Batch))
  return(table)
}
nphen<-join_phenotype_batch_info(nphen,batch.info)
bphen<-join_phenotype_batch_info(bphen,batch.info)
nphen<-nphen%>%mutate(IsBatch4 = Batch == "batch4")
bphen<-bphen%>%mutate(IsBatch4 = Batch == "batch4")
# scale the cell count information 
# Mike Love: I'd recommend however to center and scale the numeric covariates for model fitting improvement.
# scale the columns named with source.cell.log
nphen<-mutate_at(nphen,vars(all_of(source.cell.log)),scale) # scales and mutates all log-transformed cell counts 
bphen<-mutate_at(bphen,vars(all_of(source.cell.log)),scale) # scales and mutates all log-transformed cell counts 

# decide which analysis to perform, then set the phenotype data as phen

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
#set colData (phenotype data) for nasal RNAseq experiments
##########################################################

# if nasal cell analysis, phen<-nphen. 
phen<-nphen  # If nasal analysis, use this

# make df which is a list of sampleIDs which has non-NA values for each of the cell types (i.e., cell count >=0)
pi<-lapply(phen[,source.cell.log],function(data){a<-!is.na(data);return(a)})
df<-vector("list",length=10) # list of data framese used as an input for deseq2. all cell counts
names(df)<-paste(source.cell.log,"all",sep="_")
for(i in 1:10){
  df[[i]]<-phen[pi[[i]],c("SampleID",source.cell.log[i], "Batch","IsBatch4")]
}
print(sapply(df,dim)[1,]) # shows number of 

# select RNAseq counts
id<-phen$SampleID
cols<-colnames(counts)%in%id
ct<-counts[,cols] # First column is actually gene name 
genes<-counts$SampleID
rownames(ct)<-genes

###############
#
#### WGCNA ####
#
###############

enableWGCNAThreads()

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
abline(h = 2.5e6, col = "red")

# height of 2.5e6 removes sample ID N277
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 2.5e6, minSize = 10) #returns numeric vector
#Remove outlier
expression.data <- expression[cut.sampleTree==1, ]
dim(expression)
dim(expression.data)

##
# 4. network formulation
##

# Determining the Soft Power Threshold
spt <- pickSoftThreshold(expression.data) 

# Plot the R2 values as a function of the soft thresholds
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")

#Plot mean connectivity as a function of soft thresholds
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")

# spt should be 11, since after it, the power is above 0.8, for a short period lol

#calling the adjacency fx
softPower <- 11
adjacency <- adjacency(expression.data, power = softPower)

###
# 5. module Construction
### 
# convert the adjacency matrix into a TOM similarity matrix we can call the WGCNA function 
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1-TOM

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
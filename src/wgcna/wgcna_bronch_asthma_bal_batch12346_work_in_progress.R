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
phenotype<-file.path("./resources/processed_data/Nasal_Biomarkers_BAL_transformed.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
numeric_id<-phenotype$ID
phenotype$ID<-sprintf("%03d",numeric_id) # adds padded zeros in front of the subject ID numbers
phenotype<-mutate(phenotype,pos_cellcount=phenotype[,source.cell]>0)%>%arrange(ID) # check which cell counts are positive. 



###########################################################################################
## subset phenotype data for which the samples exist for nasal/bronchial RNAseq experiments   
###########################################################################################
bID<-paste0("B",phenotype$ID) # B*** indicates bronchial sample ID, sequence data is not available for all as of 2023-10-04
bexist<-bID%in%counts.ID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
bsample<-bID[bexist] # bronchial sample ID in the readcount matrix (batch 1-4,6) that has BAL phenotype data
bphen<-phenotype[phenotype$ID%in%substring(bsample,2),] # phenotype table with bsample
bphen<-mutate(bphen, SampleID=bsample)%>%relocate(SampleID, .before=1) # include sample ID for bronchial RNAseq samples

# left join batch info table with nasal/bronchial phenotype table  
## get batch information
batch<-file.path("./resources/raw_data/MS_asthma/MS_asthma_phenotype.batch12346.final.csv")
batch.info<-if(file.exists(batch)){read.csv(batch)}
bronch.batch.info<-batch.info[1:75,2:3]
bronch.batch.info$SampleID<-substr(bronch.batch.info$SampleID,1,4)

## define function join_phenotype_batch_info. p is phenotype table. b is batch info table. Factorize the batch info. 
join_phenotype_batch_info<-function(p,b){
  table<-left_join(p,b, by="SampleID")
  table$Batch<-factor(table$Batch, levels=unique(table$Batch))
  return(table)
}
bphen<-join_phenotype_batch_info(bphen,bronch.batch.info)
bphen<-bphen%>%mutate(IsBatch4 = Batch == "batch4")
# scale the cell count information 
# Mike Love: I'd recommend however to center and scale the numeric covariates for model fitting improvement.
# scale the columns named with source.cell.log

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


##########################################################
#set colData (phenotype data) for bronchial RNAseq experiments
##########################################################

phen<-bphen

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

save(spt,"./resources/processed_data/Rdata/wgcna_bronch_spt_batch12346.Rdata") # save the spt as Rdata
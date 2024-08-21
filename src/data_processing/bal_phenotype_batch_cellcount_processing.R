################
# data prep
# raw count: "./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt"
# goal: normalizes the gene counts with TMM method 
################


library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)
library(DESeq2)

# load biomarker phenotype file saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")
phen<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}

# subset phenotype data based on nasal vs bronchial samples
nphen<-phen%>%filter(grepl("^N",SampleID),!grepl("F",SampleID))
nphen_sampleID<-nphen$SampleID

bphen<-phen%>%filter(grepl("^B",SampleID))
bphen_sampleID<-bphen$SampleID

# load gene count table

countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}

# remove sampleIDs starting with T 
sampleIDs<-colnames(counts)
b_n_sampleIDs<-!grepl("^T",sampleIDs)
sampleIDs<-sampleIDs[b_n_sampleIDs]
counts<-counts[sampleIDs]
row.names(counts)<-counts[,"SampleID"]
counts<-counts[,-1]

# subset count matrix to nasal and bronchial samples 
sampleIDs<-colnames(counts)
sID_list<-list()
sID_list[[1]]<-sampleIDs[grep("^B",colnames(counts))]
sID_list[[2]]<-sampleIDs[grep("^N",colnames(counts))]
names(sID_list)<-c("bronch","nasal")
print(sID_list)

bronch.counts<-counts[,c(sID_list$bronch)]
bronch.counts<-bronch.counts[,bphen_sampleID]

nasal.counts<-counts[,c(sID_list$nasal)]
nasal.counts<-nasal.counts[,nphen_sampleID]
###################################
# custom functions for gene filter
###################################
# should load the following fx:
## filter_low_expressed_genes_method2: Filters low counts genes using TMM normalized lcpm as a cutoff point. Requires 'limma'
## rowgenes_counttable: changes the row names of the count table with gene names

source("./src/function/deg_custom_functions_v2.R")

################ prepare normalized gene count table for bronchial data######################
phen<-bphen
# filtering counts table to remove low expressed genes

## select RNAseq counts
id<-phen$SampleID
cols<-colnames(bronch.counts)%in%id
ct<-bronch.counts[,cols] # First column is actually gene name 
genes<-rownames(ct)


## Filter counts with low gene counts 
c2<-filter_low_expressed_genes_method2(ct,round(length(id)*0.1,0)) # CPM cut off is 10/M + 2/L, M is median library size, and L is mean lib size, in million. if M and L are 50 million each, cut off would be 0.2 + 0.04=0.24 CPM in 10% of samples 
ct<-rowgenes_counttable(ct,c2) # low bcounts will be filtered 

### normalize counts DESEq2
x<-as.matrix(ct)
vsd<-varianceStabilizingTransformation(x,blind=TRUE)
boxplot(x,main="")
normalized_count_table_path<-"./resources/processed_data/bronch_batch12346_normalized_ct.txt"
write.table(vsd,normalized_count_table_path,row.names = FALSE,sep = "\t")


################ prepare normalized gene count table for nasal data######################
################
# gsea data prep
################


library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)

countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
bronch.samples<-grepl("^B",colnames(counts))
bronch.samples<-which(bronch.samples==TRUE)
bronch.counts<-counts[,c(1,bronch.samples)]
colnames(bronch.counts)[bronch.samples]<-substr(colnames(bronch.counts)[bronch.samples],1,4)
head(bronch.counts)
counts.ID<-colnames(bronch.counts)
counts<-bronch.counts[,-1]
rownames(counts)<-bronch.counts$SampleID
# prepare normalized gene count table 

x<-counts
L<-mean(colSums(x))*1e-6 # mean library size
M<-median(colSums(x))*1e-6 # median library size
lcpm.cutoff <- log2(10/M + 2/L) # lcpm cutoff for filtering genes with very low counts

### normalize counts with TMM
norm.factor<-calcNormFactors(x, method = "TMM")
sample.size<-length(colnames(x))
for(i in 1:sample.size){
  x[,i]<-x[,i]/norm.factor[i]
}

normalized.count<-cbind(name=rownames(counts),counts)

write.table(normalized.count,"./resources/processed_data/gsea/bronch_batch12346_normalized_ct.txt",row.names = FALSE,sep = "\t")

### calculate lcpm based on TMM normalized counts 
lcpm.x<-cpm(x,log=TRUE)


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


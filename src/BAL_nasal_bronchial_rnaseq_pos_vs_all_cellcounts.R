library(limma)
library(edgeR)
library(dplyr)
library(DESeq2)

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

# asthma biomarker phenotype file saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/Nasal_Biomarkers_BAL_transformed.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
numeric_id<-phenotype$ID
phenotype$ID<-sprintf("%03d",numeric_id) # adds padded zeros in front of the subject ID numbers

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
deg.design<-paste("~",source.cell.log,"+ Batch")

#######################################################################
## find subject assignment ID with nasal and nasal cell RNAseq data  
#######################################################################
nID<-paste0("N",phenotype$ID) # N*** indicates nasal sample ID, sequence data is not available for all as of 2023-10-04
bID<-paste0("B",phenotype$ID) # B*** indicates bronchial sample ID, sequence data is not available for all as of 2023-10-04

nexist<-nID%in%counts.ID # find which subjects s/p nasal and had bronchial sample RNAseq completed. Nasal samples in batch 1-4 only sequenced to ID 337
bexist<-bID%in%counts.ID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
nsample<-nID[nexist] # nasal sample ID in the readcount matrix (batch 1-4) that has BAL phenotype data
bsample<-bID[bexist] # bronchial sample ID in the readcount matrix (batch 1-4) that has BAL phenotype data

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

# scale the cell count information 
# Mike Love: I'd recommend however to center and scale the numeric covariates for model fitting improvement.
# scale the columns named with source.cell.log
nphen<-mutate_at(nphen,vars(all_of(source.cell.log)),scale) # scales and mutates all log-transformed cell counts 
bphen<-mutate_at(bphen,vars(all_of(source.cell.log)),scale) # scales and mutates all log-transformed cell counts 

# decide which analysis to perform, then set the phenotype data as phen

phen<-bphen
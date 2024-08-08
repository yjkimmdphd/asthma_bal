################
# data prep
# raw count: "./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt"
# phenotype data: "./resources/processed_data/phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv" (log transformed cell counts are not scaled )
# goal1: normalizes the gene counts with TMM method 
# goal2: scale the log-transformed cell counts
################


library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)

countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}

sampleIDs<-colnames(counts)
b_n_sampleIDs<-!grepl("^T",sampleIDs)
sampleIDs<-sampleIDs[b_n_sampleIDs]
counts<-counts[sampleIDs]
row.names(counts)<-counts[,"SampleID"]
counts<-counts[,-1]


colnames(counts)<-substr(colnames(counts), 1, 4)# remove "sampleID"  # nasal samples have duplicates
sampleIDs<-colnames(counts)

sID_list<-list()
sID_list[[1]]<-sampleIDs[grep("^B",colnames(counts))]
sID_list[[2]]<-sampleIDs[grep("^N",colnames(counts))]
names(sID_list)<-c("bronch","nasal")
print(sID_list)

bronch.counts<-counts[,c(sID_list$bronch)]

# prepare normalized gene count table 

x<-bronch.counts
L<-mean(colSums(x))*1e-6 # mean library size
M<-median(colSums(x))*1e-6 # median library size
lcpm.cutoff <- log2(10/M + 2/L) # lcpm cutoff for filtering genes with very low counts

### normalize counts with TMM
norm.factor<-calcNormFactors(x, method = "TMM")
sample.size<-length(colnames(x))
for(i in 1:sample.size){
  x[,i]<-x[,i]/norm.factor[i]
}

normalized.count<-cbind(name=rownames(counts),x)

normalized_count_table_path<-"./resources/processed_data/bronch_batch12346_normalized_ct.txt"
write.table(normalized.count,normalized_count_table_path,row.names = FALSE,sep = "\t")


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
phenotype<-file.path("./resources/processed_data/phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
numeric_id<-phenotype$ID
phenotype$ID<-sprintf("%03d",numeric_id) # adds padded zeros in front of the subject ID numbers

# scale the cell count information 
# Mike Love: I'd recommend however to center and scale the numeric covariates for model fitting improvement.
# scale the columns named with source.cell.log

bphen<-mutate_at(phenotype,vars(all_of(source.cell.log)),scale) # scales and mutates all log-transformed cell counts 


bphen_path<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")
write.csv(bphen,bphen_path)

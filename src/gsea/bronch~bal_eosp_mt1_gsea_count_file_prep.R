################
# gsea data prep
################


library(tidyverse)

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

bphen_path<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")
bphen<-read.csv(bphen_path)%>%filter(grepl("^B",SampleID))


normalized_count_table_path<-"./resources/processed_data/bronch_batch12346_normalized_ct.txt"

counts<-if(file.exists(normalized_count_table_path)){read.delim(normalized_count_table_path, check.names = FALSE)}
bronch.counts<-counts
rownames(bronch.counts)<-bronch.counts$name

# subset the normalized count table that only contains genes that were identified as DEG in bronch deg ~ BAL Eos p>1 + batch12346 
bphen_more_than_1<-bphen$BAL_eos_p>1
names(bphen_more_than_1)<-bphen$SampleID
write.table(bphen_more_than_1,"./resources/processed_data/gsea/bal_eos_p_more_than_1p.txt",sep="\t")

count_bphen_more_than_1<-bronch.counts[,!is.na(match(colnames(bronch.counts),names(bphen_more_than_1)))]
deg_more_than_1<-read.csv("./reports/local_only/deg_2024-08-10_deg_bronch~bal_cell(dichot)+batch12346/deg_bronch_res_sig_5_~ bal_Eos_p_more_1 + Batch_2024-08-10_.csv")

genes_bphen_more_than_1<-row.names(count_bphen_more_than_1)%in%deg_more_than_1_down$X
count_deg_bphen_more_than_1<-count_bphen_more_than_1[genes_bphen_more_than_1,]
count_file_path<-file.path("./resources/processed_data/gsea/count_deg_bphen_more_than_1.txt")
write.table(count_deg_bphen_more_than_1,count_file_path,sep="\t",col.names = NA)

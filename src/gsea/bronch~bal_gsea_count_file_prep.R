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
bphen<-read.csv(bphen_path)

normalized_count_table_path<-"./resources/processed_data/bronch_batch12346_normalized_ct.txt"

counts<-if(file.exists(normalized_count_table_path)){read.delim(normalized_count_table_path, check.names = FALSE)}
bronch.counts<-counts
rownames(bronch.counts)<-bronch.counts$name

# subset the normalized count table that only contains genes that were identified as DEG in bronch deg ~ BAL Eos p>3 + batch12346 
bphen_more_than_3<-bphen$BAL_eos_p>3
names(bphen_more_than_3)<-bphen$SampleID
write.table(bphen_more_than_3,"./resources/processed_data/gsea/bal_eos_p_more_than_3p.txt",sep="\t")

count_bphen_more_than_3<-bronch.counts[,!is.na(match(colnames(bronch.counts),names(bphen_more_than_3)))]
deg_more_than_3<-read.csv("./reports/local_only/deg_bronch~bal_cell(dichot)+batch12346/deg_bronch_dichot_res_all_5_~ bal_Eos_p_more_3 + Batch_2024-05-29_.csv")

genes_bphen_more_than_3<-row.names(count_bphen_more_than_3)%in%deg_more_than_3$X
count_deg_bphen_more_than_3<-count_bphen_more_than_3[genes_bphen_more_than_3,]
count_file_path<-file.path("./resources/processed_data/gsea/count_deg_bphen_more_than_3.txt")
write.table(count_deg_bphen_more_than_3,count_file_path,sep="\t",col.names = NA)


# subset the normalized count table that only contains genes that were identified as DEG in bronch deg ~ BAL Eos p>1 + batch12346 
bphen_more_than_1<-bphen$BAL_eos_p>1
names(bphen_more_than_1)<-bphen$SampleID
write.table(bphen_more_than_1,"./resources/processed_data/gsea/bal_eos_p_more_than_1p.txt",sep="\t")

count_bphen_more_than_1<-bronch.counts[,!is.na(match(colnames(bronch.counts),names(bphen_more_than_1)))]
deg_more_than_1<-read.csv("./reports/local_only/deg_bronch~bal_cell(dichot)+batch12346/deg_bronch_dichot_res_all_4_~ bal_Eos_p_more_1 + Batch_2024-05-29_.csv")

genes_bphen_more_than_1<-row.names(count_bphen_more_than_1)%in%deg_more_than_1$X
count_deg_bphen_more_than_1<-count_bphen_more_than_1[genes_bphen_more_than_1,]
count_file_path<-file.path("./resources/processed_data/gsea/count_deg_bphen_more_than_1.txt")
write.table(count_deg_bphen_more_than_1,count_file_path,sep="\t",col.names = NA)


# BAL AEC > 1.15; bronch deg ~ BAL Eos aec>1.15 + batch12346 
bphen_aec_m_1.15<-bphen$BAL_eos_ct>1.15
names(bphen_aec_m_1.15)<-bphen$SampleID
write.table(bphen_aec_m_1.15,"./resources/processed_data/gsea/bal_eos_ct_m_1-15.txt",sep="\t")

count_bphen_aec_m_1.15<-bronch.counts[,!is.na(match(colnames(bronch.counts),names(bphen_aec_m_1.15)))]
deg_bal_aec_m_1.15<-read.csv("./reports/local_only/deg_bronch~bal_cell(dichot)+batch12346/deg_bronch_dichot_res_all_2_~ bal_AEC_more_1.15 + Batch_2024-05-29_.csv")

genes_bphen_bal_aec_m_1.15<-row.names(count_bphen_aec_m_1.15)%in%deg_bal_aec_m_1.15$X
count_deg_bphen_bal_aec_m_1.15<-count_bphen_aec_m_1.15[genes_bphen_bal_aec_m_1.15,]
count_file_path<-file.path("./resources/processed_data/gsea/count_deg_bphen_bal_aec_m_1.15.txt")
write.table(count_deg_bphen_bal_aec_m_1.15,count_file_path,sep="\t",col.names = NA)



# subset the normalized count table that only contains genes that has BAL ANC > 0 
bphen_anc<-bphen[!is.na(bphen$BAL_neut_ct),]%>%filter(BAL_neut_ct>0)
bphen_bal_anc<-bphen_anc$BAL_neut_ct_log

names(bphen_bal_anc)<-bphen_anc[,"SampleID"]
write.table(bphen_bal_anc,"./resources/processed_data/gsea/bal_anc.txt",sep="\t")

count_bphen_bal_anc_m_0<-bronch.counts[,!is.na(match(colnames(bronch.counts),names(bphen_bal_anc)))]
deg_bal_anc_m_0<-read.csv("./reports/local_only/deg_bronch~bal_cell(continuous)+batch12346/deg_2024-05-29/deg_bronch_continuous_poscell+batch12346_res_all_3_~ BAL_neut_ct_log + Batch_2024-05-29_.csv")

genes_bphen_bal_anc_m_0<-row.names(count_bphen_bal_anc_m_0)%in%deg_bal_anc_m_0$X
count_deg_bphen_bal_anc_m_0<-count_bphen_bal_anc_m_0[genes_bphen_bal_anc_m_0,]
count_file_path<-file.path("./resources/processed_data/gsea/count_deg_bphen_bal_anc_m_0.txt")
write.table(count_deg_bphen_bal_anc_m_0,count_file_path,sep="\t",col.names = NA)



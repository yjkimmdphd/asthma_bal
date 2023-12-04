### 
# making dataFrame containing information about 
# DEG analysis design, coldata, countdata
# DESeqDataSetFromMatrix
# for bronchial RNAseq
###

################################
## load phenotype and batch data
################################

# asthma biomarker phenotype file, nasal, saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/asthma-phenotype-Rnaseq-2023-12-01.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
org.numb<-phenotype$subjID
phenotype$subjID<-sprintf("%03d",org.numb) # adds padded zeros in front of the subject ID numbers


source.cell.log<-c("BAL_Eos_ct_log",
               "BAL_Eos_perc_log",
               "BAL_neut_ct_log10",
               "BAL_neut_perc_log",
               "BAL_WBC_log",
               "serum_Eos_log10",
               "serum_Eos_perc_log",
               "serum_Neut_log10",
               "serum_Neut_perc_log",
               "serum_WBC_log10")
source.cell<-c("BAL_Eos_ct",
               "BAL_Eos_perc",
               "BAL_neut_ct",
               "BAL_neut_perc",
               "BAL_WBC",
               "serum_Eos",
               "serum_Eos_perc",
               "serum_Neut",
               "serum_Neut_perc",
               "serum_WBC")
# following samples have no value in BAL-neut, serum-eos, or serum-Neut
bphen.na<-bphen[is.na(bphen)%>%rowSums()>0,]

# save data.frame of samples have missing cell counts
exclude.bronch<-data.frame(source_cell=c("BAL_neut_ct","BAL_neut_perc",rep(c("serum_Eos","serum_Eos_perc","serum_Neut", "serum_Neut_perc","serum_WBC"),each=6)),
                           SampleID=c("B332","B332",bphen.na$SampleID[1:6]%>%rep(5)))

write.csv(exclude.bronch,"./resources/processed_data/bronch_samples_excluded.csv",row.names = FALSE)

                       
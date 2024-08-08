# making BAL phenotype data that has asthma phenotype and sequence batch information  
processed_data_folder<-"./resources/processed_data/"
phenotype_file<-file.path(processed_data_folder,"yjk_Nasal_Biomarkers_BAL_08172023.csv")
batch_info_file<-file.path(processed_data_folder,"id_sampleID_batch_12346.csv")
phenotype<-read.csv(phenotype_file)
batch_info<-read.csv(batch_info_file)

phenotype<-inner_join(batch_info,phenotype,by="ID")
phenotype<-mutate(phenotype,BAL_eos_ct=BAL_eos_p/100*BAL_wbc,
                  BAL_neut_ct=BAL_neut_p/100*BAL_wbc)

# log transform. code from asthma_bal_phenotype_data_processing.R
phenotype<-mutate(phenotype, 
                  BAL_eos_ct_log=log(BAL_eos_ct+0.001)%>%round(2),
                  BAL_eos_p_log = log(BAL_eos_p+0.001)%>%round(2), 
                  BAL_wbc_log=log(BAL_wbc+0.001)%>%round(2), 
                  BAL_neut_p_log=log(BAL_neut_p+0.001)%>%round(2), 
                  BAL_neut_ct_log=log10(BAL_neut_ct+0.001)%>%round(2),
                  blood_eos_log=log10(blood_eos+0.001)%>%round(2),
                  blood_eos_p_log=log(blood_eos_p+0.001)%>%round(2),
                  blood_neut_log=log10(blood_neut+0.001)%>%round(2),
                  blood_neut_p_log=log(blood_neut_p+0.001)%>%round(2),
                  blood_wbc_log=log10(blood_wbc)%>%round(2))

# Desired order of specific columns
desired_order <- c("BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log", 
                   "BAL_wbc_log", "blood_eos_log", "blood_eos_p_log", "blood_neut_log", 
                   "blood_neut_p_log", "blood_wbc_log", "BAL_eos_ct", "BAL_eos_p", 
                   "BAL_neut_ct", "BAL_neut_p", "BAL_wbc", "blood_eos", "blood_eos_p", 
                   "blood_neut", "blood_neut_p", "blood_wbc")
phenotype[,desired_order]

# Reorder columns in the data frame
df <- phenotype%>%
  select(everything(),all_of(desired_order))

write.csv(df,"./resources/processed_data/phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv",row.names = FALSE)

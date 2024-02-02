library(tidyverse)

# asthma biomarker phenotype file saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/Nasal_Biomarkers_BAL.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
numeric_id<-phenotype$ID
phenotype$ID<-sprintf("%03d",numeric_id) # adds padded zeros in front of the subject ID numbers

phenotype<-mutate(phenotype, 
                  BAL_eos_ct=BAL_wbc*BAL_eos_p*0.01,
                  BAL_neut_ct=BAL_wbc*BAL_neut_p*0.01)

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
                  blood_wbc_log=log10(blood_wbc)%>%round(2))%>%
  select(ID,
         BAL_eos_ct_log,BAL_eos_p_log,BAL_neut_ct_log,BAL_neut_p_log,BAL_wbc_log,blood_eos_log,blood_eos_p_log,blood_neut_log,blood_neut_p_log,blood_wbc_log, 
         BAL_eos_ct,BAL_eos_p,BAL_neut_ct,BAL_neut_p,BAL_wbc,blood_eos,blood_eos_p,blood_neut,blood_neut_p,blood_wbc, 
         everything())

write.csv(phenotype,"./resources/processed_data/Nasal_Biomarkers_BAL_transformed.csv",row.names = FALSE)


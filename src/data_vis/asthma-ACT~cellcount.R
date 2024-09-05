########
# data visualization for asthma phenotype and bal cell count
########
library(dplyr)
# load biomarker phenotype file saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")
phenotype_big<-file.path("./resources/processed_data/Nasal_Biomarkers_BAL_transformed_with_raceinfo.csv")
phen<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
phen_big<-if(file.exists(phenotype_big)){read.csv(phenotype_big, row.names = NULL)}%>%select("ID","Age_at_visit","sex","ethnicity","Race")
phen<-left_join(phen,phen_big,by="ID")
bphen<-phen%>%filter(grepl("^B",SampleID))

# BAL Eos% 1% cutoff
## t-test
bphen_more_than_1<-bphen$BAL_eos_p>1

eos_p_1<-bphen%>%mutate(BAL_eos_p_m_1 = BAL_eos_p>1)%>%group_by(BAL_eos_p_m_1)
eos_p_1<-eos_p_1%>%arrange(by=BAL_eos_p_m_1)%>%select(SampleID,BAL_eos_p_m_1,asthma_phen_ACT.score,sex,Race)

t_result<-t.test(asthma_phen_ACT.score~BAL_eos_p_m_1, data =eos_p_1)
boxplot(asthma_phen_ACT.score~BAL_eos_p_m_1,data=eos_p_1,
        xlab="BAL Eos% >1",
        ylab="ACT score",
        main=paste("ACT ~ BAL Eos% >1; t-test p-value",signif(t_result$p.value, digits=2),sep=":"))

eos_p_1 %>%
  dplyr::summarize(median_ACT = median(asthma_phen_ACT.score, na.rm = TRUE))

## linear model
lm_result<-lm(asthma_phen_ACT.score~BAL_eos_p_m_1+sex+Race, data =eos_p_1) # lm controlling for sex and race  
summary(lm_result)

# BAL Eos% 3% cutoff
## t-test 
bphen_more_than_3<-bphen$BAL_eos_p>3

eos_p_3<-bphen%>%mutate(BAL_eos_p_m_3 = BAL_eos_p>3)%>%group_by(BAL_eos_p_m_3)
eos_p_3<-eos_p_3%>%arrange(by=BAL_eos_p_m_3)%>%select(SampleID,BAL_eos_p_m_3,asthma_phen_ACT.score,sex,Race)

t_result<-t.test(asthma_phen_ACT.score~BAL_eos_p_m_3,data=eos_p_3)
boxplot(asthma_phen_ACT.score~BAL_eos_p_m_3,data=eos_p_3,
        xlab="BAL Eos% >3",
        ylab="ACT score",
        main=paste("ACT ~ BAL Eos% >3; t-test p-value",signif(t_result$p.value, digits=2),sep=":"))

eos_p_3%>%
  dplyr::summarize(median_ACT = median(asthma_phen_ACT.score, na.rm = TRUE))

## linear model
lm_result<-lm(asthma_phen_ACT.score~BAL_eos_p_m_3+sex+Race, data =eos_p_3) # lm controlling for sex and race  
summary(lm_result)

# ACT~BAL AEC 1 above or below
## t-test
eos_aec_1<-bphen%>%mutate(BAL_eos_aec_m_1 = BAL_eos_ct>1)%>%group_by(BAL_eos_aec_m_1)
eos_aec_1<-eos_aec_1%>%arrange(by=BAL_eos_aec_m_1)%>%select(SampleID,BAL_eos_aec_m_1,asthma_phen_ACT.score,sex,Race)

t_result<-t.test(asthma_phen_ACT.score~BAL_eos_aec_m_1,data=eos_aec_1)
boxplot(asthma_phen_ACT.score~BAL_eos_aec_m_1,data=eos_aec_1,
        xlab="BAL AEC >1",
        ylab="ACT score",
        main=paste("ACT ~ BAL AEC above 1; t-test p-value",signif(t_result$p.value, digits=2),sep=":"))

eos_aec_1%>%
  dplyr::summarize(median_ACT = median(asthma_phen_ACT.score, na.rm = TRUE))


## linear model
lm_result<-lm(asthma_phen_ACT.score~BAL_eos_aec_m_1+sex+Race, data =eos_aec_1) # lm controlling for sex and race  
summary(lm_result)

# ACT~BAL ANC 1 above or below
## t-test
neut_anc_1<-bphen%>%mutate(BAL_neut_anc_m_1 = BAL_neut_ct>1)%>%group_by(BAL_neut_anc_m_1)
neut_anc_1<-neut_anc_1%>%arrange(by=BAL_neut_anc_m_1)%>%select(SampleID,BAL_neut_anc_m_1,asthma_phen_ACT.score,sex,Race)

t_result<-t.test(asthma_phen_ACT.score~BAL_neut_anc_m_1,data=neut_anc_1)
boxplot(asthma_phen_ACT.score~BAL_neut_anc_m_1,data=neut_anc_1,
        xlab="BAL ANC >1",
        ylab="ACT score",
        main=paste("ACT ~ BAL ANC above 1; t-test p-value",signif(t_result$p.value, digits=2),sep=":"))
## linear model
lm_result<-lm(asthma_phen_ACT.score~BAL_neut_anc_m_1+sex+Race, data =neut_anc_1) # lm controlling for sex and race  
summary(lm_result)

# Add the p-value to the plot
text(x = 1.5, 
     y = max(neut_anc_1$asthma_phen_ACT.score) * 0.9, 
     labels = paste("p-value =", format(t_result$p.value, digits = 2)), 
     pos = 4, 
     col = "red")
neut_anc_1%>%
  dplyr::summarize(median_ACT = median(asthma_phen_ACT.score, na.rm = TRUE))


  
  
  
  
  
  
  
  
  
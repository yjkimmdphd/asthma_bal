# find how many people had at least one ED visit or admission in the cohort
library(dplyr)
countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
rownames(counts)<-counts[,"SampleID"]

phen_org<-read.csv("./resources/processed_data/BiomarkersOfAsthma_original.csv")
phen_org<-select(phen_org,ID,Age_at_visit,Sex,Ethnicity,Race,asthma_sx_onset,asthma_dx_by_md,age_asthma_dx,asthma_ED,ED_visits,asthma_admit_hx,admit_count,reg_asthma_med)

phen_proc<-read.csv("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv")

phen<-left_join(phen_org,phen_proc,by=join_by(ID))
phen<-phen%>%filter(grepl("^B",SampleID),!is.na(BAL_wbc))
phen<-phen[phen$SampleID%in%colnames(counts),] # select for sampleIDs that underwent RNAseq 
dim(phen)

phen$exac<-rowSums(select(phen,ED_visits,admit_count),na.rm=TRUE)

phen<- phen %>%
  mutate(bal_eos_p_mt1 = factor(case_when(BAL_eos_p > 1 ~ 1 , #"high_eos"
                                          BAL_eos_p <= 1 ~ 0), # "low_eos"
                                levels = c(0, 1)))


phen_summary<-phen%>%summarize(age_visit_mean= mean(Age_at_visit),
                 age_visit_sd = sd(Age_at_visit),
                 act_mean = mean(asthma_phen_ACT.score, na.rm=TRUE),
                 act_sd = sd(asthma_phen_ACT.score, na.rm=TRUE),
                 exac_mean = mean(exac,na.rm = TRUE),
                 exac_sd = sd(exac,na.rm = TRUE),
                 exac_rate = sum(exac>=1)/sum(exac>=0))


phen_all_summary<-signif(phen_summary,2)
print(phen_all_summary)                 


phen_summary<-phen%>%group_by(bal_eos_p_mt1)%>%filter(!is.na(bal_eos_p_mt1))%>%summarize(age_visit_mean= mean(Age_at_visit),
                                     age_visit_sd = sd(Age_at_visit),
                                     sex_f_mean = mean(Sex=="Female"),
                                     sex_m_mean = mean(Sex=="Male"),
                                     sex_f_sum = sum(Sex=="Female"),
                                     sex_m_sum = sum(Sex=="Male"),
                                     act_mean = mean(asthma_phen_ACT.score, na.rm=TRUE),
                                     act_sd = sd(asthma_phen_ACT.score, na.rm=TRUE),
                                     exac_mean = mean(exac,na.rm = TRUE),
                                     exac_sd = sd(exac,na.rm = TRUE),
                                     exac_rate = sum(exac>=1)/sum(exac>=0),
                                     exac_yes = sum(exac>=1),
                                     exac_non = sum(exac==0),
                                     total = sum(!is.na(bal_eos_p_mt1)))

print(phen_summary)
view(phen_summary)

exac_ttest<-t.test(exac~ bal_eos_p_mt1, data=phen)
age_ttest<-t.test(Age_at_visit~ bal_eos_p_mt1, data=phen)
act_ttest<-t.test(asthma_phen_ACT.score~ bal_eos_p_mt1, data=phen)
exac_ttest<-t.test(exac~ bal_eos_p_mt1, data=phen)
sex_fisher<-fisher.test(as.matrix(phen_summary[,c("sex_f_sum","sex_m_sum")]))
exac_fisher<-fisher.test(as.matrix(phen_summary[,c("exac_yes","exac_non")]))

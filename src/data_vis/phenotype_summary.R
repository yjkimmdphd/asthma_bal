# find how many people had at least one ED visit or admission in the cohort
library(tidyverse)
countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
rownames(counts)<-counts[,"SampleID"]

phen<-read.csv("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2025-02-14.csv")
pft<-c("FEV1","FEV1_percent","FVC","FVC_percent","FEV1.FVC")
colnames(phen)

phen<-phen%>%filter(grepl("^B",SampleID),BAL_eos_p>=0 & BAL_neut_p>=0) # instead of !is.na(BAL_wbc) - this results in n=70 now vs n=73 previously
phen<-phen[phen$SampleID%in%colnames(counts),] # select for sampleIDs that underwent RNAseq 
dim(phen)

# fill in NA values for ED visit counts and admit counts
str(phen$ED_visits)
phen$ED_visits_new<-phen$ED_visits
no_ed_visits<-which(phen$asthma_ED=="No")
phen$ED_visits_new[no_ed_visits]<-0
if(summarise(phen,diff_ed=mean(ED_visits-ED_visits_new,na.rm=TRUE))==0){
  phen$ED_visits<-phen$ED_visits_new
}

str(phen$asthma_admit_hx)
phen$admit_count_new<-phen$admit_count
no_admits<-which(phen$asthma_admit_hx=="No")
phen$admit_count_new[no_admits]<-0
if(summarise(phen,diff_admit=mean(admit_count-admit_count_new,na.rm=TRUE))==0){
  phen$admit_count<-phen$admit_count_new
}

phen$exac<-rowSums(phen[,c("ED_visits","admit_count")])

phen$total_IgE_num<-as.numeric(phen$total_IgE_repeat)
phen$total_IgE_num[grep(">",phen$total_IgE_repeat)]<-5000
print(phen$total_IgE_num)
# add columns with categorical variables of BAL eos/neut data

phen<- phen %>%
  mutate(bal_eos_p_mt1 = factor(case_when(BAL_eos_p > 1 ~ 1 , #"high_eos"
                                          BAL_eos_p <= 1 ~ 0), # "low_eos"
                                levels = c(0, 1)),
         bal_neut_p_mt1 = factor(case_when(BAL_neut_p > 1 ~ 1 , #"high_neut"
                                          BAL_neut_p <= 1 ~ 0), # "low_neut"
                                levels = c(0, 1)),
         bal_neut_p_mt2 = factor(case_when(BAL_neut_p > 2 ~ 1 , #"high_neut"
                                          BAL_neut_p <= 2 ~ 0), # "low_neut"
                                levels = c(0, 1)),
         bal_neut_p_mt3 = factor(case_when(BAL_neut_p > 3 ~ 1 , #"high_neut"
                                          BAL_neut_p <= 3 ~ 0), # "low_neut"
                                levels = c(0, 1)),
         bal_neut_p_mt4 = factor(case_when(BAL_neut_p > 4 ~ 1 , #"high_neut"
                                          BAL_neut_p <= 4 ~ 0), # "low_neut"
                                levels = c(0, 1)),
         bal_neut_p_mt5 = factor(case_when(BAL_neut_p > 5 ~ 1 , #"high_neut"
                                          BAL_neut_p <= 5 ~ 0), # "low_neut"
                                levels = c(0, 1)),
         bal_neut_p_mt6 = factor(case_when(BAL_neut_p > 6 ~ 1 , #"high_neut"
                                          BAL_neut_p <= 6 ~ 0), # "low_neut"
                                levels = c(0, 1)),
         ED_visit_log = log(ED_visits+0.00001)%>%scale,
         admit_log = log(admit_count+0.00001)%>%scale,
         exac_log = log(exac+0.00001)%>%scale
         )

ungrouped_phen_summary<-phen%>%
  filter(!is.na(bal_eos_p_mt1))%>%
  summarise(age_visit_mean= mean(Age_at_visit),
            age_visit_sd = sd(Age_at_visit),
            sex_f_mean = mean(Sex=="Female"),
            sex_m_mean = mean(Sex=="Male"),
            sex_f_sum = sum(Sex=="Female"),
            sex_m_sum = sum(Sex=="Male"),
            race_a = sum(Race == "Asian"),
            race_black = sum(Race == "Black or African American"),
            race_mt1 = sum(Race == "More Than One Race"),
            race_unk = sum(Race == "Unknown / Not Reported"),
            race_wt = sum(Race == "White"),
            eth_L_rate = mean(Ethnicity == "Hispanic or Latino"),
            eth_L = sum(Ethnicity == "Hispanic or Latino"),
            eth_nL = sum(Ethnicity == "NOT Hispanic or Latino"),
            act_mean = mean(ACT_score, na.rm=TRUE),
            act_sd = sd(ACT_score, na.rm=TRUE),
            ED_visits_rate = mean(ED_visits>=1,na.rm=TRUE),
            ED_visits_yes = sum(ED_visits>=1,na.rm=TRUE),
            ED_visits_no = sum(ED_visits<1,na.rm=TRUE),
            admit_count_rate = mean(admit_count>=1,na.rm=TRUE),
            admit_count_yes = sum(admit_count>=1,na.rm=TRUE),
            admit_count_no = sum(admit_count<1,na.rm=TRUE),
            ED_visits_mean = mean(ED_visits,na.rm=TRUE),
            ED_visits_sd = sd(ED_visits,na.rm=TRUE),
            admit_count_mean = mean(admit_count,na.rm=TRUE),
            admit_count_sd = sd(admit_count,na.rm=TRUE),
            exac_mean = mean(exac,na.rm = TRUE),
            exac_sd = sd(exac,na.rm = TRUE),
            exac_rate = sum(exac>=1)/sum(exac>=0),
            exac_yes = sum(exac>=1),
            exac_non = sum(exac==0),
            ED_visit_log_mean = mean(ED_visit_log,na.rm = TRUE),
            ED_visit_log_sd = sd(ED_visit_log,na.rm = TRUE),
            admit_log_mean = mean(admit_log,na.rm = TRUE),
            admit_log_sd = sd(admit_log,na.rm = TRUE),
            exac_log_mean = mean(exac_log,na.rm=TRUE),
            exac_log_sd = sd(exac_log,na.rm=TRUE),
            bal_eos_p_count=sum(!is.na(BAL_eos_p)),
            bal_neut_p_count=sum(!is.na(BAL_neut_p)),
            bal_neut_p_mean=mean(BAL_neut_p,na.rm=TRUE),
            bal_neut_p_sd=sd(BAL_neut_p,na.rm=TRUE),
            bal_neut_p_log_mean=mean(BAL_neut_p_log, na.rm=TRUE),
            bal_neut_p_log_sd=sd(BAL_neut_p_log, na.rm=TRUE),
            bal_neut_p_mt1_yes=sum(BAL_neut_p>1,na.rm=TRUE),
            bal_neut_p_mt1_no=sum(BAL_neut_p<=1,na.rm=TRUE),
            
            bal_neut_p_mt2_yes=sum(BAL_neut_p>2,na.rm=TRUE),
            bal_neut_p_mt2_no=sum(BAL_neut_p<=2,na.rm=TRUE),
            
            bal_neut_p_mt3_yes=sum(BAL_neut_p>3,na.rm=TRUE),
            bal_neut_p_mt3_no=sum(BAL_neut_p<=3,na.rm=TRUE),
            
            bal_neut_p_mt4_yes=sum(BAL_neut_p>4,na.rm=TRUE),
            bal_neut_p_mt4_no=sum(BAL_neut_p<=4,na.rm=TRUE),
            
            bal_neut_p_mt5_yes=sum(BAL_neut_p>5,na.rm=TRUE),
            bal_neut_p_mt5_no=sum(BAL_neut_p<=5,na.rm=TRUE),
            
            bal_neut_p_mt6_yes=sum(BAL_neut_p>6,na.rm=TRUE),
            bal_neut_p_mt6_no=sum(BAL_neut_p<=6,na.rm=TRUE),
            FEV1_per_mean = mean(FEV1_percent,na.rm=TRUE),
            FEV1_per_sd = sd(FEV1_percent,na.rm=TRUE),
            FEV1_per_na = sum(!is.na(FEV1_percent)),
            med_ics_yes_p= mean(ICS=="Checked",na.rm=TRUE),
            med_ics_no_p= mean(ICS=="Unchecked",na.rm=TRUE),
            med_saba_yes_p= mean(SABA=="Checked",na.rm=TRUE),
            med_saba_no_p= mean(SABA=="Unchecked",na.rm=TRUE),
            med_ICS.saba_yes_p = mean(ICS.LABA=="Checked",na.rm=TRUE),
            med_ICS.saba_no_p = mean(ICS.LABA=="Unchecked",na.rm=TRUE),
            med_ltr_yes_p = mean(LTR=="Checked",na.rm=TRUE),
            med_ltr_no_p = mean(LTR=="Unchecked",na.rm=TRUE),
            med_ics_yes= sum(ICS=="Checked",na.rm=TRUE),
            med_ics_no= sum(ICS=="Unchecked",na.rm=TRUE),
            med_saba_yes= sum(SABA=="Checked",na.rm=TRUE),
            med_saba_no= sum(SABA=="Unchecked",na.rm=TRUE),
            med_ICS.saba_yes = sum(ICS.LABA=="Checked",na.rm=TRUE),
            med_ICS.saba_no = sum(ICS.LABA=="Unchecked",na.rm=TRUE),
            med_ltr_yes = sum(LTR=="Checked",na.rm=TRUE),
            med_ltr_no = sum(LTR=="Unchecked",na.rm=TRUE),
            total = sum(!is.na(bal_eos_p_mt1)))

ungrouped_phen_summary$bal_eos_p_mt1<-100

print(ungrouped_phen_summary)

# Define PFT variables
# Check which PFT variables exist in the dataset
existing_pft <- pft[pft %in% names(phen)]

# Summarize the PFT variables
if (length(existing_pft) > 0) {
  pft_summary <- phen[, existing_pft, drop = FALSE] %>%
    summarise(across(everything(), list(
      mean = ~mean(. , na.rm = TRUE),
      sd = ~sd(. , na.rm = TRUE)
    ), .names = "{.col}_{.fn}"))
  
  print(pft_summary)
} else {
  print("None of the specified PFT variables exist in the dataset.")
}

# Define the medication variables
med_vars <- c("ICS", "SABA", "ICS.LABA", "LTR")

# Check which medication variables exist in the dataset
existing_med_vars <- med_vars[med_vars %in% names(phen)]

# Summarize the medication usage
if (length(existing_med_vars) > 0) {
  med_summary <- phen[, existing_med_vars, drop = FALSE] %>%
    summarise(across(everything(), ~sum(. == "Checked", na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Medication", values_to = "Count")
  
  # Calculate the percentage of users
  total_patients <- nrow(phen)
  med_summary <- med_summary %>%
    mutate(Percentage = round((Count / total_patients) * 100, 2))
  
  print(med_summary)
} else {
  print("None of the specified medication variables exist in the dataset.")
}



# ==================
# grouped by BAL Eos % > 1%
# ==================


phen_summary<-phen%>%
  group_by(bal_eos_p_mt1, .add=TRUE)%>%
  filter(!is.na(bal_eos_p_mt1))%>%
  summarise(age_visit_mean= mean(Age_at_visit),
            age_visit_sd = sd(Age_at_visit),
            sex_f_mean = mean(Sex=="Female"),
            sex_m_mean = mean(Sex=="Male"),
            sex_f_sum = sum(Sex=="Female"),
            sex_m_sum = sum(Sex=="Male"),
            race_a = sum(Race == "Asian"),
            race_black = sum(Race == "Black or African American"),
            race_mt1 = sum(Race == "More Than One Race"),
            race_unk = sum(Race == "Unknown / Not Reported"),
            race_wt = sum(Race == "White"),
            eth_L_rate = mean(Ethnicity == "Hispanic or Latino"),
            eth_L = sum(Ethnicity == "Hispanic or Latino"),
            eth_nL = sum(Ethnicity == "NOT Hispanic or Latino"),
            act_mean = mean(ACT_score, na.rm=TRUE),
            act_sd = sd(ACT_score, na.rm=TRUE),
            ED_visits_rate = mean(ED_visits>=1,na.rm=TRUE),
            ED_visits_yes = sum(ED_visits>=1,na.rm=TRUE),
            ED_visits_no = sum(ED_visits<1,na.rm=TRUE),
            admit_count_rate = mean(admit_count>=1,na.rm=TRUE),
            admit_count_yes = sum(admit_count>=1,na.rm=TRUE),
            admit_count_no = sum(admit_count<1,na.rm=TRUE),
            ED_visits_mean = mean(ED_visits,na.rm=TRUE),
            ED_visits_sd = sd(ED_visits,na.rm=TRUE),
            admit_count_mean = mean(admit_count,na.rm=TRUE),
            admit_count_sd = sd(admit_count,na.rm=TRUE),
            exac_mean = mean(exac,na.rm = TRUE),
            exac_sd = sd(exac,na.rm = TRUE),
            exac_rate = sum(exac>=1)/sum(exac>=0),
            exac_yes = sum(exac>=1),
            exac_non = sum(exac==0),
            ED_visit_log_mean = mean(ED_visit_log,na.rm = TRUE),
            ED_visit_log_sd = sd(ED_visit_log,na.rm = TRUE),
            admit_log_mean = mean(admit_log,na.rm = TRUE),
            admit_log_sd = sd(admit_log,na.rm = TRUE),
            exac_log_mean = mean(exac_log,na.rm=TRUE),
            exac_log_sd = sd(exac_log,na.rm=TRUE),
            bal_eos_p_count=sum(!is.na(BAL_eos_p)),
            bal_neut_p_count=sum(!is.na(BAL_neut_p)),
            bal_neut_p_mean=mean(BAL_neut_p,na.rm=TRUE),
            bal_neut_p_sd=sd(BAL_neut_p,na.rm=TRUE),
            bal_neut_p_log_mean=mean(BAL_neut_p_log, na.rm=TRUE),
            bal_neut_p_log_sd=sd(BAL_neut_p_log, na.rm=TRUE),
            bal_neut_p_mt1_yes=sum(BAL_neut_p>1,na.rm=TRUE),
            bal_neut_p_mt1_no=sum(BAL_neut_p<=1,na.rm=TRUE),
            
            bal_neut_p_mt2_yes=sum(BAL_neut_p>2,na.rm=TRUE),
            bal_neut_p_mt2_no=sum(BAL_neut_p<=2,na.rm=TRUE),
            
            bal_neut_p_mt3_yes=sum(BAL_neut_p>3,na.rm=TRUE),
            bal_neut_p_mt3_no=sum(BAL_neut_p<=3,na.rm=TRUE),
            
            bal_neut_p_mt4_yes=sum(BAL_neut_p>4,na.rm=TRUE),
            bal_neut_p_mt4_no=sum(BAL_neut_p<=4,na.rm=TRUE),
            
            bal_neut_p_mt5_yes=sum(BAL_neut_p>5,na.rm=TRUE),
            bal_neut_p_mt5_no=sum(BAL_neut_p<=5,na.rm=TRUE),
            
            bal_neut_p_mt6_yes=sum(BAL_neut_p>6,na.rm=TRUE),
            bal_neut_p_mt6_no=sum(BAL_neut_p<=6,na.rm=TRUE),
            FEV1_per_mean = mean(FEV1_percent,na.rm=TRUE),
            FEV1_per_sd = sd(FEV1_percent,na.rm=TRUE),
            FEV1_per_na = sum(!is.na(FEV1_percent)),
            med_ics_yes_p= mean(ICS=="Checked",na.rm=TRUE),
            med_ics_no_p= mean(ICS=="Unchecked",na.rm=TRUE),
            med_saba_yes_p= mean(SABA=="Checked",na.rm=TRUE),
            med_saba_no_p= mean(SABA=="Unchecked",na.rm=TRUE),
            med_ICS.saba_yes_p = mean(ICS.LABA=="Checked",na.rm=TRUE),
            med_ICS.saba_no_p = mean(ICS.LABA=="Unchecked",na.rm=TRUE),
            med_ltr_yes_p = mean(LTR=="Checked",na.rm=TRUE),
            med_ltr_no_p = mean(LTR=="Unchecked",na.rm=TRUE),
            med_ics_yes= sum(ICS=="Checked",na.rm=TRUE),
            med_ics_no= sum(ICS=="Unchecked",na.rm=TRUE),
            med_saba_yes= sum(SABA=="Checked",na.rm=TRUE),
            med_saba_no= sum(SABA=="Unchecked",na.rm=TRUE),
            med_ICS.saba_yes = sum(ICS.LABA=="Checked",na.rm=TRUE),
            med_ICS.saba_no = sum(ICS.LABA=="Unchecked",na.rm=TRUE),
            med_ltr_yes = sum(LTR=="Checked",na.rm=TRUE),
            med_ltr_no = sum(LTR=="Unchecked",na.rm=TRUE),
            total_IgE_mean=mean(total_IgE_num,na.rm=TRUE),
            total_IgE_sd=sd(total_IgE_num,na.rm=TRUE),
            total = sum(!is.na(bal_eos_p_mt1)))

print(as.data.frame(phen_summary))

# write the summary
composite_summary<-rbind(ungrouped_phen_summary,phen_summary)
write.csv(composite_summary,"./reports/local_only/astham-phenotype~cellcount/phenotype_summary_composite.csv",row.names = FALSE)


# bal Neut % ~ BAL Eos% >1 not significantly different
boxplot(BAL_neut_p_log ~ bal_eos_p_mt1, data=phen)
boxplot(BAL_neut_p ~ bal_eos_p_mt1, data=phen)

# ---------------
# # of exacerbations, ED visits, Admission
# age, ACT score, FEV1 percentage, BAL_neut_p
# t-test
# ---------------

# List of variables to test
variables <- c("exac", "Age_at_visit","ACT_score", "FEV1_percent","ED_visits","admit_count", "ED_visit_log","admit_log","exac_log", "BAL_neut_p","BAL_neut_p_log", "total_IgE_num")

# Apply t.test to each variable
ttest_results <- lapply(variables, function(var) {
  formula <- as.formula(paste(var, "~ bal_eos_p_mt1"))
  t.test(formula, data = phen)
})

# Name the results for easy reference
names(ttest_results) <- variables

# Extract p-values from t-test results
p_values <- sapply(ttest_results, function(x) x$p.value)

# Print the p-values
print(p_values)


# -----------------------
# neutrophil stuff
# -----------------------



# ---------------
# sex, ethnicity
# any exacerbation (ED and admission), any ED visit, any admission
# any ICS/SABA/ICS-SABA/LTR
# Fisher's or Chi Sq
# ---------------

# Fisher's exact test 
variables <- list(c("sex_f_sum","sex_m_sum"),
                  c("eth_L","eth_nL"),
                  c("exac_yes","exac_non"),
                  c("ED_visits_yes","ED_visits_no"),
                  c("admit_count_yes","admit_count_no"),
                  c("bal_neut_p_mt1_yes","bal_neut_p_mt1_no"),
                  c("bal_neut_p_mt2_yes","bal_neut_p_mt2_no"),
                  c("bal_neut_p_mt3_yes","bal_neut_p_mt3_no"),
                  c("bal_neut_p_mt4_yes","bal_neut_p_mt4_no"),
                  c("bal_neut_p_mt5_yes","bal_neut_p_mt5_no"),
                  c("bal_neut_p_mt6_yes","bal_neut_p_mt6_no"),
                  c("med_ics_yes", "med_ics_no"),
                  c("med_saba_yes", "med_saba_no"),
                  c("med_ICS.saba_yes", "med_ICS.saba_no"),
                  c("med_ltr_yes", "med_ltr_no"))

# Apply Fisher's exact test to each pair
fisher_results <- lapply(variables, function(var_pair) {
  fisher.test(as.matrix(phen_summary[, var_pair]))
}) # fisher better because some samples are counts < 5

names(fisher_results)<-sapply(variables,function(x){x[1]})
# Extract p-values from t-test results
p_values <- sapply(fisher_results, function(x) x$p.value)

# Print the p-values
print(p_values)

#chi sq test 
variables <- list(c("exac_yes","exac_non"),
                  c("ED_visits_yes", "ED_visits_no"),
                  c("admit_count_yes", "admit_count_no"))

# Apply Fisher's exact test to each pair
chisq_results <- lapply(variables, function(var_pair) {
  chisq.test(as.matrix(phen_summary[, var_pair]))
}) # fisher better because some samples are counts < 5

names(chisq_results)<-sapply(variables,function(x){x[1]})

print(chisq_results)

#-------------
# Fisher's test for race
# -----------
# Create a matrix for the racial distribution
race_data <- matrix(c(2, 1, 1,
                      24, 13, 11,
                      9, 5, 4,
                      9, 3, 6,
                      26, 19, 7),
                    nrow = 5,
                    byrow = TRUE)

# Add row and column names
rownames(race_data) <- c("Asian", "Black", "More than 1", "Unknown", "White")
colnames(race_data) <- c("Composite", "BAL Eos %<=1%", "BAL Eos %>1%")

# Run Chi-Square Test
race_fisher_res <- fisher.test(race_data)

# Print results
print(race_fisher_res)


# ==================
# grouped by Ethnicity
# ==================
phen_summary<-phen%>%group_by(Ethnicity)%>%
  filter(!is.na(bal_eos_p_mt1))%>%
  summarise(age_visit_mean= mean(Age_at_visit),
            age_visit_sd = sd(Age_at_visit),
            sex_f_mean = mean(Sex=="Female"),
            sex_m_mean = mean(Sex=="Male"),
            sex_f_sum = sum(Sex=="Female"),
            sex_m_sum = sum(Sex=="Male"),
            race_a = sum(Race == "Asian"),
            race_black = sum(Race == "Black or African American"),
            race_mt1 = sum(Race == "More Than One Race"),
            race_unk = sum(Race == "Unknown / Not Reported"),
            race_wt = sum(Race == "White"),
            eth_L = sum(Ethnicity == "Hispanic or Latino"),
            eth_nL = sum(Ethnicity == "NOT Hispanic or Latino"),
            act_mean = mean(ACT_score, na.rm=TRUE),
            act_sd = sd(ACT_score, na.rm=TRUE),
            ED_visits_rate = mean(ED_visits>=1,na.rm=TRUE),
            ED_visits_yes = sum(ED_visits>=1,na.rm=TRUE),
            ED_visits_no = sum(ED_visits<1,na.rm=TRUE),
            admit_count_rate = mean(admit_count>=1,na.rm=TRUE),
            admit_count_yes = sum(admit_count>=1,na.rm=TRUE),
            admit_count_no = sum(admit_count<1,na.rm=TRUE),
            ED_visits_mean = mean(ED_visits,na.rm=TRUE),
            ED_visits_sd = sd(ED_visits,na.rm=TRUE),
            admit_count_mean = mean(admit_count,na.rm=TRUE),
            admit_count_sd = sd(admit_count,na.rm=TRUE),
            exac_mean = mean(exac,na.rm = TRUE),
            exac_sd = sd(exac,na.rm = TRUE),
            exac_rate = sum(exac>=1)/sum(exac>=0),
            exac_yes = sum(exac>=1),
            exac_non = sum(exac==0),
            bal_eos_p_count=sum(!is.na(BAL_eos_p)),
            bal_neut_p_count=sum(!is.na(BAL_neut_p)),
            bal_neut_p_mean=mean(BAL_neut_p,na.rm=TRUE),
            bal_neut_p_sd=sd(BAL_neut_p,na.rm=TRUE),
            bal_neut_p_log_mean=mean(BAL_neut_p_log, na.rm=TRUE),
            bal_neut_p_log_sd=sd(BAL_neut_p_log, na.rm=TRUE),
            bal_neut_p_mt1_yes=sum(BAL_neut_p>1,na.rm=TRUE),
            bal_neut_p_mt1_no=sum(BAL_neut_p<=1,na.rm=TRUE),
            
            bal_neut_p_mt2_yes=sum(BAL_neut_p>2,na.rm=TRUE),
            bal_neut_p_mt2_no=sum(BAL_neut_p<=2,na.rm=TRUE),
            
            bal_neut_p_mt3_yes=sum(BAL_neut_p>3,na.rm=TRUE),
            bal_neut_p_mt3_no=sum(BAL_neut_p<=3,na.rm=TRUE),
            
            bal_neut_p_mt4_yes=sum(BAL_neut_p>4,na.rm=TRUE),
            bal_neut_p_mt4_no=sum(BAL_neut_p<=4,na.rm=TRUE),
            
            bal_neut_p_mt5_yes=sum(BAL_neut_p>5,na.rm=TRUE),
            bal_neut_p_mt5_no=sum(BAL_neut_p<=5,na.rm=TRUE),
            
            bal_neut_p_mt6_yes=sum(BAL_neut_p>6,na.rm=TRUE),
            bal_neut_p_mt6_no=sum(BAL_neut_p<=6,na.rm=TRUE),
            FEV1_per_mean = mean(FEV1_percent,na.rm=TRUE),
            FEV1_per_sd = sd(FEV1_percent,na.rm=TRUE),
            FEV1_per_na = sum(!is.na(FEV1_percent)),
            med_ics_yes_p= mean(ICS=="Checked",na.rm=TRUE),
            med_ics_no_p= mean(ICS=="Unchecked",na.rm=TRUE),
            med_saba_yes_p= mean(SABA=="Checked",na.rm=TRUE),
            med_saba_no_p= mean(SABA=="Unchecked",na.rm=TRUE),
            med_ICS.saba_yes_p = mean(ICS.LABA=="Checked",na.rm=TRUE),
            med_ICS.saba_no_p = mean(ICS.LABA=="Unchecked",na.rm=TRUE),
            med_ltr_yes_p = mean(LTR=="Checked",na.rm=TRUE),
            med_ltr_no_p = mean(LTR=="Unchecked",na.rm=TRUE),
            med_ics_yes= sum(ICS=="Checked",na.rm=TRUE),
            med_ics_no= sum(ICS=="Unchecked",na.rm=TRUE),
            med_saba_yes= sum(SABA=="Checked",na.rm=TRUE),
            med_saba_no= sum(SABA=="Unchecked",na.rm=TRUE),
            med_ICS.saba_yes = sum(ICS.LABA=="Checked",na.rm=TRUE),
            med_ICS.saba_no = sum(ICS.LABA=="Unchecked",na.rm=TRUE),
            med_ltr_yes = sum(LTR=="Checked",na.rm=TRUE),
            med_ltr_no = sum(LTR=="Unchecked",na.rm=TRUE),
            bal_eos_p_mt1_yes=sum(bal_eos_p_mt1==1,na.rm=TRUE),
            bal_eos_p_mt1_no=sum(bal_eos_p_mt1==0,na.rm=TRUE),
            total = sum(!is.na(bal_eos_p_mt1)))

print(as.data.frame(phen_summary))

write.csv(phen_summary,"./reports/local_only/astham-phenotype~cellcount/phenotype_summary_latino-vs-non_latino.csv",row.names = FALSE)

# ---------------
# # of exacerbations, ED visits, Admission
# age, ACT score, FEV1 percentage, BAL_neut_p, total_IgE
# t-test
# ---------------

# List of variables to test
variables <- c("exac", "Age_at_visit","ACT_score", "FEV1_percent","ED_visits","admit_count", "BAL_neut_p","BAL_neut_p_log")

# Apply t.test to each variable
ttest_results <- lapply(variables, function(var) {
  formula <- as.formula(paste(var, "~ Ethnicity"))
  t.test(formula, data = phen)
})

# Name the results for easy reference
names(ttest_results) <- variables

# Extract p-values from t-test results
p_values <- sapply(ttest_results, function(x) x$p.value)

# Print the p-values
print(p_values)


# -----------------------
# neutrophil stuff
# -----------------------



# ---------------
# sex, ethnicity
# any exacerbation (ED and admission), any ED visit, any admission
# any ICS/SABA/ICS-SABA/LTR
# Fisher's or Chi Sq
# ---------------

# Fisher's exact test 
variables <- list(c("sex_f_sum","sex_m_sum"),
                  c("eth_L","eth_nL"),
                  c("exac_yes","exac_non"),
                  c("ED_visits_yes","ED_visits_no"),
                  c("admit_count_yes","admit_count_no"),
                  c("bal_eos_p_mt1_yes","bal_eos_p_mt1_no"),
                  c("bal_neut_p_mt1_yes","bal_neut_p_mt1_no"),
                  c("bal_neut_p_mt2_yes","bal_neut_p_mt2_no"),
                  c("bal_neut_p_mt3_yes","bal_neut_p_mt3_no"),
                  c("bal_neut_p_mt4_yes","bal_neut_p_mt4_no"),
                  c("bal_neut_p_mt5_yes","bal_neut_p_mt5_no"),
                  c("bal_neut_p_mt6_yes","bal_neut_p_mt6_no"),
                  c("med_ics_yes", "med_ics_no"),
                  c("med_saba_yes", "med_saba_no"),
                  c("med_ICS.saba_yes", "med_ICS.saba_no"),
                  c("med_ltr_yes", "med_ltr_no"))

# Apply Fisher's exact test to each pair
fisher_results <- lapply(variables, function(var_pair) {
  fisher.test(as.matrix(phen_summary[, var_pair]))
}) # fisher better because some samples are counts < 5

names(fisher_results)<-sapply(variables,function(x){x[1]})
# Extract p-values from t-test results
p_values <- sapply(fisher_results, function(x) x$p.value)

# Print the p-values
print(p_values)

#chi sq test 
variables <- list(c("exac_yes","exac_non"),
                  c("ED_visits_yes", "ED_visits_no"),
                  c("admit_count_yes", "admit_count_no"))

# Apply Fisher's exact test to each pair
chisq_results <- lapply(variables, function(var_pair) {
  chisq.test(as.matrix(phen_summary[, var_pair]))
}) # fisher better because some samples are counts < 5

names(chisq_results)<-sapply(variables,function(x){x[1]})

print(chisq_results)

#-------------
# Fisher's test for race
# -----------
# Create a matrix for the racial distribution
race_data <- matrix(c(2, 1, 1, 
                      25, 14, 11, 
                      9, 5, 4, 
                      10, 4, 5, 
                      27, 20, 7), 
                    nrow = 5, byrow = TRUE)

# Add row and column names
rownames(race_data) <- c("Asian", "Black", "More than 1", "Unknown", "White")
colnames(race_data) <- c("Composite", "BAL Eos %<=1%", "BAL Eos %>1%")

# Run Chi-Square Test
race_fisher_res <- fisher.test(race_data)

# Print results
print(race_fisher_res)


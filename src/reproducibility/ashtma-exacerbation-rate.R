# find how many people had at least one ED visit or admission in the cohort
library(dplyr)
phen<-read.csv("./resources/processed_data/BiomarkersOfAsthma_original.csv")
id<-read.csv("./resources/processed_data/nb_studyID_sampleID_batch.csv")
id<-id%>%filter(grepl("^B",SampleID))%>%pull(ID)
phen<-filter(phen,ID%in%id)%>%select(ED_visits,admit_count)
exacerbation_rate<-mean(apply(phen,1,function(d)sum(d,na.rm = TRUE))>0)

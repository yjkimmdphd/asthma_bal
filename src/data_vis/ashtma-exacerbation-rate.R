# find how many people had at least one ED visit or admission in the cohort
library(dplyr)
countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
rownames(counts)<-counts[,"SampleID"]

phen<-read.csv("./resources/processed_data/BiomarkersOfAsthma_original.csv")
id<-read.csv("./resources/processed_data/nb_studyID_sampleID_batch.csv")
phen<-left_join(id,phen,by=join_by(ID))
phen<-phen%>%filter(grepl("^B",SampleID),!is.na(BAL_WBC))
dim(phen)
exacerbation<-phen%>%select(ED_visits,admit_count)
exacerbation_rate<-mean(apply(exacerbation,1,function(d)sum(d,na.rm = TRUE))>0)

median(phen$Age_at_visit)
sd(phen$Age_at_visit)

mean(phen$ACT_score, na.rm=TRUE)
sd(phen$ACT_score, na.rm=TRUE)

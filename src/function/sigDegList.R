
library(tidyverse)


# function to load sig.deg table
sigDeg<- function(output.folder,file)read.csv(
  file.path(getwd(),"output",output.folder,file)
)

# summarize sig.deg
sigDegList<-function(r1,r2){
  c<-sig.deg%>%filter(results%in%c(r1,r2))%>%
    mutate(fc=case_when(log2FoldChange>0 ~ "up",log2FoldChange<0~"down"))%>%
    select(genes,fc,results)%>%ungroup%>%
    group_by(genes)
  d<-c%>%summarize(n=n())
  sdl.congruence<-left_join(c,d,by="genes")%>%ungroup%>%
    select(genes,results,fc,n)%>%
    spread(key=results,value=fc)
  return(sdl.congruence)
}


######
# check if DEG fold changes of significant DEG (nasal_bronchial_BAL_rnaseq_v2.R) between res33 and res34 are congruent
######
sig.deg<-sigDeg("DEG_2023-09-30","sig_deg_2023-10-03_rnaseq_v2.csv")
sc<-sigDegList("res33","res34")
sc<-mutate(sc,congruent=ifelse(res33==res34,1,0))
table(sc[,c("res33","res34")]) #> the significant DEG found in both results have fold change in the same direction
table(sc[which(sc$congruent>0),c("res33","res34")]) #> same as above. table function automatically removes rows with only one value


######
# check if DEG fold changes of significant DEG between res25 and res26 are congruent
######
sig.deg<-sigDeg("DEG_2023-09-30","sig_deg_2023-10-03_rnaseq_v2.csv")
sc<-sigDegList("res25","res26")
sc<-mutate(sc,congruent=ifelse(res25==res26,1,0))
table(sc[,c("res25","res26")]) #> the significant DEG found in both results have fold change in the same direction
table(sc[which(sc$congruent>0),c("res25","res26")]) #> same as above. table function automatically removes rows with only one value

######
# check if DEG fold changes of significant DEG (nasal_bronchial_BAL_rnaseq_v3.R) between res33 and res34 are congruent
######
sig.deg<-sigDeg("DEG_2023-10-01","sig_deg_2023-10-02-rnaseq-v3.csv")
sc<-sigDegList("res33","res34")
sc<-mutate(sc,congruent=ifelse(res33==res34,1,0))
table(sc[,c("res33","res34")]) #> the significant DEG found in both results have fold change in the same direction
table(sc[which(sc$congruent>0),c("res33","res34")]) #> same as above. table function automatically removes rows with only one value


######
# check if DEG fold changes of significant DEG between res25 and res26 are congruent
######
sig.deg<-sigDeg("DEG_2023-10-01","sig_deg_2023-10-02-rnaseq-v3.csv")
sc<-sigDegList("res25","res26")
sc<-mutate(sc,congruent=ifelse(res25==res26,1,0))
table(sc[,c("res25","res26")]) #> the significant DEG found in both results have fold change in the same direction
table(sc[which(sc$congruent>0),c("res25","res26")]) #> same as above. table function automatically removes rows with only one value

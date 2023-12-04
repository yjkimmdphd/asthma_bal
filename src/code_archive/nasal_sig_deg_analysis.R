library(tidyverse)
fp1<-file.path(getwd(),"output","DEG_2023-09-30","sig_deg_2023-10-03_rnaseq_v2.csv")
fp2<-file.path(getwd(),"output","DEG_2023-10-01","sig_deg_2023-10-02-rnaseq-v3.csv")

file.exists(fp1)
file.exists(fp2)
v2<-read.csv(fp1)
v3<-read.csv(fp2)

v2g<-unique(v2$genes)
v2r<-unique(v2$results)
v2m<-sapply(v2r,function(v){
  d<-v2%>%filter(results==v)%>%select(genes)%>%unlist%>%as.vector
  v2g%in%d
})
rownames(v2m)<-v2g
v2summ<-data.frame(genes=v2g,v2m)


v3g<-unique(v3$genes)
v3r<-unique(v3$results)
v3m<-sapply(v3r,function(v){
  d<-v3%>%filter(results==v)%>%select(genes)%>%unlist%>%as.vector
  v3g%in%d
})
rownames(v3m)<-v3g
v3summ<-data.frame(genes=v3g,v3m)

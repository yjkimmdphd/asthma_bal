##############
#MDS by groups for nasa_bronchial_BAL_rnaseq.R
##############

# color code covariates 
col<-data.frame(matrix(nrow=45,ncol=6))
colnames(col)<-c("eos","batch","age","sex","race","AR")
for(i in 1:6) {col[,i]<-rep("black",45)}

col[which(p.counts$BAL_Eos_perc>=3),"eos"]<-"red" # BAL Eos >= 3% are red
col[which(p.counts$BAL_Eos_perc<1),"eos"]<-"green" # BAL Eos = 0% are green
col[which(p.counts$BAL_Eos_perc>10),"eos"]<-"blue"# BAL Eos >10% are blue

col[which(p.counts$Batch=="batch1"),"batch"]<-"red" # batch1 are red
col[which(p.counts$Batch=="batch2"),"batch"]<-"green" # batch2 are green
col[which(p.counts$Batch=="batch3"),"batch"]<-"blue"# batch3 are blue. Batch4 are black

col[which(p.counts$Age<=7),"age"]<-"red" # age less than or equal to 7 are red
col[which(p.counts$Age>7 & p.counts$Age<=13),"age"]<-"green" # age 8yo to 13yo are green
col[which(p.counts$Age>13 & p.counts$Age<=19),"age"]<-"blue"# age 14yo to 19yo are blue. age >19 are black

col[which(p.counts$Gender=="Male"),"sex"]<-"blue"
col[which(p.counts$Gender=="Female"),"sex"]<-"red"

Race_corrected<-p.counts$Race_corrected%>%unique
colors<-c("red","blue","green","black","orange","purple")

race.color<-data.frame(Race_corrected,colors)
race.color<-full_join(p.counts,race.color, by="Race_corrected")
col$race<-race.color$colors
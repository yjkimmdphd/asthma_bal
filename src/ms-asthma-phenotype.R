############################################
# asthma phenotype data exploratory analysis
############################################
library(tidyverse)
# pathway for BAL data
data.folder="/sc/arion/projects/asthma-allergy/MS_asthma/Young-Jin/inputs/MS_asthma_data/" 
data.file ="BiomarkersOfAsthmaSt_DATA_LABELS_2023-08-24_revised_colnames_filterd_for_nasal_biomarker_study.csv"
filtered.col = "biomarker_asthma_columnsOfInterest_2023-08-29.csv"

# check if the file exists
data.pathway= 
  if(file.exists(paste0(data.folder,data.file))){
    print(paste0(data.folder,data.file))
  }
print(data.pathway)

filtered.col.pathway = 
  if(file.exists(paste0(data.folder,filtered.col))){
    print(paste0(data.folder,filtered.col))
  }
print(filtered.col.pathway)

# load the files into a tibble
asthma.pheno = read_csv(data.pathway)
filtered.col = read_csv(filtered.col.pathway,col_names=FALSE,col_types = "c")%>%as.data.frame
print(filtered.col)
filtered.col<-filtered.col[,1]%>%as.vector

# choose only the phenotype col of interest
apf<-asthma.pheno%>%select(any_of(filtered.col))

#check which data cols has Eos related stuff
Eos.col<-grep("Eos",colnames(asthma.pheno))

# select BAL and serum cell count data, APC = asthma phenotype cell count
BAL.serum.col<-c(grep("BAL", colnames(asthma.pheno)),grep("serum", colnames(asthma.pheno)))
apc<-asthma.pheno[,c(1,BAL.serum.col)]

# find columns with character data and turn into factors
col.class<-which("character"==sapply(apf,class))

# export filtered asthma phenotype data
# write_csv(apf, paste0(getwd(),"/outputs/asthma-phenotype-filtered.csv"))

# find which columns have <0.35, <0.1 <0.14 etc
col1<-grep("<",apf)
df.apf<-apf[,col1]%>%as.data.frame
cola<-data.frame(matrix(0,nrow=nrow(df.apf),ncol=ncol(df.apf)))
for(i in 1:ncol(df.apf)){
    cola[,i]<-grepl("<",df.apf[,i])
}
colnames(cola)<-colnames(df.apf)

# find which samples had BAL Eos% < 3%
# normal is 2% or less https://www.ncbi.nlm.nih.gov/books/NBK470600/
eos.col<-grep("Eos",colnames(apf))
head(apf[,eos.col])
bal.eos.filter<-as.vector(apf[,eos.col[1]]<3)
apf$eos.lessthan3<-bal.eos.filter

# make dataframe called 'apf.eos' that shows subject assignment ID and BAL eos perc less than 3
eos.colname<-c(colnames(apf[,eos.col]),"eos.lessthan3")
apf.eos<-as.data.frame(apf[,c("subject_assgn",eos.colname)])
write.csv(apf.eos, row.names = FALSE,"/sc/arion/projects/asthma-allergy/MS_asthma/Young-Jin/inputs/MS_asthma_data/MSasthmaBALEosPerc.csv")

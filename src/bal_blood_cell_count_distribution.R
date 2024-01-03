# exploration of cell count distribution 
# comparing distribution of raw cell counts and log-transformed counts 
library(tidyverse)
phenotype<-read.csv("./resources/processed_data/Nasal_Biomarkers_BAL_transformed.csv")
a<-colnames(phenotype)[2:21]
a<-phenotype[,a]
a<-a[,order(names(a))]
xlab<-colnames(a)
index<-!sapply(a,is.na)
par(mfcol = c(2, 5))
for(i in 1:10){
  k<-a[[i]][index[,i]]
    hist(k,
       main = xlab[i], 
       breaks=seq(min(k),max(k),length.out=25),
       xlab = xlab[i], ylab = "Frequency", col = "#69b3a2", border = "black")
}
par(mfcol = c(2, 5))
for(i in 11:20){
  k<-a[[i]][index[,i]]
  hist(k,
       main = xlab[i], 
       breaks=seq(min(k),max(k),length.out=25),
       xlab = xlab[i], ylab = "Frequency", col = "#69b3a2", border = "black")
}

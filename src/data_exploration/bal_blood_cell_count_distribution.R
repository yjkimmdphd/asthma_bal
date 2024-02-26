# exploration of cell count distribution 
# comparing distribution of raw cell counts and log-transformed counts 
library(tidyverse)
phenotype<-read.csv("./resources/processed_data/Nasal_Biomarkers_BAL_transformed.csv")
a<-colnames(phenotype)[2:21]
a<-phenotype[,a]
a<-a[,order(names(a))]
xlab<-colnames(a)
index<-!sapply(a,is.na)
log_columns <- grep("_log$", names(a), value = TRUE)
standardized_data <- a %>%
  mutate_at(log_columns, scale)

### hist of standardized data

par(mfrow = c(4,5))

#### 1. BAL cell counts
for(i in seq(from = 1, to = 10, by = 2)){
  k<-standardized_data[[i]][index[,i]]
  hist(k,
       main = xlab[i], 
       breaks=seq(min(k),max(k),length.out=25),
       xlab = xlab[i], ylab = "Frequency", col = "#69b3a2", border = "black")
}

for(i in seq(from = 2, to = 10, by = 2)){
  k<-standardized_data[[i]][index[,i]]
  hist(k,
       main = xlab[i], 
       breaks=seq(min(k),max(k),length.out=25),
       xlab = paste(xlab[i],"standardized",sep=" "), ylab = "Frequency", col = "#69b3a2", border = "black")
}

#### 2. Bld cell counts 
for(i in seq(from = 11, to = 20, by = 2)){
  k<-standardized_data[[i]][index[,i]]
  hist(k,
       main = xlab[i], 
       breaks=seq(min(k),max(k),length.out=25),
       xlab = xlab[i], ylab = "Frequency", col = "#69b3a2", border = "black")
}

for(i in seq(from = 12, to = 20, by = 2)){
  k<-standardized_data[[i]][index[,i]]
  hist(k,
       main = xlab[i], 
       breaks=seq(min(k),max(k),length.out=25),
       xlab = paste(xlab[i],"standardized",sep=" "), ylab = "Frequency", col = "#69b3a2", border = "black")
}

### hist of non-standardized data

#### 1. BAL cell counts
for(i in seq(from = 1, to = 10, by = 2)){
  k<-a[[i]][index[,i]]
  hist(k,
       main = xlab[i], 
       breaks=seq(min(k),max(k),length.out=25),
       xlab = xlab[i], ylab = "Frequency", col = "#69b3a2", border = "black")
}

for(i in seq(from = 2, to = 10, by = 2)){
  k<-a[[i]][index[,i]]
  hist(k,
       main = xlab[i], 
       breaks=seq(min(k),max(k),length.out=25),
       xlab = paste(xlab[i],"standardized",sep=" "), ylab = "Frequency", col = "#69b3a2", border = "black")
}

#### 2. Bld cell counts 
for(i in seq(from = 11, to = 20, by = 2)){
  k<-a[[i]][index[,i]]
  hist(k,
       main = xlab[i], 
       breaks=seq(min(k),max(k),length.out=25),
       xlab = xlab[i], ylab = "Frequency", col = "#69b3a2", border = "black")
}

for(i in seq(from = 12, to = 20, by = 2)){
  k<-a[[i]][index[,i]]
  hist(k,
       main = xlab[i], 
       breaks=seq(min(k),max(k),length.out=25),
       xlab = paste(xlab[i],"standardized",sep=" "), ylab = "Frequency", col = "#69b3a2", border = "black")
}

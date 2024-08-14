################
# data prep
# raw count: "./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt"
# goal: normalizes the gene counts with TMM method 
################


library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)

countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}

sampleIDs<-colnames(counts)
b_n_sampleIDs<-!grepl("^T",sampleIDs)
sampleIDs<-sampleIDs[b_n_sampleIDs]
counts<-counts[sampleIDs]
row.names(counts)<-counts[,"SampleID"]
counts<-counts[,-1]

sampleIDs<-colnames(counts)

sID_list<-list()
sID_list[[1]]<-sampleIDs[grep("^B",colnames(counts))]
sID_list[[2]]<-sampleIDs[grep("^N",colnames(counts))]
names(sID_list)<-c("bronch","nasal")
print(sID_list)

bronch.counts<-counts[,c(sID_list$bronch)]

# prepare normalized gene count table 

x<-bronch.counts
L<-mean(colSums(x))*1e-6 # mean library size
M<-median(colSums(x))*1e-6 # median library size
lcpm.cutoff <- log2(10/M + 2/L) # lcpm cutoff for filtering genes with very low counts

### normalize counts with TMM
norm.factor<-calcNormFactors(x, method = "TMM")
sample.size<-length(colnames(x))
for(i in 1:sample.size){
  x[,i]<-x[,i]/norm.factor[i]
}

normalized.count<-cbind(name=rownames(counts),x)

normalized_count_table_path<-"./resources/processed_data/bronch_batch12346_normalized_ct.txt"
write.table(normalized.count,normalized_count_table_path,row.names = FALSE,sep = "\t")

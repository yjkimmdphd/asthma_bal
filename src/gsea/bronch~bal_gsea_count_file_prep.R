################
# gsea data prep
################


library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)

countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
bronch.samples<-grepl("^B",colnames(counts))
bronch.samples<-which(bronch.samples==TRUE)
bronch.counts<-counts[,c(1,bronch.samples)]
colnames(bronch.counts)[bronch.samples]<-substr(colnames(bronch.counts)[bronch.samples],1,4)
head(bronch.counts)
counts.ID<-colnames(bronch.counts)
counts<-bronch.counts[,-1]
rownames(counts)<-bronch.counts$SampleID
# prepare normalized gene count table 

x<-counts
L<-mean(colSums(x))*1e-6 # mean library size
M<-median(colSums(x))*1e-6 # median library size
lcpm.cutoff <- log2(10/M + 2/L) # lcpm cutoff for filtering genes with very low counts

### normalize counts with TMM
norm.factor<-calcNormFactors(x, method = "TMM")
sample.size<-length(colnames(x))
for(i in 1:sample.size){
  x[,i]<-x[,i]/norm.factor[i]
}
### calculate lcpm based on TMM normalized counts 
lcpm.x<-cpm(x,log=TRUE)

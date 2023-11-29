##
# RNA seq analysis of nasal/bronchial data with edgeR and limma
# YJK local machine version
## 
library(dplyr)
library(limma)
library(Glimma)
library(edgeR)

############
# count data and phenotype
############

# reading count data for batch1,2,3,4 and saving as 'counts'
file<-c("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/data/original_data/MS_asthma/batch1234_readcount_matrix_allsamples.afterQC.txt")
file.exists(file)
counts<-read.delim(file[1])
counts.ID<-colnames(counts)

# reading count data for batch 5 as 'counts.b5'
file.b5<-c("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/data/original_data/MS_asthma/AsthmaNasal_RNAseq.batch5.geneID_readcount_rmscaffold_highqualitysamples.txt")
file.exists(file.b5)
counts.b5<-read.delim(file.b5[1])
counds.b5ID<-colnames(counts.b5)

# asthma biomarker phenotype file, nasal, saved in  'phenotype'
filename2<-("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/data/revised_data/asthma-phenotype-filtered-revised.csv")
file.exists(filename2)
phenotype<-read.csv(filename2)

# get batch information
filename3<-file.path(getwd(),"data/original_data/MS_asthma/MS_asthma_phenotype.batch1234.txt")
file.exists(filename3)
batch.info<-read.delim(filename3)

# load phenotype data by sourcing the following code 
source("./codes/phenotype_cleanup_nasal_bronchial_BAL_rnaseq.R")

# select just the nasal RNAseq counts
b<-counts.ID%in%nasal.ID
x<-counts[,b]
genes<-counts$SampleID
rownames(x)<-genes

# convert raw counts to CPM
# average coverage is ~51 million
cpm0 <-cpm(x)
lcpm0<-cpm(x,log=TRUE)
L<-mean(colSums(x))*1e-6 # mean library size
M<-median(colSums(x))*1e-6 # median library size


###########################################
# Normalising gene expression distributions
###########################################

#x3 is the TMM normalized count
# normalize counts with TMM
norm.factor<-calcNormFactors(x, method = "TMM")

x3<-x
sample.size<-length(colnames(x3))
for(i in 1:sample.size){
  x3[,i]<-x3[,i]/norm.factor[i]
}

# calculate lcpm based on TMM normalized counts 
lcpm.x3<-cpm(x3,log=TRUE)


##############
#MDS by groups
##############

# color coding based on each factor levels of sex, age, BAL Eos%, and batch sourced from the following code:
source("./codes/mds_color_code_nasal_bronchial_BAL_rnaseq.R")

# plot the MDS
par(mfrow=c(2,2))
mds.eos<-plotMDS(lcpm.x3,  col= col$eos, labels=p.counts$BAL_Eos_perc, main="Eos %")
mds.batch<-plotMDS(lcpm.x3,  col= col$batch, labels=p.counts$Batch, main="batch")
mds.age<-plotMDS(lcpm.x3,  col= col$age, labels=p.counts$Age, main = "age")
mds.race<-plotMDS(lcpm.x3,  col= col$race, labels=p.counts$Race_corrected, main = "race")

# significant batch effect



#########################################
# Removing genes that are lowly expressed
#########################################

# checking how many genes have 0 count across all samples
table(rowSums(x==0)==45)
# about 4% of the genes have 0 counts across all samples 

# setting a lcpm cutoff for filtering genes with very low counts
lcpm.cutoff <- log2(10/M + 2/L)
dropCutoff<-function(cutoff){
  which(apply(lcpm.x3, 1, max) < cutoff)
}
drop <-dropCutoff(lcpm.cutoff) 
drop2<-dropCutoff(0)
dim(x3[-drop,])
dim(x3[-drop2,])


########################
#DEG with limma and voom
########################

##### model with BAL Eos percentage 
bal.eos.log<-p.counts$BAL_Eos_perc_log
batch<-factor(p.counts$Batch, levels=unique(batch.info.BAL$Batch))
mm2<-model.matrix(~0+bal.eos.log + batch)
mm2.1<-model.matrix(~bal.eos.log + batch)
# check voom to see if different filtering needed
y<-voom(x3[-drop,],mm2.1,plot=T)

# may be needs more stringent filtering after define new cutoff point 
y2<-voom(x3[-drop2,], mm2.1, plot=T)

fit<-lmFit(y2,mm2.1)
tmp<-contrasts.fit(fit, coef=2)
tmp<-eBayes(tmp)
top.table<-topTable(tmp,coef=1,sort.by="p", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
# write.csv(top.table,file.path(getwd(),"output/",paste0("toptable_balEosPerc_",Sys.Date(),".csv")))

##### model with absolute BAL Eos count 
bal.eos.ct<-p.counts$BAL_Eos_ct_log
batch<-factor(p.counts$Batch, levels=unique(batch.info.BAL$Batch))
mm3<-model.matrix(~bal.eos.ct + batch)

# filter using cutoff 0 
drop2<-dropCutoff(0)
y3<-voom(x3[-drop2,], mm3, plot=T)

fit<-lmFit(y3,mm3)
tmp<-contrasts.fit(fit, coef=2)
tmp<-eBayes(tmp)
top.table<-topTable(tmp,coef=1,sort.by="p", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
# write.csv(top.table,file.path(getwd(),"output/",paste0("toptable_balEos_count_",Sys.Date(),".csv")))


##### model with absolute BAL WBC count 
bal.wbc<-p.counts$BAL_WBC_log
batch<-factor(p.counts$Batch, levels=unique(batch.info.BAL$Batch))
mm4<-model.matrix(~bal.wbc + batch)

# filter using cutoff 0 
drop2<-dropCutoff(0)
y4<-voom(x3[-drop2,], mm4, plot=T)

fit<-lmFit(y4,mm4)
tmp<-contrasts.fit(fit, coef=2)
tmp<-eBayes(tmp)
top.table<-topTable(tmp,coef=1,sort.by="p", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
# write.csv(top.table,file.path(getwd(),"output/",paste0("toptable_bal_WBC_",Sys.Date(),".csv")))


########################
#DEG with DEseq2
########################
library(DESeq2)

#### analyze based on log transformed BAL Eos %
dds<-DESeqDataSetFromMatrix(countData = x,colData=p.counts, design=~BAL_Eos_perc_log+Batch)
dds<-DESeq(dds)
res<-results(dds, name=c("BAL_Eos_perc_log"))
res <- res[order(res$padj),]
head(res)

# signficant results with padj<0.05
res.sig<-res[which(res$padj<0.05),]
res.sig
# write.csv(res.sig,file.path(getwd(),"output/",paste0("Deseq2_balEosPerc_",Sys.Date(),".csv")))

#### analyze based on log transformed BAL Eos Count
dds2<-DESeqDataSetFromMatrix(countData = x,colData=p.counts, design=~BAL_Eos_ct_log+Batch)
dds2<-DESeq(dds2)
res2<-results(dds2, name=c("BAL_Eos_ct_log"))
res2 <- res2[order(res2$padj),]
head(res2)

# signficant results with padj<0.05
res.sig2<-res2[which(res2$padj<0.05),]
res.sig2

# write.csv(res.sig2,file.path(getwd(),"output/",paste0("Deseq2_balEos_ct_",Sys.Date(),".csv")))

#### analyze based on log transformed BAL WBC Count
dds3<-DESeqDataSetFromMatrix(countData = x,colData=p.counts, design=~BAL_WBC_log+Batch)
dds3<-DESeq(dds3)
res3<-results(dds3, name=c("BAL_WBC_log"))
res3 <- res3[order(res3$padj),]
head(res3)

# signficant results with padj<0.05
res.sig3<-res3[which(res3$padj<0.05),]
res.sig3

# write.csv(res.sig2,file.path(getwd(),"output/",paste0("Deseq2_bal_WBC_ct_",Sys.Date(),".csv")))

# Define a function to make a basic volcano plot
asthmaVolcanoPlot<-function(res){
  # res is the result of the DEG with DESeq2
  # Add colored points: cyan if padj<0.01, red if log2FC>1 and padj<0.05)
  par(mfrow=c(1,1))
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
  with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="cyan"))
  with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  }
asthmaVolcanoPlot(res)


#######
# DEG modeled based on serum counts 


#### analyze based on serum Eos count 
# some serum Eos counts are NA because the serum Eos % is NA. 
# make phenotype and count table removign NA
hasna<-which(p.counts$serum_Eos_log%>%is.na)
nona<-p.counts[-hasna,"SampleID"]
p.counts.rmna<-p.counts[-hasna,]
x.rmna<-x[,nona]

p.counts.rmna # some 
x.rmna # count table without 

dds4<-DESeqDataSetFromMatrix(countData = x.rmna,colData=p.counts.rmna, design=~serum_Eos_log+Batch)
dds4<-DESeq(dds4)
res4<-results(dds4, name=c("serum_Eos_log"))
res4 <- res4[order(res4$padj),]
head(res4)

# signficant results with padj<0.05
res.sig4<-res4[which(res4$padj<0.05),]
res.sig4
res.sig4%>%nrow

# write.csv(res.sig2,file.path(getwd(),"output/",paste0("Deseq2_serEos_ct_",Sys.Date(),".csv")))


par(mfrow=c(1,2))

# Make a basic volcano plot
with(res4, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res4, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(res4, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

vsdata4 <- vst(dds4, blind=FALSE)
plotPCA(vsdata4, intgroup="Batch") 



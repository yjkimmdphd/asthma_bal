## 
# Bronchial RNAseq analysis with new samples:
# additional samples run in the batch 6
# alignment based on new reference genome GRCh38
# phenotype table was updated. loading was simplified
##

library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)


bphen_path<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")
bphen<-read.csv(bphen_path)

bphen<-bphen%>%filter(asthma_phen_ACT.score>0)
bphen<-mutate_at(bphen,vars(asthma_phen_ACT.score),scale) # scale ACT score for downstream analysis 
bphen<-bphen[c("SampleID","asthma_phen_ACT.score", "Batch")]
bphen$Batch<-factor(bphen$Batch,levels=unique(bphen$Batch))

# DESeq2 takes un-normalized gene counts as input 
countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
bronch.samples<-grepl("^B",colnames(counts))
bronch.samples<-which(bronch.samples==TRUE)

bronch.counts<-counts[,c(1,bronch.samples)]
colnames(bronch.counts)[bronch.samples]<-substr(colnames(bronch.counts)[bronch.samples],1,4)
head(bronch.counts)
counts.ID<-colnames(bronch.counts)

# compare the sampleID of the count matrix and the phenotype data
length(counts.ID[-1])-sum(counts.ID[-1]%in%bphen$SampleID)

source("./src/function/deg_custom_functions.R")

# filtering counts table to remove low expressed genes

## select RNAseq counts
phen<-bphen
id<-phen$SampleID
cols<-colnames(bronch.counts)%in%id
ct<-bronch.counts[,cols] # First column is actually gene name 
genes<-bronch.counts$SampleID
rownames(ct)<-genes

## Filter counts 
c2<-filter_low_expressed_genes_method2(ct,7)

# compare the sampleID of the count matrix and the phenotype data
length(colnames(c2))==sum(colnames(c2)%in%phen$SampleID)

# make colData for DESeq2
df<-list()
df[[1]]<-phen

# run the DEG for continuous predictors
deg.design<-paste("~","asthma_phen_ACT.score","+ Batch")
print(deg.design)
ct<-rowgenes_counttable(ct,c2) # low bcounts will be filtered using method 2: use TMM normalized lcpm as a cutoff point

print(deg.design)

# make a list of count table. Each element is a count table with samples for each of the experimental design. 
df.input<-df
count.table<-lapply(df.input,function(df){d<-df; ct<-ct[,colnames(ct)%in%d$SampleID]; return(ct)}) 

dds<-vector("list",length=length(df))
res<-vector("list",length=length(df))
res.sig<-vector("list",length=length(df))

names(res)<-deg.design
names(res.sig)<-deg.design

for(i in 1:length(df)){
  dds[[i]]<-run_deseq2_DEG_analysis(count.table[[i]], df.input[[i]], deg.design[i],deg.design[i])
  res[[i]]<-get_DEG_results(dds[[i]], source.cell.log[i])
  res.sig[[i]]<-res[[i]][which(res[[i]]$padj<=0.05),]
  head(res.sig[[i]])
  
}

res<-results(dds[[1]], name="asthma_phen_ACT.score",pAdjustMethod = "fdr") # use source.cell.log element for resultname
res_list<-list()
res_sig_list<-list()
res_list[[1]] <- res[order(res$log2FoldChange,decreasing = TRUE),]
res_sig_list[[1]]<-res[which(res$padj<=0.05),]


## writing the significant and all results 
deg.folder<-paste("deg",Sys.Date(),sep="_")
deg.dir<-file.path("./reports",deg.folder)
if(!dir.exists(deg.dir)){
  dir.create(deg.dir)
}

if(dir.exists(deg.dir)){
  for(i in 1:length(df)){
    a<-res_sig_list[[i]]
    b<-res_list[[i]]
    write.csv(a,row.names=TRUE,file.path(deg.dir,paste("deg","bronch","ACT-score+batch12346","res_sig",i,deg.design[[i]],Sys.Date(),".csv",sep="_")))
    write.csv(b,row.names=TRUE,file.path(deg.dir,paste("deg","bronch","ACT-score+batch12346","res_all",i,deg.design[[i]],Sys.Date(),".csv",sep="_")))
  }
}

## summarize the data input 

generate_DEG_input_summary_table<-function(original_ct,filtered_ct,dds,res,des){
  filter_method<-"TMM normalized LCPM cutoff"
  n_filtered_genes<-paste("analyzed n_genes:", nrow(filtered_ct),",","filtered n_genes:",nrow(original_ct)-nrow(filtered_ct))
  samples<-sapply(dds, function(d){colData(d)$SampleID%>%paste(collapse = ",")})
  dds<-paste("dds",1:length(dds),sep="")
  results<-paste("res",1:length(res),sep="")
  design<-des
  df<-data.frame(dds=dds,results=results,design=design,samples=samples,filter_method=filter_method,n_filtered_genes=n_filtered_genes)
  return(df)
}


if(dir.exists(deg.dir)){
  a<-generate_DEG_input_summary_table(bronch.counts[,cols],ct,dds,res,deg.design)
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("deg","bronch","continuous","analysis_input","ACT-score+Batch12346",Sys.Date(),".csv",sep="_")))
}

## summary table of the DEG analysis
generate_DEG_summary_table<-function(results_significant,deg_design,variable){
  res.sig<-results_significant
  deg.design<-deg_design
  var<-variable
  
  reslist<-paste("res.sig",1:length(res.sig),sep="")
  n_sig_deg<-unlist(sapply(res.sig,nrow))
  design<-deg.design
  
  df<-data.frame(type="bronch",results=reslist,n_sig_deg=n_sig_deg,design=design,variable=var, row.names = NULL)
  return(df)
}

if(dir.exists(deg.dir)){
  a<-generate_DEG_summary_table(res_sig_list,deg.design,"ACT-score")
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("dds","bronch","continuous","res_summary","ACT-score+Batch12346",Sys.Date(),".csv",sep="_")))
}

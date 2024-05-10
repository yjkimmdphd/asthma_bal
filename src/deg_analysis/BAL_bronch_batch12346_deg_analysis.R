## 
# Bronchial RNAseq analysis with new samples:
# additional samples run in the batch 6
# alignment based on new reference genome GRCh38
##

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

################################
## load phenotype and batch data
################################

# make vectors of variables for later use as an input for function 'run_deseq2_DEG_analysis'

source.cell.log<-c(
  "BAL_eos_ct_log",
  "BAL_eos_p_log",
  "BAL_neut_ct_log",
  "BAL_neut_p_log",
  "BAL_wbc_log",
  "blood_eos_log",
  "blood_eos_p_log",
  "blood_neut_log",
  "blood_neut_p_log",
  "blood_wbc_log")
source.cell<-c(
  "BAL_eos_ct",
  "BAL_eos_p",
  "BAL_neut_ct",
  "BAL_neut_p",
  "BAL_wbc",
  "blood_eos",
  "blood_eos_p",
  "blood_neut",
  "blood_neut_p",
  "blood_wbc")


# asthma biomarker phenotype file saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/Nasal_Biomarkers_BAL_transformed.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
numeric_id<-phenotype$ID
phenotype$ID<-sprintf("%03d",numeric_id) # adds padded zeros in front of the subject ID numbers
phenotype<-mutate(phenotype,pos_cellcount=phenotype[,source.cell]>0)%>%arrange(ID) # check which cell counts are positive. 



###########################################################################################
## subset phenotype data for which the samples exist for nasal/bronchial RNAseq experiments   
###########################################################################################
bID<-paste0("B",phenotype$ID) # B*** indicates bronchial sample ID, sequence data is not available for all as of 2023-10-04
bexist<-bID%in%counts.ID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
bsample<-bID[bexist] # bronchial sample ID in the readcount matrix (batch 1-4,6) that has BAL phenotype data
bphen<-phenotype[phenotype$ID%in%substring(bsample,2),] # phenotype table with bsample
bphen<-mutate(bphen, SampleID=bsample)%>%relocate(SampleID, .before=1) # include sample ID for bronchial RNAseq samples

# left join batch info table with nasal/bronchial phenotype table  
## get batch information
batch<-file.path("./resources/raw_data/MS_asthma/MS_asthma_phenotype.batch12346.final.csv")
batch.info<-if(file.exists(batch)){read.csv(batch)}
bronch.batch.info<-batch.info[1:75,2:3]
bronch.batch.info$SampleID<-substr(bronch.batch.info$SampleID,1,4)

## define function join_phenotype_batch_info. p is phenotype table. b is batch info table. Factorize the batch info. 
join_phenotype_batch_info<-function(p,b){
  table<-left_join(p,b, by="SampleID")
  table$Batch<-factor(table$Batch, levels=unique(table$Batch))
  return(table)
}
bphen<-join_phenotype_batch_info(bphen,bronch.batch.info)
bphen<-bphen%>%mutate(IsBatch4 = Batch == "batch4")
# scale the cell count information 
# Mike Love: I'd recommend however to center and scale the numeric covariates for model fitting improvement.
# scale the columns named with source.cell.log

bphen<-mutate_at(bphen,vars(all_of(source.cell.log)),scale) # scales and mutates all log-transformed cell counts 

# decide which analysis to perform, then set the phenotype data as phen


###################################
# custom functions for DEG analysis
###################################
# should load the following fx:
## filter_low_expressed_genes_method2: Filters low counts genes using TMM normalized lcpm as a cutoff point. Requires 'limma'
## rowgenes_counttable: changes the row names of the count table with gene names
## run_deseq2_DEG_analysis: takes countdata,coldata,design,des as variables, and runs DEG on DESeq2
## get_DEG_results: saves result of DESeq2 output, ordered in padj 
## generate_DEG_input_summary_table: makes a table of input information
## generate_DEG_summary_table: makes results summary (i.e., # of DEG for each analysis)

source("./src/function/deg_custom_functions.R")

generate_DEG_summary_table<-function(){
  reslist<-paste("res.sig",1:length(res.sig),sep="")
  n_sig_deg<-sapply(res.sig,nrow)
  design<-deg.design
  source_cell<-source.cell
  df<-data.frame(type="bronch",results=reslist,n_sig_deg,design=design,source_cell=source_cell, row.names = NULL)
  return(df)
} # needed to fix "type=unique(phen$Type)" to "type="bronch""
#-----------------------
##############################################################
#set colData (phenotype data) for bronchial RNAseq experiments
##############################################################
phen<-bphen  # If bronchial analysis, use this

# all non-NA values. cell count >=0
pi<-lapply(phen[,source.cell.log],function(data){a<-!is.na(data);return(a)})
df<-vector("list",length=10) # list of data framese used as an input for deseq2. all cell counts
names(df)<-paste(source.cell.log,"all",sep="_")
for(i in 1:10){
  df[[i]]<-phen[pi[[i]],c("SampleID",source.cell.log[i], "Batch","IsBatch4")]
}
print(sapply(df,dim)[1,])

# all cell count >0
pi.pos<-lapply(phen[,source.cell],function(data){a<-which(data>0);return(a)})
df.pos<-vector("list",length=10) # list of data framese used as an input for deseq2. subset cell count > 0
names(df.pos)<-paste(source.cell.log,"pos",sep = "_")
for(i in 1:10){
  df.pos[[i]]<-phen[pi.pos[[i]],c("SampleID",source.cell.log[i], "Batch","IsBatch4")]
}
print(sapply(df.pos,dim)[1,])
# if analyzing only cell counts >0, use df<-df.pos

#################################################################
# bronchial expression ~ log(cell count>=0) + Batch
#################################################################
# coldata for DESeq2
df.input<-df


# filtering counts table to remove low expressed genes

## select RNAseq counts
id<-phen$SampleID
cols<-colnames(bronch.counts)%in%id
ct<-bronch.counts[,cols] # First column is actually gene name 
genes<-bronch.counts$SampleID
rownames(ct)<-genes

## Filter counts (readcount table for nasal sample
c2<-filter_low_expressed_genes_method2(ct,4)


# run the DEG for continuous predictors

deg.design<-paste("~",source.cell.log,"+ Batch") # set design: nasal expression ~ log(cell count>0) + Batch 
ct<-rowgenes_counttable(ct,c2) # low bcounts will be filtered using method 2: use TMM normalized lcpm as a cutoff point

print(deg.design)
count.table<-lapply(df.input,function(df){d<-df; ct<-ct[,colnames(ct)%in%d$SampleID]; return(ct)}) # list of subsetted count table. Each element is a count table with samples for each of the experimental design. 

dds<-vector("list",length=10)
res<-vector("list",length=10)
res.sig<-vector("list",length=10)

names(res)<-deg.design
names(res.sig)<-deg.design

for(i in 1:10){
  dds[[i]]<-run_deseq2_DEG_analysis(count.table[[i]], df.input[[i]], deg.design[i],deg.design[i])
  res[[i]]<-get_DEG_results(dds[[i]], source.cell.log[i])
  res.sig[[i]]<-res[[i]][which(res[[i]]$padj<=0.05),]
  head(res.sig[[i]])
  
}

## writing the results 
deg.folder<-paste("deg",Sys.Date(),sep="_")
deg.dir<-file.path("./reports",deg.folder)
if(!dir.exists(deg.dir)){
  dir.create(deg.dir)
}

### write all results


### write only the significant results
if(dir.exists(deg.dir)){
  for(i in 1:10){
    a<-res.sig[[i]]
    write.csv(a,row.names=TRUE,file.path(deg.dir,paste("deg","bronch","allcells",deg.design[[i]],"deg_res",i,Sys.Date(),".csv",sep="_"))) } #specify allcells vs poscells
}

if(dir.exists(deg.dir)){
  for(i in 1:10){
    a<-data.frame(res[[i]])
    write.csv(a,row.names=TRUE,file.path(deg.dir,paste("deg","bronch","allcells",deg.design[[i]],"all_results_of_res",i,Sys.Date(),".csv",sep="_"))) } #specify allcells vs poscells
}

## summarize the data input 
if(dir.exists(deg.dir)){
  a<-generate_DEG_input_summary_table()
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("deg","bronch","allcells","cellcount+batch","analysis_input",Sys.Date(),".csv",sep="_"))) #specify allcells vs poscells
}

## summary table of the DEG analysis
if(dir.exists(deg.dir)){
  a<-generate_DEG_summary_table()
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("dds","bronch","allcells","cellcount+batch","res_summary",Sys.Date(),".csv",sep="_"))) #specify allcells vs poscells
}


#################################################################
# bronchial expression ~ log(cell count>0) + Batch
#################################################################
# coldata for DESeq2
df.input<-df.pos

# filtering counts table to remove low expressed genes

## select RNAseq counts
id<-phen$SampleID
cols<-colnames(bronch.counts)%in%id
ct<-bronch.counts[,cols] # First column is actually gene name 
genes<-bronch.counts$SampleID
rownames(ct)<-genes

## Filter counts (readcount table for nasal sample
c2<-filter_low_expressed_genes_method2(ct,4)


# run the DEG for continuous predictors

deg.design<-paste("~",source.cell.log,"+ Batch") # set design: nasal expression ~ log(cell count>0) + Batch 
ct<-rowgenes_counttable(ct,c2) # low bcounts will be filtered using method 2: use TMM normalized lcpm as a cutoff point

print(deg.design)
count.table<-lapply(df.input,function(df){d<-df; ct<-ct[,colnames(ct)%in%d$SampleID]; return(ct)}) # list of subsetted count table. Each element is a count table with samples for each of the experimental design. 

dds<-vector("list",length=10)
res<-vector("list",length=10)
res.sig<-vector("list",length=10)

names(res)<-deg.design
names(res.sig)<-deg.design

for(i in 1:10){
  dds[[i]]<-run_deseq2_DEG_analysis(count.table[[i]], df.input[[i]], deg.design[i],deg.design[i])
  res[[i]]<-get_DEG_results(dds[[i]], source.cell.log[i])
  res.sig[[i]]<-res[[i]][which(res[[i]]$padj<=0.05),]
  head(res.sig[[i]])
  
}

## writing the results 
deg.folder<-paste("deg",Sys.Date(),sep="_")
deg.dir<-file.path("./reports",deg.folder)
if(!dir.exists(deg.dir)){
  dir.create(deg.dir)
}

if(dir.exists(deg.dir)){
  for(i in 1:10){
    a<-res.sig[[i]]
    write.csv(a,row.names=TRUE,file.path(deg.dir,paste("deg","bronch","poscells",deg.design[[i]],"deg_res",i,Sys.Date(),".csv",sep="_"))) } #specify allcells vs poscells
}

if(dir.exists(deg.dir)){
  for(i in 1:10){
    a<-data.frame(res[[i]])
    write.csv(a,row.names=TRUE,file.path(deg.dir,paste("deg","bronch","poscells",deg.design[[i]],"all_results_of_res",i,Sys.Date(),".csv",sep="_"))) } #specify allcells vs poscells
}


## summarize the data input 
if(dir.exists(deg.dir)){
  a<-generate_DEG_input_summary_table()
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("deg","bronch","poscells","cellcount+batch","analysis_input",Sys.Date(),".csv",sep="_"))) #specify allcells vs poscells
}

## summary table of the DEG analysis
if(dir.exists(deg.dir)){
  a<-generate_DEG_summary_table()
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("dds","bronch","poscells","cellcount+batch","res_summary",Sys.Date(),".csv",sep="_"))) #specify allcells vs poscells
}
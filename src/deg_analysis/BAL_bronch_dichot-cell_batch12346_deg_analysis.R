######
# DEG Bronch mRNA ~  dichotomous cell count + batch12346
######
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
phenotype<-file.path("./resources/processed_data/nasal_biomarker_phenotype_batch12346_merged.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
numeric_id<-phenotype$ID
phenotype$ID<-sprintf("%03d",numeric_id) # adds padded zeros in front of the subject ID numbers
phenotype<-mutate(phenotype,pos_cellcount=phenotype[,source.cell]>0)%>%arrange(ID) # check which cell counts are positive. 

###########################################################################################
## subset phenotype data for which the samples exist for bronchial RNAseq experiments   
###########################################################################################
bexist<-phenotype$SampleID%in%counts.ID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
bphen<-phenotype[bexist,]
bphen<-mutate_at(bphen,vars(all_of(source.cell.log)),scale)

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

# Function to get quartiles for each column
quartile_thresholds <- function(df) {
  apply(df, 2, function(x) {
    quantile(x, probs = c(0.25, 0.5, 0.75),na.rm=TRUE)
  })
}

# Calculate quartile thresholds
quartile_results <- quartile_thresholds(bphen[,c(source.cell,source.cell.log)])
print(quartile_results)


# compare the each samples bal Eos data against the threshold 
bphen<-bphen%>%mutate(bal_AEC_more_0=BAL_eos_ct>0,
                      bal_AEC_more_1.15=BAL_eos_ct>1.15,
                      bal_Eos_p_more_0 = BAL_eos_p>0,
                      bal_Eos_p_more_1 = BAL_eos_p>1,
                      bal_Eos_p_more_3 = BAL_eos_p>3,
                      bal_ANC_more_0=BAL_neut_ct>0,
                      bal_ANC_more_5=BAL_neut_ct>5,
                      bal_ANC_more_13=BAL_neut_ct>13,
                      bal_neut_p_more_0 = BAL_neut_p>0,
                      bal_neut_p_more_2 = BAL_neut_p>2,
                      bal_neut_p_more_5 = BAL_neut_p>5,
                      bld_AEC_more_0 = blood_eos>0,
                      bld_AEC_more_100 = blood_eos>100,
                      bld_AEC_more_300 = blood_eos>300)

var_dichot<-c("bal_AEC_more_0","bal_AEC_more_1.15","bal_Eos_p_more_0",
              "bal_Eos_p_more_1","bal_Eos_p_more_3","bal_ANC_more_0",
              "bal_ANC_more_5","bal_ANC_more_13","bal_neut_p_more_0",
              "bal_neut_p_more_2","bal_neut_p_more_5","bld_AEC_more_0",
              "bld_AEC_more_100","bld_AEC_more_300")
##############################################################
#set colData (phenotype data) for bronchial RNAseq experiments
##############################################################
phen<-bphen  # If bronchial analysis, use this

# all non-NA values. cell count >=0
pi<-lapply(phen[,var_dichot],function(data){a<-!is.na(data);return(a)})
df<-vector("list",length=length(var_dichot)) # list of data framese used as an input for deseq2. all cell counts
names(df)<-var_dichot
for(i in 1:length(var_dichot)){
  df[[i]]<-phen[pi[[i]],c("SampleID",var_dichot[i], "Batch")]
}
print(sapply(df,dim)[1,])

#################################################################
# bronchial expression ~ cellcount (dichotomous) + Batch
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
c2<-filter_low_expressed_genes_method2(ct,7)

## design: gene expression ~ is_cellcount_threshold + Batch

deg.design<-paste("~",names(df),"+ Batch")
print(deg.design)
ct<-rowgenes_counttable(ct,c2) # low counts will be filtered 

print(deg.design)

count.table<-lapply(df.input,function(df){d<-df; ct<-ct[,colnames(ct)%in%d$SampleID]; return(ct)}) # list of subsetted count table. Each element is a count table with samples for each of the experimental design. 

dds<-vector("list",length=length(df))
res<-vector("list",length=length(df))
res.sig<-vector("list",length=length(df))

names(res)<-deg.design
names(res.sig)<-deg.design

for(i in 1:length(df)){
  dds[[i]]<-run_deseq2_DEG_analysis(count.table[[i]], df.input[[i]], deg.design[i],deg.design[i])
  res[[i]]<-get_DEG_results(dds[[i]], paste(names(df),"TRUE",sep="")[i])
  res.sig[[i]]<-res[[i]][which(res[[i]]$padj<=0.05),]
  head(res.sig[[i]])
}


for(i in 1:1){
  dds[[i]]<-run_deseq2_DEG_analysis(count.table[[i]], df.input[[i]], deg.design[i],deg.design[i])
  res[[i]]<-get_DEG_results(dds[[i]], paste(names(df),"TRUE",sep="")[i])
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
  for(i in 1:length(df)){
    a<-res.sig[[i]]
    write.csv(a,row.names=TRUE,file.path(deg.dir,paste("deg",unique(phen$Type),"allcells","res",i,deg.design[[i]],Sys.Date(),".csv",sep="_"))) }
}

## summarize the data input 

if(dir.exists(deg.dir)){
  a<-generate_DEG_input_summary_table()
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("deg",unique(phen$Type),"allcells","analysis_input","cellcount_thr+Batch",Sys.Date(),".csv",sep="_")))
}

## summary table of the DEG analysis

if(dir.exists(deg.dir)){
  a<-generate_DEG_summary_table()
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("dds",unique(phen$Type),"allcells","res_summary","cellcount_thr+Batch",Sys.Date(),".csv",sep="_")))
}


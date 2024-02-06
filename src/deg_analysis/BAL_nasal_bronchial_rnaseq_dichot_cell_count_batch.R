##
# RNAseq with nasal/bronchial rna expression ~ dichotomous cell count (ie percent > threshold, count > threshold) + batch
## 
library(limma)
library(edgeR)
library(dplyr)
library(DESeq2)
######################
## load readcount data
######################

# load count data from seq batch 1-4
counts<-file.path("./resources/working_data/copy_of_batch1234_readcount_matrix_allsamples.afterQC.txt") # original RNAseqs count data table is in the MS_asthma folder 
counts<-if(file.exists(counts)){read.delim(counts)}
counts.ID<-colnames(counts)

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

nID<-paste0("N",phenotype$ID) # N*** indicates nasal sample ID, sequence data is not available for all as of 2023-10-04
bID<-paste0("B",phenotype$ID) # B*** indicates bronchial sample ID, sequence data is not available for all as of 2023-10-04

nexist<-nID%in%counts.ID # find which subjects s/p nasal and had bronchial sample RNAseq completed. Nasal samples in batch 1-4 only sequenced to ID 337
bexist<-bID%in%counts.ID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
nsample<-nID[nexist] # nasal sample ID in the readcount matrix (batch 1-4) that has BAL phenotype data
bsample<-bID[bexist] # bronchial sample ID in the readcount matrix (batch 1-4) that has BAL phenotype data

nphen<-phenotype[phenotype$ID%in%substring(nsample,2),] # phenotype table with nsample
bphen<-phenotype[phenotype$ID%in%substring(bsample,2),] # phenotype table with bsample

nphen<-mutate(nphen, SampleID=nsample)%>%relocate(SampleID, .before=1) # include sample ID for nasal RNAseq samples
bphen<-mutate(bphen, SampleID=bsample)%>%relocate(SampleID, .before=1) # include sample ID for bronchial RNAseq samples

# left join batch info table with nasal/bronchial phenotype table  
## get batch information
batch<-file.path("./resources/processed_data/asthma_nasal_bronchial_batch1234_info.txt")
batch.info<-if(file.exists(batch)){read.delim(batch)}

## define function join_phenotype_batch_info. p is phenotype table. b is batch info table. Factorize the batch info. 
join_phenotype_batch_info<-function(p,b){
  table<-left_join(p,b, by="SampleID")
  table$Batch<-factor(table$Batch, levels=unique(table$Batch))
  return(table)
}
nphen<-join_phenotype_batch_info(nphen,batch.info)
bphen<-join_phenotype_batch_info(bphen,batch.info)
nphen<-nphen%>%mutate(IsBatch4 = Batch == "batch4")
bphen<-bphen%>%mutate(IsBatch4 = Batch == "batch4")
# scale the cell count information 
# Mike Love: I'd recommend however to center and scale the numeric covariates for model fitting improvement.
# scale the columns named with source.cell.log
nphen<-mutate_at(nphen,vars(all_of(source.cell.log)),scale) # scales and mutates all log-transformed cell counts 
bphen<-mutate_at(bphen,vars(all_of(source.cell.log)),scale) # scales and mutates all log-transformed cell counts 

# decide which analysis to perform, then set the phenotype data as phen

###################################
# custom functions for DEG analysis
###################################

## gene filter
filter_low_expressed_genes_method2<-function(readcounts, n_sample){ # function to filter low counts genes using TMM normalized lcpm as a cutoff point
  x<-readcounts
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
  drop.genes<-which(rowSums(lcpm.x<lcpm.cutoff)>n_sample) # which genes have normalized lcpm less than the cutoff in >10% of the samples
  x<-x[-drop.genes,] # readcounts retained 
  print(paste("# of dropped genes:",length(drop.genes))) # number of dropped genes
  return(x)
}
## count table-related
rowgenes_counttable<-function(c,cf){ # function to change rownames of the count table with gene names 
  g<-rownames(c)
  gf<-rownames(cf)
  c<-c[g%in%gf,]
  return(c)
}

## DEsesq
run_deseq2_DEG_analysis<-function(countdata,coldata,design,des){ # function to run deseq using data from matrix 
  print(paste("design:",unique(phen$Type),"expression =",des, sep = " "))
  dds<-DESeqDataSetFromMatrix(countData = countdata,colData=coldata, design=as.formula(design))
  dds<-DESeq(dds)
  print("resultnames")
  print(resultsNames(dds))
  return(dds)
}

get_DEG_results<-function(dds,resultname){ # function to call results of the deseq analysis with 'results' function
  res<-results(dds, name=resultname) # use source.cell.log element for resultname.
  res <- res[order(res$padj),]
  return(res)
}

## summarize the data input 

generate_DEG_input_summary_table<-function(){
  filter_method<-"TMM normalized LCPM cutoff"
  n_filtered_genes<-paste("analyzed n_genes:", nrow(ct),",","filtered n_genes:",nrow(counts)-nrow(ct))
  samples<-sapply(dds, function(d){colData(d)$SampleID%>%paste(collapse = ",")})
  dds<-paste("dds",1:length(dds),sep="")
  results<-paste("res",1:length(res),sep="")
  design<-deg.design
  df<-data.frame(dds=dds,results=results,design=design,samples=samples,filter_method=filter_method,n_filtered_genes=n_filtered_genes)
  return(df)
}

##  summary table of the DEG analysis
generate_DEG_summary_table<-function(){
  reslist<-paste("res.sig",1:length(res.sig),sep="")
  n_sig_deg<-sapply(res.sig,nrow)
  design<-deg.design
  source_cell<-source.cell
  df<-data.frame(type=unique(phen$Type),results=reslist,n_sig_deg,design=design,source_cell=source_cell, row.names = NULL)
  return(df)
}

##########################################################
#set colData (phenotype data) for nasal RNAseq experiments
##########################################################

# if nasal cell analysis, phen<-nphen. 
phen<-nphen  # If nasal analysis, use this
# phen<-bphen # If bronchial analysis, use this

# determine which phenotypoic thresholds to use by assessing the proportion of samples below the potential thresholds

checkThreshold<-function(data){
  a<-c(max(which(data[,2]<0.54)),  max(which(data[,2]<0.60)),  max(which(data[,2]<0.70)),  max(which(data[,2]<0.80)))
  b<-c(max(which(data[,3]<0.54)),  max(which(data[,3]<0.60)),  max(which(data[,3]<0.70)),  max(which(data[,3]<0.80)))
  return(data.frame(proportion=c("54%","60%","70%","80%"),bal_eos_p_thr=a,bal_eos_ct_thr=b))
} # function that helps to determine the thresholds


## For BAL Eos

### threshold assessment for bal Eos %
a<-numeric()
for(i in 1:40){
  a[i]<-sum(phen$BAL_eos_p<i)/nrow(phen)
}

### threshold assessment for bal Eos count 
b<-numeric()
for(i in 1:40){
  b[i]<-sum(phen$BAL_eos_ct<i)/nrow(phen)
}
### Eos summary
eos_thr<-data.frame(threshold=1:40,eos_p_fraction= round(a,2),eos_ct_fraction=round(b,2))
checkThreshold(eos_thr)


## For blood Eos

### threshold assessment for blood Eos %
c<-numeric()
na.omit(phen$blood_eos_p)%>%range
for(i in 1:18){
  c[i]<-sum(na.omit(phen$blood_eos_p)<i)/length(na.omit(phen$blood_eos_p))
}

### threshold assessment for blood Eos count 
d<-numeric()
for(i in seq(1,1701,by=100)){
  d[round(i/100+1,0)]<-sum(na.omit(phen$blood_eos)<i)/length(na.omit(phen$blood_eos))
}

### blood Eos summary
print(data.frame(blood_eos_p_thr=1:18,proportion=c, blood_eos_ct_thr= seq(1,1701,by=100),proportion=d))

# compare the each samples bal Eos data against the threshold 
phen<-phen%>%mutate(isBalEosCtNa=is.na(BAL_eos_ct),
                    isbloodEosCtNa=is.na(blood_eos),
              isBalEosCtZero=!BAL_eos_ct==0,
              isBalEosct2=!BAL_eos_ct<2,
              isBalEosP2=!BAL_eos_p<2,
              isBalEosct6=!BAL_eos_ct<6,
              isBalEosP6=!BAL_eos_p<6,
              isbldEosCtZero=!blood_eos==0,
              isbldEosp4=!blood_eos_p<4,
              isbldEosp6=!blood_eos_p<6,
              isbldEosct200=!blood_eos<200,
              isbldEosct400=!blood_eos<400
              ) # if greater than threshold, 'TRUE'; if lower than threshold, 'False'. 

df<-vector("list",length=length(49:58)) # list of data framese used as an input for deseq2. all cell counts
names(df)<-colnames(phen)[49:58]
pi<-lapply(phen[,names(df)],function(data){a<-!is.na(data);return(a)})
for(i in 1:length(df)){
  df[[i]]<-phen[pi[[i]],c("SampleID",names(df)[i], "Batch","isBalEosCtNa","isbloodEosCtNa")]
}

###########################################---------------------###########################################
# sample type: nasal
# gene filter: remove low counts
###########################################---------------------###########################################
# coldata for DESeq2
df.input<-df

# filtering counts table to remove low expressed genes

## select RNAseq counts
id<-phen$SampleID
cols<-colnames(counts)%in%id
ct<-counts[,cols] # First column is actually gene name 
genes<-counts$SampleID
rownames(ct)<-genes

## Filter counts (readcount table for nasal sample

c2<-filter_low_expressed_genes_method2(ct,4)

# run the DEG for continuous predictors


## design: gene expression ~ is_cellcount_threshold + Batch

deg.design<-paste("~",names(df),"+ Batch")
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


##########################################################
#set colData (phenotype data) for bronchial RNAseq experiments
##########################################################

# if nasal cell analysis, phen<-nphen. 
phen<-bphen  # If nasal analysis, use this
# phen<-bphen # If bronchial analysis, use this

# determine which phenotypoic thresholds to use by assessing the proportion of samples below the potential thresholds

checkThreshold<-function(data){
  a<-c(max(which(data[,2]<0.54)),  max(which(data[,2]<0.60)),  max(which(data[,2]<0.70)),  max(which(data[,2]<0.80)))
  b<-c(max(which(data[,3]<0.54)),  max(which(data[,3]<0.60)),  max(which(data[,3]<0.70)),  max(which(data[,3]<0.80)))
  return(data.frame(proportion=c("54%","60%","70%","80%"),bal_eos_p_thr=a,bal_eos_ct_thr=b))
} # function that helps to determine the thresholds


## For BAL Eos

### threshold assessment for bal Eos %
a<-numeric()
for(i in 1:40){
  a[i]<-sum(phen$BAL_eos_p<i)/nrow(phen)
}

### threshold assessment for bal Eos count 
b<-numeric()
for(i in 1:40){
  b[i]<-sum(phen$BAL_eos_ct<i)/nrow(phen)
}
### Eos summary
eos_thr<-data.frame(threshold=1:40,eos_p_fraction= round(a,2),eos_ct_fraction=round(b,2))
checkThreshold(eos_thr)


## For blood Eos

### threshold assessment for blood Eos %
c<-numeric()
na.omit(phen$blood_eos_p)%>%range
for(i in 1:18){
  c[i]<-sum(na.omit(phen$blood_eos_p)<i)/length(na.omit(phen$blood_eos_p))
}

### threshold assessment for blood Eos count 
d<-numeric()
for(i in seq(1,1701,by=100)){
  d[round(i/100+1,0)]<-sum(na.omit(phen$blood_eos)<i)/length(na.omit(phen$blood_eos))
}

### blood Eos summary
print(data.frame(blood_eos_p_thr=1:18,proportion=c, blood_eos_ct_thr= seq(1,1701,by=100),proportion=d))

# compare the each samples bal Eos data against the threshold 
phen<-phen%>%mutate(isBalEosCtNa=is.na(BAL_eos_ct),
                    isbloodEosCtNa=is.na(blood_eos),
                    isBalEosCtZero=BAL_eos_ct==0,
                    isBalEosct2=BAL_eos_ct<2,
                    isBalEosP2=BAL_eos_p<2,
                    isBalEosct6=BAL_eos_ct<6,
                    isBalEosP6=BAL_eos_p<6,
                    isbldEosCtZero=blood_eos==0,
                    isbldEosp4=blood_eos_p<4,
                    isbldEosp6=blood_eos_p<6,
                    isbldEosct200=blood_eos<200,
                    isbldEosct400=blood_eos<400
)

df<-vector("list",length=length(49:58)) # list of data framese used as an input for deseq2. all cell counts
names(df)<-colnames(phen)[49:58]
pi<-lapply(phen[,names(df)],function(data){a<-!is.na(data);return(a)})
for(i in 1:length(df)){
  df[[i]]<-phen[pi[[i]],c("SampleID",names(df)[i], "Batch","isBalEosCtNa","isbloodEosCtNa")]
}

###########################################---------------------###########################################
# sample type: bronchial
# gene filter: remove low counts
###########################################---------------------###########################################
# coldata for DESeq2
df.input<-df

# filtering counts table to remove low expressed genes

## select RNAseq counts
id<-phen$SampleID
cols<-colnames(counts)%in%id
ct<-counts[,cols] # First column is actually gene name 
genes<-counts$SampleID
rownames(ct)<-genes

## Filter counts (readcount table for nasal sample

c2<-filter_low_expressed_genes_method2(ct,4)

# run the DEG for continuous predictors


## design: gene expression ~ is_cellcount_threshold + Batch

deg.design<-paste("~",names(df),"+ Batch")
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

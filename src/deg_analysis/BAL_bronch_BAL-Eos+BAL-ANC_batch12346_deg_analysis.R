
######
# DEG Bronch mRNA ~  Eos + ANC (>0) cell count + batch12346
######
library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)

# load cell count table
countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
genes<-counts[,"SampleID"]

# select bronchial samples 
bronch.samples<-grepl("^B",colnames(counts))
bronch.counts<-counts[,bronch.samples]
rownames(bronch.counts)<-genes
head(bronch.counts)
counts.ID<-colnames(bronch.counts)

################################
## load phenotype and batch data
################################

# make vectors of variables for later use as an input for function 'run_deseq2_DEG_analysis'

source.cell.log<-c(
  "BAL_eos_ct_log",
  "BAL_eos_p_log")
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
phenotype<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
phenotype<-phenotype%>%filter(grepl("^B",SampleID))


#####################################################################################
## subset phenotype data for which the samples exist for bronchial RNAseq experiments   
#####################################################################################
bexist<-phenotype$SampleID%in%counts.ID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
bphen<-phenotype[bexist,]

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

source("./src/function/deg_custom_functions_v2.R")



bphen<-bphen%>%mutate(bal_AEC_more_0=BAL_eos_ct>0,
                      bal_AEC_more_1=BAL_eos_ct>1,
                      bal_AEC_more_1.2=BAL_eos_ct>1.2,
                      bal_AEC_more_3=BAL_eos_ct>3,
                      bal_AEC_more_5=BAL_eos_ct>5,
                      bal_Eos_p_more_0 = BAL_eos_p>0,
                      bal_Eos_p_more_1 = BAL_eos_p>1,
                      bal_Eos_p_more_3 = BAL_eos_p>3,
                      bal_Eos_p_more_5 = BAL_eos_p>5)


var_dichot_bal<-c("bal_AEC_more_0","bal_AEC_more_1","bal_AEC_more_1.2","bal_AEC_more_3",
                  "bal_AEC_more_5","bal_Eos_p_more_0","bal_Eos_p_more_1","bal_Eos_p_more_3","bal_Eos_p_more_5")



# the model will be BAL Eos % or count + BAL-ANC + batch, so filter for rows that has values for BAL ANC 
phen<-bphen%>%filter(!is.na(BAL_neut_p_log))

############ select variables to test for all non-NA values. cell count >=0
var_to_test<-c(source.cell.log,var_dichot_bal)

var_to_test_res<-c(source.cell.log,paste(c(var_dichot_bal),"TRUE",sep=""))

# make a list of the phenotype colData that will be used for DESeq2
pi<-lapply(phen[,var_to_test],function(data){a<-!is.na(data);return(a)})
df<-vector("list",length(var_to_test)) # list of data framese used as an input for deseq2. all cell counts
names(df)<-paste(var_to_test,"all",sep="_")
for(i in 1:length(var_to_test)){
  df[[i]]<-phen[pi[[i]],c("SampleID",var_to_test[i], "BAL_neut_p_log", "Batch")]
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

## Filter counts (readcount table for nasal sample
c2<-filter_low_expressed_genes_method2(ct,round(length(id)*0.1,0))

## design: gene expression ~ is_cellcount_threshold + Batch


deg.design<-paste("~",var_to_test,"+BAL_neut_p_log + Batch")
ct<-rowgenes_counttable(ct,c2) # low bcounts will be filtered 

print(deg.design)
count.table<-lapply(df.input,function(df){d<-df; ct<-ct[,colnames(ct)%in%d$SampleID]; return(ct)}) # list of subsetted count table. Each element is a count table with samples for each of the experimental design. 

dds<-vector("list",length=length(var_to_test))
res<-vector("list",length=length(var_to_test))
res.sig<-vector("list",length=length(var_to_test))

names(res)<-deg.design
names(res.sig)<-deg.design

# at this time, testing only the nasal ~ blood AEC 
assay_index<-c(1:length(deg.design))
for(i in assay_index){
  dds[[i]]<-run_deseq2_DEG_analysis(count.table[[i]], df.input[[i]], deg.design[i],deg.design[i])
  res[[i]]<-get_DEG_results(dds[[i]], var_to_test_res[i])
  res.sig[[i]]<-res[[i]][which(res[[i]]$padj<=0.05),]
  head(res.sig[[i]])
  
}


## writing the significant and all results 
deg.folder<-paste("deg",Sys.Date(),sep="_")
deg.dir<-file.path("./reports","temporary",deg.folder)
if(!dir.exists(deg.dir)){
  dir.create(deg.dir)
}

if(dir.exists(deg.dir)){
  for(i in assay_index){
    a<-res.sig[[i]]
    b<-res[[i]]
    write.csv(a,row.names=TRUE,file.path(deg.dir,paste("deg","bronch","res_sig",i,deg.design[[i]],Sys.Date(),".csv",sep="_")))
    write.csv(b,row.names=TRUE,file.path(deg.dir,paste("deg","bronch","res_all",i,deg.design[[i]],Sys.Date(),".csv",sep="_")))
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
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("deg","bronch","dichot","analysis_input","cellcount_thr+Batch",Sys.Date(),".csv",sep="_")))
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
  a<-generate_DEG_summary_table(res.sig,deg.design,var_to_test)
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("dds","bronch","dichot","res_summary","cellcount_thr+Batch",Sys.Date(),".csv",sep="_")))
}


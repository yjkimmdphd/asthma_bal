
## reanalysis of nasal ~ cell count 


library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)

# load count table 
countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}

# set rownames as gene names
rownames(counts)<-counts[,"SampleID"]

# remove f/u samples
filtered_samples<-counts[,!grepl("F",colnames(counts))]
filtered_ID<-grepl("^N",colnames(filtered_samples))
ncounts<-filtered_samples[,filtered_ID]
counts.ID<-colnames(ncounts)

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

var_dichot_bal<-c("bal_AEC_more_0","bal_AEC_more_1","bal_Eos_p_more_0",
              "bal_Eos_p_more_1","bal_Eos_p_more_3","bal_ANC_more_0",
              "bal_ANC_more_5","bal_ANC_more_13","bal_neut_p_more_0",
              "bal_neut_p_more_2","bal_neut_p_more_5")
var_dichot_blood<-c("bld_AEC_more_0",
              "bld_AEC_more_100",
              "bld_AEC_more_300",
              "bld_AEC_more_400")

# asthma biomarker phenotype file saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
phenotype<-phenotype%>%filter(grepl("^N",SampleID))%>%filter(!grepl("F",SampleID))

# to the phenotype data left join data for sampling date differences 
sampling_date_diff<-"./resources/processed_data/sampling_dates/swab-bal-cbc_differences_in_days.txt"
sampling_date_diff<-if(file.exists(sampling_date_diff)){read.table(sampling_date_diff,row.names = NULL,header = TRUE)}
sampling_date_diff<-sampling_date_diff%>%filter(Comparison=="blood_nasal")
colnames(sampling_date_diff)[1:3]<-c("ID","sampling_date_comp","sampling_date_diff_days")

phenotype<-left_join(phenotype,sampling_date_diff,by="ID")

###########################################################################################
## subset phenotype data for which the samples exist for nasal/bronchial RNAseq experiments   
###########################################################################################

nID<-phenotype$SampleID

nexist<-nID%in%counts.ID # find which subjects s/p nasal sample RNAseq completed.
nsample<-nID[nexist] 

nphen<-phenotype[phenotype$SampleID%in%nsample,]

nphen<-mutate(nphen, SampleID=nsample)%>%relocate(SampleID, .before=1)


# make new columns for categorical variable (i.e., cell count threshold)
nphen<-nphen%>%mutate(bal_AEC_more_1 = BAL_eos_ct>1,
                      bld_AEC_more_0 = blood_eos>0,
                      bld_AEC_more_100 = blood_eos>100,
                      bld_AEC_more_300 = blood_eos>300,
                      bld_AEC_more_400 = blood_eos>400)
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


phen<-nphen  # If bronchial analysis, use this



############ select variables to test for all non-NA values. cell count >=0
var_to_test<-c(source.cell.log,var_dichot_bal[2],var_dichot_blood)
var_to_test_bld<-var_to_test[c(grep("blood",var_to_test),grep("bld",var_to_test))]
var_to_test_res<-c(source.cell.log,paste(c(var_dichot_bal[2],var_dichot_blood),"TRUE",sep=""))

# make a list of the phenotype colData that will be used for DESeq2
pi<-lapply(phen[,var_to_test],function(data){a<-!is.na(data);return(a)})
df<-vector("list",length(var_to_test)) # list of data framese used as an input for deseq2. all cell counts
names(df)<-paste(var_to_test,"all",sep="_")
for(i in 1:length(var_to_test)){
  df[[i]]<-phen[pi[[i]],c("SampleID",var_to_test[i], "Batch")]
}
print(sapply(df,dim)[1,])

# identify the samples for which cbc information will be used as a variable. sampling_date_diff_days should be less than a year 
cbc_sampleID<-nphen%>%filter(!is.na(vars(var_to_test_bld)),abs(sampling_date_diff_days)<365)%>%pull(SampleID)
blood_df<-df[paste(var_to_test_bld,"all",sep="_")]
blood_df_filtered<-lapply(blood_df,function(df)filter(df,SampleID%in%cbc_sampleID))
sapply(blood_df,dim)[1,] # samples before filtering
sapply(blood_df_filtered,dim)[1,] # samples before filtering

df[paste(var_to_test_bld,"all",sep="_")]<-blood_df_filtered

#################################################################
# nasal expression ~ cellcount + Batch
#################################################################
# coldata for DESeq2
df.input<-df


# filtering counts table to remove low expressed genes

## select RNAseq counts
id<-phen$SampleID
cols<-colnames(ncounts)%in%id
ct<-ncounts[,cols] # First column is actually gene name 
genes<-rownames(ct)


## Filter counts (readcount table for nasal sample
c2<-filter_low_expressed_genes_method2(ct,round(length(id)*0.1,0))


## design: Batches. all cell counts 

deg.design<-paste("~",var_to_test,"+ Batch")
ct<-rowgenes_counttable(ct,c2) # low bcounts will be filtered 

print(deg.design)
count.table<-lapply(df.input,function(df){d<-df; ct<-ct[,colnames(ct)%in%d$SampleID]; return(ct)}) # list of subsetted count table. Each element is a count table with samples for each of the experimental design. 

dds<-vector("list",length=length(var_to_test))
res<-vector("list",length=length(var_to_test))
res.sig<-vector("list",length=length(var_to_test))

names(res)<-deg.design
names(res.sig)<-deg.design

# at this time, testing only the nasal ~ blood AEC 
assay_index<-seq_along(deg.design)
for(i in assay_index){
  dds[[i]]<-run_deseq2_DEG_analysis(count.table[[i]], df.input[[i]], deg.design[i],deg.design[i])
  res[[i]]<-get_DEG_results(dds[[i]], var_to_test_res[i])
  res.sig[[i]]<-res[[i]][which(res[[i]]$padj<=0.05),]
  head(res.sig[[i]])
  
}

## writing the significant and all results 
deg.folder<-paste("deg",Sys.Date(),sep="_")
deg.dir<-file.path("./reports",deg.folder)
if(!dir.exists(deg.dir)){
  dir.create(deg.dir)
}

if(dir.exists(deg.dir)){
  for(i in assay_index){
    a<-res.sig[[i]]
    b<-res[[i]]
    write.csv(a,row.names=TRUE,file.path(deg.dir,paste("deg","nasal","res_sig",i,deg.design[[i]],Sys.Date(),".csv",sep="_")))
    write.csv(b,row.names=TRUE,file.path(deg.dir,paste("deg","nasal","res_all",i,deg.design[[i]],Sys.Date(),".csv",sep="_")))
  }
}

## summarize the data input 

generate_DEG_input_summary_table<-function(original_ct,filtered_ct,dds,res,des){
  filter_method<-"TMM normalized LCPM cutoff"
  n_filtered_genes<-paste("analyzed n_genes:", nrow(filtered_ct),",","filtered n_genes:",nrow(original_ct)-nrow(filtered_ct))
  samples<-sapply(dds, function(d){colData(d)$SampleID%>%paste(collapse = ",")})
  dds<-paste("dds", assay_index,sep="")
  results<-paste("res",assay_index,sep="")
  design<-des
  df<-data.frame(dds=dds,results=results,design=design,samples=samples,filter_method=filter_method,n_filtered_genes=n_filtered_genes)
  return(df)
}

if(dir.exists(deg.dir)){
  a<-generate_DEG_input_summary_table(ncounts[,cols],ct,dds,res,deg.design)
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("deg","nasal","continuous_all_and_dich","analysis_input","cellcount+Batch12346",Sys.Date(),".csv",sep="_")))
}

## summary table of the DEG analysis
generate_DEG_summary_table<-function(results_significant,deg_design,variable){
  res.sig<-results_significant
  deg.design<-deg_design
  var<-variable
  
  reslist<-paste("res.sig",assay_index,sep="")
  n_sig_deg<-unlist(sapply(res.sig,nrow))
  design<-deg.design
  
  df<-data.frame(type="nasal",results=reslist,n_sig_deg=n_sig_deg,design=design,variable=var, row.names = NULL)
  return(df)
}

if(dir.exists(deg.dir)){
  a<-generate_DEG_summary_table(res.sig,deg.design,var_to_test)
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("dds","nasal","continuous_all_and_dich","res_summary","cellcount+Batch12346",Sys.Date(),".csv",sep="_")))
}

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

######################
## load phenotype data
######################

# make vectors of variables for later use as an input for function 'run_deseq2_DEG_analysis'
# load biomarker phenotype file saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}

# to the phenotype data left join data for sampling date differences 
sampling_date_diff<-"./resources/processed_data/sampling_dates/swab-bal-cbc_differences_in_days.txt"
sampling_date_diff<-if(file.exists(sampling_date_diff)){read.table(sampling_date_diff,row.names = NULL,header = TRUE)}
sampling_date_diff<-sampling_date_diff%>%filter(Comparison=="blood_bal")
colnames(sampling_date_diff)[1:3]<-c("ID","sampling_date_comp","sampling_date_diff_days")

phenotype<-left_join(phenotype,sampling_date_diff,by="ID")
#####################################################################################
## subset phenotype data for which the samples exist for bronchial RNAseq experiments   
#####################################################################################
bexist<-phenotype$SampleID%in%counts.ID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
bphen<-phenotype[bexist,]

# make categorical variables that will be used for DEG based on various thresholds
bphen <- bphen %>% filter(BAL_eos_p >= 0 & BAL_neut_p >= 0) %>%
  mutate(type1 = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 6 ~ "mixed",
                                  BAL_eos_p > 1 & BAL_neut_p <= 6 ~ "eos",
                                  BAL_eos_p <= 1 & BAL_neut_p > 6 ~ "neut",
                                  BAL_eos_p <= 1 & BAL_neut_p <= 6 ~ "pauci"), levels = c("pauci", "neut", "mixed", "eos")),
         type2 = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 6 ~ "mixed",
                                  BAL_eos_p > 1 & BAL_neut_p <= 6 ~ "eos",
                                  BAL_eos_p <= 1 & BAL_neut_p > 6 ~ "neut",
                                  BAL_eos_p <= 1 & BAL_neut_p <= 6 ~ "pauci"), levels = c("neut", "pauci", "mixed", "eos")),
         type3 = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 6 ~ "mixed",
                                  BAL_eos_p > 1 & BAL_neut_p <= 6 ~ "eos",
                                  BAL_eos_p <= 1 & BAL_neut_p > 6 ~ "neut",
                                  BAL_eos_p <= 1 & BAL_neut_p <= 6 ~ "pauci"), levels = c("mixed", "pauci", "neut", "eos")),
         type4 = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 6 ~ "mixed",
                                  BAL_eos_p > 1 & BAL_neut_p <= 6 ~ "eos",
                                  BAL_eos_p <= 1 & BAL_neut_p > 6 ~ "neut",
                                  BAL_eos_p <= 1 & BAL_neut_p <= 6 ~ "pauci"), levels = c("eos", "pauci", "neut", "mixed")),
    )
phen<-bphen

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


## source.cell.log is log-transformed, scaled, and centered cell counts. 
## source.cell is the original BAL/blood cell counts 
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


# these are categorical variables to test using BAL cell counts 

var_dichot_bal<-c("bal_Eos_p_more_1","bal_Eos_p_more_3","bal_ANC_more_0",
                  "bal_ANC_more_5","bal_ANC_more_13","bal_neut_p_more_0",
                  "bal_neut_p_more_2","bal_neut_p_more_5")





############ select variables to test for all non-NA values
var_to_test<-c("type1","type2","type3","type4") # select continuous and categorical variables 
var_to_test_res<-list()
var_to_test_res[[1]]<-c("type_mixed_vs_pauci",
                   "type_neut_vs_pauci",
                   "type_eos_vs_pauci") 
var_to_test_res[[2]]<-c("type_mixed_vs_neut",
                        "type_eos_vs_neut",
                        "type_pauci_vs_neut") 
var_to_test_res[[3]]<-c("type_eos_vs_mixed",
                        "type_neut_vs_mixed",
                        "type_pauci_vs_mixed") 
var_to_test_res[[4]]<-c("type_mixed_vs_eos",
                        "type_neut_vs_eos",
                        "type_pauci_vs_eos") 



# make a list of the phenotype colData that will be used for DESeq2
pi<-lapply(phen[,var_to_test],function(data){a<-!is.na(data);return(a)})
df<-vector("list",length(var_to_test)) # list of data framese used as an input for deseq2. all cell counts
names(df)<-paste(var_to_test,"all",sep="_")
for(i in 1:length(var_to_test)){
  df[[i]]<-phen[pi[[i]],c("SampleID",var_to_test[i], "Batch")]
}
print(sapply(df,dim)[1,])


#################################################################
# bronchial expression ~ log(cell count>=0) + Batch12346
#################################################################
# coldata for DESeq2
df.input<-df

# filtering counts table to remove low expressed genes

## select RNAseq counts
id<-phen$SampleID
cols<-colnames(bronch.counts)%in%id
ct<-bronch.counts[,cols] # First column is actually gene name 

## previous stringet filtering: 
# c2<-filter_low_expressed_genes_method2(ct,round(length(id)*0.1,0))
# ct<-rowgenes_counttable(ct,c2) # low bcounts will be filtered 

count.table<-lapply(df.input,function(df){d<-df; ct<-ct[,colnames(ct)%in%d$SampleID]; return(ct)}) # list of subsetted count table. Each element is a count table with samples for each of the experimental design.



# design: gene expression ~ is_cellcount_threshold + Batch
deg.design<-paste("~",var_to_test,"+ Batch") 
print(deg.design)

# make empty lists for deg analysis, analysis result, and significant results
dds<-vector("list",length=length(var_to_test))
res<-vector("list",length=length(unlist(var_to_test_res)))
res.sig<-vector("list",length=length(unlist(var_to_test_res)))

names(res)<-paste(deg.design,unlist(var_to_test_res),sep="-")
names(res.sig)<-paste(deg.design,unlist(var_to_test_res),sep="-")

# start running DESeq2
# filter genes that have less than 10 counts across all samples 
assay_index<-seq_along(deg.design)
for(i in assay_index){
  dds[[i]]<-run_deseq2_DEG_analysis(count.table[[i]], df.input[[i]], deg.design[i],deg.design[i])
  for(j in seq_along(var_to_test_res)){
    dds_temp<-dds[[i]]
    keep <- rowSums(counts(dds_temp)) >= 10
    dds_temp <- dds_temp[keep,]
    res[[i*j]]<-get_DEG_results(dds_temp, var_to_test_res[i*j])
    res.sig[[i*j]]<-res[[i*j]][which(res[[i*j]]$padj<=0.05),]
    head(res.sig[[i*j]])
    
  }
}


## writing the significant and all results 
deg.folder<-paste("deg","eos-neut-mixed-pauci",Sys.Date(),sep="_")
deg.dir<-file.path("./reports","temporary",deg.folder)
if(!dir.exists(deg.dir)){
  dir.create(deg.dir)
}

if(dir.exists(deg.dir)){
  for(i in seq_along(var_to_test_res)){
    a<-res.sig[[i]]
    b<-res[[i]]
    write.csv(a,row.names=TRUE,file.path(deg.dir,paste("deg","bronch","res_sig",i,var_to_test_res[[i]],Sys.Date(),".csv",sep="_")))
    write.csv(b,row.names=TRUE,file.path(deg.dir,paste("deg","bronch","res_all",i,var_to_test_res[[i]],Sys.Date(),".csv",sep="_")))
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
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("deg","bronch","analysis_input","eos-neut-mixed-pauci_ge0_+Batch",Sys.Date(),".csv",sep="_")))
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
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("dds","bronch","res_summary","cellcount_cont_ge0+Batch",Sys.Date(),".csv",sep="_")))
}
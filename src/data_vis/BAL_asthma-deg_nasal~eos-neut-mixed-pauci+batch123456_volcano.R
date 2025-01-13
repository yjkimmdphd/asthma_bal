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
######################
# load cell count table
######################
countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
genes<-counts[,"SampleID"]

nasal.samples<-phenotype$SampleID[grepl("^N",phenotype$SampleID)]
# Filter samples that do not contain "F"
filtered_nasal_samples <- nasal.samples[!grepl("F", nasal.samples)]
# Display the filtered samples
nasal.counts<-counts[,filtered_nasal_samples]
rownames(nasal.counts)<-genes

head(nasal.counts)
counts.ID<-colnames(nasal.counts)


#####################################################################################
## subset phenotype data for which the samples exist for bronchial RNAseq experiments   
#####################################################################################
nexist<-phenotype$SampleID%in%counts.ID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
nphen<-phenotype[nexist,]

# make categorical variables that will be used for DEG based on various thresholds
## there are 4 levels to the categorical variable being tested
## based on BAL Eos% threshold of 1% and Neut% threshold of 6, cells are either mixed, eos-dominant (Eos), neut-dominant (neut), or paucigranulocytic (pauci)
## the models used for the DESeq2 are ~ type + batch
## the categorical variables are factors, and the levels are set differently for each type (i.e., type1, type2, type3, type4) for different contrasts 
## comp 1 will compare vs pauci
## comp 2 will compare vs neut
## comp 3 will compare vs mixed
## comp 4 will compare vs eos 
nphen <- nphen %>% filter(BAL_eos_p >= 0 & BAL_neut_p >= 0) %>%
  mutate(comp1 = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 4 ~ "mixed",
                                  BAL_eos_p > 1 & BAL_neut_p <= 4 ~ "eos",
                                  BAL_eos_p <= 1 & BAL_neut_p > 4 ~ "neut",
                                  BAL_eos_p <= 1 & BAL_neut_p <= 4 ~ "pauci"), levels = c("pauci", "neut", "mixed", "eos")),
         comp2 = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 4 ~ "mixed",
                                  BAL_eos_p > 1 & BAL_neut_p <= 4 ~ "eos",
                                  BAL_eos_p <= 1 & BAL_neut_p > 4 ~ "neut",
                                  BAL_eos_p <= 1 & BAL_neut_p <= 4 ~ "pauci"), levels = c("neut", "pauci", "mixed", "eos")),
         comp3 = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 4 ~ "mixed",
                                  BAL_eos_p > 1 & BAL_neut_p <= 4 ~ "eos",
                                  BAL_eos_p <= 1 & BAL_neut_p > 4 ~ "neut",
                                  BAL_eos_p <= 1 & BAL_neut_p <= 4 ~ "pauci"), levels = c("mixed", "pauci", "neut", "eos")),
         comp4 = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 4 ~ "mixed",
                                  BAL_eos_p > 1 & BAL_neut_p <= 4 ~ "eos",
                                  BAL_eos_p <= 1 & BAL_neut_p > 4 ~ "neut",
                                  BAL_eos_p <= 1 & BAL_neut_p <= 4 ~ "pauci"), levels = c("eos", "pauci", "neut", "mixed")),
  )
phen<-nphen

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


############ select variables to test for all non-NA values
var_to_test<-c("comp1","comp2","comp3","comp4") # select continuous and categorical variables 
var_to_test_res<-list()
var_to_test_res[[1]]<-c("comp1_mixed_vs_pauci",
                        "comp1_neut_vs_pauci",
                        "comp1_eos_vs_pauci") 
var_to_test_res[[2]]<-c("comp2_eos_vs_neut",
                        "comp2_mixed_vs_neut",
                        "comp2_pauci_vs_neut") 
var_to_test_res[[3]]<-c("comp3_eos_vs_mixed",
                        "comp3_neut_vs_mixed",
                        "comp3_pauci_vs_mixed") 
var_to_test_res[[4]]<-c("comp4_mixed_vs_eos",
                        "comp4_neut_vs_eos",
                        "comp4_pauci_vs_eos") 



# make a list of the phenotype colData that will be used for DESeq2
pi<-lapply(phen[,var_to_test],function(data){a<-!is.na(data);return(a)})
df<-vector("list",length(var_to_test)) # list of data framese used as an input for deseq2. all cell counts
names(df)<-paste(var_to_test)
for(i in 1:length(var_to_test)){
  df[[i]]<-phen[pi[[i]],c("SampleID",var_to_test[i], "Batch")]
}
print(sapply(df,dim)[1,])


#################################################################
# nasal expression ~ log(cell count>=0) + Batch12346
#################################################################
# coldata for DESeq2
df.input<-df

# select RNAseq counts

id<-phen$SampleID
cols<-colnames(nasal.counts)%in%id
ct<-nasal.counts[,cols] # First column is actually gene name 

count.table<-lapply(df.input,function(df){d<-df; ct<-ct[,colnames(ct)%in%d$SampleID]; return(ct)}) # list of subsetted count table. Each element is a count table containing samples for each part of the experimental design.



# design: gene expression ~ is_cellcount_threshold + Batch
deg.design<-paste("~",var_to_test,"+ Batch") 
print(deg.design)

# make empty lists for deg analysis, analysis result, and significant results
dds<-vector("list",length=length(var_to_test))
# Initialize the list
res <- list()

assay_index<-seq_along(deg.design)

# Loop over the indices in assay_index

for(i in assay_index){
  # Initialize an inner list for res[[i]]
  res[[i]] <- list()
  
  # Loop to create 3 elements within each res[[i]]
  for(j in 1:3){
    res[[i]][[j]] <- NA
  }
}

res.sig<-list()
for(i in assay_index){
  # Initialize an inner list for res[[i]]
  res.sig[[i]] <- list()
  
  # Loop to create 3 elements within each res[[i]]
  for(j in 1:3){
    res.sig[[i]][[j]] <- NA
  }
}


# start running DESeq2
# filter genes that have less than 10 counts across all samples 
assay_index<-seq_along(deg.design)
for(i in assay_index){
  dds[[i]]<-run_deseq2_DEG_analysis(count.table[[i]], df.input[[i]], deg.design[i],deg.design[i])
}


for (i in assay_index) {
  for (j in 1:3) {
    dds_temp <- dds[[i]]
    keep <- rowSums(counts(dds_temp)) >= 10
    dds_temp <- dds_temp[keep, ]
    res[[i]][[j]] <- get_DEG_results(dds_temp, var_to_test_res[[i]][j])
    res.sig[[i]][[j]] <- res[[i]][[j]][which(res[[i]][[j]]$padj <= 0.05), ]
    head(res.sig[[i]][[j]])
    
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
    for( j in 1:3){
      a<-res.sig[[i]][[j]]
      b<-res[[i]][[j]]
      write.csv(a,row.names=TRUE,file.path(deg.dir,paste("deg","nasal","res_sig",i,j,var_to_test_res[[i]][j],Sys.Date(),".csv",sep="_")))
      write.csv(b,row.names=TRUE,file.path(deg.dir,paste("deg","nasal","res_all",i,j,var_to_test_res[[i]][j],Sys.Date(),".csv",sep="_")))
    }
  }
}
## summarize the data input 

generate_DEG_input_summary_table<-function(original_ct,filtered_ct,dds,res,des){
  filter_method<-"at least 10 counts across samples, then DESeq2 automatic filter"
  n_filtered_genes<-paste("analyzed n_genes:", mean(sapply(filtered_ct,nrow)),",","filtered n_genes:",mean(sapply(original_ct,nrow))-mean(sapply(filtered_ct,nrow)))
  samples<-sapply(dds, function(d){colData(d)$SampleID%>%paste(collapse = ",")})
  dds<-paste("dds",1:length(dds),sep="")
  results<-paste("res",1:length(res),sep="")
  design<-des
  df<-data.frame(dds=dds,results=results,design=design,samples=samples,filter_method=filter_method,n_filtered_genes=n_filtered_genes)
  return(df)
}


if(dir.exists(deg.dir)){
  a<-generate_DEG_input_summary_table(nasal.counts[,cols],ct,dds,res,deg.design)
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("deg","nasal","analysis_input","eos-neut-mixed-pauci_ge0_+Batch",Sys.Date(),".csv",sep="_")))
}

## summary table of the DEG analysis
generate_DEG_summary_table<-function(results_significant,deg_design,variable){
  res.sig<-results_significant
  deg.design<-deg_design
  var<-variable
  
  reslist<-paste("res.sig",1:length(res.sig),sep="")
  n_sig_deg<-unlist(sapply(res.sig,nrow))
  design<-deg.design
  
  df<-data.frame(type="nasal",results=reslist,n_sig_deg=n_sig_deg,design=design,variable=var, row.names = NULL)
  return(df)
}

if(dir.exists(deg.dir)){
  a<-generate_DEG_summary_table(res.sig,deg.design,var_to_test)
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("dds","nasal","res_summary","cellcount_cont_ge0+Batch",Sys.Date(),".csv",sep="_")))
}

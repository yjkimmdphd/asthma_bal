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

# asthma biomarker phenotype file saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/Nasal_Biomarkers_BAL_transformed.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
numeric_id<-phenotype$ID
phenotype$ID<-sprintf("%03d",numeric_id) # adds padded zeros in front of the subject ID numbers
phenotype<-phenotype[,-28]
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
source.cell<-c("BAL_eos_ct",
               "BAL_eos_p",
               "BAL_neut_ct",
               "BAL_neut_p",
               "BAL_wbc",
               "blood_eos",
               "blood_eos_p",
               "blood_neut",
               "blood_neut_p",
               "blood_wbc")
deg.design<-paste("~",source.cell.log,"+ Batch")

#######################################################################
## find subject assignment ID with nasal and bronchial cell RNAseq data  
#######################################################################
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

## define function join_phenotype_batch_info. 
join_phenotype_batch_info<-function(p,b){ 
  #p is phenotype table. b is batch info table. Factorize the batch info. 
  table<-left_join(p,b, by="SampleID")
  table$Batch<-factor(table$Batch, levels=unique(table$Batch))
  return(table)
}
nphen<-join_phenotype_batch_info(nphen,batch.info)
bphen<-join_phenotype_batch_info(bphen,batch.info)

# scale the cell count information 
# Mike Love: I'd recommend however to center and scale the numeric covariates for model fitting improvement.
# scale the columns named with source.cell.log
nphen<-mutate_at(nphen,vars(all_of(source.cell.log)),scale) # scales and mutates all log-transformed cell counts 
bphen<-mutate_at(bphen,vars(all_of(source.cell.log)),scale) # scales and mutates all log-transformed cell counts 


# decide which analysis to perform, then set the phenotype data as phen
# will perform rnaseq on bronchial samples so 'bphen' should be 'phen'
phen<-bphen

# make a list in which each element is a numeric vector of row numbers corresponding to rows of interests in phenotype data matching the following condition: 
## nasal RNAseq samples with cell counts at least 0
source.cell.all<-sapply(bphen[,source.cell],function(p){a<-p; b<-which(a>=0); return(b)})

# define function 'subset_phenotype' 
subset_phenotype<-function(phenotype_data, source_cell_list, source_cell, sample_id, v1, v2){## make function that subsets the phenotype data referencing the rownames saved in the list data made above using the following function
  ## the inputs are the phenotype_data to subset, the list file (source_cell_df), source (BAL/serum) and cell type (eos,neut,wbc) to subset with, sample ID, and independent variables that will be used for downstream DEG analysis (v1, v2) 
  l<-source_cell_list
  sc<-source_cell
  p<-phenotype_data #phenotype data for nasal/bronchial RNAseq samples 
  rn<-as.numeric(unlist(l[sc]))
  sp<-p[rn,c(sample_id,v1,v2)]
  return(sp)
}

##df with bal eos>0
bal.all.e.c<-subset_phenotype(phen,source.cell.all,"BAL_eos_ct","SampleID","BAL_eos_ct_log", "Batch")
bal.all.e.p<-subset_phenotype(phen,source.cell.all,"BAL_eos_p","SampleID","BAL_eos_p_log", "Batch")
##df with bal neut>0
bal.all.n.c<-subset_phenotype(phen,source.cell.all,"BAL_neut_ct","SampleID","BAL_neut_ct_log","Batch")
bal.all.n.p<-subset_phenotype(phen,source.cell.all,"BAL_neut_p","SampleID","BAL_neut_p_log","Batch")
##df with bal WBC>0
bal.all.w.c<-subset_phenotype(phen,source.cell.all,"BAL_wbc","SampleID","BAL_wbc_log","Batch")

##df with serum_Eos>0
ser.all.e.c<-subset_phenotype(phen,source.cell.all,"blood_eos","SampleID","blood_eos_log", "Batch")
ser.all.e.p<-subset_phenotype(phen,source.cell.all,"blood_eos_p","SampleID","blood_eos_p_log","Batch")
##df with serum_Neut>0
ser.all.n.c<-subset_phenotype(phen,source.cell.all,"blood_neut","SampleID","blood_neut_log","Batch")
ser.all.n.p<-subset_phenotype(phen,source.cell.all,"blood_neut_p","SampleID","blood_neut_p_log","Batch")
##df with serum WBC>0
ser.all.w.c<-subset_phenotype(phen,source.cell.all,"blood_wbc","SampleID","blood_wbc_log","Batch")


## make dictionary table of df names and corresponding source cell counts
df.all.cell<-data.frame(source_cell = paste(source.cell,">=0"),
                        df_name=c("bal.all.e.c",
                                  "bal.all.e.p",
                                  "bal.all.n.c",
                                  "bal.all.n.p",
                                  "bal.all.w.c",
                                  "ser.all.e.c",
                                  "ser.all.e.p",
                                  "ser.all.n.c",
                                  "ser.all.n.p",
                                  "ser.all.w.c"),
                        row.names = NULL)

## check if any df has NA 
sapply(df.all.cell$df_name,function(df){a<-get(df);is.na(a[,2])})%>%lapply(sum)

## check # of samples for each bronchial RNAseq sample phenotype 
n_sample_pos<-unlist(sapply(df.all.cell$df_name,function(df){a<-get(df);nrow(a)})%>%lapply(sum))
df.all.cell<-cbind(df.all.cell,n_sample_pos)


###########################################
## filtering counts table to remove low expressed genes
###########################################

# select just the bronchial RNAseq counts
id<-counts.ID%in%bID
cols<-colnames(counts)[id]
bcounts<-counts[,cols] # First column is actually gene name 
genes<-counts$SampleID
rownames(bcounts)<-genes

# Filter bcounts (readcount table for bronchial sample
## filter method 1:
### Remove genes with ≤ 2 counts per million in >10% of samples (>4 samples) were removed to reduce noise on low counts and low abundance genes (PMID: 37730635)
filter_low_expressed_genes_method1<-function(readcounts, cutoff, n_sample){
  x<-readcounts
  cpm0 <-cpm(readcounts)
  drop.genes<-which(rowSums(cpm0<=cutoff)>n_sample)
  x<-x[-drop.genes,]
  print(paste("# of dropped genes:",length(drop.genes)))
  return(x)
}
c1<-filter_low_expressed_genes_method1(bcounts, 2, 4)

## filter method 2: use TMM normalized lcpm as a cutoff point
### x will be the TMM normalized count
filter_low_expressed_genes_method2<-function(readcounts, n_sample){
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
c2<-filter_low_expressed_genes_method2(bcounts,4)

#######################################
# run the DEG for continuous predictors
#######################################

# count data = ncounts (readcount matrix)
# coldata = bphen (phenotype data)
# design = statistical design eq

run_deseq2_DEG_analysis<-function(countdata,coldata,design,resultname){
  print(paste("design",resultname, sep = ": "))
  dds<-DESeqDataSetFromMatrix(countData = countdata,colData=coldata, design=as.formula(design))
  dds<-DESeq(dds)
  print("resultnames")
  print(resultsNames(dds))
  return(dds)
}

get_DEG_results<-function(dds,resultname){
  res<-results(dds, name=resultname) # use source.cell.log element for resultname.
  res <- res[order(res$padj),]
  return(res)
}

filter_readcounts<-function(nc,ncf){
  n<-nc
  gn<-rownames(nc)
  gnf<-rownames(ncf)
  c<-n[gn%in%gnf,]
  return(c)
}

bcounts<-filter_readcounts(bcounts,c2) # low bcounts will be filtered using method 2: use TMM normalized lcpm as a cutoff point

subset_readcounts<-function(nid,np){ # function to subset bronchial RNAseq readcounts with bronchial RNAseq samples for each analysis
  np<-np
  nID<-nid
  
  ns<-nID%in%np$SampleID
  ns<-nID[ns]
  nc<-select(bcounts,all_of(ns)) #'bcount' is a hard-coded dataframe variable 
  return(nc)
}

# DEG with design "~  BAL_Eos_ct_log  + Batch"
df<-get(df.all.cell[1,2])
count.table<-subset_readcounts(bID,df)

dds1<-run_deseq2_DEG_analysis(count.table, df, deg.design[1],source.cell.log[1])
res1<-get_DEG_results(dds1, source.cell.log[1])
res1.sig<-res1[which(res1$padj<=0.05),]

# DEG with design "~  BAL_Eos_perc_log  + Batch"
df<-get(df.all.cell[2,2])
count.table<-subset_readcounts(bID,df)

dds2<-run_deseq2_DEG_analysis(count.table, df, deg.design[2],source.cell.log[2])
res2<-get_DEG_results(dds2, source.cell.log[2])
res2.sig<-res2[which(res2$padj<=0.05),]

# DEG with design "~  BAL_neut_ct_log10  + Batch"  
df<-get(df.all.cell[3,2])
count.table<-subset_readcounts(bID,df)

dds3<-run_deseq2_DEG_analysis(count.table, df, deg.design[3],source.cell.log[3])
res3<-get_DEG_results(dds3, source.cell.log[3])
res3.sig<-res3[which(res3$padj<=0.05),]

# DEG with design "~  BAL_neut_perc_log  + Batch"   
df<-get(df.all.cell[4,2])
count.table<-subset_readcounts(bID,df)

dds4<-run_deseq2_DEG_analysis(count.table, df, deg.design[4],source.cell.log[4])
res4<-get_DEG_results(dds4, source.cell.log[4])
res4.sig<-res4[which(res4$padj<=0.05),]

# DEG with design "~  BAL_WBC_log  + Batch"       
df<-get(df.all.cell[5,2])
count.table<-subset_readcounts(bID,df)

dds5<-run_deseq2_DEG_analysis(count.table, df, deg.design[5],source.cell.log[5])
res5<-get_DEG_results(dds5, source.cell.log[5])
res5.sig<-res5[which(res5$padj<=0.05),]

# DEG with design "~  serum_Eos_log10  + Batch"    
df<-get(df.all.cell[6,2])
count.table<-subset_readcounts(bID,df)

dds6<-run_deseq2_DEG_analysis(count.table, df, deg.design[6],source.cell.log[6])
res6<-get_DEG_results(dds6, source.cell.log[6])
res6.sig<-res6[which(res6$padj<=0.05),]


# DEG with design "~  serum_Eos_perc_log  + Batch"  
df<-get(df.all.cell[7,2])
count.table<-subset_readcounts(bID,df)

dds7<-run_deseq2_DEG_analysis(count.table, df, deg.design[7],source.cell.log[7])
res7<-get_DEG_results(dds7, source.cell.log[7])
res7.sig<-res7[which(res7$padj<=0.05),]

# DEG with design "~  serum_Neut_log10  + Batch"   
df<-get(df.all.cell[8,2])
count.table<-subset_readcounts(bID,df)

dds8<-run_deseq2_DEG_analysis(count.table, df, deg.design[8],source.cell.log[8])
res8<-get_DEG_results(dds8, source.cell.log[8])
res8.sig<-res8[which(res8$padj<=0.05),]

# DEG with design "~  serum_Neut_perc_log  + Batch"
df<-get(df.all.cell[9,2])
count.table<-subset_readcounts(bID,df)

dds9<-run_deseq2_DEG_analysis(count.table, df, deg.design[9],source.cell.log[9])
res9<-get_DEG_results(dds9, source.cell.log[9])
res9.sig<-res9[which(res9$padj<=0.05),]
# DEG with design "~  serum_WBC_log10  + Batch"    
df<-get(df.all.cell[10,2])
count.table<-subset_readcounts(bID,df)

dds10<-run_deseq2_DEG_analysis(count.table, df, deg.design[10],source.cell.log[10])
res10<-get_DEG_results(dds10, source.cell.log[10])
res10.sig<-res10[which(res10$padj<=0.05),]

## writing the results 
deg.folder<-paste("deg",Sys.Date(),sep="_")
deg.dir<-file.path("./reports",deg.folder)
dir.create(deg.dir)
if(dir.exists(deg.dir)){
  for(i in 1:10){
    a<-get(paste0("res",i,".sig"))
    write.csv(a,row.names=TRUE,file.path(deg.dir,paste0("deg_","bronch_allcell_","res",i,"_","design_",Sys.Date(),".csv"))) }
}

# summarize the data input 

generate_DEG_input_summary_table<-function(){
  filter_method<-"TMM normalized LCPM cutoff"
  n_filtered_genes<-paste("analyzed n_genes:", nrow(bcounts),",","filtered n_genes:",nrow(counts)-nrow(bcounts))
  samples<-c(colData(dds1)$SampleID%>%paste(collapse = ","),
             colData(dds2)$SampleID%>%paste(collapse = ","),
             colData(dds3)$SampleID%>%paste(collapse = ","),
             colData(dds4)$SampleID%>%paste(collapse = ","),
             colData(dds5)$SampleID%>%paste(collapse = ","),
             colData(dds6)$SampleID%>%paste(collapse = ","),
             colData(dds7)$SampleID%>%paste(collapse = ","),
             colData(dds8)$SampleID%>%paste(collapse = ","),
             colData(dds9)$SampleID%>%paste(collapse = ","),
             colData(dds10)$SampleID%>%paste(collapse = ","))
  dds<-paste("dds",1:10,sep="")
  results<-paste("res",1:10,sep="")
  design<-deg.design
  df<-data.frame(dds=dds,results=results,design=design,samples=samples,filter_method=filter_method,n_filtered_genes=n_filtered_genes)
  df<-cbind(df,df.all.cell)
  return(df)
}

if(dir.exists(deg.dir)){
  a<-generate_DEG_input_summary_table()
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("deg_bronch_allcell_analysis_input_",Sys.Date(),".csv",sep="")))
}

# summary table of the DEG analysis
generate_DEG_summary_table<-function(){
  reslist<-paste("res",1:10,".sig",sep="")
  n_sig_deg<-sapply(reslist,function(a){nrow(get(a))})
  design<-deg.design
  source_cell<-source.cell
  df<-data.frame(results=reslist,n_sig_deg,design=design,source_cell=df.all.cell$source_cell,row.names = NULL)
  return(df)
}

if(dir.exists(deg.dir)){
  a<-generate_DEG_summary_table()
  write.csv(a,row.names=FALSE,file.path(deg.dir,paste("dds_bronch_allcell_res_summary_",Sys.Date(),".csv",sep="")))
}


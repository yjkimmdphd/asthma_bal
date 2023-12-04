#########################
# cleaning up patient phenotype data table for bronchial_BAL_rnaseq.R
#########################
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

# asthma biomarker phenotype file, nasal, saved in  'phenotype'
phenotype<-file.path("./resources/processed_data/asthma-phenotype-Rnaseq-2023-12-01.csv")
phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
org.numb<-phenotype$subjID
phenotype$subjID<-sprintf("%03d",org.numb) # adds padded zeros in front of the subject ID numbers

# make vectors of variables for later use 
source.cell.log<-c("BAL_Eos_ct_log",
                   "BAL_Eos_perc_log",
                   "BAL_neut_ct_log10",
                   "BAL_neut_perc_log",
                   "BAL_WBC_log",
                   "serum_Eos_log10",
                   "serum_Eos_perc_log",
                   "serum_Neut_log10",
                   "serum_Neut_perc_log",
                   "serum_WBC_log10")
source.cell<-c("BAL_Eos_ct",
               "BAL_Eos_perc",
               "BAL_neut_ct",
               "BAL_neut_perc",
               "BAL_WBC",
               "serum_Eos",
               "serum_Eos_perc",
               "serum_Neut",
               "serum_Neut_perc",
               "serum_WBC")
deg.design<-paste("~",source.cell.log,"+ Batch")
# samples to exclude in certain analysis due to absent cell count values
# used for filtering phenotype df and readcount matrix in DEseq
exclude.bronch<-read.csv("./resources/processed_data/bronch_samples_excluded.csv")

# scale the cell count information 
# Mike Love: I'd recommend however to center and scale the numeric covariates for model fitting improvement.
# scale the columns named with source.cell.log
phenotype.scale<-phenotype[,c("subjID", source.cell.log)] # subsets subject ID and the source.cell that will be used for DEG analysis
phenotype.scale<-mutate_at(phenotype.scale,vars(all_of(source.cell.log)),scale) # scales and mutates all cell counts 

#######################################################################
## find subject assignment ID with nasal and bronchial cell RNAseq data  
#######################################################################
nID<-paste0("N",phenotype$subjID) # N*** indicates bronchial sample ID, sequence data is not available for all as of 2023-10-04
bID<-paste0("B",phenotype$subjID) # B*** indicates bronchial sample ID, sequence data is not available for all as of 2023-10-04
nexist<-counts.ID%in%nID # find which subjects s/p nasal and had bronchial sample RNAseq completed. Nasal samples in batch 1-4 only sequenced to subjID 337
bexist<-counts.ID%in%bID # find which subjects s/p BAL and had bronchial sample RNAseq completed 
nsample<-counts.ID[nexist] # nasal sample ID in the readcount matrix (batch 1-4) that has BAL phenotype data
bsample<-counts.ID[bexist] # bronchial sample ID in the readcount matrix (batch 1-4) that has BAL phenotype data

nphen<-phenotype.scale[phenotype$subjID%in%substring(nsample,2),] # phenotype table with nsample
bphen<-phenotype.scale[phenotype$subjID%in%substring(bsample,2),] # phenotype table with bsample

nphen<-mutate(nphen, SampleID=nsample)%>%relocate(SampleID, .before=1)
bphen<-mutate(bphen, SampleID=bsample)%>%relocate(SampleID, .before=1)

# left join batch info table with nasal/bronchial phenotype table  

# get batch information
batch<-file.path("./resources/processed_data/asthma_nasal_bronchial_batch1234_info.txt")
batch.info<-if(file.exists(batch)){read.delim(batch)}

## define function joinBatch. p is phenotype table. b is batch info table. Factorize the batch info. 
joinBatch<-function(p,b){
  table<-left_join(p,b, by="SampleID")
  table$Batch<-factor(table$Batch, levels=unique(table$Batch))
  return(table)
}

nphen<-joinBatch(nphen,batch.info)
bphen<-joinBatch(bphen,batch.info)


###########################################
## filtering counts table to remove low expressed genes
###########################################

# select just the bronchial RNAseq counts
b<-counts.ID%in%bID
bcols<-colnames(counts)[b]
bcounts<-counts[,bcols] # First column is actually gene name 
genes<-counts$SampleID
rownames(bcounts)<-genes

# Filter bcounts (readcount table for bronchial sample
## filter method 1:
### Remove genes with ≤ 2 counts per million in >10% of samples (>4 samples) were removed to reduce noise on low counts and low abundance genes (PMID: 37730635)
lowReadFilter1<-function(readcounts, cutoff, n_sample){
  x<-readcounts
  cpm0 <-cpm(readcounts)
  drop.genes<-which(rowSums(cpm0<=cutoff)>n_sample)
  x<-x[-drop.genes,]
  print(paste("# of dropped genes:",length(drop.genes)))
  return(x)
}
bc1<-lowReadFilter1(bcounts, 2, 4)

## filter method 2: use TMM normalized lcpm as a cutoff point
### x will be the TMM normalized count
lowReadFilter2<-function(readcounts, n_sample){
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
bc2<-lowReadFilter2(bcounts,4)

#######################################
# run the DEG for continuous predictors
#######################################

# count data = bcounts (readcount matrix)
# coldata = bphen (phenotype data)
# design = statistical design eq

deseq2DEG<-function(countdata,coldata,design,resultname){
  print(paste("design",resultname, collapse = ":"))
  dds<-DESeqDataSetFromMatrix(countData = countdata,colData=coldata, design=as.formula(design))
  dds<-DESeq(dds)
  print("resultnames")
  print(resultsNames(dds))
  return(dds)
}

degRes<-function(dds,resultname){
  res<-results(dds, name=resultname) # use source.cell.log element for resultname.
  res <- res[order(res$padj),]
  return(res)
}

bcounts<-bc1 # bcounts will be using readcount method 1: cpm <=2 in 10% of samples filtered

# DEG with design "~  BAL_Eos_ct_log  + Batch"
# no excluded samples
dds1<-deseq2DEG(bcounts, bphen, deg.design[1],source.cell.log[1])
res1<-degRes(dds1, source.cell.log[1])
res1.sig<-res1[which(res1$padj<=0.05),]

# DEG with design "~  BAL_Eos_perc_log  + Batch"
# no excluded samples
dds2<-deseq2DEG(bcounts, bphen, deg.design[2],source.cell.log[2])
res2<-degRes(dds2, source.cell.log[2])
res2.sig<-res2[which(res2$padj<=0.05),]

# DEG with design "~  BAL_neut_ct_log10  + Batch"  
# Excluded sample B332 due to missing BAL neut info
dds3<-deseq2DEG(select(bcounts,-"B332"), filter(bphen,SampleID!="B332"), deg.design[3],source.cell.log[3])
res3<-degRes(dds3, source.cell.log[3])
res3.sig<-res3[which(res3$padj<=0.05),]

# DEG with design "~  BAL_neut_perc_log  + Batch"   
# Excluded sample B332 due to missing BAL neut info
dds4<-deseq2DEG(select(bcounts,-"B332"), filter(bphen,SampleID!="B332"), deg.design[4],source.cell.log[4])
res4<-degRes(dds4, source.cell.log[4])
res4.sig<-res4[which(res4$padj<=0.05),]

# DEG with design "~  BAL_WBC_log  + Batch"       
# no excluded samples
dds5<-deseq2DEG(bcounts, bphen, deg.design[5],source.cell.log[5])
res5<-degRes(dds5, source.cell.log[5])
res5.sig<-res5[which(res5$padj<=0.05),]

# DEG with design "~  serum_Eos_log10  + Batch"    
esc1<-filter(exclude.bronch,source_cell=="serum_Eos")%>%select(SampleID)%>%unlist%>%as.character # samples to exclude due to absent source cell info 
k<-!(bphen$SampleID%in%esc1) 
dds6<-deseq2DEG(bcounts[,k], bphen[k,], deg.design[6],source.cell.log[6])
res6<-degRes(dds6, source.cell.log[6])
res6.sig<-res6[which(res6$padj<=0.05),]

# DEG with design "~  serum_Eos_perc_log  + Batch"  
esc1<-filter(exclude.bronch,source_cell=="serum_Eos")%>%select(SampleID)%>%unlist%>%as.character # samples to exclude due to absent source cell info 
k<-!(bphen$SampleID%in%esc1) 
dds7<-deseq2DEG(bcounts[,k], bphen[k,], deg.design[7],source.cell.log[7])
res7<-degRes(dds7, source.cell.log[7])
res7.sig<-res7[which(res7$padj<=0.05),]

# DEG with design "~  serum_Neut_log10  + Batch"   
esc2<-filter(exclude.bronch,source_cell=="serum_Neut")%>%select(SampleID)%>%unlist%>%as.character # samples to exclude due to absent source cell info 
k<-!(bphen$SampleID%in%esc2) 
dds8<-deseq2DEG(bcounts[,k], bphen[k,], deg.design[8],source.cell.log[8])
res8<-degRes(dds8, source.cell.log[8])
res8.sig<-res8[which(res8$padj<=0.05),]

# DEG with design "~  serum_Neut_perc_log  + Batch"
esc2<-filter(exclude.bronch,source_cell=="serum_Neut")%>%select(SampleID)%>%unlist%>%as.character # samples to exclude due to absent source cell info 
k<-!(bphen$SampleID%in%esc2) 
dds9<-deseq2DEG(bcounts[,k], bphen[k,], deg.design[9],source.cell.log[9])
res9<-degRes(dds9, source.cell.log[9])
res9.sig<-res9[which(res9$padj<=0.05),]

# DEG with design "~  serum_WBC_log10  + Batch"    
esc3<-filter(exclude.bronch,source_cell=="serum_WBC")%>%select(SampleID)%>%unlist%>%as.character # samples to exclude due to absent source cell info 
k<-!(bphen$SampleID%in%esc3)
dds10<-deseq2DEG(bcounts[,k], bphen[k,], deg.design[10],source.cell.log[10])
res10<-degRes(dds10, source.cell.log[10])
res10.sig<-res10[which(res10$padj<=0.05),]


## writing the results *** in the future
# for(i in c(cont.var,zero.var)){
#  a<-get(paste0("res",i))
#  write.csv(a[[1]],row.names=TRUE,file.path(getwd(),"output",paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")))
#  res.table[i,"output"]<-paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")
# }



for(i in cont.var){
  print(df.deseq2input[i,4])
  assign(paste0("res",i),deseq2DEG(df.deseq2input[i,1],df.deseq2input[i,2],df.deseq2input[i,3],df.deseq2input[i,4]))
}


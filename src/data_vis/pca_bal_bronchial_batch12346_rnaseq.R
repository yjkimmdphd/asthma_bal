##
# perform PCA analysis with normalized lcpm RNA seq gene counts of nasal/bronchial data with edgeR and limma
## 
library(dplyr)
library(limma)
library(Glimma)
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
phenotype<-file.path("./resources/processed_data/Nasal_Biomarkers_BAL_transformed_with_raceinfo.csv")
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

calculate_Cell_Count_Quantile<-function(normalized_cellcounts){
  p<-normalized_cellcounts
  q<-quantile(p, probs = c(0.25, 0.5, 0.75), na.rm=TRUE)
  return(q)
} # find quantile cutoffs of normalized cell counts 
cellCountRange1<-function(cell,q1,q2){
  p<-phen[,cell]
  q1p<-cell.count.quantiles[q1/25,cell]
  q2p<-cell.count.quantiles[q2/25,cell]
  range<-which(p>q1p & p<=q2p)
  return(range)
} # function to set range of cell count between 25th-75th percentile 
cellCountRange2<-function(cell,q1){
  p<-phen[,cell]
  q1p<-cell.count.quantiles[q1/25,cell]
  range<-which(p>q1p)
  return(range)
} # function to set range of cell count for PCA >75th percentile
cellCountRange3<-function(cell,q1){
  p<-phen[,cell]
  q1p<-cell.count.quantiles[q1/25,cell]
  range<-which(p<=q1p)
  return(range)
} # function to set range of cell count for PCA >75th percentile

###########################################---------------------###########################################
# sample type: nasal-bronchial
# gene filter: all genes
###########################################---------------------###########################################
phen<-if(exists("nphen")){
  rbind(nphen,bphen)
  }else{
    bphen
  }

##
# Normalising gene expression distributions
##

# select just the nasal RNAseq counts
sid<-phen$SampleID

if(exists("nphen")){
  print("check counts matrix")
}else{
  counts<-bronch.counts
}

ct<-counts[,sid] # First column is actually gene name 
genes<-counts$SampleID
rownames(ct)<-genes

## use TMM normalized lcpm as a cutoff point ()
ct_filtered<-filter_low_expressed_genes_method2(ct,9) 


# normalize counts with TMM
norm.factor<-calcNormFactors(ct, method = "TMM")

sample.size<-length(colnames(ct))
for(i in 1:sample.size){
  ct[,i]<-ct[,i]/norm.factor[i]
}

# calculate lcpm based on TMM normalized counts 
lcpm.ct<-cpm(ct,log=TRUE)


###
#MDS by groups
###

#MDS color scheme for nasal 
## color coding based on each factor levels of sex, age, cell counts, and batch sourced from the following code:
## color code covariates 

# make color matrix. all entry is black at baseline 
col<-data.frame(matrix(nrow=sample.size,ncol=16))
colnames(col)<-c(source.cell.log,"nasal_bronch","batch","age","sex","Race","ethnicity")
rownames(col)<-phen$SampleID
for(i in 1:ncol(col)) {col[,i]<-rep("black",sample.size)}

# color for cell counts. 
## grey if cell count = 0
## red, green, blue depending on quantile cutoff among >0 cell counts
### find sample ID of samples with zero cell counts 
source.cell.zero<-sapply(phen[,source.cell],function(p){a<-p; b<-which(a==0); return(b)})
zerosample<-vector("list", length = 10)
names(zerosample)<-source.cell.log

for(i in 1:10){
  zerosample[[i]]<-sid[source.cell.zero[[i]]]
}

### find quantiles for each of the normalized cell counts
source.cell.pos<-sapply(phen[,source.cell],function(p){a<-p; b<-which(a>0); return(b)})

### scaled cell counts list of each cell count > 0
posphen<-vector("list",length=10) 
names(posphen)<-source.cell.log
for(i in 1:10){
  posphen[[i]]<- phen[source.cell.pos[[i]],source.cell.log[i]]
}

### quantile cutoff for each scaled cell count > 0
cell.count.quantiles<-sapply(posphen,calculate_Cell_Count_Quantile) 

## set color 

for(i in 1:10){
  col[cellCountRange1(source.cell.log[i],25,50),source.cell.log[i]]<-"red" # red if pos cell count between 25th - 50th percentile
  col[cellCountRange1(source.cell.log[i],50,75),source.cell.log[i]]<-"green" # green if pos cell count between 50th - 75th percentile
  col[cellCountRange2(source.cell.log[i],75),source.cell.log[i]]<-"blue" # blue if pos cell count above 75th percentile
  col[cellCountRange3(source.cell.log[i],25),source.cell.log[i]]<-"black" # black if pos cell count below 25th percentile or unavailable
  col[source.cell.zero[[i]],source.cell.log[i]]<-"grey" # grey if cell count is 0
}

# colors by batch
# batch1/2/3/4/6 is red/gree/blue/orange/black, respectively
col[which(phen$Batch=="batch1"),"batch"]<-"red" 
col[which(phen$Batch=="batch2"),"batch"]<-"green"
col[which(phen$Batch=="batch3"),"batch"]<-"blue"
col[which(phen$Batch=="batch4"),"batch"]<-"orange"

# color by sample type nasal_bronch
col[which(phen$Type=="Nasal"),"nasal_bronch"]<-"red" 
col[which(phen$Type=="Bronchial"),"nasal_bronch"]<-"blue"

# colors by age
col[which(phen$Age<=7),"age"]<-"red" # age less than or equal to 7 are red
col[which(phen$Age>7 & phen$Age<=13),"age"]<-"green" # age 8yo to 13yo are green
col[which(phen$Age>13 & phen$Age<=19),"age"]<-"blue"# age 14yo to 19yo are blue. age >19 are black

# color by sex
col[which(phen$sex=="Male"),"sex"]<-"blue"
col[which(phen$sex=="Female"),"sex"]<-"red"

# color by race 
Race<-phen$Race%>%unique
colors<-c("red","blue","green","black","orange")
race.color<-data.frame(Race,colors)
race.color<-full_join(phen,race.color, by="Race")
col$Race<-race.color$colors

# color by latino vs non-latino
col[which(phen$ethnicity=="latino"),"ethnicity"]<-"red" # latino = red
col[which(phen$ethnicity=="not-latino"),"ethnicity"]<-"green" # not latino = green
col[which(phen$ethnicity=="unk"),"ethnicity"] <-"black" # black if unknown ethnicity

# plot the MDS
par(mfrow=c(3,5),xpd=TRUE)

par(mfrow=c(3,5),mar=c(5,5,2,2),mai=c(0,0.5,1,0),xpd=TRUE)
color_factor<- factor(c("25-50%tile","50-75%tile",">75%tile","<25%tile","cell count=0"),
                      levels= c("25-50%tile","50-75%tile",">75%tile","<25%tile","cell count=0"))


library(ggplot2)
library(gridExtra)
library(patchwork)
library(gridGraphics)
library(grid)

mds.ct1<-plotMDS(lcpm.ct,  col= col$BAL_eos_ct_log, labels=phen$BAL_eos_ct, main=source.cell[1])
mds.ct2<-plotMDS(lcpm.ct,  col= col$BAL_eos_p_log, labels=phen$BAL_eos_p, main=source.cell[2])
mds.ct3<-plotMDS(lcpm.ct,  col= col$BAL_neut_ct_log, labels=phen$BAL_neut_ct, main=source.cell[3])
mds.ct4<-plotMDS(lcpm.ct,  col= col$BAL_neut_p_log, labels=phen$BAL_neut_p, main=source.cell[4])
mds.ct5<-plotMDS(lcpm.ct,  col= col$BAL_wbc_log, labels=phen$BAL_wbc, main=source.cell[5])

mds.ct6<-plotMDS(lcpm.ct,  col= col$blood_eos_log, labels=phen$blood_eos, main=source.cell[6])
mds.ct7<-plotMDS(lcpm.ct,  col= col$blood_eos_p_log, labels=phen$blood_eos_p, main=source.cell[7])
mds.ct8<-plotMDS(lcpm.ct,  col= col$blood_neut_log, labels=phen$blood_neut, main=source.cell[8])
mds.ct9<-plotMDS(lcpm.ct,  col= col$blood_neut_p_log, labels=phen$blood_neut_p, main=source.cell[9])
mds.ct10<-plotMDS(lcpm.ct,  col= col$blood_wbc_log, labels=phen$blood_wbc, main=source.cell[10])

mds.batch<-plotMDS(lcpm.ct,  col= col$batch, labels=phen$Batch, main="batch")
mds.nasal_bronch<-plotMDS(lcpm.ct,  col= col$nasal_bronch, labels=phen$Type, main="sample type (nasal/bronch)")
mds.age<-plotMDS(lcpm.ct,  col= col$age, labels=phen$Age, main = "age")
mds.race<-plotMDS(lcpm.ct,  col= col$Race, labels=phen$Race, main = "Race")
mds.eth<-plotMDS(lcpm.ct,  col= col$ethnicity, labels=phen$ethnicity, main = "ethnicity")

# Convert the last plot to a grob
grid.echo()

legend(x=-7, y=5, cex=0.8,col=c("red","green","blue","black","grey"), horiz=TRUE, pch=c(10,10), legend=color_factor, title="Quantiles,ct>0")
legend("bottomleft", cex=0.8,inset=c(-0.2,0),col=c("red","green","blue"), pch=c(10,10), legend=c("<=7yo","8-13yo","14-19yo"), title="Age")


grid.arrange(mds.ct1, mds.ct2,mds.ct3,mds.ct4,mds.ct5,mds.ct6,mds.ct7,mds.ct8,mds.ct9,mds.ct10,
             mds.batch, mds.nasal_bronch,mds.age,mds.race,mds.eth, ncol = 3)
  # significant batch effect



###########################################---------------------###########################################
# sample type: nasal-bronchial
# gene filter: low counts 
###########################################---------------------###########################################
phen<-rbind(nphen,bphen)

##
# Normalising gene expression distributions
##

# select just the nasal RNAseq counts
sid<-phen$SampleID
ct<-counts[,sid] # First column is actually gene name 
genes<-counts$SampleID
rownames(ct)<-genes

## use TMM normalized lcpm as a cutoff point ()
ct_filtered<-filter_low_expressed_genes_method2(ct,9) 
ct<-ct_filtered 

# normalize counts with TMM
norm.factor<-calcNormFactors(ct, method = "TMM")

sample.size<-length(colnames(ct))
for(i in 1:sample.size){
  ct[,i]<-ct[,i]/norm.factor[i]
}

# calculate lcpm based on TMM normalized counts 
lcpm.ct<-cpm(ct,log=TRUE)


##
#MDS by groups
##

#MDS color scheme for nasal 
## color coding based on each factor levels of sex, age, cell counts, and batch sourced from the following code:
## color code covariates 

# make color matrix. all entry is black at baseline 
col<-data.frame(matrix(nrow=sample.size,ncol=16))
colnames(col)<-c(source.cell.log,"nasal_bronch","batch","age","sex","Race","ethnicity")
rownames(col)<-phen$SampleID
for(i in 1:ncol(col)) {col[,i]<-rep("black",sample.size)}

# color for cell counts. 
## grey if cell count = 0
## red, green, blue depending on quantile cutoff among >0 cell counts
### find sample ID of samples with zero cell counts 
source.cell.zero<-sapply(phen[,source.cell],function(p){a<-p; b<-which(a==0); return(b)})
zerosample<-vector("list", length = 10)
names(zerosample)<-source.cell.log

for(i in 1:10){
  zerosample[[i]]<-sid[source.cell.zero[[i]]]
}

### find quantiles for each of the normalized cell counts
source.cell.pos<-sapply(phen[,source.cell],function(p){a<-p; b<-which(a>0); return(b)})

### scaled cell counts list of each cell count > 0
posphen<-vector("list",length=10) 
names(posphen)<-source.cell.log
for(i in 1:10){
  posphen[[i]]<- phen[source.cell.pos[[i]],source.cell.log[i]]
}

### quantile cutoff for each scaled cell count > 0
cell.count.quantiles<-sapply(posphen,calculate_Cell_Count_Quantile) 

## set color 

for(i in 1:10){
  col[cellCountRange1(source.cell.log[i],25,50),source.cell.log[i]]<-"red" # red if pos cell count between 25th - 50th percentile
  col[cellCountRange1(source.cell.log[i],50,75),source.cell.log[i]]<-"green" # green if pos cell count between 50th - 75th percentile
  col[cellCountRange2(source.cell.log[i],75),source.cell.log[i]]<-"blue" # blue if pos cell count above 75th percentile
  col[cellCountRange3(source.cell.log[i],25),source.cell.log[i]]<-"black" # black if pos cell count below 25th percentile or unavailable
  col[source.cell.zero[[i]],source.cell.log[i]]<-"grey" # grey if cell count is 0
}

# colors by batch
# batch1/2/3/4 is red/gree/blue/black, respectively
col[which(phen$Batch=="batch1"),"batch"]<-"red" 
col[which(phen$Batch=="batch2"),"batch"]<-"green"
col[which(phen$Batch=="batch3"),"batch"]<-"blue"

# color by sample type nasal_bronch
col[which(phen$Type=="Nasal"),"nasal_bronch"]<-"red" 
col[which(phen$Type=="Bronchial"),"nasal_bronch"]<-"blue"

# colors by age
col[which(phen$Age<=7),"age"]<-"red" # age less than or equal to 7 are red
col[which(phen$Age>7 & phen$Age<=13),"age"]<-"green" # age 8yo to 13yo are green
col[which(phen$Age>13 & phen$Age<=19),"age"]<-"blue"# age 14yo to 19yo are blue. age >19 are black

# color by sex
col[which(phen$sex=="Male"),"sex"]<-"blue"
col[which(phen$sex=="Female"),"sex"]<-"red"

# color by race 
Race<-phen$Race%>%unique
colors<-c("red","blue","green","black","orange")
race.color<-data.frame(Race,colors)
race.color<-full_join(phen,race.color, by="Race")
col$Race<-race.color$colors

# color by latino vs non-latino
col[which(phen$ethnicity=="latino"),"ethnicity"]<-"red" # latino = red
col[which(phen$ethnicity=="not-latino"),"ethnicity"]<-"green" # not latino = green
col[which(phen$ethnicity=="unk"),"ethnicity"] <-"black" # black if unknown ethnicity

# plot the MDS

par(mfrow=c(2,5),mar=c(5,5,2,2),xpd=TRUE)
par(mfrow=c(2,3),mar=c(5,5,2,2),xpd=TRUE)
color_factor<- factor(c("25-50%tile","50-75%tile",">75%tile","<25%tile","cell count=0"),
                      levels= c("25-50%tile","50-75%tile",">75%tile","<25%tile","cell count=0"))
mds.ct1<-plotMDS(lcpm.ct,  col= col$BAL_eos_ct_log, labels=phen$BAL_eos_ct, main=source.cell[1])
mds.ct2<-plotMDS(lcpm.ct,  col= col$BAL_eos_p_log, labels=phen$BAL_eos_p, main=source.cell[2])
mds.ct3<-plotMDS(lcpm.ct,  col= col$BAL_neut_ct_log, labels=phen$BAL_neut_ct, main=source.cell[3])
mds.ct4<-plotMDS(lcpm.ct,  col= col$BAL_neut_p_log, labels=phen$BAL_neut_p, main=source.cell[4])
mds.ct5<-plotMDS(lcpm.ct,  col= col$BAL_wbc_log, labels=phen$BAL_wbc, main=source.cell[5])
mds.ct6<-plotMDS(lcpm.ct,  col= col$blood_eos_log, labels=phen$blood_eos, main=source.cell[6])
mds.ct7<-plotMDS(lcpm.ct,  col= col$blood_eos_p_log, labels=phen$blood_eos_p, main=source.cell[7])
mds.ct8<-plotMDS(lcpm.ct,  col= col$blood_neut_log, labels=phen$blood_neut, main=source.cell[8])
mds.ct9<-plotMDS(lcpm.ct,  col= col$blood_neut_p_log, labels=phen$blood_neut_p, main=source.cell[9])
mds.ct10<-plotMDS(lcpm.ct,  col= col$blood_wbc_log, labels=phen$blood_wbc, main=source.cell[10])
legend("bottomleft", cex=0.8,inset=c(-0.2,0),col=c("red","green","blue","black","grey"), pch=c(10,10), legend=color_factor, title="Quantiles,ct>0")

mds.batch<-plotMDS(lcpm.ct,  col= col$batch, labels=phen$Batch, main="batch")
mds.nasal_bronch<-plotMDS(lcpm.ct,  col= col$nasal_bronch, labels=phen$Type, main="sample type (nasal/bronch)")
mds.age<-plotMDS(lcpm.ct,  col= col$age, labels=phen$Age, main = "age")
legend("bottomleft", cex=0.8,inset=c(-0.2,0),col=c("red","green","blue"), pch=c(10,10), legend=c("<=7yo","8-13yo","14-19yo"), title="Age")
mds.race<-plotMDS(lcpm.ct,  col= col$Race, labels=phen$Race, main = "Race")
mds.race<-plotMDS(lcpm.ct,  col= col$ethnicity, labels=phen$ethnicity, main = "ethnicity")
mds.sex<-plotMDS(lcpm.ct,  col= col$sex, labels=phen$sex, main = "Sex")
# significant batch effect

library(ggpubr) #load in library for multi-panel figures
#put all three plots together into one multipanel plot
multi_plot<- ggarrange(mds.ct1,mds.ct2,mds.ct3,mds.ct4,mds.ct5,
                       mds.ct6,mds.ct7,mds.ct8,mds.ct9,mds.ct10,#plots that are going to be included in this multipanel figure
                       labels = c(letters[1:10]), #labels given each panel 
                       ncol = 5, nrow = 2, #adjust plot space 
                       common.legend = T) #does the plot have a common legend
#add titles and labels to the multi-panel graph
multi_plot <- annotate_figure(multi_plot,
                              top = text_grob("PCA by cell counts", color = "black", face = "bold", size = 11))
multi_plot 

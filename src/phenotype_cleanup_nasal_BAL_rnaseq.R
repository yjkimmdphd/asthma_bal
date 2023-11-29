#########################
# cleaning up patient phenotype data table for nasal_BAL_rnaseq.R
#########################
library(limma)
library(edgeR)
library(dplyr)
############################################
## set working directory and load count data
############################################
# if working on HPC, set WD as YJ's asthma_transcriptome folder
# if working on local PC, set WD as YJ's asthma-allergy-bunyavanich folder on windows
ifelse(getwd()=="/sc/arion/projects/asthma-allergy/MS_asthma/Young-Jin/asthma_transcriptome",
       setwd("/sc/arion/projects/asthma-allergy/MS_asthma/Young-Jin/asthma_transcriptome"),
       setwd("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich"))

ifelse(getwd()=="/sc/arion/projects/asthma-allergy/MS_asthma/Young-Jin/asthma_transcriptome",WDlocation<-"Minerva",WDlocation<-"localPC")

# setting directory for count data when working on local PC
file.pc<-c("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/input/original_data/MS_asthma/batch1234_readcount_matrix_allsamples.afterQC.txt")
file.b5.pc<-c("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/input/original_data/MS_asthma/AsthmaNasal_RNAseq.batch5.geneID_readcount_rmscaffold_highqualitysamples.txt")

# setting directory for count data when working on HPC minerva
file.minerva<-c("/sc/arion/projects/asthma-allergy/MS_asthma/outputs/batch1234_combined/read_count/batch1234_readcount_matrix_allsamples.afterQC.txt")
file.b5.minerva<-c("/sc/arion/projects/asthma-allergy/MS_asthma/outputs/batch5/read_count/AsthmaNasal_RNAseq.batch5.geneID_readcount_rmscaffold_highqualitysamples.txt")



# reading count data for batch1,2,3,4 and saving as 'counts'
file<-ifelse(WDlocation=="Minerva",file.minerva,file.pc)
file.exists(file)
counts<-read.delim(file[1])
counts.ID<-colnames(counts)

# reading count data for batch 5 as 'counts.b5'
file.b5<-ifelse(WDlocation=="Minerva",file.b5.minerva,file.b5.pc)
file.exists(file.b5)
counts.b5<-read.delim(file.b5[1])
counds.b5ID<-colnames(counts.b5)


######################
## load phenotype data
######################

# asthma biomarker phenotype file, nasal, saved in  'phenotype'
filename2<-file.path(getwd(),"input/asthma-phenotype-filtered-revised-2023-10-04.csv")
file.exists(filename2)
phenotype<-read.csv(filename2)

# get batch information
filename3<-file.path(getwd(),"input/original_data/MS_asthma/MS_asthma_phenotype.batch1234.txt")
file.exists(filename3)
batch.info<-read.delim(filename3)

####################################################################
## transform BAL and CBC information in the phenotype data as needed 
####################################################################

# calculate absolute Eos from total BAL WBC and BAL Eos %
phenotype<-phenotype%>%mutate(
  bal_Eos_ct=BAL_Eos_perc*0.01*BAL_WBC,
  BAL_neut_ct=BAL_neut_perc*0.01*BAL_WBC)

# log transform serum and BAL cell counts
phenotype<-mutate(phenotype, 
       BAL_Eos_ct_log=log(bal_Eos_ct+0.001)%>%round(2),
       BAL_Eos_perc_log = log(BAL_Eos_perc+0.001)%>%round(2),
       BAL_WBC_log=log(BAL_WBC+0.00001)%>%round(2),
       BAL_neut_perc_log=log(BAL_neut_perc+0.001)%>%round(2),
       BAL_neut_ct_log10=log10(BAL_neut_ct+0.00001)%>%round(2),
       serum_Eos_log10=log10(serum_Eos+0.001)%>%round(2),
       serum_Eos_perc_log=log(serum_Eos_perc+0.001)%>%round(2),
       serum_Neut_log10=log10(serum_Neut+0.001)%>%round(2),
       serum_Neut_perc_log=log(serum_neut_perc+0.001)%>%round(2),
       serum_WBC_log10=log10(serum_WBC)%>%round(2)
       )
# find which cell count data is zero and nonzero 
phenotype<-mutate(phenotype,
       zero.BALEos=ifelse(BAL_Eos_perc>0,"nonZero","zero"),
       zero.BALWbc=ifelse(BAL_WBC>0,"nonZero","zero"),
       zero.BALNeut=ifelse(BAL_neut_perc>0,"nonZero","zero"),
       zero.serEos=ifelse(serum_Eos_perc>0,"nonZero","zero"),
       zero.serNeut=ifelse(serum_neut_perc>0,"nonZero","zero"))
z<-lapply(phenotype[,c("zero.BALEos","zero.BALWbc", "zero.BALNeut", "zero.serEos", "zero.serNeut")],
       function(d){factor(d,levels=c("zero","nonZero"))})%>%data.frame()
phenotype[,c("zero.BALEos","zero.BALWbc", "zero.BALNeut", "zero.serEos", "zero.serNeut")]<-z

# subset  'phenotype' with BAL and blood cell count information into 'phenotype.counts'
p.col<-colnames(phenotype)
count.col<-sapply(c("Eos","WBC","serum","BAL","bal","Neut"), function(keyword){grep(keyword,p.col)})%>%unlist%>%unique()
phenotype.counts<-phenotype[,c("subject_assgn",p.col[count.col],"Blood_draw_date")]
phenotype.counts<-phenotype.counts[,c("subject_assgn",colnames(phenotype.counts)[-1]%>%sort)]
colnames(phenotype.counts)[1]<-"SampleID"

#######################################################################
## find subject assignment ID with nasal and bronchial cell RNAseq data  
#######################################################################
nasal.ID<-phenotype[,"subject_assgn"] # nasal study ID that have been collected
bronch.ID<-sub("N","B",nasal.ID) # bronchial subject ID 
nb.ID<-c(nasal.ID,bronch.ID) # nasal and bronchial subject ID
a<-counts.ID%in%nb.ID  
nb.exist.ID<-counts.ID[a] # subject ID nasal and bronchial samples that are in the count matrix

# subset 'phenotype.counts' for which RNA counts data exist to 'p.counts'
sample.ID<-counts.ID[counts.ID%in%nasal.ID]
p.counts<-phenotype.counts[nasal.ID%in%sample.ID,]

# subset batch info filtered based on nb.exist.ID.
batch.info.ID<-which(batch.info$SampleID%in%nb.exist.ID)
batch.info.BAL<-batch.info[batch.info.ID,]
(batch.info.BAL$Type=="Nasal")%>%table

# left join batch info with 'p.counts'
p.counts<-left_join(p.counts,batch.info.BAL, by="SampleID")
p.counts$Batch<-factor(p.counts$Batch, levels=unique(batch.info.BAL$Batch))
p.counts$Race_corrected<-factor(p.counts$Race_corrected, levels=p.counts$Race_corrected%>%unique)

```
################################################
## calculate time difference between BAL and CBC
################################################
``
p.counts<-mutate(p.counts,BAL_Date=BAL_Date%>%as.Date(format="%m/%d/%Y"), 
	Blood_draw_date=Blood_draw_date%>%as.Date(format="%m/%d/%Y"), 
	BAL_CBC_delay=as.Date(BAL_Date,format="%m/%d/%Y")-as.Date(Blood_draw_date,format="%m/%d/%Y"))
```
## some serum or Eos counts are NA of missing values
# make phenotype and count table removing NA
# could conveniently use na.omit(), but that would remove too many data points and reduce power 
``

a<-sapply(p.counts,is.na)
a1<- unique(which(a==TRUE)%%45)
p.counts[a1,]
p.count.BalNeut<-p.counts[-c(41,42),]
p.count.SerCt<-p.counts[-c(13,22,28,37,38),]

sapply(p.counts,is.na)%>%colSums() 
sapply(p.count.BalNeut,is.na)%>%colSums() # check if BAL cols with Neut of the new dataframe has 'NA'
sapply(p.count.SerCt,is.na)%>%colSums() # check if cols with serum of the new dataframe has 'NA'

# should use 'p.count.BalNeut' for DEG using BAL Neut Ct or %
# should use 'p.count.SerCt' for DEG using Ser counts 
# summary statistics of each cell count

pc.cn<-colnames(p.counts)
pc.cn<-pc.cn[-grep("log",pc.cn)]
pc.cn<-pc.cn[c(3:7,9:13)]
pcrange=sapply(p.counts[,pc.cn],function(d){range(d,na.rm=TRUE)})
pc.cn.df<-data.frame(mean=sapply(p.counts[,pc.cn],function(d){mean(d,na.rm=TRUE)%>%round(2)}),
                     sd=sapply(p.counts[,pc.cn],function(d){sd(d,na.rm=TRUE)%>%round(2)}),
                     range=apply(pcrange,2,function(d){paste(d[1],d[2],sep="-")}))

delay.input<-p.counts
delay.input<-delay.input%>%filter(BAL_CBC_delay<3000)
cbc.bal.delay=data.frame(
  mean.delay=delay.input$BAL_CBC_delay%>%mean(na.rm=TRUE)%>%round(2),
  sd.delay=delay.input$BAL_CBC_delay%>%sd(na.rm=TRUE)%>%round(2),
  range=paste(range(delay.input$BAL_CBC_delay,na.rm=TRUE)[1],
              range(delay.input$BAL_CBC_delay,na.rm=TRUE)[2],
              sep="-"))
cbc.bal.delay.abs=data.frame(
  mean.delay=delay.input$BAL_CBC_delay%>%abs%>%mean(na.rm=TRUE)%>%round(2),
  sd.delay=delay.input$BAL_CBC_delay%>%abs%>%sd(na.rm=TRUE)%>%round(2),
  range=paste(range(delay.input$BAL_CBC_delay%>%abs,na.rm=TRUE)[1],
              range(delay.input$BAL_CBC_delay%>%abs,na.rm=TRUE)[2],
              sep="-"))

print(pc.cn.df)    
print(cbc.bal.delay)
print(p.counts%>%select(SampleID,BAL_CBC_delay)%>%arrange(desc(BAL_CBC_delay)))
# bal.cbc.hist<-hist(delay.input$BAL_CBC_delay%>%as.numeric, 
#     breaks=seq(-1000,4000,100), 
#     xlab="BAL date - CBC date (days)", 
#     ylab="number of samples")


###########################################
## filtering counts table to remove low expressed genes
###########################################

# select just the nasal RNAseq counts
b<-counts.ID%in%nasal.ID
x<-counts[,b]
genes<-counts$SampleID
rownames(x)<-genes

# convert raw counts to CPM
# average coverage is ~51 million
cpm0 <-cpm(x)
lcpm0<-cpm(x,log=TRUE)
L<-mean(colSums(x))*1e-6 # mean library size
M<-median(colSums(x))*1e-6 # median library size


#x3 is the TMM normalized count
# normalize counts with TMM
norm.factor<-calcNormFactors(x, method = "TMM")

x3<-x
sample.size<-length(colnames(x3))
for(i in 1:sample.size){
  x3[,i]<-x3[,i]/norm.factor[i]
}

# calculate lcpm based on TMM normalized counts 
lcpm.x3<-cpm(x3,log=TRUE)

#### Removing genes that are lowly expressed

# checking how many genes have 0 count across all samples
table(rowSums(x==0)==45)
# about 4% of the genes have 0 counts across all samples 

# setting a lcpm cutoff for filtering genes with very low counts
lcpm.cutoff <- log2(10/M + 2/L) # M is median. L is mean. library size
dropCutoff<-function(cutoff){
  which(apply(lcpm.x3, 1, max) < cutoff)
}
d1<-dropCutoff(0) 
d2<-dropCutoff(lcpm.cutoff)
dim(x3[-d1,])
dim(x3[-d2,])

################################################################################################
## subsetting counts table based on the BAL and CBC data availability and gene expression filter
################################################################################################
x.BalNeut<-x[,c( p.count.BalNeut$SampleID)] # count table for DEG using BAL Neut information as predictor
x.SerCt<-x[,c(p.count.SerCt$SampleID)]# count table for DEG using serum cell counts information as predictor

x2<-x[-d1,]

x2.BalNeut<-x2[,c( p.count.BalNeut$SampleID)] # count table for DEG using BAL Neut information as predictor
x2.SerCt<-x2[,c(p.count.SerCt$SampleID)]# count table for DEG using serum cell counts information as predictor

#### filtering phenotype table based on cell counts
p.BalEos.pos<-p.counts%>%filter(bal_Eos_ct>0)
p.BalNeut.pos<-p.counts%>%filter(BAL_neut_ct>0)
p.serEos.pos<-p.counts%>%filter(serum_Eos>0)
p.serNeut.pos<-p.counts%>%filter(serum_Neut>0)

#### make new count tables with sampleID filtered for cell counts
x.BalEos.pos<-x[,p.BalEos.pos$SampleID]
x.BalNeut.pos<-x[,p.BalNeut.pos$SampleID]
x.serEos.pos<-x[,p.serEos.pos$SampleID]
x.serNeut.pos<-x[,p.serNeut.pos$SampleID]

x2.BalEos.pos<-x2[,p.BalEos.pos$SampleID]
x2.BalNeut.pos<-x2[,p.BalNeut.pos$SampleID]
x2.serEos.pos<-x2[,p.serEos.pos$SampleID]
x2.serNeut.pos<-x2[,p.serNeut.pos$SampleID]

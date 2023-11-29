##
# RNA seq analysis of nasal data with edgeR and limma
# 
## 
library(dplyr)
library(limma)
library(Glimma)
library(edgeR)
library(gplots)
library(DESeq2)

# load phenotype data by sourcing the following code 
source("./codes/phenotype_cleanup_bronchial_BAL_rnaseq.R")


# print dims of count data
print(list("unfiltered countdata",`dim table`=dim(x)))
print(list("filtered countdata",`dim table`=dim(x2)))
write.csv(
  list("unfiltered countdata",`dim table`=dim(x),"filtered countdata",`dim table`=dim(x2)),file.path(getwd(),"output",paste0("count_data_dims_bronch_bal_rnaseq",Sys.Date(),".csv")),row.names=FALSE)



##########################################
# make input dataframe for DEG with DESeq2
##########################################


#### find which are going to be predictor variables for gene expressions, i.e., names of the log transformed cell counts
# save as 'p.counts.var'
p.counts.var<-
  colnames(p.counts)[grep("log",colnames(p.counts))]

##### make dataframe of input values that will be used for DEG with DESeq2
print("make dataframe of input values that will be used for DEG with DESeq2")
# save as 'df.deseq2input'
# x is original count table
# x2 is x filtered for low gene counts
df.deseq2input<-data.frame(count.data="x",
                           col.data="p.counts", 
                           design = paste0("~ ",p.counts.var," + Batch"), 
                           resoutput = p.counts.var)
df.deseq2input[c(3,4),1:2]<-data.frame(rep("x.BalNeut",2),rep("p.count.BalNeut",2))
df.deseq2input[6:10,1:2]<-data.frame(rep("x.SerCt",5),rep("p.count.SerCt",5))
df.input2<-df.deseq2input%>%mutate(count.data=sub("x","x2",count.data))
df.deseq2input<-rbind(df.deseq2input,df.input2)
print(df.deseq2input)


# update 'df.deseq2input'
input.pos<-data.frame(count.data=rep(c("x.BalEos.pos", "x.BalNeut.pos", "x.serEos.pos", "x.serNeut.pos",
                                       "x2.BalEos.pos","x2.BalNeut.pos","x2.serEos.pos","x2.serNeut.pos"),
                                     each=2))
input.pos$col.data<-rep(rep(c("p.BalEos.pos","p.BalNeut.pos","p.serEos.pos","p.serNeut.pos"),each=2),2)
input.pos$design<-rep((df.deseq2input[c(1:4,6:9),"design"]),2)
input.pos$resoutput<-rep((df.deseq2input[c(1:4,6:9),"resoutput"]),2)
df.deseq2input<-rbind(df.deseq2input,input.pos)
print(df.deseq2input)

input<-df.deseq2input
z.counts<-c("zero.BALEosnonZero","zero.BALNeutnonZero","zero.serEosnonZero")
z.input<-data.frame(count.data=c("x","x.BalNeut","x.SerCt" ,
                                 "x2","x2.BalNeut","x2.SerCt" ),
                    col.data=rep(c("p.counts","p.count.BalNeut", "p.count.SerCt"),2),
                    design = paste0("~ ",z.counts," + Batch")%>%rep(2),
                    resoutput = z.counts%>%rep(2))
df.deseq2input<-rbind(input,z.input)
df.deseq2input$sample<-"bronch"
print(df.deseq2input)


#### summary of analysis is saved in 'res.table'
fluid.cell.sample<-sapply(df.deseq2input$resoutput,function(d){substring(d,1,7)})
perc<-grep("perc",df.deseq2input$resoutput)
fcp<-data.frame(fluid.cell=fluid.cell.sample)
fcp[perc,"unit"]<-"%"
fcp[-perc,"unit"]<-"abs count"
res.table<-data.frame(software="DESEq2",
                      fluid_cell=fcp[,1],
                      units=fcp[,2],
                      model=df.deseq2input$design,
                      sig.genes=0,
                      count_data=df.deseq2input$count.data,
                      sample="bronch")
print(res.table)

#######################################
# run the DEG for continuous predictors
#######################################

#### define function deseq2DEG 
# saves DEG results as a list consisting of the results and significant DEG results
deseq2DEG<-function(countdata,coldata,design,resultname){
  print(resultname)
  dds<-DESeqDataSetFromMatrix(countData = get(countdata),colData=get(coldata), design=as.formula(design))
  dds<-DESeq(dds)
  res<-results(dds, name=resultname)
  res <- res[order(res$padj),]
  # signficant results with padj<0.05
  res.sig<-res[which(res$padj<0.05),]
  return(list(res,res.sig))
}


cont.var<-grep("zero",res.table$fluid_cell, invert=TRUE)

print(cont.var)

for(i in cont.var){
  print(df.deseq2input[i,4])
  assign(paste0("res",i),deseq2DEG(df.deseq2input[i,1],df.deseq2input[i,2],df.deseq2input[i,3],df.deseq2input[i,4]))
}

########################################################
# DEG comparing zero vs nonZero Eos/Neut in BAL or serum
########################################################
#### define function deseq2DE.disc for ananlysis using model matrix containing discrete variable
# saves DEG results as a list consisting of the results and significant DEG results
deseq2DEG.disc<-function(countdata,coldata,des,resultname){
  print(resultname)
  dds<-DESeqDataSetFromMatrix(countData = get(countdata),colData=get(coldata), design=des)
  dds<-DESeq(dds)
  res<-results(dds, name=resultname)
  res <- res[order(res$padj),]
  # signficant results with padj<0.05
  res.sig<-res[which(res$padj<0.05),]
  return(list(res,res.sig))
}

### run the DEG on discrete predictors 

# make model
ml<-list(zeroBalEos=model.matrix(~ zero.BALEos + Batch, p.counts),
         zero.BALNeut=model.matrix(~ zero.BALNeut + Batch, p.count.BalNeut),
         zero.serEos=model.matrix(~ zero.serEos + Batch, p.count.SerCt),
         zeroBalEos=model.matrix(~ zero.BALEos + Batch, p.counts),
         zero.BALNeut=model.matrix(~ zero.BALNeut + Batch, p.count.BalNeut),
         zero.serEos=model.matrix(~ zero.serEos + Batch, p.count.SerCt))

zero.var<-grep("zero",res.table$fluid_cell)
print(zero.var)

for(i in zero.var){
  print(df.deseq2input[i,4])
  assign(paste0("res",i),deseq2DEG.disc(df.deseq2input[i,1],df.deseq2input[i,2],ml[[i-length(cont.var)]],df.deseq2input[i,4]))
}

##########################################################
# summarize number of significant genes in the 'res.table'
##########################################################
for(i in c(cont.var,zero.var)){
  res.table[i,"sig.genes"]<-get(paste0("res",i))[[2]]%>%nrow
  res.table[i,"results"]<-paste0("res",i)
}
print(res.table)


# write results

for(i in c(cont.var,zero.var)){
  a<-get(paste0("res",i))
  write.csv(a[[1]],row.names=TRUE,file.path(getwd(),"output",paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")))
  res.table[i,"output"]<-paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")
}

print(res.table)

write.csv(res.table,file.path(getwd(),"output",paste0("res.table_bronch_rnaseq",Sys.Date(),".csv")))


#Make a basic volcano plot

#with(res4, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
#with(subset(res4, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#with(subset(res4, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


######
## write DEG table
table(deg.tab$genes)%>%as.data.frame()%>%arrange(desc(Freq))

# save results as CSV files
write.csv(deg.tab,file.path(getwd(),"output",paste0("DEG_table_bronch_rnaseq",Sys.Date(),".csv")))
write.csv(res.table,file.path(getwd(),"output",paste0("res.table_bronch_rnaseq",Sys.Date(),".csv")))
for(i in zero.var){
  a<-get(paste0("res",i))
  write.csv(a[[1]],file.path(getwd(),"output",paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")))
  res.table[i,"output"]<-paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")
}

print(sessionInfo())
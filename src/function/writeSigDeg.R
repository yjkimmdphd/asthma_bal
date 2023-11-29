##############################################
## make data frame saving all significant DEG
##############################################

# define function that 
# 1) load all DEG results in a DEG results folder
# 2) extract significant DEG results 
# 3) combine all analysis to one data.frame
# fp = file path for the folder containing the DEG results
# rt = name for the file containing results table 
# example: fp<-file.path(getwd(),"output/DEG_2023-09-30")
# example: rt<-"res.table_2023-09-30.csv"
writeSigDeg<-function(fp,rt){
  library(tidyverse)
  library(data.table)
  if(file.exists(fp)&file.exists(file.path(fp,rt))){
    print(paste(as.character(file.path(fp,rt)),"exists"))
    flist<-list.files(fp)
    # DEG output file name convention is 'DEG_res#_fluid-cell_date.csv'
    fl<-flist[grep("DEG_res",flist)] # filter for DEG files
    flist<-fl%>%unlist%>%as.vector
    print(c("flist:",flist))
    fl<-sapply(fl,function(d)str_split(d,"DEG_res"))%>%unlist%>%as.vector
    fl<-fl[!is.na(fl)&fl!=""]%>%as.data.frame() # select for non-empty strings 
    print(cbind(flist,fl))
    colnames(fl)<-"fl" 
    fl<-separate(fl, col=fl, into=c('a', 'b'), sep='_', extra="merge")
    fl<-fl$a%>%as.numeric
    f.df<-data.frame(fl,flist)%>%arrange(fl)
    print(f.df)
    for(i in 1:nrow(f.df)){
      assign(paste0("res",i),
             read.csv(file.path(fp,f.df[i,"flist"]),row.names = 1)%>%mutate(results=paste0("res",i))%>%rownames_to_column("genes"),
             envir=.GlobalEnv)
    }
    
    # summarize number of significant genes in the 'res.table'
    res.table<-read.csv(file.path(fp,rt))
    cont.var<-grep("zero",res.table$fluid_cell, invert=TRUE)
    zero.var<-grep("zero",res.table$fluid_cell)
    ##############################
    # make list of significant DEG 
    ##############################
    res.list=c(paste0("res",c(cont.var,zero.var)))
    dt<-rbindlist((lapply(res.list,get)))
    sig.deg<-dt%>%filter(padj<0.05)%>%group_by(results)
    write.csv(
      sig.deg,file.path(getwd(),"output",paste0("sig_deg_",Sys.Date(),".csv")),row.names=FALSE
    )
    print(
      paste("list of significant DEG saved in",
                file.path(getwd(),"output",paste0("sig_deg_",Sys.Date(),".csv"))%>%as.character
            )
      )
  } else{
    print(paste("fix input file path"))
  }
  
}

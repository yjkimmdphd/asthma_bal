###################################################################################
# analyzing overlap of DEG in different analysis in nasal_bronchial_BAL_rnaseq_v2.R
###################################################################################
library(dplyr)
library(ggvenn)
library(gplots)
sd.path<-file.path(getwd(),"output/sig_deg_2023-09-30.csv")
rt.path<-file.path(getwd(),"output/DEG_2023-09-30/res.table_2023-09-30.csv")

if(file.exists(sd.path)){ print(paste("sd.path exists?",file.exists(sd.path)))
  sig.deg<-read.csv(sd.path)}
if(file.exists(rt.path)){print(paste("rt.path exists?",file.exists(rt.path)))
  res.table<-read.csv(rt.path)}
sig.deg<-left_join(sig.deg,res.table)%>%group_by(fluid_cell,results,units,count_data)
pred<-sig.deg$fluid_cell%>%unique

# vendiagram of DEG of serum Eos abs count and % as predictors, 

# listing of of fluid, cell count type, results, and count data from which we are analyzing 
summarize(sig.deg,)%>%write.csv(file.path(getwd(),"output/sigdeg_summary_2023_10-02.csv"))
print(list(paste("number of predictors =",length(pred)),
           paste("predictors:", paste(pred,collapse=", ")))
      )

# significant DEG subset by predictor (BAL, serum, Eos, Neut)
deg.bn<-sig.deg%>%filter(fluid_cell==pred[1])
deg.bw<-sig.deg%>%filter(fluid_cell==pred[2])
deg.se<-sig.deg%>%filter(fluid_cell==pred[3])
deg.sn<-sig.deg%>%filter(fluid_cell==pred[4])
deg.be<-sig.deg%>%filter(fluid_cell==pred[5])

# make list of data for venndiagram of DEG 
bn.gl<-ungroup(deg.bn)%>%as.data.frame()%>%split(deg.bn$results)%>%lapply(function(d)select(d,genes)) #for BAL Neut
bw.gl<-ungroup(deg.bw)%>%as.data.frame()%>%split(deg.bw$results)%>%lapply(function(d)select(d,genes)) #for BAL WBC
se.gl<-ungroup(deg.se)%>%as.data.frame()%>%split(deg.se$results)%>%lapply(function(d)select(d,genes)) #for ser Eos
sn.gl<-ungroup(deg.sn)%>%as.data.frame()%>%split(deg.sn$results)%>%lapply(function(d)select(d,genes)) #for ser Neut
be.gl<-ungroup(deg.be)%>%as.data.frame()%>%split(deg.be$results)%>%lapply(function(d)select(d,genes)) #for BAL Eos

# make venndiagram
customVenn<-function(res,genelist){
  #res is a character vector of results that you want to make VD with
  a<-genelist
  r<-res
  res.names<-sapply(r,
                    function(d)filter(res.table,results==d)%>%
                      select(count_data)%>%
                      unlist%>%as.vector%>%unname)
  venn(a[r],names=paste0(res.names,"(count) ",names(res.names)),simplify=TRUE)
}


coda<-sapply(c("res4","res14"),
                   function(d)filter(res.table,results==d)%>%
                     select(count_data)%>%
                     unlist%>%as.character)
coda



bn.v.all.bn<-customVenn(c("res4","res14"),bn.gl)

bn.v.pos.bn<-customVenn(c("res23","res24","res31","res32"),bn.gl)
bw.v<-bw.gl%>%venn
se.v.all.se<-customVenn(c("res6","res7","res16","res17"),se.gl)
se.v.pos.se<-customVenn(c("res25","res26","res33","res34"),se.gl)
se.v.pos.lowR.se<-customVenn(c("res16","res17","res33","res34"),se.gl) # 
sn.

sn.v<-sn.gl[1:4]%>%venn
be.v<-be.gl%>%venn
attr(bn.v,"intersections")
deg.bn%>%summarise()%>%print #> count data used were when BAL Neut>0

# venndiagram of DEG for BAL Neut
bn.gl<-ungroup(deg.bn)%>%as.data.frame()%>%split(deg.bn$results)%>%lapply(function(d)select(d,genes))
bn.v<-bn.gl[2:5]%>%venn
attr(bn.v,"intersections")
attr(bw.v,"intersections")
attr(se.v,"intersections")
deg.bn%>%summarise()%>%print #> count data used were when BAL Neut>0
deg.bw%>%summarise()%>%print #> count data used were when BAL Neut>0
deg.se%>%summarise()%>%print #> count data used were when BAL Neut>0
deg.sn%>%summarise()%>%print #> count data used were when BAL Neut>0

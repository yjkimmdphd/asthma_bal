library(dplyr)
library(pheatmap)
file_path<-"./reports/local_only/nasal_comparison_01-09-24_vs_08-07-24/nasal_comp.txt"
file.exists(file_path)
ncomp<-read.table(file_path,row.names=NULL,header=TRUE)
res<-unique(ncomp$res)
# ncomp%>%filter(res==res[1])%>%mutate_at(vars(log2FoldChange),rank)
# # 
# ncomp<-ncomp%>%group_by(res)%>%mutate_at(vars(log2FoldChange),rank)
ncomp_pivot<-pivot_wider(ncomp%>%select(res,genes,log2FoldChange),values_from = log2FoldChange,names_from = res,values_fill = 0)

ncomp_matrix<-ncomp_pivot[,-1]%>%as.matrix()
rownames(ncomp_matrix)<-ncomp_pivot$genes
plot(x=ncomp_pivot$`n~bld_eos_p_cont_all`,y=ncomp_pivot$`n~bld_aec_cont_mt0`)
model<-lm(data=ncomp_pivot,`n~bld_aec_cont_mt0`~`n~bld_eos_p_cont_all`)
abline(model, col='red',lwd=2)

ncomp_up_overlap<-ncomp_pivot%>%filter(`n~bld_eos_p_cont_all`>0,`n~bld_aec_cont_mt0`>0)
ncomp_down_overlap<-ncomp_pivot%>%filter(`n~bld_eos_p_cont_all`<0,`n~bld_aec_cont_mt0`<0)

write.table(ncomp_up_overlap,"./reports/local_only/nasal_comparison_01-09-24_vs_08-07-24/nasal_comp_overlap_up.txt",row.names = FALSE,col.names = TRUE, sep="\t")
write.table(ncomp_down_overlap,"./reports/local_only/nasal_comparison_01-09-24_vs_08-07-24/nasal_comp_overlap_down.txt",row.names = FALSE,col.names = TRUE, sep="\t")

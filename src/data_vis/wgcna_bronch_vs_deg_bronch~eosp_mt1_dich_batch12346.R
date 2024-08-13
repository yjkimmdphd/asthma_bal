
# DEG analysis was previously performed on bronch ~ BAL Eos% mt1 (dich) + batch12346 to find several unique genes associated with inflammation

library(dplyr)
wgcna_folder<-"./reports/local_only/wgcna/bronch"

deg_folder<-file.path("./resources/processed_data/focused_analysis")
deg_file<-list.files(deg_folder)[grep(".csv",list.files(deg_folder))]

deg<-read.csv(file.path(deg_folder,"deg_bronch_res_sig_5_~ bal_Eos_p_more_1 + Batch_2024-08-10_.csv"), row.names = 1)
eos_down_genes<-deg%>%filter(log2FoldChange<(-0.585))%>%rownames

list.files(wgcna_folder)
module_list<-c("bluebatch12346.txt",
               "brownbatch12346.txt",
               "greenyellowbatch12346.txt",
               "magentabatch12346.txt",
               "pinkbatch12346.txt",
               "purplebatch12346.txt",
               "tanbatch12346.txt",
               "turquoisebatch12346.txt")
modules<-file.path(wgcna_folder,module_list)

m_gene_list<-lapply(modules,read.table)

names(m_gene_list)<-module_list

for(i in seq_along(m_gene_list)){
  print(mean(unlist(m_gene_list[[i]])%in%eos_down_genes))
} # no notable overlap between bronch down gene with WGCNA module genes 



m_gene_list[[4]][unlist(m_gene_list[[4]])%in%eos_down_genes,]
m_gene_list[[7]][unlist(m_gene_list[[7]])%in%eos_down_genes,]

names()

write.table(data.frame(feature="bronch Tan module-n-bronch~bal-ANC DEG",genes=tan_anc_overlap),"./reports/local_only/wgcna/bronch/bronch_tan_anc_overlap.txt",sep="\t",row.names =FALSE)

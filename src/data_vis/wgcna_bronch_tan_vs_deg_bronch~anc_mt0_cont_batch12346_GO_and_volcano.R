# wgcna of Bronchial genes was previously performed to find Tan module 
# DEG analysis was previously performed on bronch ~ BAL ANC mt0 (continuous) + batch12346 to find several unique genes associated with inflammation
# looking at GO terms and volcano plot 

library(dplyr)
wgcna_folder<-"C:/Users/kimyo/Dropbox/Research/asthma_bal/reports/local_only/wgcna/bronch"
deg_list<-"C:/Users/kimyo/Dropbox/Research/asthma_bal/reports/local_only/deg_result_clustering/deg_cluster_utd_08-03-24.csv"
deg<-read.csv(deg_list, row.names = 1)
anc_genes<-deg%>%filter(analysis=="Bronchial . BAL ANC, >0 (continuous) + batch(12346)",direction=="up")%>%pull(gene)

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

for(i in 1:length(m_gene_list)){
  print(mean(unlist(m_gene_list[[i]])%in%anc_genes))
} # 71% of the genes in the Tan module overlap with bronchial DEG found in clustering 



tan_anc_overlap<-m_gene_list[[7]][unlist(m_gene_list[[7]])%in%anc_genes,]

write.table(data.frame(feature="bronch Tan module-n-bronch~bal-ANC DEG",genes=tan_anc_overlap),"./reports/local_only/wgcna/bronch/bronch_tan_anc_overlap.txt",sep="\t",row.names =FALSE)

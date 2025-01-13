# wgcna of Bronchial genes was previously performed to find Tan module 
# DEG analysis was previously performed on bronch ~ BAL ANC mt0 (continuous) + batch12346 to find several unique genes associated with inflammation
# looking at GO terms and volcano plot 

library(dplyr)
wgcna_folder<-"./reports/local_only/wgcna/bronch"
deg_folder<-"./reports/local_only/deg_bronch~eos-neut-mixed-pauci+batch123456"
deg_1<-file.path(deg_folder,"deg_bronch_res_sig_3_2_comp3_neut_vs_mixed_2024-11-05_.csv")
deg_2<-file.path(deg_folder,"deg_bronch_res_sig_4_2_comp4_neut_vs_eos_2024-11-05_.csv")
deg_3<-file.path(deg_folder,"deg_bronch_res_sig_1_1_comp1_mixed_vs_pauci_2024-11-05_.csv")

ifelse(file.exists(deg_1),deg_1_read<-read.csv(deg_1),print("can't"))
ifelse(file.exists(deg_2),deg_2_read<-read.csv(deg_2),print("can't"))
ifelse(file.exists(deg_3),deg_3_read<-read.csv(deg_3),print("can't"))

deg_1_genes<-deg_1_read%>%filter(abs(log2FoldChange)>=1)%>%pull(X)
deg_2_genes<-deg_2_read%>%filter(abs(log2FoldChange)>=1)%>%pull(X)
deg_3_genes<-deg_3_read%>%filter(abs(log2FoldChange)>=1)%>%pull(X)

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
  print(mean(unlist(m_gene_list[[i]])%in%deg_1_genes))
}

for(i in 1:length(m_gene_list)){
  print(mean(unlist(m_gene_list[[i]])%in%deg_2_genes))
}

for(i in 1:length(m_gene_list)){
  print(mean(unlist(m_gene_list[[i]])%in%deg_3_genes))
}


deg_1_overlap<-m_gene_list[[4]][unlist(m_gene_list[[4]])%in%deg_1_genes,]
deg_2_overlap<-m_gene_list[[4]][unlist(m_gene_list[[4]])%in%deg_2_genes,]

write.table(data.frame(feature="bronch Tan module-n-bronch~bal-ANC DEG",genes=tan_anc_overlap),"./reports/local_only/wgcna/bronch/bronch_tan_anc_overlap.txt",sep="\t",row.names =FALSE)


library(tidyverse)
#===============================================================================
# Setup Project Paths
#===============================================================================
project_base <- "."
resources_dir <- file.path(project_base, "resources", "processed_data")
reports_dir <- file.path(project_base, "reports", "local_only")
wgcna_dir <- file.path(reports_dir, "wgcna", "bronch")
output_dir <- file.path(wgcna_dir, "output_2025-05-07") # to analyze wGNCA from january use 'output'

#=============
# load relevant data
#=============

# load wgcna network nodes and edge information  
wgcna_edge<-read.table(file.path(output_dir,"wgcna_bronch_edge_list.txt"),header = TRUE)
wgcna_nodes<-read.table(file.path(output_dir,"wgcna_bronch_node_attributes.txt"),header = TRUE)

head(wgcna_edge)
head(wgcna_nodes)

# load list of modules of interest and hubgenes, and intra/inter modular connectivity information 
top_kwithin<-read.table(file.path(output_dir,"top_kWithin_by_module.txt"),header = TRUE)
top_kwithin <- top_kwithin %>%
  group_by(module) %>%
  mutate(hubgene = kWithin == max(kWithin)) %>%
  ungroup()
top_kwithin <- top_kwithin %>% rename(Gene = gene)

top_kwithin_deg<-read.table(file.path(output_dir,"top_kWithin_by_module_deg_overlap.txt"),header = TRUE)
modules_interest<-unique(top_kwithin_deg$module)

#load table of DEGs
DEG_eos_mt1<-read.csv(file.path(output_dir,"deg_bronch_res_sig_16_~ bal_Eos_p_more_1 + Batch_2025-01-03_.csv"),header = TRUE)
colnames(DEG_eos_mt1)[1]<-"Gene"

#=============
# subset the network attributes
#=============
# subset nodes
wgcna_nodes_OI<-wgcna_nodes[wgcna_nodes$Module%in%modules_interest,]
wgcna_nodes_OI<-left_join(wgcna_nodes_OI,DEG_eos_mt1,by="Gene")
wgcna_nodes_OI<-left_join(wgcna_nodes_OI,top_kwithin,by="Gene")
wgcna_nodes_OI <- wgcna_nodes_OI %>% select(-module)
nodes_oi<-wgcna_nodes_OI$Gene
# subset edges
wgcna_edge<-wgcna_edge%>%mutate(source_oi=Source%in%nodes_oi,target_oi=Target%in%nodes_oi)
edges_oi<-which(rowSums(data.frame(wgcna_edge$source_oi,wgcna_edge$target_oi))/2>=1)
wgcna_edge_oi<-wgcna_edge[edges_oi,][,1:3]

#================
# write the files
#================
write.table(wgcna_nodes_OI,file.path(output_dir,"wgcna_nodes_OI.txt"),quote=FALSE,sep="\t", col.names = TRUE,row.names = FALSE)
write.table(wgcna_edge_oi,file.path(output_dir,"wgcna_edge_oi.txt"),quote=FALSE,sep="\t", col.names = TRUE,row.names = FALSE)

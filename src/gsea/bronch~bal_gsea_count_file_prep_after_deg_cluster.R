library(dplyr)

deg_list<-read.csv("./reports/local_only/deg_result_clustering/deg_cluster_utd_08-03-24.csv", row.names = 1)

analysis_list<-deg_list%>%group_by(analysis)%>%select(analysis,direction)%>%summarize(direction=unique(direction))%>%as.data.frame()

bphen_path<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")
bphen<-read.csv(bphen_path)

pos_sampleID<-bphen%>%filter(BAL_neut_ct>0)%>%pull(SampleID)

normalized_count_table_path<-"./resources/processed_data/bronch_batch12346_normalized_ct.txt"

counts<-if(file.exists(normalized_count_table_path)){read.delim(normalized_count_table_path, check.names = FALSE)}
bronch.counts<-counts
rownames(bronch.counts)<-bronch.counts$name

anc_up<-deg_list%>%filter(analysis=="Bronchial . BAL ANC, >0 (continuous) + batch(12346)",direction=="up")

counts_anc_up<-bronch.counts[rownames(bronch.counts)%in%anc_up$gene,pos_sampleID]

dim(counts_anc_up)

write.csv(counts_anc_up,"./resources/processed_data/gsea/count_deg_bphen_bal_anc_m_0_clustered.txt")

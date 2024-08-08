library(tidyverse)
library(pheatmap)
library(RColorBrewer)

d1<-"./reports/local_only/deg_bronch~bal_cell(continuous)+batch12346/gene_list_sig_results/deg_bronch_continuous_allcell+batch12346_res_sig_2_~ BAL_eos_p_log + Batch_2024-08-02_.csv"

file.exists(d1)
d2<-"./reports/local_only/deg_bronch~bal_cell(continuous)+batch12346/gene_list_sig_results/deg_bronch_continuous_allcell+batch12346_res_sig_5_~ BAL_wbc_log + Batch_2024-08-02_.csv"

d3<-"./reports/local_only/deg_bronch~bal_cell(continuous)+batch12346/gene_list_sig_results/deg_bronch_continuous_poscell+batch12346_res_sig_3_~ BAL_neut_ct_log + Batch_2024-05-29_.csv"

d4<-"./reports/local_only/deg_bronch~bal_cell(dichot)+batch12346/deg_bronch_dichot_res_sig_5_~ bal_Eos_p_more_3 + Batch_2024-05-29_.csv"

d5<-"./reports/local_only/deg_bronch~bal_cell(dichot)+batch12346/deg_bronch_dichot_res_sig_14_~ bld_AEC_more_300 + Batch_2024-05-29_.csv"

d6<-"./reports/local_only/deg_nasal-bronch~bal-blood_cell(continuous)+batch/deg_gene_list/deg_Nasal_poscells_~ blood_eos_log + Batch_res_6_2024-01-09_.csv"

d7<-"./reports/local_only/deg_bronch~bal_cell(dichot)+batch12346/deg_bronch_dichot_res_sig_4_~ bal_Eos_p_more_1 + Batch_2024-05-29_.csv"

d8<-"./reports/local_only/deg_bronch~bal_cell(dichot)+batch12346/deg_bronch_dichot_res_sig_2_~ bal_AEC_more_1 + Batch_2024-08-05_.csv"

sapply(c(d1,d2,d3,d4,d5,d6,d7,d8),file.exists)
d1<-read.csv(d1)
d2<-read.csv(d2)
d3<-read.csv(d3)
d4<-read.csv(d4)
d5<-read.csv(d5)
d6<-read.csv(d6)
d7<-read.csv(d7)
d8<-read.csv(d8)

d1<-d1%>%mutate(
  analysis = rep("bronch. BAL_eos_p_all_cont + Batch12346", nrow(d1)),
  anatomy = rep("bronch", nrow(d1)),
  `cont-dich` = rep("cont", nrow(d1))
)

d2<-d2%>%mutate(
  analysis = rep("bronch. BAL_wbc_all_cont + Batch12346", nrow(d2)),
  anatomy = rep("bronch", nrow(d2)),
  `cont-dich` = rep("cont", nrow(d2))
)

d3<-d3%>%mutate(
  analysis = rep("bronch.bal_anc_m_0_cont+batch12346", nrow(d3)),
  anatomy = rep("bronch", nrow(d3)),
  `cont-dich` = rep("cont", nrow(d3))
)

d4<-d4%>%mutate(
  analysis = rep("bronch.bal_Eos_p_m_3_dich + Batch12346", nrow(d4)),
  anatomy = rep("bronch", nrow(d4)),
  `cont-dich` = rep("dich", nrow(d4))
)

d5<-d5%>%mutate(
  analysis = rep("bronch.bld_AEC_more_300_dich + Batch12346", nrow(d5)),
  anatomy = rep("bronch", nrow(d5)),
  `cont-dich` = rep("dich", nrow(d5))
)

d6<-d6%>%mutate(
  analysis = rep("nasal.blood_aec_m_0_cont+batch1234", nrow(d6)),
  anatomy = rep("nasal", nrow(d6)),
  `cont-dich` = rep("cont", nrow(d6))
)

d7<-d7%>%mutate(
  analysis = rep("bronch.bal_Eos_p_m_1_dich + Batch12346", nrow(d7)),
  anatomy = rep("bronch", nrow(d7)),
  `cont-dich` = rep("dich", nrow(d7))
)
d8<-d8%>%mutate(
  analysis = rep("bronch.bal_aec_1.15_dich + Batch12346", nrow(d8)),
  anatomy = rep("bronch", nrow(d8)),
  `cont-dich` = rep("dich", nrow(d8))
)

deg_data<-rbind(d1,d2,d3,d4,d5,d6,d7,d8)
deg_data_noWBC<-rbind(d1,d3,d4,d5,d6,d7,d8)

deg_file_path<-"./resources/processed_data/deg_list_as_of_08-03-24.csv"
deg_file_path_noWBC<-"./resources/processed_data/deg_list_as_of_08-03-24_noWBC.csv"

# write.csv(deg_data_noWBC,deg_file_path_noWBC)

deg_data<-read.csv(deg_file_path_noWBC,row.names = NULL)

# Select relevant columns
clustering_data <- deg_data %>%
  select(analysis, X,log2FoldChange) %>%
  na.omit()  # Remove any rows with missing values


# Spread the data to transform it into the desired format
spread_data <- clustering_data %>%
  select(analysis, X,log2FoldChange)%>%
  pivot_wider(values_from =log2FoldChange, names_from = analysis,values_fill = 0 )


# design_list<-c("Bronchial . BAL Eos%, all (continuous) + batch(12346)",
#                "Bronchial . BAL WBC, all (continuous) + batch(12346)",
#                "Bronchial . BAL ANC, >0 (continuous) + batch(12346)",
#                "Bronchial . BAL Eos% > 3% vs <3% + batch(12346)",
#                "Bronchial . Blood AEC > 300 vs <300 + batch(12346)",
#                "Nasal . Blood AEC, >0 (continuous) + batch(1234)",
#                "Bronchial . BAL Eos% > 1% vs <1% + batch(12346)",
#                "Bronchial . BAL AEC >1.15 vs <1.15 + batch(12346)")

design_list<-c("Bronchial . BAL Eos%, all (continuous) + batch(12346)",
               "Bronchial . BAL ANC, >0 (continuous) + batch(12346)",
               "Bronchial . BAL Eos% > 3% vs <3% + batch(12346)",
               "Bronchial . Blood AEC > 300 vs <300 + batch(12346)",
               "Nasal . Blood AEC, >0 (continuous) + batch(1234)",
               "Bronchial . BAL Eos% > 1% vs <1% + batch(12346)",
               "Bronchial . BAL AEC >1.15 vs <1.15 + batch(12346)")

# Prepare the matrix for pheatmap
log2fc_matrix <- as.matrix(spread_data %>% select(-X))
rownames(log2fc_matrix) <- spread_data$X
colnames(log2fc_matrix)
colnames(log2fc_matrix)<-design_list

# export the matrix
deg_matrix_asof080324<-write.table(log2fc_matrix,"./reports/local_only/deg_result_clustering/deg_matrix_asof080324.txt",sep="\t",row.names = TRUE,col.names = NA)

# Identify columns where the maximum absolute value is greater than 2
max_abs_values <- apply(log2fc_matrix, 1, function(x) max(abs(x)))

# Set the maximum and minimum values for the color scale
max_value <- 4
min_value <- -4

# Define a custom color palette that includes grey at the midpoint
custom_colors <- colorRampPalette(c("blue", "dodgerblue2","white" ,"pink", "red"))(100)

# Create the heatmap with custom color scale limits, colors, and column label angle

# pheatmap(log2fc_matrix, 
#          scale = "row",
#          clustering_distance_rows = "correlation", 
#          clustering_method = "complete", 
#          main = "Heatmap of log2FoldChange by Sample",
#          breaks = seq(min_value, max_value, length.out = 100),
#          color = custom_colors,
#          angle_col = 315)

set.seed(100)

out <- pheatmap(log2fc_matrix, 
                            scale = "row",
                            kmeans_k = 10,
                            clustering_distance_rows = "correlation", 
                            clustering_method = "ward.D2", 
                            main = "Heatmap of log2FoldChange by Sample",
                            breaks = seq(min_value, max_value, length.out = 100),
                            color = custom_colors,
                            angle_col = 315)
geneclusters <- out[["kmeans"]][["cluster"]]

Nasal_up_cluster<-names(geneclusters[geneclusters==6])
Bronch_._blood_AEC_down_cluster<-names(geneclusters[geneclusters==2])
Bronch_._blood_AEC_up_cluster<-names(geneclusters[geneclusters==5])
Bronch_._ANC_up_cluster<-names(geneclusters[geneclusters==1])
Bronch_._ANC_down_cluster<-names(geneclusters[geneclusters==4])
Bronch_._High.mod_Eos_cluster<-names(geneclusters[geneclusters==10])
Bronch_._High_Eos_p_up_cluster<-names(geneclusters[geneclusters==9])
Bronch_._High_Eos_p_down_cluster<-names(geneclusters[geneclusters==7])
Nasal_down_cluster<-names(geneclusters[geneclusters==3|geneclusters==8])

nasal_up_genes<-data.frame(analysis=design_list[5],direction="up",gene=Nasal_up_cluster)
nasal_down_genes<-data.frame(analysis=design_list[5],direction="down",gene=Nasal_down_cluster)
br_blood_aec_up<-data.frame(analysis=design_list[4],direction="up",gene=Bronch_._blood_AEC_up_cluster)
br_blood_aec_down<-data.frame(analysis=design_list[4],direction="down",gene=Bronch_._blood_AEC_down_cluster)
br_bal_ANC_up<-data.frame(analysis=design_list[2],direction="up",gene=Bronch_._ANC_up_cluster)
br_bal_ANC_down<-data.frame(analysis=design_list[2],direction="down",gene=Bronch_._ANC_down_cluster)
br_bal_eos_hm_up<-data.frame(analysis=paste(design_list[3],design_list[6],design_list[7],sep="_"),
                          direction="up",gene=Bronch_._High.mod_Eos_cluster)
br_bal_eos_p_up<-data.frame(analysis=design_list[3],direction="up",gene=Bronch_._High_Eos_p_up_cluster)
br_bal_eos_p_down<-data.frame(analysis=design_list[3],direction="down",gene=Bronch_._High_Eos_p_down_cluster)




cluster_list<-rbind(nasal_up_genes,
                    nasal_down_genes,
                    br_blood_aec_up,
                    br_blood_aec_down,
                    br_bal_ANC_up,
                    br_bal_ANC_down,
                    br_bal_eos_hm_up,
                    br_bal_eos_p_up,
                    br_bal_eos_p_down)

write.csv(cluster_list,"./reports/local_only/deg_result_clustering/deg_cluster_utd_08-03-24.csv")


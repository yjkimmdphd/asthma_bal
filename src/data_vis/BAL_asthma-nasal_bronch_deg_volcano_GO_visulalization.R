library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(gridExtra)
library(grid)
data.folder<-file.path("./reports/local_only/deg~bal-blood_cell(continuous)+batch/deg_gene_list")
####################################################################
# data exploration of DEG analysis using nasal/bronchial rnaseq data
# model: "Nasal DEG ~ blood AEC (>0) + Batch
# model: "Bronch DEG ~ BAL ANC(>0) + Batch
####################################################################

res_nasal<-read.csv(file.path(data.folder,"nasal_res6.csv"),row.names = 1)
res_bronch<-read.csv(file.path(data.folder,"bronch_res3.csv"),row.names = 1)
par(mfrow=c(1,2))

nasal_kda_up_genes<-c("FAM92B", "AK7", "DNAH9", "TSNAXIP1", "LRRC23", "RSPH14","RSPH4A","SPA17", "SAXO2","NME9","CCDC108","DZIP1L")
nasal_kda_down_genes<-c("MIER3","FGL2","VSIG10L","FAM49B","HOPX","NAMPT","IL1RAP","TGFBI","PHF20L1","PLXNC1")
anh_do_nasal_kda_cilia_module<-c("FOXJ1","CAPS","TEKT2","TPPP3","C19orf51","C2orf70","LRRC10B","C2orf81","LRRC23","TMEM231")
anh_do_nasal_kda_inf_module<-c("FGL2","PTPRC","FYB","TLR8","SNX10")
nasal_genes<-rownames(res_nasal)

bronch_genes<-rownames(res_bronch)

# create custom 'nasal' gene key-value pairs for color based on whethe rsomething is KDA DEG or not
keyvals_n <- ifelse(
  nasal_genes%in%nasal_kda_down_genes, 'royalblue',
  ifelse(nasal_genes%in%nasal_kda_up_genes, 'gold',
         'grey'))
keyvals_n[which(nasal_genes%in%anh_do_nasal_kda_cilia_module)]<-"cyan"
keyvals_n[which(nasal_genes%in%anh_do_nasal_kda_inf_module)]<-"magenta"


keyvals_n[is.na(keyvals)] <- 'grey'
names(keyvals_n)[keyvals_n == 'gold'] <- 'nasal_KDA_up'
names(keyvals_n)[keyvals_n == 'grey'] <- 'nasal_deg'
names(keyvals_n)[keyvals_n == 'royalblue'] <- 'nasal_KDA_down'
names(keyvals_n)[keyvals_n == 'cyan']<-'nasal_anhdo_kda_cilia'
names(keyvals_n)[keyvals_n == 'magenta']<-'nasal_anhdo_kda_inf'

# create custom 'bronch' gene key-value pairs for color based on whethe rsomething is KDA DEG or not
keyvals_b <- ifelse(
  bronch_genes%in%nasal_kda_down_genes, '#648FFF',
  ifelse(bronch_genes%in%nasal_kda_up_genes, '#FFB000',
         '#DC267F'))

keyvals_b[is.na(keyvals_b)] <- 'grey'
names(keyvals_b)[keyvals_b == 'gold'] <- 'nasal_KDA_up'
names(keyvals_b)[keyvals_b == '#DC267F'] <- 'bronch_deg'
names(keyvals_b)[keyvals_b == 'royalblue'] <- 'nasal_KDA_down'

p1<-EnhancedVolcano(res_nasal,
                lab = rownames(res_nasal),
                title='Nasal DEG',
                subtitle = '~ blood AEC + batch',
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                xlim=c(-10,8),
                ylim=c(0,8),
                pCutoff = 5e-3,
                FCcutoff = 1,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 4.0,
                labSize = 3,
                selectLab = rownames(res_nasal)[which(names(keyvals_n)%in%c('nasal_KDA_up','nasal_KDA_down','nasal_anhdo_kda_cilia','nasal_anhdo_kda_inf'))],
                colAlpha = 0.4,
                colCustom = keyvals_n,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-adj (<0.005) & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 5.0,    
                drawConnectors = TRUE,
                widthConnectors = 0.75)
p1.2<-EnhancedVolcano(res_nasal,
                    lab = rownames(res_nasal),
                    title='Nasal DEG',
                    subtitle = '~ blood AEC + batch',
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    xlab = bquote(~Log[2]~ 'fold change'),
                    xlim=c(-10,8),
                    ylim=c(0,8),
                    pCutoff = 5e-3,
                    FCcutoff = 1,
                    cutoffLineType = 'twodash',
                    cutoffLineWidth = 0.8,
                    pointSize = 4.0,
                    labSize = 3,
                    colAlpha = 0.4,
                    legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                   'p-adj (<0.005) & Log (base 2) FC'),
                    legendPosition = 'right',
                    legendLabSize = 10,
                    legendIconSize = 5.0,    
                    drawConnectors = TRUE,
                    widthConnectors = 0.75)

p1
p1.2

p2<-EnhancedVolcano(res_bronch,
                lab = rownames(res_bronch),
                title = "Bronch DEG",
                subtitle = "~ BAL ANC + batch",
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                xlim=c(-2.5,2.5),
                ylim=c(0,8),
                pCutoff = 5e-3,
                FCcutoff = 0.58,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 4.0,
                labSize = 3,
                selectLab = rownames(res_bronch)[which(names(keyvals_b)%in%c('nasal_KDA_up','nasal_KDA_down'))],
                colAlpha = 0.4,
                colCustom = keyvals_b,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-adj (<0.005) & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 5.0,    
                drawConnectors = TRUE,
                widthConnectors = 0.75)
p2.2<-EnhancedVolcano(res_bronch,
                    lab = rownames(res_bronch),
                    title = "Bronch DEG",
                    subtitle = "~ BAL ANC + batch",
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    xlab = bquote(~Log[2]~ 'fold change'),
                    xlim=c(-2.5,2.5),
                    ylim=c(0,8),
                    pCutoff = 5e-3,
                    FCcutoff = 0.58,
                    cutoffLineType = 'twodash',
                    cutoffLineWidth = 0.8,
                    pointSize = 4.0,
                    labSize = 3,
                    colAlpha = 0.4,
                    legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                   'p-value & Log (base 2) FC (>0.58)'),
                    legendPosition = 'right',
                    legendLabSize = 10,
                    legendIconSize = 5.0,    
                    drawConnectors = TRUE,
                    widthConnectors = 0.75)
p2
p2.2

###########################
## GO analysis  
###########################
go_folder<-file.path("./reports/local_only/deg~bal-blood_cell(continuous)+batch/GO")

# 1. GO term for Nasal DEG ~ blood AEC_+batch 
nasal_go<-read.csv(file.path(go_folder,"GO_david_Nasal_poscells_~ blood_eos_log + Batch_pos_cells_res_6_up_down.csv"))
nasal_go<-nasal_go%>%
  filter(FDR<0.05)%>%
  mutate(Fold.Enrichment=ifelse(deg_lfc2=="down",-Fold.Enrichment,Fold.Enrichment))%>%
  arrange(desc(deg_lfc2),desc(Fold.Enrichment))
nasal_go<-nasal_go[c(1:10,129:138),]
nasal_go$term<-factor(nasal_go$term, levels=unique(nasal_go$term))
library(ggplot2)

# Plot
ggplot(nasal_go, aes(x = Fold.Enrichment, y = term, fill = -log10(FDR))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(-log10(FDR), 1)), hjust = -0.1, size = 3, color = "black") +  # Add text labels for log10(FDR)
  scale_fill_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
  labs(x = "Fold Enrichment", y = "Gene Ontology Term", fill = "-log10(FDR)") +
  theme(axis.text.y = element_text(size = 10))  # Adjust y-axis label size for better readability

# 2. GO term for Bronch DEG ~ BAL ANC +batch
bronch_go<-read.csv(file.path(go_folder,"GO_david_Bronchial_poscells_~ BAL_neut_ct_log + Batch_pos_cells_res_3_min_abs_lfc_0.58_top100.csv"))
bronch_go<-bronch_go%>%
  filter(FDR<0.05)%>%
  mutate(Fold.Enrichment=ifelse(deg_lfc2=="down",-Fold.Enrichment,Fold.Enrichment))%>%
  arrange(desc(deg_lfc2),desc(Fold.Enrichment))
bronch_go<-bronch_go[c(1:10,41:50),]
bronch_go$term<-factor(bronch_go$term, levels=unique(bronch_go$term))
library(ggplot2)
ggplot(bronch_go, aes(x = Fold.Enrichment, y = term, fill = -log10(FDR))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(-log10(FDR), 1)), hjust = -0.1, size = 3, color = "black") +  # Add text labels for log10(FDR)
  scale_fill_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
  labs(x = "Fold Enrichment", y = "Gene Ontology Term", fill = "-log10(FDR)") +
  theme(axis.text.y = element_text(size = 10))  # Adjust y-axis label size for better readability

bronch_go%>%filter(term=="serine protease inhibitor complex")

library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(gridExtra)
library(grid)
data.folder<-file.path("./reports/local_only/deg~bal-blood_cell(continuous)+batch/deg_gene_list")
####################################################################
# data exploration of DEG analysis using bronchial rnaseq data 
# new reference genome and includes batch 12346
# model: Bronch DEG ~ BAL cells + Batch
####################################################################
data.folder<-file.path("./reports/local_only/deg_bronch~bal_cell(continuous)+batch12346/all_results")
filelist<-list.files(data.folder)
res_bronch<-vector("list",length=20)
for(i in 1:20){
  res_bronch[[i]]<-read.csv(file.path(data.folder,filelist[[i]]),row.names = 1)
  }

# padj are BH adjusted pvalues 

par(mfrow=c(2,3))
for(i in 1:6){
  EnhancedVolcano(res_bronch[[i]],
                  lab = rownames(res_bronch[[i]]),
                  title='Bronch DEG',
                  subtitle = '~ bronch cell + batch12346',
                  x = 'log2FoldChange',
                  y = 'padj',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  xlim=c(-2,2),
                  ylim=c(0,6),
                  pCutoff = 5e-2,
                  FCcutoff = 1,
                  cutoffLineType = 'twodash',
                  cutoffLineWidth = 0.8,
                  pointSize = 4.0,
                  labSize = 3,
                  colAlpha = 0.4,
                  legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                 'p-adj (<0.05) & Log (base 2) FC'),
                  legendPosition = 'right',
                  legendLabSize = 10,
                  legendIconSize = 5.0,    
                  drawConnectors = TRUE,
                  widthConnectors = 0.75)
}

EnhancedVolcano(res_bronch[[14]],
                lab = rownames(res_bronch[[2]]),
                title='Bronch DEG',
                subtitle = '~ bronch cell + batch12346',
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                xlim=c(-2,2),
                ylim=c(0,6),
                pCutoff = 5e-2,
                FCcutoff = 1,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 4.0,
                labSize = 3,
                colAlpha = 0.4,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-adj (<0.05) & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 5.0,    
                drawConnectors = TRUE,
                widthConnectors = 0.75)
par(mfrow = c(2, 3))


p1.2<-EnhancedVolcano(res_nasal,
                    lab = rownames(res_nasal),
                    title='Nasal DEG',
                    subtitle = '~ blood AEC + batch',
                    x = 'log2FoldChange',
                    y = 'padj',
                    xlab = bquote(~Log[2]~ 'fold change'),
                    xlim=c(-10,8),
                    ylim=c(0,8),
                    pCutoff = 5e-2,
                    FCcutoff = 1,
                    cutoffLineType = 'twodash',
                    cutoffLineWidth = 0.8,
                    pointSize = 4.0,
                    labSize = 3,
                    colAlpha = 0.4,
                    legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                   'p-adj (<0.05) & Log (base 2) FC'),
                    legendPosition = 'right',
                    legendLabSize = 10,
                    legendIconSize = 5.0,    
                    drawConnectors = TRUE,
                    widthConnectors = 0.75)
res_nasal<-res_nasal%>%mutate(deg_sig=ifelse(padj<0.05&log2FoldChange>1,'cyan',
                                  ifelse(padj<0.05&log2FoldChange< -1,'magenta',
                                         'grey')))


EnhancedVolcano(res_nasal,
                      lab = rownames(res_nasal),
                      title='Nasal DEG',
                      subtitle = '~ blood AEC + batch, p-adj (BH)',
                      x = 'log2FoldChange',
                      y = 'padj',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      xlim=c(-10,8),
                      ylim=c(0,3.5),
                      pCutoff = 5e-2,
                      FCcutoff = 1,
                      cutoffLineType = 'twodash',
                      cutoffLineWidth = 0.8,
                      pointSize = 4.0,
                      labSize = 3,
                      colAlpha = 0.2,
                      colCustom = keyvals_n,
                      selectLab = NA,
                      legendLabSize = 10,
                      legendIconSize = 5.0,    
                      drawConnectors = TRUE,
                      widthConnectors = 0.75)

p1
p1.3

p2<-EnhancedVolcano(res_bronch,
                lab = rownames(res_bronch),
                title = "Bronch DEG",
                subtitle = "~ BAL ANC + batch",
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                xlim=c(-2.5,2.5),
                ylim=c(0,8),
                pCutoff = 5e-2,
                FCcutoff = 0.58,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 4.0,
                labSize = 3,
                selectLab = rownames(res_bronch)[which(names(keyvals_b)%in%c('nasal_KDA_up','nasal_KDA_down'))],
                colAlpha = 0.4,
                colCustom = keyvals_b,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-adj (<0.05) & Log (base 2) FC'),
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
                    y = 'padj',
                    xlab = bquote(~Log[2]~ 'fold change'),
                    xlim=c(-2.5,2.5),
                    ylim=c(0,8),
                    pCutoff = 5e-2,
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
### 
# needs revision below here 
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
# focused GO term plots
focused_nasal_GO_blood_AEC<-c("epithelial cilium movement","cilium movement involved in cell motility","cilium-dependent cell motility","inflammasome complex","negative regulation of T cell mediated cytotoxicity","positive regulation of regulatory T cell differentiation","B cell homeostasis","inflammatory response")
focused_nasal_GO_blood_AEC<-nasal_go[nasal_go$term%in%focused_nasal_GO_blood_AEC,]
focused_nasal_GO_blood_AEC$term<-factor(focused_nasal_GO_blood_AEC$term, levels=unique(focused_nasal_GO_blood_AEC$term))

nasal_go_top10<-nasal_go[c(1:10,129:138),]
nasal_go_top10$term<-factor(nasal_go_top10$term, levels=unique(nasal_go_top10$term))
library(ggplot2)
library(stringr)
wrapped_label<-str_wrap(focused_nasal_GO_blood_AEC$term,width = 20)
# Plot
ggplot(focused_nasal_GO_blood_AEC, aes(x = Fold.Enrichment, y = term, fill = -log10(FDR))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = round(-log10(FDR), 1)), fill="white",nudge_y=0.3, hjust = -0.1, size = 3, color = "black") +  # Add text labels for log10(FDR)
  scale_fill_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
  labs(x = "Fold Enrichment", y = "Gene Ontology Term", fill = "-log10(FDR)") +
  theme(axis.text.x = element_text(size = 15),  # Change size of x-axis labels
        axis.text.y = element_text(size = 17),
        axis.title= element_text(size=20))+  # Adjust y-axis label size for better readability
  scale_y_discrete(labels=wrapped_label)

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
  geom_label(aes(label = round(-log10(FDR), 1)), fill="white",nudge_y=0.3, hjust = -0.1, size = 3, color = "black") +  # Add text labels for log10(FDR)
  scale_fill_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
  labs(x = "Fold Enrichment", y = "Gene Ontology Term", fill = "-log10(FDR)") +
  theme(axis.text.x = element_text(size = 15),  # Change size of x-axis labels
        axis.text.y = element_text(size = 15),
        axis.title= element_text(size=20))  # Change size of y-axis labels  # Adjust y-axis label size for better readability

bronch_go%>%filter(term=="serine protease inhibitor complex")


######
# volcano plot of bronch DE genes assoc with BAL neutrophilia 
# plot the genes assoc with the GO terms
######

rownames(bronch_go)<-c(1:nrow(bronch_go))
sapply(bronch_go$Genes[1:16], function(genes){strsplit(genes, ", ")})
  
    unlist()
# create custom 'bronch' gene key-value pairs for color based on whethe rsomething is KDA DEG or not
keyvals_b <- ifelse(
  bronch_genes%in%nasal_kda_down_genes, '#648FFF',
  ifelse(bronch_genes%in%nasal_kda_up_genes, '#FFB000',
         '#DC267F'))

keyvals_b[is.na(keyvals_b)] <- 'grey'
names(keyvals_b)[keyvals_b == 'gold'] <- 'nasal_KDA_up'
names(keyvals_b)[keyvals_b == '#DC267F'] <- 'bronch_deg'
names(keyvals_b)[keyvals_b == 'royalblue'] <- 'nasal_KDA_down'

p3<-EnhancedVolcano(res_bronch,
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
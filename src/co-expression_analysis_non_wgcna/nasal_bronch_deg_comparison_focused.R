# DEG and GO term visualization

focused_list<-c("deg_bronch_res_sig_3_",
  "deg_bronch_res_sig_5_",
  "deg_nasal_res_sig_13_",
  "deg_nasal_res_sig_7_")

focused_list_all<-c("deg_bronch_res_all_3_",
                "deg_bronch_res_all_5_",
                "deg_nasal_res_all_13_",
                "deg_nasal_res_all_7_")

library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(gridExtra)
library(grid)

deg_folder<-file.path("./resources/processed_data/focused_analysis")
deg_file<-list.files(deg_folder)
matching_files <- deg_file[sapply(deg_file, function(file) any(sapply(focused_list, grepl, file)))]
deg_file<-matching_files

all_deg_folder<-file.path(deg_folder,"deg_all")
all_deg_file<-list.files(all_deg_folder)
matching_files <- all_deg_file[sapply(all_deg_file, function(file) any(sapply(focused_list_all, grepl, file)))]
all_deg_file<-matching_files


analysis_list<-c("br.bal.aec.mt1.2",
                 "br.bal.eosp.mt1",
                 "n.bld.aec.mt100",
                 "n.bld.eosp.cont")

# Tables containing the gene names, l2fc, padj of only DEGs in each analysis
deg_file_call_df<-data.frame(ananlysis=analysis_list,file=deg_file)
deg_list<-lapply(file.path(deg_folder,deg_file),function(d)read.csv(d,header=TRUE,row.names = 1))
names(deg_list)<-deg_file_call_df$ananlysis

# Find # of DEGs abs(log2FC) > 0.585 (1.5fold increase) and Padj < 0.05
print(lapply(deg_list,function(d)filter(d,abs(log2FoldChange)>0.585)%>%dim))

# Tables containing the gene names, l2fc, padj of all genes analyzed in each analysis
all_deg_file_call_df<-data.frame(ananlysis=analysis_list,file=all_deg_file)
all_deg_list<-lapply(file.path(all_deg_folder,all_deg_file),function(d)read.csv(d,header=TRUE,row.names = 1))
names(all_deg_list)<-all_deg_file_call_df$ananlysis

# find pairwise overlapping DEGs between each analysis 
  # Define the threshold for down-regulated genes
  threshold <- -0.585
  
  # Filter each data frame for down-regulated genes and extract their row names
  downregulated_genes_list <- lapply(deg_list, function(df) {
    rownames(df[df$log2FoldChange < threshold, ])
  })
  
  # Function to perform pairwise comparisons
  pairwise_comparisons <- function(genes_list) {
    n <- length(genes_list)
    comparisons <- list()
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        pair_name <- paste0("Comparison_", i, "_vs_", j)
        shared_genes <- intersect(genes_list[[i]], genes_list[[j]])
        comparisons[[pair_name]] <- shared_genes
      }
    }
    
    return(comparisons)
  }
  
  # Perform pairwise comparisons
  pairwise_results <- pairwise_comparisons(downregulated_genes_list)
  
  # Print the results
  print(pairwise_results)
  
  # $Comparison_2_vs_3 is essentially comparing down-regulated genes in bronchial ~ BAL-Eos%>1% (True/false, br.bal.eosp.mt1) and nasal ~ blood AEC > 100 (True/false, n.bld.aec.mt100)


# volcano plots
for(i in seq_along(all_deg_list)){
  res<-all_deg_list[[i]]
  res<-res%>%mutate(deg_sig=ifelse(padj<0.05&log2FoldChange>0.585,'cyan',
                                   ifelse(padj<0.05&log2FoldChange< -0.585,'magenta',
                                          'grey')))
  
  keyvals <- ifelse(res$padj<0.05&res$log2FoldChange>0.585,'cyan',
                    ifelse(res$padj<0.05&res$log2FoldChange< -0.585,'magenta',
                           'grey'))
  
  names(keyvals)[keyvals == 'cyan'] <- 'up'
  names(keyvals)[keyvals == 'magenta'] <- 'down'
  names(keyvals)[keyvals == 'grey'] <- 'nonsig'
  
  title= names(all_deg_list[i])
  labs_up= rownames(filter(res,padj<0.05,log2FoldChange>0)%>%arrange(desc(log2FoldChange)))[1:10]
  labs_down= rownames(filter(res,padj<0.05,log2FoldChange<0)%>%arrange(desc(abs(log2FoldChange))))[1:10]
  labs=c(labs_up,labs_down,"CXCR1","CXCR2","CXCR4")
  p<-EnhancedVolcano(res,
                     lab = rownames(res),
                     selectLab = labs,
                     title=title,
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlab = bquote(~Log[2]~ 'fold change'),
                     xlim=c(min(res$log2FoldChange),max(res$log2FoldChange)),
                     ylim=c(0,-log(min(res$padj),10)),
                     pCutoff = 5e-2,
                     FCcutoff = 0.585,
                     cutoffLineType = 'twodash',
                     cutoffLineWidth = 0.8,
                     pointSize = 4.0,
                     labSize = 3,
                     colAlpha = 0.4,
                     colCustom = keyvals,
                     legendPosition = 'right',
                     legendLabSize = 10,
                     legendIconSize = 5.0,    
                     drawConnectors = TRUE,
                     widthConnectors = 0.75)
  assign(x = paste0("p",i),p)
}
  
  
vp<-lapply(paste0("p",seq_along(all_deg_list)),get)
names(vp)<-names(all_deg_list)
vp_folder<-file.path(deg_folder,"volcano_2")
# Check if the directory exists, and create it if it doesn't
if (!dir.exists(vp_folder)) {
  dir.create(vp_folder, recursive = TRUE)
}

# Save each grob as a PNG file
for (i in seq_along(vp)) {
  # Create a filename for each grob
  filename <- file.path(vp_folder,paste0("volcano", names(vp)[i], ".png"))
  
  # Save the grob to a PNG file
  ggsave(
    filename = filename,
    plot = grid.draw(vp[[i]]),
    device = "png",
    width = 6, height = 6, units = "in"
  )
} 



# GO term visualization
# for AAAAI abstract, focus on the following:
### "GO-05_sig-bronch-deg~bal_Eos_p_more_1_down",
### "GO-06_sig-bronch-deg~bal_Eos_p_more_1_up",
### "GO-11_sig-nasal-deg~blood_Eos_p_cont_down",
### "GO-12_sig-nasal-deg~blood_Eos_p_cont_up",
### "GO-17_nasal~bld_AEC_mt100_down",
### "GO-18_nasal~bld_AEC_mt100_up",
### "GO-21_overlap-sig-nasal~bld_AEC_mt100_down-bronch~bal_eos_p_mt1_down"

go_folder<-file.path(deg_folder,"go")
go_folder_files<-list.files(go_folder)
go_files<-file.path(go_folder,go_folder_files[grep(".txt",go_folder_files)])
go_deg_terms<-lapply(go_files,read.table,sep="\t",header=TRUE)

# Filter each element of go_deg_terms to remove rows where Category is "GOTERM_CC_DIRECT"
go_deg_terms_filtered <- lapply(go_deg_terms, function(df) {
  df %>% filter(Category != "GOTERM_CC_DIRECT")
})

lapply(go_deg_terms,head)
lapply(go_deg_terms,colnames)

go_analysis<-c(
  "GO-01_sig-bronch-deg~bal_AEC_mt_1.2_down",
  "GO-02_sig-bronch-deg~bal_AEC_mt_1.2_up",
  "GO-03_sig-bronch-deg~bal_AEC_mt_1_down",
  "GO-04_sig-bronch-deg~bal_AEC_mt_1_up",
  "GO-05_sig-bronch-deg~bal_Eos_p_more_1_down",
  "GO-06_sig-bronch-deg~bal_Eos_p_more_1_up",
  "GO-07_sig-bronch-deg~bal_Eos_p_more_3_down",
  "GO-08_sig-bronch-deg~bal_Eos_p_more_3_up",
  "GO-09_sig-bronch-deg~bld_AEC_more_300_down",
  "GO-10_sig-bronch-deg~bld_AEC_more_300_up",
  "GO-11_sig-nasal-deg~blood_Eos_p_cont_down",
  "GO-12_sig-nasal-deg~blood_Eos_p_cont_up",
  "GO-13_sig_bronch~deg_bal_ANCmt0_cont_down",
  "GO-14_sig_bronch~deg_bal_ANCmt0_cont_up",
  "GO-15_wgcna-bronch_magenta_module-overlap_br-bal-eosp-mt1",
  "GO-16_wgcna-bronch_magenta_module",
  "GO-17_nasal~bld_AEC_mt100_down",
  "GO-18_nasal~bld_AEC_mt100_up",
  "GO-19_sig_bronch_deg~bal_anc_mt13_down_none",
  "GO-20_sig_bronch_deg~bal_anc_mt13_up_none",
  "GO-21_overlap-sig-nasal~bld_AEC_mt100_down-bronch~bal_eos_p_mt1_down"
)
names(go_deg_terms)<-go_analysis

# select relevant columns column 'X' contains the descriptive GO terms. Then filter and arranged based on FDR
go_deg_terms<-lapply(go_deg_terms,function(d){
  d[,c("X","Genes","Fold.Enrichment","FDR")]%>%
    filter(FDR<0.05)%>%
    arrange(FDR)})

# factorize and set levels of descriptive GO terms for the plot
gt<-lapply(go_deg_terms,
           function(data){
             gt<-factor(data$X,levels=data$X)
           })

# replace the go term table with descriptive GO terms factorized 
for(i in 1:length(go_deg_terms)){
  go_deg_terms[[i]]$X<-gt[[i]]
}

# y axis label for GO terms
wrapped_label<-lapply(go_deg_terms,
                      function(data){
                        go_term<-data$X
                        str_wrap(go_term, width=40)})

# plot FDR vs ontology (fill by FC)
for(i in seq_along(go_deg_terms)){
  p<-ggplot(go_deg_terms[[i]][1:10,], aes(fill = Fold.Enrichment, y = X, x = -log10(FDR))) + # show only top 15
    geom_bar(stat = "identity") +
    #geom_label(aes(label = round(Fold.Enrichment, 1)), fill="white",nudge_y=0.3, hjust = -0.1, size = 3, color = "black") +  # Add text labels for log10(FDR)
    scale_fill_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
    labs(fill = "Fold Enrichment", y = "Gene Ontology Term", x = "-log10(FDR)", title = names(go_deg_terms[i])) +
    theme(axis.text.x = element_text(size = 10),  # Change size of x-axis labels
          axis.text.y = element_text(size=12),
          axis.title= element_text(size=12),
          legend.title=element_text(size=10),
          title = element_text(size=12))+  # Adjust y-axis label size for better readability
    scale_y_discrete(labels=wrapped_label[[i]])+
    xlim(0,max(-log(go_deg_terms[[i]]$FDR,10),na.rm=TRUE))
  
  assign(paste("GO", i,"by_FDR", sep = "_"), p)
}

# Arrange the plots in a 2x1 grid
grobs<-lapply(paste("GO", seq_along(go_deg_terms), "by_FDR",sep = "_"),get)

gp_3<-grid.arrange(grobs = grobs[5:6], nrow = 2)
gp_6<-grid.arrange(grobs = grobs[11:12], nrow = 2)
gp_9<-grid.arrange(grobs = grobs[17:18], nrow = 2)
gp_11<-grid.arrange(grobs = grobs[21], nrow = 2)

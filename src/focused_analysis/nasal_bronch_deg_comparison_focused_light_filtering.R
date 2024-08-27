# DEG volcano plot and GO term visualization

focused_list_all<-c("bronch_res_all_16_~ bal_Eos_p_more_1",
  "nasal_res_all_7_~ blood_eos_p_log",
  "nasal_res_all_15_~ bld_AEC_more_500")

focused_list<-c("bronch_res_sig_16_~ bal_Eos_p_more_1",
                "nasal_res_sig_7_~ blood_eos_p_log",
                "nasal_res_sig_15_~ bld_AEC_more_500")

library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(gridExtra)
library(grid)

deg_folder<-file.path("./resources/processed_data/focused_analysis_light_filter")
deg_file<-list.files(deg_folder)
matching_files <- deg_file[sapply(deg_file, function(file) any(sapply(focused_list, grepl, file)))]
deg_file<-matching_files

all_deg_folder<-file.path(deg_folder,"deg_all")
all_deg_file<-list.files(all_deg_folder)
matching_files <- all_deg_file[sapply(all_deg_file, function(file) any(sapply(focused_list_all, grepl, file)))]
all_deg_file<-matching_files


analysis_list<-c("br.bal.eosp.mt1",
                 "n.bld.aec.mt500",
                 "n.bld.eosp.cont")

# Tables containing the gene names, l2fc, padj of only DEGs in each analysis
deg_file_call_df<-data.frame(ananlysis=analysis_list,file=deg_file)
deg_list<-lapply(file.path(deg_folder,deg_file),function(d)read.csv(d,header=TRUE,row.names = 1))
names(deg_list)<-deg_file_call_df$ananlysis


# Find # of DEGs abs(log2FC) > 1 (2 fold difference) and Padj < 0.05
print(lapply(deg_list,function(d)filter(d,abs(log2FoldChange)>1)%>%dim))

# Tables containing the gene names, l2fc, padj of all genes analyzed in each analysis
all_deg_file_call_df<-data.frame(ananlysis=analysis_list,file=all_deg_file)
all_deg_list<-lapply(file.path(all_deg_folder,all_deg_file),function(d)read.csv(d,header=TRUE,row.names = 1))
names(all_deg_list)<-all_deg_file_call_df$ananlysis

# find pairwise overlapping DEGs between each analysis 
  # Define the threshold for down-regulated genes
  threshold <- -1
  
  # Filter each data frame for down-regulated genes and extract their row names
  downregulated_genes_list <- lapply(deg_list, function(df) {
    df[df$log2FoldChange < threshold, "X" ]
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

FC_thr<-1
# volcano plots
for(i in seq_along(all_deg_list)){
  res<-all_deg_list[[i]]
  res<-res%>%mutate(deg_sig=ifelse(padj<0.05&log2FoldChange>FC_thr,'cyan',
                                   ifelse(padj<0.05&log2FoldChange< -FC_thr,'magenta',
                                          'grey')))
  res_sig<-res%>%filter(padj<0.05)
  keyvals <- ifelse(res$padj<0.05&res$log2FoldChange>FC_thr,'cyan',
                    ifelse(res$padj<0.05&res$log2FoldChange< -FC_thr,'magenta',
                           'grey'))
  
  names(keyvals)[keyvals == 'cyan'] <- 'up'
  names(keyvals)[keyvals == 'magenta'] <- 'down'
  names(keyvals)[keyvals == 'grey'] <- 'nonsig'
  
  title= names(all_deg_list[i])
  labs_up= rownames(filter(res,padj<0.05,log2FoldChange>0)%>%arrange(desc(log2FoldChange)))[1:10]
  labs_down= rownames(filter(res,padj<0.05,log2FoldChange<0)%>%arrange(desc(abs(log2FoldChange))))[1:10]
  labs=c(labs_up,labs_down)
  x_lim<-c(min(res_sig$log2FoldChange),max(res_sig$log2FoldChange))
  p<-EnhancedVolcano(res,
                     lab = rownames(res),
                     selectLab = labs,
                     title=title,
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlab = bquote(~Log[2]~ 'fold change'),
                     xlim=x_lim,
                     ylim=c(0,-log(min(res_sig$padj),10)),
                     pCutoff = 5e-2,
                     FCcutoff = FC_thr,
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
vp_folder<-file.path(deg_folder,"volcano")
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

go_folder<-file.path(deg_folder,"go")
go_folder_files<-list.files(go_folder)
go_files<-file.path(go_folder,go_folder_files[grep(".txt",go_folder_files)])
go_deg_terms<-lapply(go_files,read.table,sep="\t",header=TRUE)

go_deg_terms<-lapply(go_deg_terms,function(d)separate(d,Term, into = c("Term", "Term_description"), sep = "~"))


# Extract the file names from the paths
file_names <- basename(go_files)

# Remove the '.txt' extension
go_analysis <- sub("\\.txt$", "", file_names)
# Print the results
print(go_analysis)

# Filter each element of go_deg_terms to remove rows where Category is "GOTERM_CC_DIRECT"
go_deg_terms_filtered <- lapply(go_deg_terms, function(df) {
  df %>% filter(Category != "GOTERM_CC_DIRECT")
})

lapply(go_deg_terms,head)
lapply(go_deg_terms,colnames)

names(go_deg_terms)<-go_analysis

# select relevant columns column 'Term_description' contains the descriptive GO terms. Then filter and arranged based on FDR
go_deg_terms<-lapply(go_deg_terms,function(d){
  d[,c("Term_description","Genes","Fold.Enrichment","FDR")]%>%
    filter(FDR<0.05)%>%
    arrange(FDR)})

# factorize and set levels of descriptive GO terms for the plot
gt<-lapply(go_deg_terms,
           function(data){
             gt<-factor(data$Term_description,levels=data$Term_description)
           })

# replace the go term table with descriptive GO terms factorized 
for(i in 1:length(go_deg_terms)){
  go_deg_terms[[i]]$Term_description<-gt[[i]]
}

# y axis label for GO terms
wrapped_label<-lapply(go_deg_terms,
                      function(data){
                        go_term<-data$Term_description
                        str_wrap(go_term, width=40)})

# plot FDR vs ontology (fill by FC)
for(i in seq_along(go_deg_terms)){
  p<-ggplot(go_deg_terms[[i]][1:10,], aes(fill = Fold.Enrichment, y = Term_description, x = -log10(FDR))) + # show only top 15
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
names(grobs)<- names(go_deg_terms)
# examine the go plots of interest
grid.arrange(grobs = grobs[1:2], nrow = 2)
grid.arrange(grobs = grobs[3:4], nrow = 2)


# select GO plot index 
go_plot_index<-c(1:4)

# make dir for GO plots
go_plot_folder<-file.path(deg_folder,"go_plot")
# Check if the directory exists, and create it if it doesn't
if (!dir.exists(go_plot_folder)) {
  dir.create(go_plot_folder, recursive = TRUE)
}

# Save each grob as a PNG file
for (i in go_plot_index) {
  # Create a filename for each grob
  filename <- file.path(go_plot_folder,paste0("GO_plot", names(grobs)[i], ".png"))
  
  # Save the grob to a PNG file
  ggsave(
    filename = filename,
    plot = grid.draw(grobs[[i]]),
    device = "png",
    width = 10, height = 6, units = "in"
  )
} 

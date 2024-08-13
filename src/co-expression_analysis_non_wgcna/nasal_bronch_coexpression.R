####
# Bronch = BAL-ANC (cont var, ANC>0) ->   Increased Neutrophilic chemotaxis, response to LPS, complement
# Bronch = BAL-AEC (dich var, AEC>1.2) ->   Decreased Neutrophilic chemotaxis, response to LPS
# Bronch = BAL-Eos%(dich var, Eos%>3%) -> Decreased response to LPS, dendritic chemotaxis
# Bronch = BLD-AEC (dich var, AEC>300)-> Increased O-glycan, Deuterosome 
# Nasal =  BLD-AEC (cont var, AEC>0)(old)  ->  Increased Cilia function genes 
# Nasal =  BLD-Eos% (cont var, all) ->  Increased Cilia functino genes, Decreased Neutrophil genes 
library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(gridExtra)
library(grid)

deg_folder<-file.path("./resources/processed_data/focused_analysis")
all_deg_folder<-file.path(deg_folder,"deg_all")
deg_file<-list.files(deg_folder)[grep(".csv",list.files(deg_folder))]
all_deg_file<-list.files(all_deg_folder)[grep(".csv",list.files(all_deg_folder))]

analysis_list<-c("br.bal.anc.cont",
                 "br.bld.aec.mt300",
                 "br.bal.aec.mt1",
                 "br.bal.aec.mt1.2",
                 "br.bal.eosp.mt1",
                 "br.bal.eosp.mt3",
                 "br.bal.anc.mt13",
                 "n.bld.aec.mt100",
                 "n.bld.eosp.cont")

deg_file_call_df<-data.frame(ananlysis=analysis_list,file=deg_file)
all_deg_file_call_df<-data.frame(ananlysis=analysis_list,file=all_deg_file)

deg_list<-lapply(file.path(deg_folder,deg_file),function(d)read.csv(d,header=TRUE,row.names = 1))
all_deg_list<-lapply(file.path(all_deg_folder,all_deg_file),function(d)read.csv(d,header=TRUE,row.names = 1))
names(deg_list)<-deg_file_call_df$ananlysis
names(all_deg_list)<-all_deg_file_call_df$ananlysis




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
} # too small for vp[[8]] to display CXCR2 in it. consider increasing it. 

lapply(vp,print)
grid.arrange(vp[[1]],vp[[2]],vp[[3]],nrow=3,ncol=1)




# GO term 
go_folder<-file.path(deg_folder,"go")
go_folder_files<-list.files(go_folder)
go_files<-file.path(go_folder,go_folder_files[grep(".txt",go_folder_files)])
go_deg_terms<-lapply(go_files,read.table,sep="\t",header=TRUE)
lapply(go_deg_terms,head)
lapply(go_deg_terms,colnames)

go_analysis<-c("GO-sig-bronch-deg~bal_AEC_mt_1.2_down",
               "GO-sig-bronch-deg~bal_AEC_mt_1.2_up",
               "GO-sig-bronch-deg~bal_AEC_mt_1_down",
               "GO-sig-bronch-deg~bal_AEC_mt_1_up",
               "GO-sig-bronch-deg~bal_Eos_p_more_1_down",
               "GO-sig-bronch-deg~bal_Eos_p_more_1_up(1)",
               "GO-sig-bronch-deg~bal_Eos_p_more_3_down",
               "GO-sig-bronch-deg~bal_Eos_p_more_3_up(1)",
               "GO-sig-bronch-deg~bld_AEC_more_300_down",
               "GO-sig-bronch-deg~bld_AEC_more_300_up",
               "GO-sig-nasal-deg~blood_Eos_p_cont_down",
               "GO-sig-nasal-deg~blood_Eos_p_cont_up",
               "GO-sig_bronch~deg_bal_ANCmt0_cont_down",
               "GO-sig_bronch~deg_bal_ANCmt0_cont_up",
               "GO-wgcna-bronch_magenta_module-overlap_br-bal-eosp-mt1",
               "GO-wgcna-bronch_magenta_module",
               "go_nasal~bld_AEC_mt100_down",
               "go_nasal~bld_AEC_mt100_up",  
               "go_sig_bronch_deg~bal_anc_mt13_down_none",
               "go_sig_bronch_deg~bal_anc_mt13_up_none")


names(go_deg_terms)<-go_analysis

# analysis-Terms-Genes
  Go_Term_Gene_Table<-function(go_analysis){
    # Step 1: Load the data (assuming the file is a tab-separated text file)
    data <- go_deg_terms[[go_analysis]]
    data <- filter(data,FDR<0.05)
    
    # Check if the data has zero rows
    if (nrow(data) == 0) {
      # Return a data frame with NULL values
      return(data.frame(analysis = NULL, Term = NULL, Genes = NULL))
    }
    
    # Step 2: Split the "Genes" column by commas
    genes_split <- strsplit(as.character(data$Genes), ",")
    
    # Step 3: Repeat the "Term" column to match the number of genes
    terms_repeated <- rep(data$X, sapply(genes_split, length))
    
    # Step 4: Flatten the list of genes into a single vector
    genes_vector <- unlist(genes_split)
    
    # Step 5: Combine the vectors into a new data frame
    result_df <- data.frame(analysis=go_analysis, Term = terms_repeated, Genes = genes_vector)
    
    # Remove spaces from the "Genes" column
    result_df$Genes <- gsub(" ", "", result_df$Genes)
    
    # Display the result
    print(result_df)
    
    return(result_df)
  }
  
  # Apply the function to each element of go_analysis and combine results
  go_result_list <- lapply(go_analysis, Go_Term_Gene_Table)
  
  # Combine all the resulting data frames into one
  final_result_df <- do.call(rbind, go_result_list)
  
  # View the combined result
  print(final_result_df)
  
  write.table(final_result_df,file.path(go_folder,"final_go_df","go_genes_final_df.txt"),sep="\t",col.names=TRUE,row.names=FALSE)
  
  final_result_df%>%filter(grepl("CXCR1",Genes))%>%reframe(unique(analysis))
  final_result_df%>%filter(grepl("CXCR2",Genes))%>%reframe(unique(analysis))
  final_result_df%>%filter(grepl("CXCR4",Genes))%>%reframe(unique(analysis))
  
# select relevant columns
  go_deg_terms<-lapply(go_deg_terms,function(d){
    d[,c("X","Genes","Fold.Enrichment","FDR")]%>%
      filter(FDR<0.05)%>%
      arrange(desc(Fold.Enrichment))})
  
  gt<-lapply(go_deg_terms,
             function(data){
               gt<-factor(data$X,levels=data$X)
             })

  for(i in 1:length(go_deg_terms)){
    go_deg_terms[[i]]$X<-gt[[i]]
  }

  wrapped_label<-lapply(go_deg_terms,
                        function(data){
                          go_term<-data$X
                          str_wrap(go_term, width=40)})
  # Plot
  
  for(i in seq_along(go_deg_terms)){
    p<-ggplot(go_deg_terms[[i]][1:10,], aes(x = Fold.Enrichment, y = X, fill = -log10(FDR))) + # show only top 15
      geom_bar(stat = "identity") +
      geom_label(aes(label = round(-log10(FDR), 1)), fill="white",nudge_y=0.3, hjust = -0.1, size = 3, color = "black") +  # Add text labels for log10(FDR)
      scale_fill_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
      labs(x = "Fold Enrichment", y = "Gene Ontology Term", fill = "-log10(FDR)", title = names(go_deg_terms[i])) +
      theme(axis.text.x = element_text(size = 10),  # Change size of x-axis labels
            axis.text.y = element_text(size=12),
            axis.title= element_text(size=12),
            legend.title=element_text(size=10),
            title = element_text(size=12))+  # Adjust y-axis label size for better readability
      scale_y_discrete(labels=wrapped_label[[i]])+
      xlim(0,max(go_deg_terms[[i]]$Fold.Enrichment,na.rm=TRUE)+5)
    
    assign(paste("GO", i, sep = "_"), p)
  }
  grobs<-lapply(paste("GO", seq_along(go_deg_terms), sep = "_"),get)
  # Arrange the plots in a 2x4 grid
  gp_1<-grid.arrange(grobs = grobs[1:2], nrow = 2)
  gp_2<-grid.arrange(grobs = grobs[3:4], nrow = 2)
  gp_3<-grid.arrange(grobs = grobs[5:6], nrow = 2)
  gp_4<-grid.arrange(grobs = grobs[7:8], nrow = 2)
  gp_5<-grid.arrange(grobs = grobs[9:10], nrow = 2)
  gp_6<-grid.arrange(grobs = grobs[11:12], nrow = 2)
  gp_7<-grid.arrange(grobs = grobs[13:14], nrow = 2)
  gp_8<-grid.arrange(grobs = grobs[15:16], nrow = 2)
  gp_8<-grid.arrange(grobs = grobs[17:18], nrow = 2)
  gp_9<-grid.arrange(grobs = grobs[19:20], nrow = 2)

  go_enrich_plot<-lapply(paste("gp", 1:8, sep = "_"),get)
  
  go_plot_folder<-file.path(go_folder,"go_plots")
  # Check if the directory exists, and create it if it doesn't
  if (!dir.exists(go_plot_folder)) {
    dir.create(go_plot_folder, recursive = TRUE)
  }
  
  # Save each grob as a PNG file
  for (i in seq_along(go_enrich_plot)) {
    # Create a filename for each grob
    filename <- file.path(go_plot_folder,paste0("go_enrichment_plot",i, ".png"))
    
    # Save the grob to a PNG file
    ggsave(
      filename = filename,
      plot = grid.draw(go_enrich_plot[[i]]),
      device = "png",
      width = 6, height = 6, units = "in"
    )
  } 
  

############
source("./src/function/deg_custom_functions.R")

countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
rownames(counts)<-counts[,"SampleID"]
bronch.samples<-grepl("^B",colnames(counts))
bronch.samples<-which(bronch.samples==TRUE)
bronch.counts<-counts[,c(1,bronch.samples)]
colnames(bronch.counts)[bronch.samples]<-substr(colnames(bronch.counts)[bronch.samples],1,4)
head(bronch.counts)

filtered_samples<-counts[,!nchar(colnames(counts))>4]
filtered_ID<-grepl("^N",colnames(filtered_samples))
ncounts<-filtered_samples[,filtered_ID]

bronch_sid<-colnames(bronch.counts)[-1]
nasal_sid<-colnames(ncounts)

bronch_sid_filtered<-bronch_sid[substr(bronch_sid,2,4)%in%substr(nasal_sid,2,4)]
nasal_sid_filtered<-nasal_sid[substr(nasal_sid,2,4)%in%substr(bronch_sid,2,4)]
BAL_samples<-data.frame(bronch_sid_filtered,nasal_sid_filtered)

ct_n<-ncounts[,BAL_samples$nasal_sid_filtered]
ct_b<-bronch.counts[,BAL_samples$bronch_sid_filtered]

# normalize counts

norm_ct<-function(readcounts){
  # normalize counts LCPM
  x<-readcounts
  ### normalize counts with TMM
  norm.factor<-calcNormFactors(x, method = "TMM")
  sample.size<-length(colnames(x))
  for(i in 1:sample.size){
    x[,i]<-x[,i]/norm.factor[i]
  }
  ### calculate lcpm based on TMM normalized counts 
  lcpm.x<-cpm(x,log=TRUE)
  return(lcpm.x)
}
ct_n<-norm_ct(ct_n)
ct_b<-norm_ct(ct_b)

ct_n_cxcr<-ct_n[c(which(rownames(ct_n)=="CXCR1"),which(rownames(ct_n)=="CXCR2"),which(rownames(ct_n)=="CXCR4")),]
ct_b_cxcr<-ct_b[c(which(rownames(ct_b)=="CXCR1"),which(rownames(ct_b)=="CXCR2"),which(rownames(ct_b)=="CXCR4")),]




cor(unlist(ct_n_cxcr[1,]),unlist(ct_b_cxcr[1,]))
cor(unlist(ct_n_cxcr[2,]),unlist(ct_b_cxcr[2,]))
cor(unlist(ct_n_cxcr[3,]),unlist(ct_b_cxcr[3,]))

plot(unlist(ct_n_cxcr[1,]),unlist(ct_b_cxcr[1,]))

plot(unlist(ct_n_cxcr[2,]),unlist(ct_b_cxcr[2,]))

plot(unlist(ct_n_cxcr[3,]),unlist(ct_b_cxcr[3,]))

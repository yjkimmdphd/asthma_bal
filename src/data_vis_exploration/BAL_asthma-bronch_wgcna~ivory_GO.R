###########################
## GO analysis of genes that overlap between Bronchial DEG ~ BAL Eos% >1 && Bronchial WGCNA Ivory module  
## be sure to have run ./src/data_processing/go_david_filtering.py before
###########################
library(tidyverse)
go_deg_folder<-file.path("./reports/local_only/wgcna/bronch/output/module_GO_filtered")

go_deg_filelist<-if(file.exists(go_deg_folder)){list.files(path=go_deg_folder)}

# Subset elements that end with ".txt"
txt_files <- go_deg_filelist[grep("\\.txt$", go_deg_filelist)]
txt_files<-txt_files[grepl("ivory",txt_files)]
print(txt_files)

go_deg_filelist<-txt_files
go_deg_filelist<-go_deg_filelist[1:length(go_deg_filelist)]
print(go_deg_filelist)
print(file.exists(file.path(go_deg_folder,go_deg_filelist)))
###
# 1. GO term for  DEG ~ BAL
###

go_deg_terms<-lapply(file.path(go_deg_folder,go_deg_filelist),read.table,sep="\t",header=TRUE)

lapply(go_deg_terms,colnames)

# select relevant columns
go_deg_terms<-lapply(go_deg_terms,
                     function(d){
                       d[,colnames(d)%in%c("Category","GO_ID","Description","Genes","Fold.Enrichment","FDR")]%>%
                         filter(FDR<0.05)%>%
                         arrange(Fold.Enrichment)})

names(go_deg_terms)<-gsub("\\+batch12346", "", substr(go_deg_filelist, 1, nchar(go_deg_filelist) - 4))

gt<-lapply(go_deg_terms,
           function(data){
             gt<-factor(data$Description,levels=data$Description)
           })

for(i in 1:length(go_deg_terms)){
  go_deg_terms[[i]]$Description<-gt[[i]]
}

wrapped_label<-lapply(go_deg_terms,
                      function(data){
                        go_term<-data$Description
                        str_wrap(go_term, width=40)})
# Plot
library(grid)
library(gridExtra)
for(i in 1:length(go_deg_terms)){
  if(nrow(go_deg_terms[[i]])>=10){
    n<-4
  }else{
    n<-nrow(go_deg_terms[[i]])
  }
  p<-ggplot(go_deg_terms[[i]][1:n,], aes(x = Fold.Enrichment, y = Description, fill = -log10(FDR))) + # show only top 10
    geom_bar(stat = "identity") +
    # geom_label(aes(label = round(-log10(FDR), 1)), fill="white",nudge_y=0.3, hjust = -0.1, size = 3, color = "black") +  # Add text labels for log10(FDR)
    scale_fill_gradient(low = "blue", high = "blue") +  # Adjust color gradient as needed
    labs(x = "Fold Enrichment", y = "Gene Ontology Term", fill = "-log10(FDR)", title = names(go_deg_terms[i])) +
    theme(axis.text.x = element_text(size = 10),  # Change size of x-axis labels
          axis.text.y = element_text(size = 12),
          axis.title= element_text(size=12),
          legend.title=element_text(size=10),
          title = element_text(size=12))+  # Adjust y-axis label size for better readability
    scale_y_discrete(labels=wrapped_label[[i]])
  assign(paste("a", i, sep = ""), p)
}
grobs_names<-paste0("a",seq(1,length(go_deg_terms),by=1))
grobs<-lapply(grobs_names,get)

lapply(grobs,print)





###########################
## GO analysis  
# bronch expression ~ BAL cell counts (dichot) + batch12346
###########################
library(tidyverse)
library(ggplot2)
library(stringr)

go_deg_folder<-file.path("./reports/local_only/deg_bronch~bal_cell(dichot)+batch12346/GO")
go_deg_filelist<-if(file.exists(go_deg_folder)){list.files(path=go_deg_folder)}

# Subset elements that end with ".txt"
txt_files <- go_deg_filelist[grep("\\.txt$", go_deg_filelist)]
print(txt_files)

go_deg_filelist<-txt_files
print(file.exists(file.path(go_deg_folder,go_deg_filelist)))

###
# 1. GO term for  DEG ~ BAL+batch12346 
###

go_deg_terms<-lapply(file.path(go_deg_folder,go_deg_filelist),read.table,sep="\t",header=TRUE)

lapply(go_deg_terms,colnames)

# select relevant columns
go_deg_terms<-lapply(go_deg_terms,
                     function(d){
                       d[,c(2,3,6,10,13)]%>%
                         filter(FDR<0.05)%>%
                         arrange(desc(Fold.Enrichment))})

names(go_deg_terms)<-gsub("\\+batch12346", "", substr(go_deg_filelist, 1, nchar(go_deg_filelist) - 4))

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
library(grid)
library(gridExtra)
for(i in 1:length(go_deg_terms)){
  p<-ggplot(go_deg_terms[[i]][1:15,], aes(x = Fold.Enrichment, y = X, fill = -log10(FDR))) + # show only top 15
    geom_bar(stat = "identity") +
    geom_label(aes(label = round(-log10(FDR), 1)), fill="white",nudge_y=0.3, hjust = -0.1, size = 3, color = "black") +  # Add text labels for log10(FDR)
    scale_fill_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
    labs(x = "Fold Enrichment", y = "Gene Ontology Term", fill = "-log10(FDR)", title = names(go_deg_terms[i])) +
    theme(axis.text.x = element_text(size = 10),  # Change size of x-axis labels
        axis.text.y = element_text(size = 10),
        axis.title= element_text(size=15),
        legend.title=element_text(size=10),
        title = element_text(size=15))+  # Adjust y-axis label size for better readability
    scale_y_discrete(labels=wrapped_label[[i]])
  assign(paste("a", i, sep = ""), p)
}
grobs<-list(a1,a2,a3,a4,a5,a6,a7,a8)
# Arrange the plots in a 2x4 grid
grid.arrange(grobs = grobs[1:2], nrow = 2)
grid.arrange(grobs = grobs[3:4], nrow = 2)
grid.arrange(grobs = grobs[5:6], nrow = 2)
grid.arrange(grobs = grobs[7:8], nrow = 2)
grid.arrange(grobs = grobs[8], nrow = 1)



###
# 2. GO term for WGCNA ~ bronch batch12346
###

go_wgcna_terms<-lapply(file.path(go_wgcna_folder,go_wgcna_filelist),read.table,sep="\t",header=TRUE)

lapply(go_wgcna_terms,colnames)%>%unlist%>%matrix(nrow=13)

goSubset<-function(data){
  data[,c(2,3,6,10,13)]%>% # 
    filter(FDR<0.05)%>%
    arrange(desc(Fold.Enrichment))}

go_wgcna_terms<-lapply(go_wgcna_terms,goSubset)
names(go_wgcna_terms)<-substr(go_wgcna_filelist, 1, nchar(go_wgcna_filelist) - 4)

gt<-lapply(go_wgcna_terms,
           function(data){
             gt<-factor(data$X,levels=data$X)
           })

for(i in 1:length(go_wgcna_terms)){
  print(paste("success",i))
  go_wgcna_terms[[i]]$X<-gt[[i]]
  
}

wrapped_label<-lapply(go_wgcna_terms,
                      function(data){
                        go_term<-data$X
                        str_wrap(go_term, width=40)})
# Plot

library(grid)
library(gridExtra)
for(i in 1:length(gt)){
  p<-ggplot(go_wgcna_terms[[i]], aes(x = Fold.Enrichment, y = X, fill = -log10(FDR))) +
    geom_bar(stat = "identity") +
    geom_label(aes(label = round(-log10(FDR), 1)), fill="white",nudge_y=0.3, hjust = -0.1, size = 3, color = "black") +  # Add text labels for log10(FDR)
    scale_fill_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
    labs(x = "Fold Enrichment", y = "Gene Ontology Term", fill = "-log10(FDR)", title = names(go_wgcna_terms[i])) +
    theme(axis.text.x = element_text(size = 10),  # Change size of x-axis labels
          axis.text.y = element_text(size = 10),
          axis.title= element_text(size=15),
          legend.title=element_text(size=10),
          title = element_text(size=15))+  # Adjust y-axis label size for better readability
    scale_y_discrete(labels=wrapped_label[[i]])
  assign(paste("a", i, sep = ""), p)
}
grobs<-list(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
print(a1)
# Arrange the plots in a 2x4 grid
grid.arrange(grobs = grobs[1:3], nrow = 3)
grid.arrange(grobs = grobs[7], nrow = 1)
grid.arrange(grobs = grobs[8], nrow = 1)


###
# 3. GO term for  WGCNA ~ bronch batch12346 top 15 for each
# focusing on GO terms for DE genes from Bronch~BAL_ANC+batch12346 that have overlap with
# WGCNA ~ bronch batch12346 Magenta module and Tan module
###
go_wgcna_terms_top15<-go_wgcna_terms
for(i in 1:length(go_wgcna_terms_top15)){
  k<-go_wgcna_terms_top15[[i]]
  go_wgcna_terms_top15[[i]]<-k[1:15,]
}
wrapped_label<-lapply(go_wgcna_terms_top15,
                      function(data){
                        go_term<-data$X
                        str_wrap(go_term, width=40)})
# Plot

library(grid)
library(gridExtra)
for(i in 1:length(go_wgcna_terms_top15)){
  p<-ggplot(go_wgcna_terms_top15[[i]], aes(x = Fold.Enrichment, y = X, fill = -log10(FDR))) +
    geom_bar(stat = "identity") +
    geom_label(aes(label = round(-log10(FDR), 1)), fill="white",nudge_y=0.3, hjust = -0.1, size = 3, color = "black") +  # Add text labels for log10(FDR)
    scale_fill_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
    labs(x = "Fold Enrichment", y = "Gene Ontology Term", fill = "-log10(FDR)", title = names(go_wgcna_terms[i])) +
    theme(axis.text.x = element_text(size = 10),  # Change size of x-axis labels
          axis.text.y = element_text(size = 8),
          axis.title= element_text(size=15),
          legend.title=element_text(size=10),
          title = element_text(size=15))+  # Adjust y-axis label size for better readability
    scale_y_discrete(labels=wrapped_label[[i]])
  assign(paste("a", i, sep = ""), p)
}
grobs<-list(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)

# Arrange the plots in a 2x6 grid
grid.arrange(grobs = grobs[1:2], nrow = 2) 

grid.arrange(grobs = grobs[3:4], nrow = 2)

grid.arrange(grobs = grobs[5:6], nrow = 2)

grid.arrange(grobs = grobs[7:8], nrow = 2)

grid.arrange(grobs = grobs[9:10], nrow = 2)

grid.arrange(grobs = grobs[11:12], nrow = 2)

ggsave(filename = "./reports/figures/GO plot/GO_term_figure.pdf", plot = grid.arrange(grobs = grobs[1:12], nrow = 12) , width = 10, height = 40)

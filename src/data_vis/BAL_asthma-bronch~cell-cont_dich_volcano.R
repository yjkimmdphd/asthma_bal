library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(gridExtra)
library(grid)

# load deg results (all results FDR > 0) bronch ~ cell count, continuous
cont_bronch<-file.path("./reports/local_only/deg_bal_bronch~cell2025-01-03")
file_path<-cont_bronch
file_names<-list.files(file_path)
file_names<-file_names[grep("deg_bronch_res_all",file_names)]
deg_results<-lapply(file.path(file_path,file_names),function(d)read.csv(d,row.names = 1))

# Extract the string after '~' and remove the '.csv' extension to name each list element
extracted_strings <- sapply(file_names, function(x) {
  string_after_tilde <- trimws(strsplit(x, "~")[[1]][2])
  string_without_csv <- sub("\\+ Batch_2025-01-03_.csv$", "", string_after_tilde)
  return(string_without_csv)
})
names(deg_results)<-extracted_strings
lapply(deg_results,head)
####################################################################
# data exploration of DEG analysis using bronchial rnaseq data
# model: "Bronch DEG ~ cell count continuous + Batch
####################################################################

deg_results<-lapply(deg_results, function(d)
  mutate(d, deg_sig = ifelse(
    padj < 0.05 & log2FoldChange > 1,
    'cyan',
    ifelse(padj < 0.05 &
             log2FoldChange < -1, 'magenta',
           'grey')
  )))

keyvals <- lapply(deg_results,function(res)ifelse(res$padj<0.05&res$log2FoldChange>1,'cyan',
                  ifelse(res$padj<0.05&res$log2FoldChange< -1,'magenta',
                         'grey')))
# Iterate over each element in the list
for (name in names(keyvals)) {
  # Assign names to the character vector elements within each list item
  names(keyvals[[name]]) <- sapply(keyvals[[name]], function(x) {
    if (is.na(x)) {
      return(NA)
    } else if (x == "magenta") {
      return("down")
    } else if (x == "cyan") {
      return("up")
    } else if (x == "grey") {
      return("nonsig")
    } else {
      return(NA)  # Handle any other unexpected values
    }
  }, USE.NAMES = FALSE)  # Ensure that sapply doesn't auto-assign names
}

# To check the result
print(keyvals)

# how many have NA?
lapply(lapply(keyvals,is.na),sum)

all_deg_list<-deg_results

for(i in seq_along(all_deg_list)){
  res<-all_deg_list[[i]]
  kv<-keyvals[[i]]
  title= names(all_deg_list[i])
  labs_up= rownames(filter(res,padj<0.05,log2FoldChange>0)%>%arrange(desc(log2FoldChange)))[1:10]
  labs_down= rownames(filter(res,padj<0.05,log2FoldChange<0)%>%arrange(desc(abs(log2FoldChange))))[1:10]
  labs=c(labs_up,labs_down)
  p<-EnhancedVolcano(res,
                     lab = rownames(res),
                     selectLab = labs,
                     title=title,
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlab = bquote(~Log[2]~ 'fold change'),
                     xlim=c(min(res$log2FoldChange)-0.5,max(res$log2FoldChange)+0.5),
                     ylim=c(0,-log(min(res$padj),10)),
                     pCutoff = 5e-2,
                     FCcutoff = 1,
                     cutoffLineType = 'twodash',
                     cutoffLineWidth = 0.8,
                     pointSize = 4.0,
                     labSize = 3,
                     colAlpha = 0.4,
                     colCustom = kv,
                     legendPosition = 'right',
                     legendLabSize = 10,
                     legendIconSize = 5.0,    
                     drawConnectors = TRUE,
                     widthConnectors = 0.75)
  assign(x = paste0("p",i),p)
}

names(all_deg_list)


# Create a sequence of plot names
plot_names <- paste0("p", seq_along(all_deg_list))

# Retrieve the ggplot objects and save them into a list
plot_list <- mget(plot_names)


# Directory where the plots will be saved
output_dir <- "./reports/figures/deg/volcano_plot/bronch~cell(cont_or_dich)"

# Loop through the list and save each plot
for (i in 1:length(plot_list)) {
  # Construct the file name
  file_name <- paste("volcano", names(all_deg_list)[i], sep = "_")
  
  # Save the plot
  ggsave(filename = paste0(output_dir, "/", file_name, ".png"), plot = plot_list[[i]], width = 8, height = 6, dpi = 300)
}

###########################
## GO analysis  
###########################
go_deg_folder<-file.path("./reports/local_only/deg_bal_bronch~cell2025-01-03/GO_output")

go_deg_filelist<-if(file.exists(go_deg_folder)){list.files(path=go_deg_folder)}

# Subset elements that end with ".txt"
txt_files <- go_deg_filelist[grep("\\.txt$", go_deg_filelist)]
txt_files<-txt_files[grepl("filtered",txt_files)]
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
                         arrange(desc(Fold.Enrichment))})

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
  p<-ggplot(go_deg_terms[[i]][1:10,], aes(x = Fold.Enrichment, y = Description, fill = -log10(FDR))) + # show only top 10
    geom_bar(stat = "identity") +
    geom_label(aes(label = round(-log10(FDR), 1)), fill="white",nudge_y=0.3, hjust = -0.1, size = 3, color = "black") +  # Add text labels for log10(FDR)
    scale_fill_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
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


plot_list <- list()

# Loop to create and store each plot
for (i in 1:length(go_deg_terms)) {
  # Create the plot
  p <- ggplot(go_deg_terms[[i]][1:10,], aes(x = Fold.Enrichment, y = Description, fill = -log10(FDR))) +
    geom_bar(stat = "identity") +
    geom_label(aes(label = round(-log10(FDR), 1)), fill = "white", nudge_y = 0.3, hjust = -0.1, size = 3, color = "black") +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(x = "Fold Enrichment", y = "Gene Ontology Term", fill = "-log10(FDR)", title = names(go_deg_terms[i])) +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 9.5),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 10),
          plot.title = element_text(size = 12)) +
    scale_y_discrete(labels = wrapped_label[[i]])
  
  # Add the plot to the list
  plot_list[[i]] <- p
}
# Determine how many plots go into each figure
third_length <- ceiling(length(plot_list) / 3)
plot_list1 <- plot_list[1:third_length]
plot_list2 <- plot_list[(third_length + 1):(2 * third_length)]
plot_list3 <- plot_list[(2 * third_length + 1):length(plot_list)]

# Arrange each list in a separate figure with a single column layout
grid.arrange(grobs = plot_list1, ncol = 1)  # First figure
grid.arrange(grobs = plot_list2, ncol = 1)  # Second figure
grid.arrange(grobs = plot_list3, ncol = 1)  # Third figure


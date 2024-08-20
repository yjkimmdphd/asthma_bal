library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(gridExtra)
library(grid)

# load deg results (all results FDR > 0) nasal ~ cell count, continuous and dichotomous
nasal<-file.path("./reports/local_only/deg_nasal~cell(cont_and_dich)+batch12346/ge0_deg_2024-08-20")
file_path<-nasal
file_names<-list.files(file_path)
file_names<-file_names[grep("res_all",file_names)]
deg_results<-lapply(file.path(file_path,file_names),function(d)read.csv(d,row.names = 1))
nasal<-deg_results
head(nasal[1])

# Extract the string after '~' and remove the '.csv' extension to name each list element
extracted_strings <- sapply(file_names, function(x) {
  string_after_tilde <- trimws(strsplit(x, "~")[[1]][2])
  string_without_csv <- sub("\\_.csv$", "", string_after_tilde)
  return(string_without_csv)
})
names(nasal)<-extracted_strings
names(nasal)
####################################################################
# data exploration of DEG analysis using bronchial rnaseq data
# model: "nasal DEG ~ cell count cont or dichotomous + Batch
####################################################################



nasal<-lapply(nasal, function(d)
  mutate(d, deg_sig = ifelse(
    padj < 0.05 & log2FoldChange > 1,
    'cyan',
    ifelse(padj < 0.05 &
             log2FoldChange < -1, 'magenta',
           'grey')
  )))

keyvals <- lapply(nasal,function(res)ifelse(res$padj<0.05&res$log2FoldChange>1,'cyan',
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

all_deg_list<-nasal

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
  assign(x = paste0("pn",i),p)
}

# Assuming p1 to p10 are your ggplot objects stored in a list
plots <- list(pn1, pn2, pn3, pn4, pn5, pn6, pn7, pn8, pn9, pn10, pn11,pn12,pn13,pn14,pn15,pn16)

# Names vector corresponding to dich_bronch
names_nasal <- names(nasal)

# Directory where the plots will be saved
output_dir <- "./reports/figures/deg/volcano_plot/nasal~cell(cont_or_dich)"


# Loop through the list and save each plot
for (i in 1:length(plots)) {
  # Construct the file name
  file_name <- paste("volcano","nasal", names_nasal[i], sep = "_")
  
  # Save the plot
  ggsave(filename = paste0(output_dir, "/", file_name, ".png"), plot = plots[[i]], width = 8, height = 6, dpi = 300)
}

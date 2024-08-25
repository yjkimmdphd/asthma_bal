library(dplyr)

# set DEG results file folder 
res_path<-file.path("./reports/temporary/deg_nasal_2024-08-24") 
res_file_name<-list.files(res_path)
res_index<-grep("res_sig",res_file_name)
res_file<-file.path(res_path,res_file_name[res_index])

# read the DEG results 
deg_res<-lapply(res_file,read.csv)
names(deg_res)<-res_file_name[res_index]

# set DEG filters
fdr <- 0.05  # Adjusted significance threshold
l2fc <- 1    # Log2 fold change threshold

filtered_res <- lapply(deg_res, function(res) {
  res_filtered <- res %>%
    dplyr::filter(padj < fdr, abs(log2FoldChange) > l2fc)
  return(res_filtered)
})

# save the results as new csv files 
filtered_res_file_name<-paste("l2fc_mt1",res_file_name[res_index],sep="_")
for(i in seq_along(filtered_res)){
  write.csv(filtered_res[[i]],file.path(res_path,filtered_res_file_name[i]))
}
  



# 
# subsetting the network edge list by modules
# 
library(tidyverse)
network_dir<-file.path("./reports/local_only/wgcna/bronch/output")
edge_list<-read.delim(file.path(network_dir,"wgcna_bronch_edge_list.txt"),sep="\t",header=TRUE)
node_att<-read.delim(file.path(network_dir,"wgcna_bronch_node_attributes.txt"),sep="\t",header=TRUE)
head(node_att)
colnames(node_att)[1]<-"Source"
left_join(edge_list,node_att, by = "Source")

# Split the data frame by 'Module'
module_list <- split(new_edge_list, new_edge_list$Module)

# Directory to save the files (change as needed)
output_dir <-"./reports/local_only/wgcna/bronch/output/modules_output"
if(!dir.exists(output_dir)){dir.create(output_dir, showWarnings = FALSE)}

# Loop through each subset and save it as a tab-delimited file
for (module_name in names(module_list)) {
  file_name <- paste0(output_dir,"/", module_name, ".txt")  # Create file name

  write.table(module_list[[module_name]], file = file_name, sep = "\t",
              row.names = FALSE, quote = FALSE)
}

cat("Files saved in", output_dir, "\n")


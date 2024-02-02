##################################################################
# exploring the DEG results of BAL project with volcano plot
##################################################################
library(dplyr)
library(DESeq2)
#######################
# loading DEG results #
#######################
# Directory path where you want to search
directory_path <- "./reports/deg_result"


sort_filename_bynumbers<-function(file_names){
  # Extracting the numbers following 'res'
  extracted_numbers <- sapply(file_names, function(file_name) {
    extracted <- regmatches(file_name, regexpr("(?<=res)\\d+", file_name, perl = TRUE))
    as.numeric(extracted)
  })
  
  # Sorting file names based on the extracted numbers
  sorted_file_names <- file_names[order(extracted_numbers)]
  return(sorted_file_names)
}


# List all files starting with "deg_res" in the results folder and subdirectories
deg_list1 <- list.files(
  path = directory_path,
  pattern = "^deg_nasal_allcell_res",  # Pattern to match filenames starting with "deg_nasal_allcell_res"
  full.names = TRUE      # Return full file paths
)%>%sort_filename_bynumbers()

deg_list2 <- list.files(
  path = directory_path,
  pattern = "^deg_nasal_poscell_res",  # Pattern to match filenames starting with "deg_nasal_poscell_res"
  full.names = TRUE      # Return full file paths
)%>%sort_filename_bynumbers()

deg_list3 <- list.files(
  path = directory_path,
  pattern = "^deg_bronch_allcell_res",  # Pattern to match filenames starting with "deg_bronch_allcell_res"
  full.names = TRUE      # Return full file paths
)%>%sort_filename_bynumbers()

deg_list4 <- list.files(
  path = directory_path,
  pattern = "^deg_bronch_poscell_res",  # Pattern to match filenames starting with "deg_bronch_poscell_res"
  full.names = TRUE      # Return full file paths
)%>%sort_filename_bynumbers()

############
# List all files starting with "deg_analysis_input_" in the results folder and subdirectories
deg_analysis_input1 <- list.files(
  path = directory_path,
  pattern = "^deg_nasal_allcell_analysis",  # Pattern to match filenames starting with "deg_nasal_allcell_analysis"
  full.names = TRUE      # Return full file paths
)

deg_analysis_input2 <- list.files(
  path = directory_path,
  pattern = "^deg_nasal_poscell_analysis",  # Pattern to match filenames starting with "deg_nasal_poscell_analysis"
  full.names = TRUE      # Return full file paths
)

deg_analysis_input3 <- list.files(
  path = directory_path,
  pattern = "^deg_bronch_allcell_analysis",  # Pattern to match filenames starting with "deg_bronch_allcell_analysis"
  full.names = TRUE      # Return full file paths
)

deg_analysis_input4 <- list.files(
  path = directory_path,
  pattern = "^deg_bronch_poscell_analysis",  # Pattern to match filenames starting with "deg_bronch_poscell_analysis"
  full.names = TRUE      # Return full file paths
)

# Display the list of files found
print(deg_list1); print(deg_list2); print(deg_list3); print(deg_list4)
print(deg_analysis_input1); print(deg_analysis_input2); print(deg_analysis_input3); print(deg_analysis_input4)

# read the analysis input files and extract results, design, samples, source_cell
analysis_input1<-read.csv(deg_analysis_input1)[,c("results","design", "source_cell")]
analysis_input2<-read.csv(deg_analysis_input2)[,c("results","design", "source_cell")]
analysis_input3<-read.csv(deg_analysis_input3)[,c("results","design", "source_cell")]
analysis_input4<-read.csv(deg_analysis_input4)[,c("results","design", "source_cell")]

source.cell<-c("BAL_eos_ct",
               "BAL_eos_p",
               "BAL_neut_ct",
               "BAL_neut_p",
               "BAL_wbc",
               "blood_eos",
               "blood_eos_p",
               "blood_neut",
               "blood_neut_p",
               "blood_wbc")
############

# deg_list1= nasal, all cells
# deg_list2= nasal, >0 cells
# deg_list3= bronchial, all cells
# deg_list4= bronchial, >0 cells 

# read the results list 
x<-lapply(deg_list1,read.csv)
genes1 <- lapply(x, function(df)genes=df$X)
genes1<-lapply(genes1,function(list){if(class(list)=="logical")list<-NA else(list)})
genes1_names<-paste("nasal","all",source.cell,sep="_")
names(genes1)<-genes1_names


x<-lapply(deg_list2,read.csv)
genes2 <- lapply(x, function(df)genes=df$X)
genes2<-lapply(genes2,function(list){if(class(list)=="logical")list<-NA else(list)})
genes2_names<-paste("nasal","pos",source.cell,sep="_")
names(genes2)<-genes2_names

x<-lapply(deg_list3,read.csv)
genes3 <- lapply(x, function(df)genes=df$X)
genes3<-lapply(genes3,function(list){if(class(list)=="logical")list<-NA else(list)})

genes3_names<-paste("bronchial","all",source.cell,sep="_")
names(genes3)<-genes3_names

x<-lapply(deg_list4,read.csv)
genes4 <- lapply(x, function(df)genes=df$X)
genes4<-lapply(genes4,function(list){if(class(list)=="logical")list<-NA else(list)})

genes4_names<-paste("bronchial","pos",source.cell,sep="_")
names(genes4)<-genes4_names
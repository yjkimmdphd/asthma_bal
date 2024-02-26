##
# demonstrate the distinct gene expression profile in nasal and bronchial 
## 
library(limma)
library(edgeR)
library(dplyr)
library(DESeq2)
library(pheatmap)
######################
## load readcount data
data.folder<-file.path("./reports/local_only/deg~bal-blood_cell(continuous)+batch/deg_gene_list")

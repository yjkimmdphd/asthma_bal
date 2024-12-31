################
# gsea data prep
################
library(tidyverse)

# -----------------------------------------------------------------------------
# 1) Load phenotype data
# -----------------------------------------------------------------------------
phen_path <- file.path(
  "./resources/processed_data",
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv"
)
phen <- read.csv(phen_path)

phen<-phen %>% filter(BAL_eos_p >= 0 & BAL_neut_p >= 0) %>%
  mutate(comp1 = factor(case_when(BAL_eos_p > 1 & BAL_neut_p > 4 ~ "mixed",
                                  BAL_eos_p > 1 & BAL_neut_p <= 4 ~ "eos",
                                  BAL_eos_p <= 1 & BAL_neut_p > 4 ~ "neut",
                                  BAL_eos_p <= 1 & BAL_neut_p <= 4 ~ "pauci"), levels = c("pauci", "neut", "mixed", "eos")),
         comp2 = factor(case_when(BAL_eos_p > 1 ~ "high_eos",
                                  BAL_eos_p <= 1 ~ "low_eos"), levels = c("high_eos", "low_eos")))
phen$Batch<-factor(phen$Batch,levels=unique(phen$Batch))

phen_bronch<-phen[grepl("^B",phen$SampleID),]
phen_nasal<-phen[grepl("^N",phen$SampleID),]
phen_nasal<-phen_nasal[-grep("F", phen_nasal$SampleID),]

# -----------------------------------------------------------------------------
# 2) Normalized nasal count table with vsd
# -----------------------------------------------------------------------------

# load count table
countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
genes<-counts[,"SampleID"]

phen_input<-phen_nasal

# select bronchial samples 
sample_id<-phen_input$SampleID
counts_selected<-counts[,sample_id]
rownames(counts_selected)<-genes
head(counts_selected)

# Assuming your DESeq2 object is called 'dds'
library(DESeq2)
library(vsn)
countdata<-counts_selected
coldata_cols<-c("comp2","Batch")
coldata<-phen_input[,coldata_cols]
rownames(coldata)<-sample_id

dds<-DESeqDataSetFromMatrix(countData = countdata,colData=coldata, design= ~ comp2+Batch)

# prefilter low count genes
smallestGroupSize <- min(apply(table(coldata),1,sum)) 
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

vsd <- vst(dds, blind=FALSE) 

meanSdPlot(assay(vsd))

## vsd with vs without batch effect removal 

### vsd without batch effect removal 
normalized_counts <- assay(vsd)  # This is now your transformed expression matrix
sampleDists <- dist(t(assay(vsd)))

library(pheatmap)
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Batch
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(vsd, intgroup=c("Batch"))
plotPCA(vsd, intgroup=c("comp2"))

### vsd with batch effect removed 
library(limma)
mat <- assay(vsd)
mm <- model.matrix(~comp2, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
assay(vsd) <- mat
meanSdPlot(assay(vsd))

normalized_counts <- assay(vsd)  # This is now your transformed expression matrix
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Batch
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(vsd, intgroup=c("Batch"))
plotPCA(vsd, intgroup=c("comp2"))


# -----------------------------------------------------------------------------
# 3) export normalized bronchial gene count table
# -----------------------------------------------------------------------------

# write.table(
#   assay(vsd),
#   "./resources/processed_data/gsea/gsea_normalized_gene_count_nasal_vsd_batch-corrected.txt",
#   sep = "\t",
#   row.names = TRUE, 
#   col.names = NA
# )
# 


### will only export vsd with batch effect removal
# -----------------------------------------------------------------------------
# 4) Extract class labels from 'coldata'
#    Ensure 'comp1' is a factor or character. If a factor, convert to character.
# -----------------------------------------------------------------------------
labels <- as.character(coldata$comp2)


# -----------------------------------------------------------------------------
# 5) Create the lines for the .cls file
# -----------------------------------------------------------------------------
num_samples <- length(labels)                     # e.g., 70
unique_classes <- unique(labels)                  # e.g., c("high_eos", "low_eos")
num_classes <- length(unique_classes)             # e.g., 4

# Line 1: "<num_samples> <num_classes> 1"
line1 <- paste(num_samples, num_classes, "1")

# Line 2: "# <class1> <class2> ... <classN>"
line2 <- paste("#", paste(unique_classes, collapse = " "))

# Line 3: actual class labels (for each sample) in order
line3 <- paste(labels, collapse = " ")

# -----------------------------------------------------------------------------
# 6) Write the .cls file
# -----------------------------------------------------------------------------
cls_content <- c(line1, line2, line3)
output_cls_file <- "./resources/processed_data/gsea/gsea_nasal_eos-p_mt1.cls"
writeLines(cls_content, con = output_cls_file)


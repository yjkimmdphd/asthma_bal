# --------------------------------------------------------------------
# generation of scaled gene counts from RNAseq experiment of Batch 1-6
# --------------------------------------------------------------------

library(tidyverse)

# -----------------------------------------------------------------------------
# 1) Load phenotype data
# -----------------------------------------------------------------------------
phen_path <- file.path(
  "./resources/processed_data",
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2025-02-14.csv"
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
# 2) Normalized bronchial count table with vsd
# -----------------------------------------------------------------------------

# load count table
countdata<-file.path("./resources/raw_data/MS_asthma/MS_asthma.batch12346.GRCh38.geneID_readcount.all_samples.QCed_final.txt")
counts<-if(file.exists(countdata)){read.delim(countdata, check.names = FALSE)}
genes<-counts[,"SampleID"]

# count normalization without batch correction
normalizeCounts <- function(input_df, var1="comp2", var2="Batch", formula) {
  
  phen_input <- input_df
  
  # select bronchial samples 
  sample_id <- phen_input$SampleID
  counts_selected <- counts[,sample_id]
  rownames(counts_selected) <- genes
  head(counts_selected)
  
  # Assuming your DESeq2 object is called 'dds'
  library(DESeq2)
  library(vsn)
  countdata <- counts_selected
  coldata_cols <- c(var1, var2)
  coldata <- phen_input[,coldata_cols]
  rownames(coldata) <- sample_id
  
  dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata, design= as.formula(formula))
  
  # prefilter low count genes
  smallestGroupSize <- min(apply(table(coldata),1,sum)) 
  keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
  dds <- dds[keep,]
  
  vsd <- vst(dds, blind=FALSE) 
  
  # Plot the meanSD plot as a side effect
  meanSdPlot(assay(vsd))
  
  ## vsd with vs without batch effect removal 
  
  ### vsd without batch effect removal 
  normalized_counts <- assay(vsd)  # This is now your transformed expression matrix
  
  sampleDists <- dist(t(assay(vsd)))
  
  library(pheatmap)
  library("RColorBrewer")
  library(gridExtra)
  
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- vsd$Batch
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  # For pheatmap to display in RStudio's plot pane, use grid.arrange or direct plotting
  
  dev.new()
  p_heatmap <- pheatmap(sampleDistMatrix,
                        clustering_distance_rows=sampleDists,
                        clustering_distance_cols=sampleDists,
                        col=colors)
  
  # Display PCA plots
  # These should display fine as they use ggplot2 under the hood
  pca_plot1 <- plotPCA(vsd, intgroup=c(var2))
  print(pca_plot1)
  
  pca_plot2 <- plotPCA(vsd, intgroup=c(var1))
  print(pca_plot2)
  
  # Return normalized counts
  return(normalized_counts)
}

### vsd with batch effect removed, count normalization without batch correction
normalizeCounts_batch_removed <- function(input_df, var1="comp2", var2="Batch", formula) {
  
  phen_input <- input_df
  
  # select bronchial samples 
  sample_id <- phen_input$SampleID
  counts_selected <- counts[,sample_id]
  rownames(counts_selected) <- genes
  head(counts_selected)
  
  # Assuming your DESeq2 object is called 'dds'
  library(DESeq2)
  library(vsn)
  countdata <- counts_selected
  coldata_cols <- c(var1, var2)
  coldata <- phen_input[,coldata_cols]
  rownames(coldata) <- sample_id
  
  dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata, design= as.formula(formula))
  
  # prefilter low count genes
  smallestGroupSize <- min(apply(table(coldata),1,sum)) 
  keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
  dds <- dds[keep,]
  
  vsd <- vst(dds, blind=FALSE) 
  
  library(limma)
  mat <- assay(vsd)
  mm <- model.matrix(~comp2, colData(vsd))
  mat <- limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
  assay(vsd) <- mat
  meanSdPlot(assay(vsd))

  ## vsd with vs without batch effect removal 
  
  ### vsd without batch effect removal 
  normalized_counts <- assay(vsd)  # This is now your transformed expression matrix
  
  sampleDists <- dist(t(assay(vsd)))
  
  library(pheatmap)
  library("RColorBrewer")
  library(gridExtra)
  
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- vsd$Batch
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  # For pheatmap to display in RStudio's plot pane, use grid.arrange or direct plotting
  
  dev.new()
  p_heatmap <- pheatmap(sampleDistMatrix,
                        clustering_distance_rows=sampleDists,
                        clustering_distance_cols=sampleDists,
                        col=colors)
  
  # Display PCA plots
  # These should display fine as they use ggplot2 under the hood
  pca_plot1 <- plotPCA(vsd, intgroup=c(var2))
  print(pca_plot1)
  
  pca_plot2 <- plotPCA(vsd, intgroup=c(var1))
  print(pca_plot2)
  
  # Return normalized counts
  return(normalized_counts)
}

normalized_bronch_count<-normalizeCounts(phen_bronch,"comp2","Batch", "~ comp2 + Batch")
normalized_nasal_count<-normalizeCounts(phen_nasal,"comp2","Batch", "~ comp2 + Batch")

normalized_bronch_count_batch_removed<-normalizeCounts_batch_removed(phen_bronch,"comp2","Batch", "~ comp2 + Batch")
normalized_nasal_count_batch_removed<-normalizeCounts_batch_removed(phen_nasal,"comp2","Batch", "~ comp2 + Batch")

### will only export vsd with batch effect removal

# -----------------------------------------------------------------------------
# 3) export normalized bronchial gene count table
# -----------------------------------------------------------------------------

write.table(
  normalized_bronch_count,
  "./resources/processed_data/normalized_gene_count/normalized_gene_count_bronch_vsd_no-batch-corrected.txt",
  sep = "\t",
  row.names = TRUE, 
  col.names = NA
)


write.table(
  normalized_bronch_count_batch_removed,
  "./resources/processed_data/normalized_gene_count/normalized_gene_count_bronch_vsd_batch-corrected.txt",
  sep = "\t",
  row.names = TRUE, 
  col.names = NA
)



# -----------------------------------------------------------------------------
# 4) export normalized nasal gene count table
# -----------------------------------------------------------------------------


write.table(
  normalized_nasal_count,
  "./resources/processed_data/normalized_gene_count/normalized_gene_count_nasal_vsd_no-batch-corrected.txt",
  sep = "\t",
  row.names = TRUE, 
  col.names = NA
)


write.table(
  normalized_nasal_count_batch_removed,
  "./resources/processed_data/normalized_gene_count/normalized_gene_count_nasal_vsd_batch-corrected.txt",
  sep = "\t",
  row.names = TRUE, 
  col.names = NA
)


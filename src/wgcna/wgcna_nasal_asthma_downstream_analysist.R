library(tidyverse)
library(WGCNA)

########################
## Load Readcount Data ##
########################

# Load cell count table
normalized_count_table_path <- "./resources/processed_data/normalized_gene_count/normalized_gene_count_nasal_vsd_batch-corrected.txt"
if (file.exists(normalized_count_table_path)) {
  counts <- read.table(normalized_count_table_path, 
                       header = TRUE, 
                       row.names = 1, 
                       sep = "\t")
}
genes <- rownames(counts)

# Define the path for the output folder
output_folder <- file.path("./reports/local_only/wgcna/nasal/output")

# Check if the folder exists; if not, create it
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE) # recursive = TRUE ensures all intermediate directories are created
}

# Now the folder is guaranteed to exist, and the path is set
output_folder

###
# 1. Create a new format for expression data
###
expression <- as.data.frame(t(counts))
colnames(expression) <- genes

# Check if genes/samples are good
gsg <- goodSamplesGenes(expression)

###
# 2. Remove outlier genes (and samples if necessary)
###

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(
      paste("Removing genes:", 
            paste(names(expression)[!gsg$goodGenes], collapse = ", "))
    )
  }
  
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(
      paste("Removing samples:", 
            paste(rownames(expression)[!gsg$goodSamples], collapse = ", "))
    )
  }
  
  expression <- expression[gsg$goodSamples, gsg$goodGenes]
}
###
# 3. Identify outlier samples using a dendrogram
###
sampleTree <- hclust(dist(expression), method = "average")

# png(file.path(output_folder,"sampleTree_nasal.png"), width = 800, height = 600)
# par(cex = 0.6)
# par(mar = c(0, 4, 2, 0))
# plot(
#   sampleTree,
#   main = "Sample clustering to detect outliers",
#   sub = "",
#   xlab = "",
#   cex.lab = 1.5,
#   cex.axis = 1.5,
#   cex.main = 2
# )
# abline(h = 250, col = "red")
# dev.off()

# Remove outliers
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 250, minSize = 7)
expression.data <- expression[cut.sampleTree == 1, ]

###
# 4. load merged colors
### 

mergedColors<-if(file.exists(file.path(output_folder,"mergedColors.txt"))){read.delim(file.path(output_folder,"mergedColors.txt"), sep="\t",header=TRUE,row.names = 1)}

mergedMEs<-if(file.exists(file.path(output_folder,"mergedMEs.txt"))){read.delim(file.path(output_folder,"mergedMEs.txt"), sep="\t",header=TRUE, row.names=1)}


# -----------------------------------------------------------------------------
# Load phenotype data
# -----------------------------------------------------------------------------
phen_path <- file.path(
  "./resources/processed_data",
  "scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_2024-12-26.csv"
)
phen <- read.csv(phen_path)

# load data for sampling date differences 
sampling_date_diff<-"./resources/processed_data/sampling_dates/swab-bal-cbc_differences_in_days.txt"
sampling_date_diff<-if(file.exists(sampling_date_diff)){read.table(sampling_date_diff,row.names = NULL,header = TRUE)}else{print("sampling date data doesn't exist")}
sampling_date_diff<-sampling_date_diff%>%filter(Comparison=="blood_bal")
colnames(sampling_date_diff)[1:3]<-c("ID","sampling_date_comp","sampling_date_diff_days")

# left join sampling date data and phenotype data 
phen<-left_join(phen,sampling_date_diff,by="ID")

phen$Batch<-factor(phen$Batch,levels=unique(phen$Batch))

phen_nasal<-phen[grepl("^N",phen$SampleID),]
phen_nasal<-phen_nasal[-grep("F", phen_nasal$SampleID),]
phen_input<-phen_nasal
phen_input$SampleID <- gsub("-", ".", phen_input$SampleID)

# Create continuous and categorical variables based on various thresholds
## Log-transformed/scaled/centered
source.cell.log <- c(
  "BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log",
  "BAL_wbc_log",    "blood_eos_log", "blood_eos_p_log", "blood_neut_log",
  "blood_neut_p_log", "blood_wbc_log"
)

## Categorical variables to test (BAL)
var_dichot_bal <- c(
  "bal_AEC_more_0", "bal_AEC_more_1", "bal_AEC_more_3", "bal_AEC_more_5",
  "bal_Eos_p_more_0", "bal_Eos_p_more_1", "bal_Eos_p_more_3",
  "bal_ANC_more_0", "bal_ANC_more_5", "bal_ANC_more_13",
  "bal_neut_p_more_0", "bal_neut_p_more_2", "bal_neut_p_more_5"
)

## Categorical variables to test (Blood)
var_dichot_blood <- c(
  "bld_AEC_more_0",   "bld_AEC_more_100",
  "bld_AEC_more_300", "bld_AEC_more_500"
)

phen_input <- phen_input %>%
  mutate(
    bal_AEC_more_0   = as.numeric(BAL_eos_ct > 0),
    bal_AEC_more_1   = as.numeric(BAL_eos_ct > 1),
    bal_AEC_more_1.2 = as.numeric(BAL_eos_ct > 1.2),
    bal_AEC_more_3   = as.numeric(BAL_eos_ct > 3),
    bal_AEC_more_5   = as.numeric(BAL_eos_ct > 5),
    
    bal_Eos_p_more_0 = as.numeric(BAL_eos_p > 0),
    bal_Eos_p_more_1 = as.numeric(BAL_eos_p > 1),
    bal_Eos_p_more_3 = as.numeric(BAL_eos_p > 3),
    
    bal_ANC_more_0   = as.numeric(BAL_neut_ct > 0),
    bal_ANC_more_5   = as.numeric(BAL_neut_ct > 5),
    bal_ANC_more_13  = as.numeric(BAL_neut_ct > 13),
    
    bal_neut_p_more_0 = as.numeric(BAL_neut_p > 0),
    bal_neut_p_more_2 = as.numeric(BAL_neut_p > 2),
    bal_neut_p_more_5 = as.numeric(BAL_neut_p > 5),
    
    bld_AEC_more_0   = as.numeric(blood_eos > 0),
    bld_AEC_more_100 = as.numeric(blood_eos > 100),
    bld_AEC_more_300 = as.numeric(blood_eos > 300),
    bld_AEC_more_500 = as.numeric(blood_eos > 500)
  ) %>%
  mutate(across(c(var_dichot_bal,var_dichot_blood), ~ factor(., levels = c(0, 1))))


#===============================================================================
# Identify Cell-Count Variables
#===============================================================================

# Combine continuous (log counts) and categorical variables
var_to_test <- c(source.cell.log, var_dichot_bal, var_dichot_blood)

# identify blood cell count phenotype variables
var_to_test_bld<-var_to_test[c(grep("blood",var_to_test),grep("bld",var_to_test))]

# Identify samples with CBC data and sampling date difference < 1 year
cbc_sampleID <- phen_input %>%
  filter(!is.na(vars(var_to_test_bld)), abs(sampling_date_diff_days) < 365) %>%
  pull(SampleID)

# subset alltraits for variables pertaining to BAL phenotype
alltraits<-phen_input[phen_input$SampleID%in%rownames(mergedMEs),]
rownames(alltraits)<-alltraits$SampleID

alltraits<-alltraits[,c("asthma_phen_ACT.score","BAL_eos_ct_log", "BAL_eos_p_log", "BAL_neut_ct_log", "BAL_neut_p_log",
                        "BAL_wbc_log", var_dichot_bal)] # only select cell counts in the phenotype data
good<-!(sapply(alltraits,is.na)%>%rowSums()>0) # remove rows with at least one NA
datTraits<-alltraits[good,]

# subset alltraits for variables pertaining to BAL phenotype
alltraits<-phen_input[phen_input$SampleID%in%rownames(mergedMEs),]
alltraits<-alltraits[alltraits$SampleID%in%cbc_sampleID,]
rownames(alltraits)<-alltraits$SampleID

alltraits_bld<-alltraits[,c("asthma_phen_ACT.score",var_to_test_bld)] # only select cell counts in the phenotype data
good_bld<-!(sapply(alltraits,is.na)%>%rowSums()>0) # remove rows with at least one NA
datTraits_bld<-alltraits_bld[good_bld,]

###
# 5. Module-Trait associations
###

gene.module.table<-data.frame(genes=colnames(expression.data),modules=mergedColors)

# ==================
# quantify the association between the expression profile and ACT score and BAL cell count phenotypes
# ==================
# calculates the correlation of the trait with previously identified module eigengenes. 
# This pairwise correlation is known as the eigengene gene significance
expression.data_bal<-expression.data[rownames(expression.data)%in%rownames(datTraits),]
mergedMEs_bal<-mergedMEs[rownames(mergedMEs)%in%rownames(datTraits),]
nSamples_bal = nrow(expression.data_bal)
module.trait.correlation_bal = cor(mergedMEs_bal, datTraits_bal, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue_bal = corPvalueStudent(module.trait.correlation_bal, nSamples_bal) #calculate the p-value associated with the correlation

# Will display correlations and their p-values
textMatrix_bal = paste(signif(module.trait.correlation_bal, 2), "\n(",
                   signif(module.trait.Pvalue_bal, 1), ")", sep = "");
dim(textMatrix_bal) = dim(module.trait.correlation_bal)
par(mar = c(6, 6, 3, 1))
# Display the correlation values within a heatmap plot
png(file.path(output_folder,"module-trait-correlation_bal.png"), width = 1200, height = 1200)
labeledHeatmap(Matrix = module.trait.correlation_bal,
               xLabels = names(datTraits_bal),
               yLabels = names(mergedMEs_bal),
               ySymbols = names(mergedMEs_bal),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_bal,
               setStdMargins = FALSE,
               cex.text = 1.1,
               zlim = c(-1,1),
               main = paste("Module-trait correlation"))
dev.off()


# ==================
# quantify the association between the expression profile and ACT score and blood cell count profile phenotypes
# ==================
# calculates the correlation of the trait with previously identified module eigengenes. 
# This pairwise correlation is known as the eigengene gene significance

expression.data_bld<-expression.data[rownames(expression.data)%in%rownames(datTraits_bld),]
mergedMEs_bld<-mergedMEs[rownames(mergedMEs)%in%rownames(datTraits_bld),]
nSamples_bld = nrow(expression.data_bld)
module.trait.correlation_bld = cor(mergedMEs_bld, datTraits_bld) #p for pearson correlation coefficient 
module.trait.Pvalue_bld = corPvalueStudent(module.trait.correlation_bld, nSamples_bld) #calculate the p-value associated with the correlation

# Will display correlations and their p-values
textMatrix_bld = paste(signif(module.trait.correlation_bld, 2), "\n(",
                   signif(module.trait.Pvalue_bld, 1), ")", sep = "");
dim(textMatrix_bld) = dim(module.trait.correlation_bld)
par(mar = c(6, 6, 3, 1))
# Display the correlation values within a heatmap plot
png(file.path(output_folder,"module-trait-correlation_bld.png"), width = 1200, height = 1200)
labeledHeatmap(Matrix = module.trait.correlation_bld,
               xLabels = names(datTraits_bld),
               yLabels = names(mergedMEs_bld),
               ySymbols = names(mergedMEs_bld),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_bld,
               setStdMargins = FALSE,
               cex.text = 1.1,
               zlim = c(-1,1),
               main = paste("Module-trait correlation"))
dev.off()


# ----------------
# proportion of module gene that overlaps with nasal DEG for 
# blood Eos %, AEC, Neut %, ANC
# ----------------
deg_folder<-file.path("./reports/local_only/deg_bal_nasal~cell2025-01-21")
deg_file<-if(file.exists(deg_folder)){list.files(deg_folder)[grep(".csv",list.files(deg_folder))]}else{print("folder doesn't exist")}

# identify results of DEG associated with blood cell counts 
deg_file_bld<-deg_file[c(grep("blood",deg_file),grep("bld",deg_file))]
deg_file_bld<-deg_file_bld[c(grep("sig",deg_file_bld))]

deg<-lapply(deg_file_bld,function(x){read.csv(file.path(deg_folder,x), row.names = 1)})

deg_abs_lfc<-lapply(deg,function(x){deg_genes<-filter(x,abs(log2FoldChange)>1); return(rownames(deg_genes))})

module_output_folder <- file.path("./reports/local_only/wgcna/nasal/output/module-gene_list")

wgcna_folder<-module_output_folder

module_list<-list.files(wgcna_folder)[grep(".txt",list.files(wgcna_folder))]

modules_files<-file.path(wgcna_folder,module_list)

m_gene_list<-lapply(modules_files,read.table)

names(m_gene_list)<-sub("batch12346.txt","",module_list)

overlap_proportion <- sapply(m_gene_list, function(x) {
  round(mean(unlist(x) %in% deg_abs_lfc),2)
})# no notable overlap between nasal down gene with WGCNA module genes 

overlap_sum <- sapply(m_gene_list, function(x) {
  sum(unlist(x) %in% deg_abs_lfc)
})

overlap_df<-data.frame(overlap_proportion,overlap_sum)
overlap_df<-overlap_df[order(-(overlap_df$overlap_sum)),]
print(overlap_df)

write.table(overlap_df,file.path(output_folder,"module-gene_deg-gene_overlap.txt"), sep="\t",quote=FALSE, row.names=TRUE, col.names=NA)

### 
# 9. Target gene identification
### 

# Whole Network connectivity - a measure for how well the node is connected throughout the entire system
# Intramodular connectivity - a measure for how well the node is connected within its assigned module. Also an indicator for how well that node belongs to its module. This is also known as module membership.

modNames = substring(names(mergedMEs), 3) #extract module names

#Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(WGCNA::cor(expression.data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# ----------------
# for log(BAL Eos %)
# ----------------

# Define variable weight containing the weight column of datTrait
BAL_eos_p_log = as.data.frame(datTraits$BAL_eos_p_log)
names(BAL_eos_p_log) = "BAL_eos_p_log"

phen_of_interest<-BAL_eos_p_log

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(WGCNA::cor(expression.data, phen_of_interest, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(phen_of_interest), sep="")
names(GSPvalue) = paste("p.GS.", names(phen_of_interest), sep="")
head(GSPvalue)

# scatter plot of gene significance vs. module membership in all the module

# sort the modules based on the cor coefficient

gene_sig_v_mod_memb<-data.frame(cor=rep(NA,length(modNames)),p_val=rep(NA,length(modNames)))
rownames(gene_sig_v_mod_memb)<-modNames

for (mod in modNames) {
  module<-mod
  column<-match(module,modNames)
  moduleGenes <- mergedColors == module
  print(paste("Module Membership in", module, "module","vs Gene significance for eos"))
  cor_test<-cor.test(abs(geneModuleMembership[moduleGenes,column]),abs(geneTraitSignificance[moduleGenes,1]))
  p_val<-WGCNA::corPvalueStudent(
    round(cor(abs(geneModuleMembership[moduleGenes,column]),abs(geneTraitSignificance[moduleGenes,1])),2),sum(moduleGenes)
    )
  gene_sig_v_mod_memb[module,"cor"]<-cor_test$estimate
  gene_sig_v_mod_memb[module,"p_val"]<-p_val
} # if cor used for corPvalueStudent fx is rounded to second decimal place, the associated p-value is essentially the same as when using verboseScatterplot to calculate r and p-val

gene_sig_v_mod_memb<-arrange(gene_sig_v_mod_memb,desc(cor))

# Initialize an empty list to store the plots
plot_list <- list()

# Loop through the module names
for (mod in modNames) {
  par(mar = c(4,4,4,4)) 
  module <- mod
  column <- match(module, modNames)
  moduleGenes <- mergedColors == module
  
  # Generate the plot and capture it
  verboseScatterplot(
    abs(geneModuleMembership[moduleGenes, column]),
    abs(geneTraitSignificance[moduleGenes, 1]),
    xlab = paste("Module Membership in", module, "module"),
    ylab = "Gene significance for eos",
    main = paste("Module membership vs. gene significance\n"),
    use = "complete.obs",
    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module
  )
  
  # Save the plot in the list using recordPlot
  plot_list[[module]] <- recordPlot()
}


# Now you can replay any plot from the list
# For example, replay the plot for a specific module
replayPlot(plot_list[["blue"]])  # Replace "blue" with an actual module name

# ----------------
# for BAL Eos % > 1 vs <=1
# ----------------

# Define variable weight containing the weight column of datTrait
phen_of_interest<-as.data.frame(datTraits$comp2_eos_p_mt1)

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(WGCNA::cor(expression.data, phen_of_interest, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(phen_of_interest), sep="")
names(GSPvalue) = paste("p.GS.", names(phen_of_interest), sep="")
head(GSPvalue)

# scatter plot of gene significance vs. module membership in all the module
# DEG are colored in red
par(mar = c(4,4,4,4),mfrow=c(6,4))

png(file.path(output_folder,"Module_membership_vs._gene_significance_eos-p_mt1.png"), width = 1200, height = 1200)
for(mod in modNames){
  module = mod
  column = match(module, modNames)
  moduleGenes = mergedColors==module
  
  module_color<-rep(mod,length(geneModuleMembership[moduleGenes,column]))
  module_color[which(rownames(geneModuleMembership)[moduleGenes]%in%deg_abs_lfc)]<-"red"
  
  verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                     abs(geneTraitSignificance[moduleGenes,1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for eos mt1",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module_color)
}
dev.off()
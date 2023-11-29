########################
# GO analysis with DEG post DESeq2
##################################

###################################################
### code chunk number 1: load_library
###################################################
library(goseq)


###################################################
### code chunk number 2: set_width
###################################################
options(width=84)

################
## Load DEG data
################
res<-read.csv(file.path(getwd(),"output/res.table_2023-09-24.csv"))
deg<-list(res25=read.csv(file.path(getwd(),"output/DEG_2023-09-19/DEG_res25_serum_E_2023-09-19.csv")),
          res26=read.csv(file.path(getwd(),"output/DEG_2023-09-19/DEG_res26_serum_E_2023-09-19.csv")))

genes.25=as.integer(p.adjust(deg$res25$padj<.05))
names(genes.25)=deg$res25[deg$res25$log2FoldChange!=0,"X"]
table(genes.25)

genes.26=as.integer(p.adjust(deg$res26$padj<.05))
names(genes.26)=deg$res25[deg$res26$log2FoldChange!=0,"X"]
table(genes.26)


###################################################
### code chunk number 9: head_organisms
###################################################
head(supportedOrganisms())


###################################################
### code chunk number 10: head_geneids
###################################################
supportedOrganisms()[supportedOrganisms()$Genome=="hg19",]


###################################################
### code chunk number 11: pwf
###################################################
pwf.25=nullp(genes.25,"hg38","knownGene")
head(pwf.25)


###################################################
### code chunk number 12: GO.wall
###################################################
GO.wall=goseq(pwf,"hg19","ensGene")
head(GO.wall)


###################################################
### code chunk number 13: GO.samp
###################################################
GO.samp=goseq(pwf,"hg19","ensGene",method="Sampling",repcnt=1000)


###################################################
### code chunk number 14: head_samp
###################################################
head(GO.samp)


###################################################
### code chunk number 15: plot_wal_v_samp
###################################################
plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]), 
     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)", 
     xlim=c(-3,0))
abline(0,1,col=3,lty=2)


###################################################
### code chunk number 16: GO.nobias
###################################################
GO.nobias=goseq(pwf,"hg19","ensGene",method="Hypergeometric")
head(GO.nobias)


###################################################
### code chunk number 17: plot_wal_v_hyper
###################################################
plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.wall[,1],GO.nobias[,1]),2]), 
     xlab="log10(Wallenius p-values)", ylab="log10(Hypergeometric p-values)", 
     xlim=c(-3,0), ylim=c(-3,0))
abline(0,1,col=3,lty=2)


###################################################
### code chunk number 18: GO.limited
###################################################
GO.MF=goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF"))
head(GO.MF)


###################################################
### code chunk number 19: enriched_GO
###################################################
enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,
                                      method="BH")<.05]
head(enriched.GO)


###################################################
### code chunk number 20: GO_explained
###################################################
library(GO.db)
for(go in enriched.GO[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}


###################################################
### code chunk number 21: getlength
###################################################
len=getlength(names(genes),"hg19","ensGene")
length(len)
length(genes)
head(len)


###################################################
### code chunk number 22: getgo
###################################################
go=getgo(names(genes),"hg19","ensGene")
length(go)
length(genes)
head(go)


###################################################
### code chunk number 23: conv_table
###################################################
goseq:::.ID_MAP
goseq:::.ORG_PACKAGES


###################################################
### code chunk number 24: norm_analysis (eval = FALSE)
###################################################
## pwf=nullp(genes,"hg19","ensGene")
## go=goseq(pwf,"hg19","ensGene")


###################################################
### code chunk number 25: verbose_analysis (eval = FALSE)
###################################################
## gene_lengths=getlength(names(genes),"hg19","ensGene")
## pwf=nullp(genes,bias.data=gene_lengths)
## go_map=getgo(names(genes),"hg19","ensGene")
## go=goseq(pwf,"hg19","ensGene",gene2cat=go_map)


###################################################
### code chunk number 26: KEGG_mappings (eval = FALSE)
###################################################
## # Get the mapping from ENSEMBL 2 Entrez
## en2eg=as.list(org.Hs.egENSEMBL2EG)
## # Get the mapping from Entrez 2 KEGG
## eg2kegg=as.list(org.Hs.egPATH)
## # Define a function which gets all unique KEGG IDs 
## # associated with a set of Entrez IDs
## grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
## # Apply this function to every entry in the mapping from 
## # ENSEMBL 2 Entrez to combine the two maps
## kegg=lapply(en2eg,grepKEGG,eg2kegg)
## head(kegg)


###################################################
### code chunk number 27: KEGG (eval = FALSE)
###################################################
## pwf=nullp(genes,"hg19","ensGene")
## KEGG=goseq(pwf,gene2cat=kegg)
## head(KEGG)


###################################################
### code chunk number 28: KEGG_goseq
###################################################
pwf=nullp(genes,'hg19','ensGene')
KEGG=goseq(pwf,'hg19','ensGene',test.cats="KEGG")
head(KEGG)


###################################################
### code chunk number 29: KEGG_from_db
###################################################
kegg=as.list(org.Hs.egPATH)
head(kegg)


###################################################
### code chunk number 30: countbias
###################################################
countbias=rowSums(counts)[rowSums(counts)!=0]
length(countbias)
length(genes)


###################################################
### code chunk number 31: GO.counts
###################################################
pwf.counts=nullp(genes,bias.data=countbias)
GO.counts=goseq(pwf.counts,"hg19","ensGene")
head(GO.counts)


###################################################
### code chunk number 32: setup
###################################################
sessionInfo()

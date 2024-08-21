
###################################
# custom functions for DEG analysis
###################################

# filter_low_expressed_genes_method2: Filters low counts genes using TMM normalized lcpm as a cutoff point. Requires 'limma'
# rowgenes_counttable: changes the row names of the count table with gene names
# run_deseq2_DEG_analysis: takes countdata,coldata,design,des as variables, and runs DEG on DESeq2
# get_DEG_results: saves result of DESeq2 output, ordered in padj 
# generate_DEG_input_summary_table: makes a table of input information
# generate_DEG_summary_table: makes results summary (i.e., # of DEG for each analysis)


## gene filter
filter_low_expressed_genes_method2<-function(readcounts, n_sample){
  x<-readcounts
  L<-mean(colSums(x))*1e-6 # mean library size
  M<-median(colSums(x))*1e-6 # median library size
  lcpm.cutoff <- log2(1/M + 1/L) # lcpm cutoff for filtering genes with very low counts
  
  ### normalize counts with TMM
  lib_size<-colSums(x)*1e-6
  norm.factor<-calcNormFactors(x, method = "TMM",lib.size = lib_size)
  sweep(x, 2, norm.factor, FUN = "/")
  
  ### calculate lcpm based on TMM normalized counts 
  lcpm.x<-cpm(x,log=TRUE)
  drop.genes<-which(rowSums(lcpm.x<lcpm.cutoff)>n_sample) # which genes have normalized lcpm less than the cutoff in >10% of the samples
  x<-x[-drop.genes,] # readcounts retained 
  print(paste("# of dropped genes:",length(drop.genes))) # number of dropped genes
  return(x)
}  # function to filter low counts genes using TMM normalized lcpm as a cutoff point


## count table-related
rowgenes_counttable<-function(c,cf){ 
  g<-rownames(c)
  gf<-rownames(cf)
  c<-c[g%in%gf,]
  return(c)
}# function to change rownames of the count table with gene names 

## DEsesq
run_deseq2_DEG_analysis<-function(countdata,coldata,design,des){ # function to run deseq using data from matrix 
  print(paste("design:",unique(phen$Type),"expression =",des, sep = " "))
  dds<-DESeqDataSetFromMatrix(countData = countdata,colData=coldata, design=as.formula(design))
  dds<-DESeq(dds)
  print("resultnames")
  print(resultsNames(dds))
  return(dds)
}

get_DEG_results<-function(dds,resultname){ # function to call results of the deseq analysis with 'results' function
  res<-results(dds, name=resultname) # use source.cell.log element for resultname.
  res <- res[order(res$padj),]
  return(res)
}

## summarize the data input 

generate_DEG_input_summary_table<-function(){
  filter_method<-"TMM normalized LCPM cutoff"
  n_filtered_genes<-paste("analyzed n_genes:", nrow(ct),",","filtered n_genes:",nrow(counts)-nrow(ct))
  samples<-sapply(dds, function(d){colData(d)$SampleID%>%paste(collapse = ",")})
  dds<-paste("dds",1:length(dds),sep="")
  results<-paste("res",1:length(res),sep="")
  design<-deg.design
  df<-data.frame(dds=dds,results=results,design=design,samples=samples,filter_method=filter_method,n_filtered_genes=n_filtered_genes)
  return(df)
}

##  summary table of the DEG analysis
generate_DEG_summary_table<-function(){
  reslist<-paste("res.sig",1:length(res.sig),sep="")
  n_sig_deg<-sapply(res.sig,nrow)
  design<-deg.design
  source_cell<-source.cell
  df<-data.frame(type=unique(phen$Type),results=reslist,n_sig_deg,design=design,source_cell=source_cell, row.names = NULL)
  return(df)
}

# BAL Asthma Bronchial Cell Analysis - DEG and GO Analysis
# Author: [Your Name]
# Date: [Date]
# Description: Differential gene expression analysis of bronchial RNA-seq data
#              Model: Bronch DEG ~ cell count continuous + Batch

# Load required libraries ----
library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(gridExtra)
library(grid)

# Set file paths and parameters ----
BASE_DIR <- "./reports/local_only/deg_bal_bronch~cell2025-01-03"
OUTPUT_DIR <- "./reports/figures/deg/volcano_plot/bronch~cell(cont_or_dich)"
GO_DIR <- file.path(BASE_DIR, "GO_output")

P_CUTOFF <- 0.05
FC_CUTOFF <- 1
TOP_GENES_DISPLAY <- 10

# Helper functions ----
extract_analysis_names <- function(file_names) {
  sapply(file_names, function(x) {
    string_after_tilde <- trimws(strsplit(x, "~")[[1]][2])
    string_without_csv <- sub("\\+ Batch_2025-01-03_.csv$", "", string_after_tilde)
    return(string_without_csv)
  })
}

create_color_mapping <- function(deg_results) {
  keyvals <- lapply(deg_results, function(res) {
    ifelse(res$padj < P_CUTOFF & res$log2FoldChange > FC_CUTOFF, 'cyan',
           ifelse(res$padj < P_CUTOFF & res$log2FoldChange < -FC_CUTOFF, 'magenta', 'grey'))
  })
  
  # Assign meaningful names to colors
  for (name in names(keyvals)) {
    names(keyvals[[name]]) <- sapply(keyvals[[name]], function(x) {
      switch(x,
             "magenta" = "down",
             "cyan" = "up",
             "grey" = "nonsig",
             NA)
    }, USE.NAMES = FALSE)
  }
  
  return(keyvals)
}

create_volcano_plot <- function(results, keyvals, title, custom_labs = NULL, 
                                xlim_custom = NULL, top_n = 10) {
  # Get top up/down regulated genes for labeling
  if (is.null(custom_labs)) {
    labs_up <- rownames(
      filter(results, padj < P_CUTOFF, log2FoldChange > 0) %>%
        arrange(desc(log2FoldChange))
    )[1:top_n]
    
    labs_down <- rownames(
      filter(results, padj < P_CUTOFF, log2FoldChange < 0) %>%
        arrange(desc(abs(log2FoldChange)))
    )[1:top_n]
    
    labs <- c(labs_up, labs_down)
  } else {
    labs <- custom_labs
  }
  
  # Set x-axis limits
  if (is.null(xlim_custom)) {
    xlim_range <- c(min(results$log2FoldChange) - 0.5, 
                    max(results$log2FoldChange) + 0.5)
  } else {
    xlim_range <- xlim_custom
  }
  
  # Create volcano plot
  EnhancedVolcano(results,
                  lab = rownames(results),
                  selectLab = labs,
                  title = title,
                  x = 'log2FoldChange',
                  y = 'padj',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  xlim = xlim_range,
                  ylab = bquote(~-Log[10]~ 'FDR'),
                  ylim = c(0, -log(min(results$padj), 10)),
                  pCutoff = P_CUTOFF,
                  FCcutoff = FC_CUTOFF,
                  cutoffLineType = 'twodash',
                  cutoffLineWidth = 0.8,
                  pointSize = 4.0,
                  labSize = 3,
                  colAlpha = 0.4,
                  colCustom = keyvals,
                  legendPosition = 'right',
                  legendLabSize = 10,
                  legendIconSize = 5.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.75)
}

create_go_barplot <- function(go_data, title, top_n = 10) {
  # Limit to top N terms
  n_terms <- min(nrow(go_data), top_n)
  plot_data <- go_data[1:n_terms, ]
  
  # Wrap labels for better display
  wrapped_labels <- str_wrap(plot_data$Description, width = 30)
  
  ggplot(plot_data, aes(y = Fold.Enrichment, x = Description, fill = -log10(FDR))) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "blue", high = "yellow") +
    labs(y = "Fold Enrichment", 
         x = "Gene Ontology Term", 
         fill = "-log10(FDR)", 
         title = title) +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 9.5),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 10),
          plot.title = element_text(size = 12)) +
    scale_x_discrete(labels = wrapped_labels) +
    coord_flip()
}

# Load and process DEG results ----
print("Loading DEG results...")

# Get DEG result files
file_names <- list.files(BASE_DIR)
deg_files <- file_names[grep("deg_bronch_res_all", file_names)]

# Load DEG results
deg_results <- lapply(file.path(BASE_DIR, deg_files), function(d) {
  read.csv(d, row.names = 1)
})

# Extract and assign meaningful names
extracted_names <- extract_analysis_names(deg_files)
names(deg_results) <- extracted_names

print("Loaded DEG results for:")
print(names(deg_results))

# Add significance classification ----
print("Adding significance classifications...")

deg_results <- lapply(deg_results, function(d) {
  mutate(d, deg_sig = ifelse(
    padj < P_CUTOFF & log2FoldChange > FC_CUTOFF, 'cyan',
    ifelse(padj < P_CUTOFF & log2FoldChange < -FC_CUTOFF, 'magenta', 'grey')
  ))
})

# Create color mapping for plots
keyvals <- create_color_mapping(deg_results)

# Check for NA values in color mapping
na_counts <- lapply(lapply(keyvals, is.na), sum)
print("NA counts in color mapping:")
print(na_counts)

# Generate volcano plots ----
print("Generating volcano plots...")

plot_list <- list()
for (i in seq_along(deg_results)) {
  title <- names(deg_results)[i]
  plot_list[[i]] <- create_volcano_plot(
    results = deg_results[[i]],
    keyvals = keyvals[[i]],
    title = title
  )
}

# Save volcano plots
print("Saving volcano plots...")
for (i in 1:length(plot_list)) {
  file_name <- paste("volcano", names(deg_results)[i], sep = "_")
  ggsave(filename = file.path(OUTPUT_DIR, paste0(file_name, ".png")), 
         plot = plot_list[[i]], 
         width = 8, height = 6, dpi = 300)
}

# Create custom BAL Eos % > 1% volcano plot ----
print("Creating custom BAL Eos % > 1% volcano plot...")

bal_eos_key <- "bal_Eos_p_more_1 "
if (bal_eos_key %in% names(deg_results)) {
  custom_genes <- c("POSTN", "SERPINB2", "CLCA1", "CST1", "CST2")
  
  p_bal_eos_custom <- create_volcano_plot(
    results = deg_results[[bal_eos_key]],
    keyvals = keyvals[[bal_eos_key]],
    title = names(deg_results)[8],
    custom_labs = custom_genes,
    xlim_custom = c(-4, max(deg_results[[bal_eos_key]]$log2FoldChange) + 0.5)
  )
  
  # Customize x-axis and print
  p_bal_eos_final <- p_bal_eos_custom + 
    scale_x_continuous(breaks = seq(-4, 5, 1)) + 
    coord_cartesian(xlim = c(-4, 5))
  
  print(p_bal_eos_final)
}

# GO Analysis ----
print("Processing GO analysis results...")

if (!file.exists(GO_DIR)) {
  stop("GO output directory does not exist: ", GO_DIR)
}

# Load GO results
go_files <- list.files(GO_DIR)
txt_files <- go_files[grep("\\.txt$", go_files)]
filtered_files <- txt_files[grepl("filtered", txt_files)]

print("Found GO files:")
print(filtered_files)

# Check file existence
file_paths <- file.path(GO_DIR, filtered_files)
print("File existence check:")
print(file.exists(file_paths))

# Load GO term data
go_results <- lapply(file_paths, function(file) {
  read.table(file, sep = "\t", header = TRUE)
})

# Process GO results
go_results <- lapply(go_results, function(d) {
  d[, colnames(d) %in% c("Category", "GO_ID", "Description", "Genes", 
                         "Fold.Enrichment", "FDR")] %>%
    filter(FDR < P_CUTOFF) %>%
    arrange(desc(Fold.Enrichment))
})

# Assign names to GO results
names(go_results) <- gsub("\\+batch12346", "", 
                          substr(filtered_files, 1, nchar(filtered_files) - 4))

# Create factor levels for proper ordering
for (i in 1:length(go_results)) {
  go_results[[i]]$Description <- factor(go_results[[i]]$Description, 
                                        levels = go_results[[i]]$Description)
}

# Generate GO plots ----
print("Generating GO term plots...")

go_plot_list <- list()
for (i in 1:length(go_results)) {
  if (nrow(go_results[[i]]) > 0) {
    go_plot_list[[i]] <- create_go_barplot(
      go_data = go_results[[i]],
      title = names(go_results)[i]
    )
  }
}

# Display first plot as example
if (length(go_plot_list) > 0) {
  print(go_plot_list[[1]])
}

# Create arranged plots
if (length(go_plot_list) >= 2) {
  print("Arranging GO plots...")
  grid.arrange(grobs = go_plot_list[1:2], ncol = 1)
}

# Special analysis for BAL Eos % > 1% GO terms ----
print("Processing BAL Eos % > 1% GO analysis...")

bal_eos_pattern <- "GO_deg_bronch_res_sig_16_~ bal_Eos_p_more_1"
bal_eos_go <- go_results[grep(bal_eos_pattern, names(go_results))]

if (length(bal_eos_go) > 0) {
  # Select top 5 terms for each direction
  bal_eos_go_top <- lapply(bal_eos_go, function(df) {
    n_terms <- min(nrow(df), 5)
    return(df[1:n_terms, ])
  })
  
  # Process positive and negative regulation separately
  if (length(bal_eos_go_top) >= 2) {
    # Flip sign for negatively regulated genes
    bal_eos_go_top[[1]] <- bal_eos_go_top[[1]] %>% 
      mutate(Fold.Enrichment = -Fold.Enrichment)
    
    # Combine data
    combined_go_data <- rbind(bal_eos_go_top[[1]], bal_eos_go_top[[2]]) %>%
      arrange(Fold.Enrichment)
    
    # Check for duplicates
    duplicated_terms <- which(duplicated(combined_go_data$Description) | 
                                duplicated(combined_go_data$Description, fromLast = TRUE))
    
    if (length(duplicated_terms) == 0) {
      print("No duplicate GO terms found")
    } else {
      print("Duplicate GO terms:")
      print(combined_go_data[duplicated_terms, ])
    }
    
    # Create factor levels and wrapped labels
    combined_go_data$Description <- factor(combined_go_data$Description, 
                                           levels = unique(combined_go_data$Description))
    wrapped_labels_bal <- str_to_title(str_wrap(levels(combined_go_data$Description), 
                                                width = 30))
    
    # Create final BAL Eos plot
    p_bal_eos_go <- ggplot(combined_go_data, 
                           aes(y = Fold.Enrichment, x = Description, fill = -log10(FDR))) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "blue", high = "yellow") +
      labs(y = "Fold Enrichment", 
           x = "Gene Ontology Term", 
           fill = "-log10(FDR)", 
           title = "bal_Eos_p_more_1") +
      theme(axis.text.x = element_text(size = 15, angle = -45, hjust = 0, vjust = 1),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 20),
            legend.title = element_text(size = 10),
            plot.title = element_text(size = 12)) +
      scale_x_discrete(labels = wrapped_labels_bal) +
      scale_y_continuous(position = "left") +
      coord_flip()
    
    print(p_bal_eos_go)
  }
}

print("Analysis complete!")
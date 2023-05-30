# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("AnnotationDbi")

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)





setwd("~/My Drive/Bench/Opossum/scRNAseq/plots/Pyramid_Subset/Species_Differences_by_Celltype/GO_by_celltype/")

# Load the table with differentially expressed genes from multiple cell types
# Assuming 'sc_data' is a data.frame with gene IDs in the first column, log2 fold changes in the second column, and cell types in the third column
sc_data.all <- read.table("~/My Drive/Bench/Opossum/scRNAseq/plots/Pyramid_Subset/Species_Differences_by_Celltype/PN_subset_Species_specific_differences_per_cluster.csv", 
                      sep = ",", header = TRUE)

# because we want this with reference to what is higher in mouse, 
# original statistical test was opossum vs mouse, so the log2FC value need to be inverted
sc_data.all$avg_log2FC <- sc_data.all$avg_log2FC*-1

# Remove those that are 100% mouse specific (ie MALAT1 etc)
sc_data <- subset(sc_data.all, pct.1 > 0)

# Function for Capitalization
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

#Fix Gene names so they match with keytype SYMBOL
sc_data$gene <- capwords(tolower(sc_data$gene))



# Function to perform GO analysis and GSEA plots for a specific celltype
perform_go_analysis <- function(data, celltype, ontology) {
  # Filter data for the current celltype
  cell_data <- data[data$celltype == celltype, ]
  
  cat("Celltype:", celltype, "\n")
  cat("Number of unique genes:", length(unique(cell_data$gene)), "\n")
  cat("Total number of genes:", nrow(cell_data), "\n")
  cat("Testing ontology ", ontology, "\n")
  
  # Convert gene symbols to Entrez IDs
  #entrez_ids <- mapIds(org.Mm.eg.db, keys = cell_data$gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  log2_fc <- cell_data$avg_log2FC 
  
  # Create a named vector of log2 fold changes with Entrez IDs as names
  names(log2_fc) <- cell_data$gene
  
  # Sort the log2_fc vector in decreasing order
  log2_fc <- sort(log2_fc, decreasing = TRUE)
  

  
  #remove any NA gene 
  log2_fc <- na.omit(log2_fc)
  
  # Perform GSEA analysis
  gsea_results <- gseGO(geneList = log2_fc,
                        OrgDb = org.Mm.eg.db,
                        keyType = "SYMBOL",
                        ont = ontology,
                        minGSSize = 20,
                        maxGSSize = 1000,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        verbose = FALSE)
  
  # Check if gsea_results has any results
  if (length(gsea_results) == 0 || is.null(gsea_results)) {
    cat("No significant enrichment found for celltype:", celltype, "\n")
    return(NULL)
  }
  
  # Visualize GSEA plots for the top 10 categories
  top_n <- 10
  top_sets <- head(as.data.frame(gsea_results)$Description, top_n)
  
  for (i in 1:top_n) {
    gseaplot(gsea_results, geneSetID = i, title = paste(celltype, "-", top_sets[i]), legendPos = "bottomright")
  }
  # Return the gsea_results object
  return(gsea_results)
}

#test <- perform_go_analysis(data = sc_data, celltype = "CA3")

# Get unique cell types
celltypes <- unique(sc_data$celltype)

# Initialize an empty list to store the gsea_results for each celltype
gsea_results_list <- list()


# Perform GO analysis and GSEA plots for each celltype
for (celltype in celltypes) {
  # Call perform_go_analysis with the desired ontology (e.g., "BP" for Biological Process or "MF" for Molecular Function)
  ontology <- "BP"
  gsea_results <- perform_go_analysis(sc_data, celltype, ontology)
  
  if (!is.null(gsea_results)) {
    gsea_results_list[[celltype]] <- gsea_results
    
    # Create a folder named 'gsea' if it doesn't exist
    if (!dir.exists("gsea")) {
      dir.create("gsea")
    }
    
    # Replace slashes with hyphens for creating directory later
    clean_celltype <- gsub("/", "-", celltype)
    
    # Create a subfolder named after the celltype if it doesn't exist
    celltype_folder <- file.path("gsea", clean_celltype)
    if (!dir.exists(celltype_folder)) {
      dir.create(celltype_folder)
    }
    
    # Generate and save GSEA plots for the top 10 categories
    top_n <- 20
    top_sets <- head(as.data.frame(gsea_results)$Description, top_n)
    
    for (i in 1:top_n) {
      pdf_filename <- file.path(celltype_folder, paste(clean_celltype, "- Rank", i, "-", ontology, "-", top_sets[i], ".pdf", sep = ""))
      pdf(pdf_filename)
      plot <- gseaplot(gsea_results, geneSetID = i, title = paste(celltype, "- Rank", i, "-", ontology, "-", top_sets[i]), legendPos = "bottomright")
      print(plot)
      dev.off()
    }
  }
}


# Dotplot of gene ontologies accross celltypes

library(ggplot2)

# Function to extract the top 20 gene sets for each cell type
extract_top_gene_sets <- function(gsea_results, n = 20) {
  gene_sets <- as.data.frame(gsea_results)[1:n, ]
  gene_sets$Description <- as.character(gene_sets$Description)
  return(gene_sets)
}

# Extract the top 20 gene sets for each cell type
top_gene_sets <- lapply(gsea_results_list, extract_top_gene_sets, n = 50)

# Combine the gene sets and add a cell type column
top_gene_sets_df <- do.call(rbind, lapply(names(top_gene_sets), function(x) {
  data.frame(CellType = x, top_gene_sets[[x]])
}))

# Replace slashes with hyphens in the cell type names
top_gene_sets_df$CellType <- gsub("/", "-", top_gene_sets_df$CellType)

# Calculate the median NES for each gene set
median_NES <- aggregate(NES ~ Description, data = top_gene_sets_df, FUN = median)

# Reorder the levels of the Description factor based on the median NES
top_gene_sets_df$Description <- factor(top_gene_sets_df$Description, levels = unique(median_NES$Description[order(median_NES$NES, decreasing = FALSE)]))


library(forcats)
ggplot(top_gene_sets_df, aes(x =  fct_rev(Description), y = CellType, size = abs(-log10(p.adjust)), color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top 20 gene sets by cell type",
       x = "Cell Type",
       y = "Gene Set",
       size = "log10(p-value)",
       color = "NES") +
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))


# UPSET PLOT
library(UpSetR)

# 

# Create an empty list for storing named vectors
list_input <- list()

# Iterate over the gsea_results_list and create named vectors for each cell type
for (celltype in names(gsea_results_list)) {
  gsea_results <- gsea_results_list[[celltype]]
  categories <- as.data.frame(gsea_results)$Description
  celltype <- gsub("/", "-", celltype) # Replace slashes with hyphens in the celltype names
  list_input[[celltype]] <- categories
}

# Print the list_input to check the structure
print(list_input)


# Create an UpSet plot using the list_input
upset_plot <- upset(fromList(list_input),empty.intersections = "on", order.by = "freq", 
                    nsets = 14
                    )
upset_plot


# Heatmap of differentially expressed genes by celltype

top_genes_list <- subset(sc_data.all, pct.2 > 0) # we will only plot mouse values, so remove Opossum genes where mouse expresison is 0
top_genes_list <- unique(top_genes_list$gene)

library(Seurat)
sc_seurat.sct <- readRDS("~/My Drive/Bench/Opossum/scRNAseq/sc_seurat.PN.sct.rds")

# To make a heatmap,
# lets subset the seurat object to mouse cells and genes that show species specific changes
sc_seurat.sct.mouse <- subset(sc_seurat.sct, Species == "Mouse")
sc_seurat.sct.mouse <- subset(sc_seurat.sct, Stage %in% c("Mouse E14","Mouse E16"))

# Scale the data using the genes that were previously determined to show species specific differences by celltype
#  Must re-scale, using all the rows, not just previously selected Variable features
sc_seurat.sct.mouse <- ScaleData(sc_seurat.sct.mouse, features=rownames(sc_seurat.sct.mouse)) 

library(dplyr)
top100 <- subset(sc_data.all, pct.2 > 0)  %>% 
  group_by(celltype) %>% 
  slice_max(n = 1000, order_by = avg_log2FC)
top100 <- unique(top100$gene)

library(viridis)
# Make a gene x cell heatmap of genes showing species specific changes
DoHeatmap(sc_seurat.sct.mouse,  
                               assay = "SCT",
                               features = top100 ) + scale_fill_viridis() + NoLegend()

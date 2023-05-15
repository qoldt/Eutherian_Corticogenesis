rm(list=ls())
library(GenomicRanges)
library(rtracklayer)


setwd("~/My Drive/Bench/Opossum/ChIPseq/")

# Load the ChIP-seq peaks BED file
chipseq_peaks <- read.table("peak_files/mouse_specific_H3K27ac_peaks_mm10.bed")
colnames(chipseq_peaks) <- c("chr","start","end","name","score","strand","signalValue", "pValue","qValue","peak")
standard_chrs <- paste0("chr", c(1:19, "X", "Y"))
chipseq_peaks <- chipseq_peaks[chipseq_peaks$chr %in% standard_chrs,]

# Make a GRanges object
chipseq_peaks.gR <- makeGRangesFromDataFrame(chipseq_peaks, keep.extra.columns = TRUE)

# Annotate ChIP peaks using encode enhancer promoter interactions

# Load bed files
enhancer_promoter_mapping <- read.table("Annotations/enc_EPD_enhc-gene.bed")
colnames(enhancer_promoter_mapping) <- c("chrom","chromStart","chromEnd","name","score","SCC","exp","color",
                                         "enhancerChrom","enhancerStart","enhancerEnd","enhancerName",
                                         "enhancerStrand","geneChrom","promoterStart","promoterEnd","geneName","geneStrand", "pvalue")


# Create enhancer and promoter GRanges objects
enhancer_ranges <- GRanges(
  seqnames = enhancer_promoter_mapping$enhancerChrom,
  ranges = IRanges(start = enhancer_promoter_mapping$enhancerStart,
                   end = enhancer_promoter_mapping$enhancerEnd),
  gene_name = enhancer_promoter_mapping$geneName,
  type = "enhancer"
)

promoter_ranges <- GRanges(
  seqnames = enhancer_promoter_mapping$geneChrom,
  ranges = IRanges(start = enhancer_promoter_mapping$promoterStart,
                   end = enhancer_promoter_mapping$promoterEnd),
  gene_name = enhancer_promoter_mapping$geneName,
  type = "promoter"
)

# Annotate function
annotate_peak <- function(peak, enhancer_ranges, promoter_ranges) {
  # Some very clear annotations were being missed, so extend the peak by 300 base pairs on both sides
  extended_peak <- resize(peak, width = width(peak) + 600, fix = "center") # 
  
  overlaps_enhancer <- findOverlaps(extended_peak, enhancer_ranges)
  overlaps_promoter <- findOverlaps(extended_peak, promoter_ranges)
  
  if (length(overlaps_enhancer) > 0) {
    idx <- subjectHits(overlaps_enhancer)[1]
    return(c(gene_name = enhancer_ranges$gene_name[idx],
             type = enhancer_ranges$type[idx]))
  } else if (length(overlaps_promoter) > 0) {
    idx <- subjectHits(overlaps_promoter)[1]
    return(c(gene_name = promoter_ranges$gene_name[idx],
             type = promoter_ranges$type[idx]))
  } else {
    return(c(gene_name = "no_gene", type = "undetermined"))
  }
}

# Convert chipseq_peaks to a list of GRanges objects with one range per element
chipseq_peaks_list <- lapply(seq_len(length(chipseq_peaks.gR)), function(i) chipseq_peaks.gR[i])

library(parallel)
# Determine the number of cores to use
num_cores <- detectCores() - 1 # Use all available cores minus one

# Annotate all peaks in parallel
annotated_peaks <- mclapply(chipseq_peaks_list, function(peak) {
  annotate_peak(peak, enhancer_ranges, promoter_ranges)
}, mc.cores = num_cores)

# Combine results into a data frame
annotated_peaks_df <- do.call(rbind.data.frame, annotated_peaks)
colnames(annotated_peaks_df) <- c("Associated_gene", "type")

# Add the new columns to the chipseq_peaks data.frame
chipseq_peaks <- cbind(chipseq_peaks, annotated_peaks_df)

# lets check the annotation:
library(ggplot2)

ggplot(chipseq_peaks, aes(x = type, fill = type))+
  geom_bar()+
  theme_minimal()


# List of ChromHMM BED files
chromHMM_bed_files <- c(
  "Annotations/E12_Ctx_ChromHMM_18_state_mm10_ENCFF242PDA.bed",
  "Annotations/E14_Ctx_ChromHMM_18_state_mm10_ENCFF562HWK.bed",
  "Annotations/E16_Ctx_ChromHMM_18_state_mm10_ENCFF499HWY.bed",
  "Annotations/E15_Ctx_ChromHMM_18_state_mm10_ENCFF163AVC.bed",
  "Annotations/E11_Ctx_ChromHMM_18_state_mm10_ENCFF749VMS.bed",
  "Annotations/E13_Ctx_ChromHMM_18_state_mm10_ENCFF197PCI.bed",
  "Annotations/P0_ctx_ChromHMM_18_state_mm10_ENCFF758EGD.bed"
)


# Function to read, process, and name data frames from a list of files
read_and_process_files <- function(file_list) {
  # Define standard chromosomes
  standard_chrs <- paste0("chr", c(1:19, "X", "Y"))
  
  # Create an empty list to store the data frames
  df_list <- list()
  
  # Iterate through the file list
  for (file_path in file_list) {
    # Read the file
    df <- read.table(file_path)
    
    # Process the data frame
    df <- df[, 1:4]
    colnames(df) <- c("chr", "start", "end", "state")
    df <- df[df$chr %in% standard_chrs, ]
    
    # Extract the first three characters of the file name as the data frame name
    df_name <- substring(basename(file_path), 1, 3)
    
    # Add the data frame to the list with its name as the list element name
    df_list[[df_name]] <- df
  }
  
  # Return the list of named data frames
  return(df_list)
}

# Call the function with your list of ChromHMM BED files
chromHMM_bed_files_list <- read_and_process_files(chromHMM_bed_files)


#library(BiocParallel)
#
# # Set up the parallel backend
# register(MulticoreParam(workers = 12))
# 
# # # Function to annotate the ChIP-seq peaks with the ChromHMM state from a given BED file
# annotate_peaks_with_chromHMM <- function(chipseq_peaks_df, chromHMM_bed_df) {
#   # Filter out non-standard chromosomes
#   standard_chrs <- paste0("chr", c(1:19, "X", "Y"))
#   chipseq_peaks_df <- chipseq_peaks_df[chipseq_peaks_df$chr %in% standard_chrs,]
#   chromHMM_bed_df <- chromHMM_bed_df[chromHMM_bed_df$chr %in% standard_chrs,]
# 
#   # Convert data.frames to GRanges objects
#   chipseq_peaks <- makeGRangesFromDataFrame(chipseq_peaks_df, keep.extra.columns = TRUE)
#   chromHMM_bed <- makeGRangesFromDataFrame(chromHMM_bed_df, keep.extra.columns = TRUE)
# 
#   # Find the overlaps between the ChIP-seq peaks and the ChromHMM BED file
#   overlaps <- findOverlaps(chipseq_peaks, chromHMM_bed)
# 
#   # Initialize a new column for the ChromHMM state annotations
#   chromHMM_state <- rep(NA, length(chipseq_peaks))
#   max_overlap_width <- rep(0, length(chipseq_peaks))
# 
#   # Handle cases where ChIP-seq peaks overlap with more than one ChromHMM annotation
#   for (i in 1:length(overlaps)) {
#     query_idx <- queryHits(overlaps)[i]
#     subject_idx <- subjectHits(overlaps)[i]
#     current_overlap_width <- width(intersect(chipseq_peaks[query_idx], chromHMM_bed[subject_idx]))
# 
#     if (is.na(chromHMM_state[query_idx]) || current_overlap_width > max_overlap_width[query_idx]) {
#       chromHMM_state[query_idx] <- mcols(chromHMM_bed)[subject_idx, "state"]
#       max_overlap_width[query_idx] <- current_overlap_width
#     }
#   }
# 
#   # Return the annotated ChromHMM state column
#   return(chromHMM_state)
# }
# 
# 
# # Iterate through the ChromHMM data frames and annotate the ChIP-seq peaks
# for (i in seq_along(chromHMM_bed_files_list)) {
#   # Get the current ChromHMM data frame
#   chromHMM_bed_df <- chromHMM_bed_files_list[[i]]
# 
#   # Generate a new column name based on the name of the current list element
#   column_name <- paste0(names(chromHMM_bed_files_list)[i], "_chromHMM_state")
# 
#   # Add a new column to the annotated_peaks data frame with the ChromHMM state annotations
#   chipseq_peaks[[column_name]] <- annotate_peaks_with_chromHMM(chipseq_peaks, chromHMM_bed_df)
# }



# Load necessary libraries
library(GenomicRanges)
library(BiocParallel)
library(foreach)
library(doParallel)

# Function to annotate the ChIP-seq peaks with the ChromHMM state from a given BED file
annotate_peaks_with_chromHMM_parallel <- function(chipseq_peaks_df, chromHMM_bed_df, num_cores = 10) {
  # Filter out non-standard chromosomes
  standard_chrs <- paste0("chr", c(1:19, "X", "Y"))
  chipseq_peaks_df <- chipseq_peaks_df[chipseq_peaks_df$chr %in% standard_chrs,]
  chromHMM_bed_df <- chromHMM_bed_df[chromHMM_bed_df$chr %in% standard_chrs,]
  
  # Convert data.frames to GRanges objects
  chipseq_peaks <- makeGRangesFromDataFrame(chipseq_peaks_df, keep.extra.columns = TRUE)
  chromHMM_bed <- makeGRangesFromDataFrame(chromHMM_bed_df, keep.extra.columns = TRUE)
  
  # Parallel processing
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  chromHMM_state <- foreach(i = seq_along(chipseq_peaks), .packages = c("GenomicRanges")) %dopar% {
    query_range = chipseq_peaks[i]
    overlaps <- findOverlaps(query_range, chromHMM_bed)
    
    if (length(overlaps) == 0) {
      NA
    } else {
      subject_indices <- subjectHits(overlaps)
      overlapping_ranges <- intersect(query_range, chromHMM_bed[subject_indices])
      
      if (length(overlapping_ranges) == 0) {
        NA
      } else {
        overlap_widths <- width(overlapping_ranges)
        largest_overlap_idx <- which.max(overlap_widths)
        result <- mcols(chromHMM_bed)[subject_indices[largest_overlap_idx], "state"]
        if (is.null(result)) {
          message("Error in process_peak for query_idx: ", i)
        }
        result
      }
    }
  }
  
  stopCluster(cl)
  
  # Return the annotated ChromHMM state column
  return(chromHMM_state)
}

# Iterate through the ChromHMM data frames and annotate the ChIP-seq peaks
for (i in seq_along(chromHMM_bed_files_list)) {
  # Get the current ChromHMM data frame
  chromHMM_bed_df <- chromHMM_bed_files_list[[i]]

  # Generate a new column name based on the name of the current list element
  column_name <- paste0(names(chromHMM_bed_files_list)[i], "_chromHMM_state")

  # Add a new column to the annotated_peaks data frame with the ChromHMM state annotations
  chipseq_peaks[[column_name]] <- annotate_peaks_with_chromHMM_parallel(chipseq_peaks, chromHMM_bed_df)
  
  # convert the list column to a vector
  chipseq_peaks[[column_name]] <- unlist(chipseq_peaks[[column_name]])
  
  #print to console
  print(paste("Finished Intersecting ", names(chromHMM_bed_files_list)[i], " ChromHMM with ChIPseq peaks"))
  
}


# Save the annotated_peaks data frame to a new file
write.csv(chipseq_peaks, "annotated_H3K27ac_peaks_mm10_multichromHMM_parallel.csv", quote = FALSE, row.names = FALSE)

key_table <- read.csv("Annotations/mm10_ChromHMM_18_state_key.csv")


# Function to convert ChromHMM states to numbers using the key table
convert_chromHMM_state_to_number <- function(chromHMM_state_col, key_table) {
  chromHMM_state_number_col <- vapply(chromHMM_state_col, function(x) {
    if (is.na(x)) {
      return(NA)
    } else {
      result <- key_table[key_table$chromhmm_state == x, "chromhmm_number"]
      if (length(result) == 0) {
        return(NA)
      } else {
        return(result)
      }
    }
  }, numeric(1))
  return(chromHMM_state_number_col)
}


# Iterate through the ChromHMM data frames and convert ChromHMM state columns to number columns
for (i in seq_along(chromHMM_bed_files_list)) {
  # Get the current ChromHMM state column name
  state_column_name <- paste0(names(chromHMM_bed_files_list)[i], "_chromHMM_state")
  
  # Generate a new column name for the ChromHMM state numbers
  number_column_name <- paste0(names(chromHMM_bed_files_list)[i], "_chromHMM_state_number")
  
  # Add a new column to the annotated_peaks data frame with the ChromHMM state numbers
  chipseq_peaks[[number_column_name]] <- convert_chromHMM_state_to_number(chipseq_peaks[[state_column_name]], key_table)
}



library(pheatmap)
# Create a matrix from the chipseq_peaks
heatmap_data <- as.matrix(chipseq_peaks[, c("E11_chromHMM_state_number",
                                                       "E12_chromHMM_state_number",
                                                       "E13_chromHMM_state_number",
                                                       "E14_chromHMM_state_number",
                                                       "E15_chromHMM_state_number",
                                                       "E16_chromHMM_state_number",
                                                       "P0__chromHMM_state_number")])

# Create a named vector with ChromHMM states and their corresponding colors
chromhmm_states <- sort(unique(unlist(chipseq_peaks[,c("E11_chromHMM_state_number",
                                                          "E12_chromHMM_state_number",
                                                          "E13_chromHMM_state_number",
                                                          "E14_chromHMM_state_number",
                                                          "E15_chromHMM_state_number",
                                                          "E16_chromHMM_state_number",
                                                          "P0__chromHMM_state_number")])))

num_states <- length(chromhmm_states) # 1-17, all but heterochromatin present

# Generate color palette
color_palette <- c("#ff0000", #Active TSS (1)
                   "#ff4500", #Flanking TSS (2)
                   "#008000", #Transcription (3)
                   "#006400", #Weak Transcription (4)
                   "#c2ff05", #Enhancer in Gene. (5)
                   "#ffc34d", #Enhancer (6)
                   "#ffff00", #Weak Enhancer (7)
                   "#bdb76b", #Poised Enhancer (8)
                   "#bdc86b", #Primed Enhancer (9)
                   "#cd5c5c", #Bivalent TSS (10)
                   "#808080", #Repressed by Polycomb (11)
                   "#c0c0c0", #Repressed by Polycomb (weak) (12)
                   "#3cf0f0", #Quiescent Gene (13)
                   "#ffffff", #Quiescent1 (14)
                   "#ffffff", #Quiescent2 (15)
                   "#ffffff", #Quiescent3 (16)
                   "#ffffff" #Quiescent4 (17) 
                   )

# Create named vector
chromhmm_colors <- setNames(color_palette, c("Active TSS (1)", 
                                             "Flanking TSS (2)",
                                             "Transcription (3)",
                                             "Weak Transcription (4)",
                                             "Enhancer in Gene. (5)",
                                             "Enhancer (6)",
                                             "Weak Enhancer (7)",
                                             "Poised Enhancer (8)",
                                             "Primed Enhancer (9)",
                                             "Bivalent TSS (10)",
                                             "Repressed by Polycomb (11)",
                                             "Repressed by Polycomb (weak) (12)",
                                             "Quiescent Gene (13)",
                                             "Quiescent1 (14)",
                                             "Quiescent2 (15)",
                                             "Quiescent3 (16)",
                                             "Quiescent4 (17)" ))

rownames(heatmap_data) <- paste(chipseq_peaks$name, "_", chipseq_peaks$Associated_gene)
heatmap_data_clean <- as.matrix(heatmap_data[!apply(is.na(heatmap_data), 1, any),])


# a lot of entries are TSS (number code (1)) or Bivalent TSS (10) that don't change, so lets remove those that don't change.
# Remove rows where all chromHMM states are equal to 1

# Define a vector of values to filter
values_to_filter <- c(1, 2, 3, 4, 9, 10, 13, 14, 15, 16, 17)

# Function to check if all values in a row are within the values_to_filter vector
is_row_filtered <- function(row) {
  all(sapply(row, function(x) x %in% values_to_filter))
}

# Remove rows where all chromHMM states are in the values_to_filter vector
heatmap_data_filtered <- heatmap_data_clean[!apply(heatmap_data_clean, 1, is_row_filtered),]

# Get the unique values from the heatmap data and sort them
unique_values <- sort(unique(as.vector(heatmap_data_filtered)))

# Generate a named vector of breaks for the color scale
breaks <- setNames(c(unique_values, max(unique_values) + 1), c(unique_values, "Max"))


library(ggplot2)

# Generate the heatmap with a custom ChromHMM legend
da.heat <- pheatmap(heatmap_data_filtered, cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = TRUE,
         main = "ChIP-seq Peaks vs. ChromHMM States",
         color = chromhmm_colors)
da.heat

# Get the dendrogram of the clustered rows
row_hclust <- da.heat$tree_row

# Cut the tree into clusters and get the cluster assignments
row_clusters <- cutree(row_hclust, k = length(row_hclust$order))

# Get the actual rownames from the original data, ordered by cluster assignments
dynamic_enhancers <- rownames(heatmap_data_filtered)[order(row_clusters, row_hclust$order)]

# Extract the text before the gene name using regular expression
dynamic_enhancers <- gsub("(.*) _ .*", "\\1", dynamic_enhancers)

# Add a new column 'change' to chipseq_peaks
chipseq_peaks$change <- ifelse(chipseq_peaks$name %in% dynamic_enhancers, "dynamic", "static")


library(ggplot2)
library(ggwordcloud)

chipseq_peaks.known.enhancers <- subset(chipseq_peaks, Associated_gene != "no_gene" & type == "enhancer" )
chipseq_peaks.known.promoters <- subset(chipseq_peaks, Associated_gene != "no_gene" & type == "promoter")

# Create a word plot using chipseq_peaks_annotated
ggplot(chipseq_peaks.known.enhancers, aes(label = Associated_gene, size = signalValue, color = change)) +
  geom_text_wordcloud_area() +
  scale_size_continuous(range = c(0.5, 10)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Mouse Specific H3K27ac Peaks attributed to Bonafide Enhancers",
       subtitle = "Size based on signal value, color based on type",
       x = NULL, y = NULL,
       color = "change", size = "Signal Value")

ggplot(chipseq_peaks.known.promoters, aes(label = Associated_gene, size = signalValue, color = change)) +
  geom_text_wordcloud() +
  scale_size_continuous(range = c(0.5, 5)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Mouse Specific H3K27ac Peaks attributed to Bonafide Promoters",
       subtitle = "Size based on signal value, color based on type",
       x = NULL, y = NULL,
       color = "change", size = "Signal Value")


sc_data.all <- read.table("~/My Drive/Bench/Opossum/scRNAseq/plots/Pyramid_Subset/Species_Differences_by_Celltype/PN_subset_Species_specific_differences_per_cluster.csv", 
                          sep = ",", header = TRUE)
sc_data.all <- subset(sc_data.all, pct.2 > 0) # remove opossum specific genes

# Function for Capitalization
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

#Fix Gene names so they match with keytype SYMBOL
sc_data.all$gene <- capwords(tolower(sc_data.all$gene))


Species_specific_DEG <- unique(sc_data.all$gene)
write.table(Species_specific_DEG, file="~/Desktop/Species_specific_expression_DEG.csv", quote = FALSE, row.names = FALSE)

chipseq_peaks.promoter <- unique(chipseq_peaks.known.promoters$Associated_gene)
write.table(chipseq_peaks.promoter, file="~/Desktop/CHIP_genes_promoters.csv", quote = FALSE, row.names = FALSE)
chipseq_peaks.enhancer <- unique(chipseq_peaks.known.enhancers$Associated_gene)
write.table(chipseq_peaks.enhancer, file="~/Desktop/CHIP_genes_enhancer.csv", quote = FALSE, row.names = FALSE)




# Add a new column 'expression' to chipseq_peaks
chipseq_peaks$expression <- ifelse(chipseq_peaks$Associated_gene %in% Species_specific_DEG, "yes", "no")

chipseq_peaks.known.enhancers <- subset(chipseq_peaks, Associated_gene != "no_gene" & type == "enhancer" & expression == "yes")
ggplot(chipseq_peaks.known.enhancers, aes(label = Associated_gene, size = signalValue, color = change)) +
  geom_text_wordcloud_area() +
  scale_size_continuous(range = c(0.5, 10)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Mouse Specific H3K27ac Peaks attributed to Bonafide Enhancers",
       subtitle = "Size based on signal value, color based on type",
       x = NULL, y = NULL,
       color = "change", size = "Signal Value")

chipseq_peaks.known.promoters <- subset(chipseq_peaks, Associated_gene != "no_gene" & type == "promoter" & expression == "yes")
ggplot(chipseq_peaks.known.promoters, aes(label = Associated_gene, size = signalValue, color = change)) +
  geom_text_wordcloud() +
  scale_size_continuous(range = c(0.5, 5)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Mouse Specific H3K27ac Peaks attributed to Bonafide Promoters",
       subtitle = "Size based on signal value, color based on type",
       x = NULL, y = NULL,
       color = "change", size = "Signal Value")


setwd(
  "~/Charite Thesis/Intersect analysis/Eutherian_Corticogenesis/Intersect_mm10_specific_H3K27ac_mm10_SATB2"
)

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("EnsDb.Mmusculus.v79")
BiocManager::install("org.Mm.eg.db") # mouse organism annotation
# BiocManager::install("org.Md.eg.db") # opossum organism annotation
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("biomaRt")
BiocManager::install("HDO.db")
BiocManager::install("regioneR")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10") # used by RegioneR
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10.masked")
BiocManager::install("rtracklayer")

library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
library(biomaRt)
library(dplyr)
library(ggimage)
library(regioneR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Mmusculus.UCSC.mm10.masked)
library(rtracklayer)
library(RColorBrewer)
library(ggrepel)

samplefiles <-
  list.files("peak_files", pattern = ".bed", full.names = TRUE)
samplefiles <- as.list(samplefiles)
names(samplefiles) <-
  c(
    "Mouse_SATB2",
    "Mouse_H3K27ac",
    "Mouse_Specific_H3K27ac",
    "Opossum_H3K27ac",
    "Intersect_Mouse_Specific_H3K27ac_Mouse_SATB2"
  )

# Sort in a desired order for the plots

samplefiles <-
  samplefiles[c(
    "Intersect_Mouse_Specific_H3K27ac_Mouse_SATB2",
    "Mouse_SATB2",
    "Mouse_Specific_H3K27ac",
    "Mouse_H3K27ac",
    "Opossum_H3K27ac"
  )]

# Peak annotation

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peakAnnoList <-
  lapply(samplefiles,
         annotatePeak,
         TxDb = txdb,
         tssRegion = c(-1000, 1000))

# Create a list with genes from each sample

genes <- lapply(peakAnnoList, function(i) {
  as.data.frame(i)$geneId
})

compKEGG <- compareCluster(
  geneCluster = genes,
  fun = "enrichKEGG",
  organism = "mouse",
  pvalueCutoff = 0.1,
  #OrgDb = org.Mm.eg.db,
  pAdjustMethod = "BH"
)

dotplot(
  compKEGG,
  showCategory = 7,
  font.size = 8,
  title = "KEGG Enrichment Analysis"
)

# Extract genes ID of the desired/dominant enrichment class

cluster_result_df <- compKEGG@compareClusterResult
intersect_result_df <-
  subset(
    cluster_result_df,
    Cluster == "Intersect_Mouse_Specific_H3K27ac_Mouse_SATB2" &
      Description == "Hippo signaling pathway - Mus musculus (house mouse)"
  )
# genes_in_major_class = intersect_result_df[which.min(intersect_result_df$p.adjust),]$geneID ; genes_in_major_class
genes_in_major_class <-
  intersect_result_df$geneID
genes_in_major_class
genes_in_major_class <- strsplit(genes_in_major_class, split = "/")

# Get annotation dataframe for each cluster

H3K27ac_SATB2_intersect <-
  data.frame(peakAnnoList[["Intersect_Mouse_Specific_H3K27ac_Mouse_SATB2"]]@anno)
SATB2_mouse <- data.frame(peakAnnoList[["Mouse_SATB2"]]@anno)
H3K27ac_mouse_spc <-
  data.frame(peakAnnoList[["Mouse_Specific_H3K27ac"]]@anno)
H3K27ac_mouse_all <-
  data.frame(peakAnnoList[["Mouse_H3K27ac"]]@anno)
H3K27ac_opossum_mm10 <-
  data.frame(peakAnnoList[["Opossum_H3K27ac"]]@anno)

# Get gene names from gene ids

mart <-
  useEnsembl(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")
tmp1 <- data.frame("entrezgene_id" = H3K27ac_SATB2_intersect$geneId)

ids <- getBM(
  attributes = c("entrezgene_id", "external_gene_name"),
  filters = "entrezgene_id",
  values = tmp1[, "entrezgene_id"],
  mart = mart
)

# Merge Gene Annotations

names(ids)[names(ids) == "entrezgene_id"] <- "geneId"
ids$geneId <- as.character(ids$geneId)
H3K27ac_SATB2_intersect <- left_join(
  x = H3K27ac_SATB2_intersect,
  y = ids,
  by = "geneId",
  relationship = "many-to-many"
)

H3K27ac_SATB2_intersect_major_class <-
  subset(H3K27ac_SATB2_intersect, geneId %in% genes_in_major_class[[1]])

# To export
# write.table(H3K27ac_SATB2_intersect_major_class, "GO_transcription_coregulator_genes.csv", row.names=FALSE, sep=";")

# Plot genes in class of interest

H3K27ac_SATB2_intersect_major_class$geneChr <- as.factor(H3K27ac_SATB2_intersect_major_class$geneChr)
c25 <- c(
  "dodgerblue2",
  "#E31A1C",
  # red
  "green4",
  "#6A3D9A",
  # purple
  "#FF7F00",
  # orange
  "grey70",
  "gold1",
  "skyblue2",
  "#FB9A99",
  # lt pink
  "palegreen2",
  "#CAB2D6",
  # lt purple
  "yellow3",
  # lt orange
  "green1",
  "khaki2",
  "maroon",
  "orchid1",
  "deeppink1",
  "blue1",
  "steelblue4",
  "darkturquoise",
  "black",
  "yellow4",
  "darkorange4",
  "#FDBF6F",
  "brown"
)

p <- ggplot(data = H3K27ac_SATB2_intersect_major_class,
            aes(
              x = geneChr,
              y = geneLength,
              col = geneChr,
              label = external_gene_name
            )) +
  ggtitle("Hippo Signaling Pathway Genes") +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c25) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")
p

# Extract genes for enrichGO

intersect.list <- H3K27ac_SATB2_intersect$geneId
SATB2.mouse.list <- SATB2_mouse$geneId
H3K27ac.mouse.spc.list <- H3K27ac_mouse_spc$geneId
H3K27ac.mouse.all.list <- H3K27ac_mouse_all$geneId
H3K27ac.opossum.list <- H3K27ac_opossum_mm10$geneId

# Check that all intersect genes are in the original subsets SATB2 and Mouse specific

H3K27ac_SATB2_intersect_unique_genes <- subset(
  H3K27ac_SATB2_intersect,!(geneId %in% SATB2.mouse.list) &
    !(geneId %in% H3K27ac.mouse.spc.list)
)

# Individual sample enrichment from GO

intersect.ego <- enrichGO(
  gene = intersect.list,
  keyType = "ENTREZID",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Trying with KEGG

intersect.kegg <- enrichKEGG(
  gene = intersect.list,
  keyType = "ncbi-geneid",
  organism = "mouse",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)


upsetplot(intersect.kegg)

###################

SATB2.mouse.ego <- enrichGO(
  gene = SATB2.mouse.list,
  keyType = "ENTREZID",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

H3K27ac.mouse.spc.ego <- enrichGO(
  gene = H3K27ac.mouse.spc.list,
  keyType = "ENTREZID",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

H3K27ac.mouse.all.ego <- enrichGO(
  gene = H3K27ac.mouse.all.list,
  keyType = "ENTREZID",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

H3K27ac.opossum.ego <- enrichGO(
  gene = H3K27ac.opossum.list,
  keyType = "ENTREZID",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Dotplot visualization

intersect.dot <-
  dotplot(intersect.ego, showCategory = 10)
intersect.dot
SATB2.mouse.dot <-
  dotplot(SATB2.mouse.ego, showCategory = 10)
SATB2.mouse.dot
H3K27ac.mouse.spc.dot <-
  dotplot(H3K27ac.mouse.spc.ego, showCategory = 10)
H3K27ac.mouse.spc.dot
H3K27ac.mouse.all.dot <-
  dotplot(H3K27ac.mouse.all.ego, showCategory = 10)
H3K27ac.mouse.all.dot
H3K27ac.opossum.dot <-
  dotplot(H3K27ac.opossum.ego, showCategory = 10)
H3K27ac.opossum.dot


###################### Some plots #########################################


plotAnnoBar(peakAnnoList)
upsetplot(peakAnnoList[["Intersect_Mouse_Specific_H3K27ac_Mouse_SATB2"]], vennpie = TRUE)
vennplot(genes) #error shows genes specific to the intersect
plotDistToTSS(peakAnnoList, title = "Distribution of transcription factor-binding loci \n relative to TSS")

peak <-
  readPeakFile("peak_files/mouse_specific_H3K27ac_peaks_mm10.bed")
covplot(peak, weightCol = "V5", lower = 200)


##################### Permutation test ###################################


getGenomeAndMask("mm10")

SATB2_peaks <- import("peak_files/GSE222608_Peaks_SATB2_sorted.bed", format =
                        "BED")

# not a normal bed file because of extra data so select the 3 first useful columns
tmp <- read.table("peak_files/mouse_specific_H3K27ac_peaks_mm10.bed")
Mouse_specific_peaks <- GRanges(tmp[, 1], IRanges(tmp[, 2], tmp[, 3]))

# Filter to have only canonical
SATB2_peaks <- filterChromosomes(SATB2_peaks, organism = "mm10", chr.type =
                                   "canonical")
Mouse_specific_peaks <- filterChromosomes(Mouse_specific_peaks,
                                          organism = "mm10",
                                          chr.type = "canonical")

len_SATB2_peaks <- length(SATB2_peaks)
len_SATB2_peaks
len_Mouse_specific_peaks <- length(Mouse_specific_peaks)
len_Mouse_specific_peaks

# Permutation test
intersect_analysis <- overlapPermTest(
  A = Mouse_specific_peaks,
  B = SATB2_peaks,
  ntimes = 500,
  genome = "mm10",
  count.once = TRUE,
  mask = NA
)
intersect_analysis
plot(intersect_analysis, ylim = c(0, 0.08))
intersect_analysis$numOverlaps$observed / len_Mouse_specific_peaks # proportion of overlapping peaks
pvalue_specific <- intersect_analysis$numOverlaps$pval

# Odds ratio

random_overlap_nb_specific <- mean(intersect_analysis$numOverlaps$permuted)
observed_overlap_nb_specific <- intersect_analysis$numOverlaps$observed
random_nonoverlap_nb_specific <- len_Mouse_specific_peaks - random_overlap_nb_specific
observed_nonoverlap_nb_specific <- len_Mouse_specific_peaks - observed_overlap_nb_specific

OR_specific = (observed_overlap_nb_specific / random_overlap_nb_specific) /
  (observed_nonoverlap_nb_specific / random_nonoverlap_nb_specific)
log(OR_specific)

# Zscore when Mouse specific peaks are shifted
lz_intersect_analysis_local <- localZScore(
  A = Mouse_specific_peaks,
  B = SATB2_peaks,
  pt = intersect_analysis,
  window = 1000,
  step = 500,
  count.once = TRUE
)
plot(lz_intersect_analysis_local)

# We can see a broader effect which might not be correlated to SATB2
lz_intersect_analysis_regional <- localZScore(
  A = Mouse_specific_peaks,
  B = SATB2_peaks,
  pt = intersect_analysis,
  window = 10000,
  step = 5000,
  count.once = TRUE
)
plot(lz_intersect_analysis_regional)

# SATB2 overlap analysis with all mouse and opossum peaks
tmp <- read.table("peak_files/mouse.H3K27ac_peaks.bed")
Mouse_H3K27ac_peaks <- GRanges(tmp[, 1], IRanges(tmp[, 2], tmp[, 3]))
Opossum_H3K27ac_peaks <- import("peak_files/opossum_H3K27ac_peaks.merged.bed", format =
                                  "BED")

Mouse_H3K27ac_peaks <- filterChromosomes(Mouse_H3K27ac_peaks,
                                         organism = "mm10",
                                         chr.type = "canonical")
Opossum_H3K27ac_peaks <- filterChromosomes(Opossum_H3K27ac_peaks,
                                           organism = "mm10",
                                           chr.type = "canonical")

len_Mouse_H3K27ac_peaks <- length(Mouse_H3K27ac_peaks)
len_Opossum_H3K27ac_peaks <- length(Opossum_H3K27ac_peaks)

Mouse_H3K27ac_SATB2_analysis <- overlapPermTest(
  A = Mouse_H3K27ac_peaks,
  B = SATB2_peaks,
  ntimes = 500,
  genome = "mm10",
  count.once = TRUE,
  mask = NA
)
Mouse_H3K27ac_SATB2_analysis
plot(Mouse_H3K27ac_SATB2_analysis, ylim = c(0, 0.08))
Mouse_H3K27ac_SATB2_analysis$numOverlaps$observed / len_Mouse_H3K27ac_peaks
pvalue_all <- Mouse_H3K27ac_SATB2_analysis$numOverlaps$pval

random_overlap_nb_all <- mean(Mouse_H3K27ac_SATB2_analysis$numOverlaps$permuted)
observed_overlap_nb_all <- Mouse_H3K27ac_SATB2_analysis$numOverlaps$observed
random_nonoverlap_nb_all <- len_Mouse_H3K27ac_peaks - random_overlap_nb_all
observed_nonoverlap_nb_all <- len_Mouse_H3K27ac_peaks - observed_overlap_nb_all

OR_all = (observed_overlap_nb_all / random_overlap_nb_all) / (observed_nonoverlap_nb_all /
                                                                random_nonoverlap_nb_all)
log(OR_all)

Opossum_H3K27ac_SATB2_analysis <- overlapPermTest(
  A = Opossum_H3K27ac_peaks,
  B = SATB2_peaks,
  ntimes = 500,
  genome = "mm10",
  count.once = TRUE,
  mask = NA
)
Opossum_H3K27ac_SATB2_analysis
plot(Opossum_H3K27ac_SATB2_analysis, ylim = c(0, 0.08))
Opossum_H3K27ac_SATB2_analysis$numOverlaps$observed / len_Opossum_H3K27ac_peaks
pvalue_opossum <- Opossum_H3K27ac_SATB2_analysis$numOverlaps$pval

random_overlap_nb_opossum <- mean(Opossum_H3K27ac_SATB2_analysis$numOverlaps$permuted)
observed_overlap_nb_opossum <- Opossum_H3K27ac_SATB2_analysis$numOverlaps$observed
random_nonoverlap_nb_opossum <- len_Opossum_H3K27ac_peaks - random_overlap_nb_opossum
observed_nonoverlap_nb_opossum <- len_Opossum_H3K27ac_peaks - observed_overlap_nb_opossum

OR_opossum = (observed_overlap_nb_opossum / random_overlap_nb_opossum) /
  (observed_nonoverlap_nb_opossum / random_nonoverlap_nb_opossum)
log(OR_opossum)

odds_ratio <- c(OR_specific, OR_all, OR_opossum)
feature <- c('Mouse_Specific_H3K27ac', 'Mouse_H3K27ac', 'Opossum_H3K27ac')
pval <- c(pvalue_specific, pvalue_all, pvalue_opossum)
df <- data.frame(feature, odds_ratio, pval)

##### Plot

x_max <- max(abs(log(df$odds_ratio)), na.rm = TRUE)

ggplot(df, aes(log(odds_ratio), feature, size = pval)) +
  geom_point(aes(color = feature)) +
  labs(y = 'Tested dataset', x = 'log(odds ratio)') +
  geom_vline(xintercept = 0) +
  xlim(c(-x_max - 0.5, x_max + 0.5)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype = 3, linewidth = 1.2)
  )


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

library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
library(biomaRt)
library(dplyr)
library(ggimage)


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
# samplefiles <- samplefiles[c("Mouse_SATB2", "Intersect_Mouse_Specific_H3K27ac_Mouse_SATB2")]

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
  # OrgDb = org.Mm.eg.db,
  pAdjustMethod = "BH"
)

dotplot(
  compKEGG,
  showCategory = 7,
  font.size = 8,
  title = "GO Enrichment Analysis"
)

# If use fun = "enrichGO"
# ego2 <- simplify(compKEGG)
# cnetplot(ego2)

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
# write.csv(H3K27ac_SATB2_intersect_major_class, "hippo_pathway_genes.csv", row.names=FALSE)


# Extract genes for enrichGO

intersect.list <- H3K27ac_SATB2_intersect$geneId
SATB2.mouse.list <- SATB2_mouse$geneId
H3K27ac.mouse.spc.list <- H3K27ac_mouse_spc$geneId
H3K27ac.mouse.all.list <- H3K27ac_mouse_all$geneId
H3K27ac.opossum.list <- H3K27ac_opossum_mm10$geneId

# Check that all intersect genes are in the original subsets SATB2 and Mouse specific

H3K27ac_SATB2_intersect_unique_genes <- subset(
  H3K27ac_SATB2_intersect,
  !(geneId %in% SATB2.mouse.list) &
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

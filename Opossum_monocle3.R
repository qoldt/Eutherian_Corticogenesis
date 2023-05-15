# https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p
# installation of Monocle3: https://cole-trapnell-lab.github.io/monocle3/docs/installation/
# also need the seurat wrapper: devtools::install_github("satijalab/seurat-wrappers")

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

setwd("~/My Drive/Bench/Opossum/scRNAseq/")
sc_seurat.sct <- readRDS("sc_seurat.PN.sct.rds")

DefaultAssay(sc_seurat.sct) <- "SCT"

# mouse_seurat <- subset(sc_seurat.sct, Species == "Mouse")
# mouse_monocle <- as.cell_data_set(mouse_seurat) # Seurat wrapper library needed for this
# cds <- as.cell_data_set(sc_seurat.sct)
# # 
# cds <- mouse_monocle
# # using monocle3 preprocessing
# cds <- preprocess_cds(mouse_monocle, num_dim = 50)
# cds <- align_cds(cds, alignment_group = "Stage", model_formula_str = "~Stage")
# cds <- reduce_dimension(cds)
# cds <- cluster_cells(cds)
# plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype")
# cds <- learn_graph(cds)
# plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype")
# cds <- order_cells(cds)
# plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
#            label_branch_points = T, label_roots = T, label_leaves = F,
#            group_label_size = 5)
# 
# 
#using seurat wrapper
mouse_seurat <- subset(sc_seurat.sct, Species == "Mouse")
mouse_monocle <- as.cell_data_set(mouse_seurat) # Seurat wrapper library needed for this
fData(mouse_monocle)$gene_short_name <- rownames(fData(mouse_monocle))
opossum_seurat <- subset(sc_seurat.sct, Species == "Opossum")
opossum_monocle <- as.cell_data_set(opossum_seurat) # Seurat wrapper library needed for this
fData(opossum_monocle)$gene_short_name <- rownames(fData(opossum_monocle))

cds <- mouse_monocle
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions
list.cluster <- mouse_seurat@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- mouse_seurat@reductions$umap@cell.embeddings

cds <- learn_graph(cds, use_partition = F)
cds.mouse <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) %in%
                                                                              c("Apical Progenitors I",
                                                                                "Apical Progenitors II",
                                                                                 "Apical Progenitors III")]))
rm(cds)
cds.mouse$monocle3_pseudotime <- pseudotime(cds.mouse)
data.pseudo.mouse <- as.data.frame(colData(cds.mouse))


# now for opossum
cds <- opossum_monocle
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions
list.cluster <- opossum_seurat@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- opossum_seurat@reductions$umap@cell.embeddings

cds <- learn_graph(cds, use_partition = F)
cds.opossum <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) %in%
                                                                               c("Apical Progenitors I",
                                                                                 "Apical Progenitors II",
                                                                                 "Apical Progenitors III")]))
rm(cds)
cds.opossum$monocle3_pseudotime <- pseudotime(cds.opossum)
data.pseudo.opossum <- as.data.frame(colData(cds.opossum))

data.pseudo <- rbind(data.pseudo.mouse, data.pseudo.opossum)

data.pseudo$Species <- factor(data.pseudo$Species, levels = c("Mouse","Opossum"))



M <- plot_cells(cds.mouse, color_cells_by = "pseudotime", show_trajectory_graph = FALSE,
                label_groups_by_cluster = T,
                label_branch_points = T, label_roots = F, label_leaves = F,
                alpha = 0.3) + ggtitle("Mouse")+ coord_fixed()

O <- plot_cells(cds.opossum,  color_cells_by = "pseudotime", show_trajectory_graph = FALSE,
                label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F,
           alpha = 0.3) + ggtitle("Opossum") +coord_fixed()

M+O



ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(celltype.species, monocle3_pseudotime), colour = Stage)) + 
  geom_boxplot(inherit.aes = FALSE, 
               data = data.pseudo, 
               aes(monocle3_pseudotime, reorder(celltype.species, monocle3_pseudotime)),
               size = 0.2,
               outlier.color = 0.1)+
  geom_jitter(size = 0.2, alpha = 0.1) +
  #facet_grid( Species~.)+
  scale_colour_manual(values = c("yellow","orange","red","green","blue","purple"))+
  theme_minimal()+
  xlab("pseudotime")+
  ylab("Celltype")


ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(celltype, monocle3_pseudotime), colour = Stage)) + 
  geom_boxplot(inherit.aes = FALSE, 
               data = data.pseudo, 
               aes(monocle3_pseudotime, reorder(celltype, monocle3_pseudotime)),
               size = 0.2,
               outlier.color = 0.1)+
  geom_jitter(size = 0.2, alpha = 0.1) +
  facet_grid(Species~.)+
  theme_minimal()+
  scale_colour_manual(values = c("red","orange","yellow","green","blue","purple"))+
  xlab("pseudotime")+
  ylab("Celltype")


# Find genes that change as a function of pseudotime
# these tests take a long time
# Mouse
deg.m <- graph_test(cds.mouse, neighbor_graph = "principal_graph")
deg.m %>% arrange(q_value) %>% filter(status == "OK") %>% head()
# Opossum
deg.o <- graph_test(cds.opossum, neighbor_graph = "principal_graph")
deg.o %>% arrange(q_value) %>% filter(status == "OK") %>% head()

gene.list <- c("PLAGL1","EOMES", "NEUROD1", "BTG2", "HS6ST2")

my_genes.m <- row.names(subset(fData(cds.mouse), gene_short_name %in% gene.list)) 
mouse_subset <- cds.mouse[my_genes.m,]
M.pt <- plot_genes_in_pseudotime(mouse_subset, color_cells_by = "monocle3_pseudotime" )

my_genes.o <- row.names(subset(fData(cds.opossum), gene_short_name %in% gene.list)) 
opossum_subset <- cds.opossum[my_genes.o,]
O.pt <- plot_genes_in_pseudotime(opossum_subset, color_cells_by = "monocle3_pseudotime" )

M.pt + O.pt

plot_genes_in_pseudotime(subset(cds.mouse, gene_short_name == "PAX6" ),
                         color_cells_by="Stage",
                         min_expr=0.5)


# For mouse
plot_cells(mouse_monocle, color_by = "pseudotime")

# For opossum
plot_cells(opossum_monocle, color_by = "pseudotime")



rm(list=ls()) ; gc();

set.seed(42)

library(Seurat)
library(ggplot2)
library(readr)
library(grid)
library(gridExtra)
library(edgeR)
library(FactoMineR)
library(dplyr) 
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(EnhancedVolcano)
library(org.Hs.eg.db)

public <- readRDS(
  "/storage/zhangyanxiaoLab/niuyuxiao/projects/AD/data/zhang_CR/GSE187332_DH_final_seurat_object.RDS"
)
Idents(public) <- "cell_type"
DimPlot(public, group.by = "cell_type", reduction = "umap", label = TRUE,  label.size = 5 ) 

setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/doublet/brainHip")


mouse <- readRDS("data/mouse_all.rds")

gila <- subset(
  mouse,
  subset = seurat_clusters %in%  c(0, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 15, 18, 19, 20, 24, 28, 29, 41)
)

DefaultAssay(object = gila) <- "RNA"
gila = NormalizeData(gila, normalization.method = "LogNormalize", scale.factor = 10000)
gila = FindVariableFeatures(gila, selection.method = "vst", nfeatures = 3000)  
gila = ScaleData(gila, vars.to.regress = c("nFeature_RNA", "percent_mt"))
gila <- RunPCA(gila, features = VariableFeatures(object = gila))
reduction = "pca"
ndims = 1:20 # make it smaller when like like a star
gila <- FindNeighbors(gila, dims = ndims, reduction=reduction, k.param=15)

gila <- RunUMAP(gila, dims = ndims, reduction=reduction)
gila <- FindClusters(gila, resolution = 0.99)



p1 <- DimPlot( gila, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4 )
p2 <- DimPlot( gila, reduction = "umap", group.by = "cell_anno", label = TRUE, label.size = 4 ) + ggtitle("cell_anno_GY")
p3 <- DimPlot( gila, reduction = "umap", group.by = "main_anno", label = TRUE, label.size = 4 ) + ggtitle("main_anno_GY")

p1 + (p2/p3)
table(gila$seurat_clusters, gila$cell_anno)

saveRDS(gila, file = "data/mouse_ExN_high_resolution.rds")
#########

p1 <- DimPlot(gila, reduction = "umap", group.by="seurat_clusters",  label = T, label.size = 5) + ggtitle("brianHip ExN re-cluster")
p2 <- DimPlot(gila, reduction = "umap", group.by="scDblFinder.class",  label = T, label.size = 5) + ggtitle("brainFC scDblFinder.class")
p2 <- FeaturePlot(gila, features = "scDblFinder.score", reduction = "umap", order=T, label=T)

# 查看每个cluster中双细胞分类的分布
cluster_doublet_table <- table(gila$seurat_clusters, gila$scDblFinder.class)
singlet_proportions <- prop.table(cluster_doublet_table, margin = 1)[, "singlet"]
singlet_summary <- data.frame(
  cluster = rownames(cluster_doublet_table),
  total_cells = rowSums(cluster_doublet_table),
  singlet_count = cluster_doublet_table[, "singlet"],
  doublet_count = cluster_doublet_table[, "doublet"],
  singlet_proportion = round(singlet_proportions * 100, 2)
)

singlet_summary$cluster <- factor(
  singlet_summary$cluster,
  levels = sort(as.numeric(unique(singlet_summary$cluster)))
)

p3 <- ggplot(singlet_summary, aes(x = cluster, y = singlet_proportion)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(singlet_proportion, "%")), vjust = -0.5) +
  ylim(0, 100) +
  labs(title = "Singlet Proportion by Cluster",
       x = "Seurat Cluster",
       y = "Singlet Proportion (%)") +
  theme_minimal()

(p1 + p2) / p3

all_markers <- c("Slc17a7","Gad2",	
                 "Apoe", # Asc
                 "C1qb",
                 "Mog", # C1qb
                 "Pdgfra", "Tie1","Foxd1")

Idents(gila) <- "seurat_clusters"
FeaturePlot(gila, features = all_markers, reduction = "umap", ncol = 4, order=T, label=T)
DotPlot(gila, features = all_markers) + RotatedAxis()

############# cluster 14 ##########  
cluster_id <- "14"
table(subset(gila, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(gila, seurat_clusters == cluster_id), features = all_markers, reduction = "umap", ncol = 4, order=T)

cells_use <- WhichCells(gila, idents = cluster_id)
n_cells <- length(cells_use)
p1 <- DimPlot( gila, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( gila, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id, "  (n = ", n_cells, ")"))

p2 <- DotPlot( subset( gila, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(gila) <- "seurat_clusters"

gila$one_vs_others <- ifelse(
  gila$seurat_clusters == cluster_id, cluster_id,
  ifelse(gila$seurat_clusters %in% c(14), "exclude", "others")
)
Idents(gila) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(gila, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(gila) <- "seurat_clusters"
FeaturePlot(gila, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 


p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label=T,
               cells.highlight = WhichCells(gila, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("Sub-cluster ", cluster_id, " highlighted in previous all-cell clusters"))

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p2 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainCB"))
p1 + p2


gila <- AddModuleScore( object = gila, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(gila, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainCB"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainCB\n source:CR"))
p1 + p2


VlnPlot( gila, features = c("Fam19a1","Rasgrf2"), ncol = 1,pt.size = 0 )

############# cluster 17 ##########  
cluster_id <- "17"
table(subset(gila, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(gila, seurat_clusters == cluster_id), features = all_markers, reduction = "umap", ncol = 4, order=T)

cells_use <- WhichCells(gila, idents = cluster_id)
n_cells <- length(cells_use)
p1 <- DimPlot( gila, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( gila, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id, "  (n = ", n_cells, ")"))

p2 <- DotPlot( subset( gila, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(gila) <- "seurat_clusters"

gila$one_vs_others <- ifelse(
  gila$seurat_clusters == cluster_id, cluster_id,
  ifelse(gila$seurat_clusters %in% c(17), "exclude", "others")
)
Idents(gila) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(gila, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(gila) <- "seurat_clusters"
FeaturePlot(gila, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 


p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label=T,
               cells.highlight = WhichCells(gila, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("Sub-cluster ", cluster_id, " highlighted in previous all-cell clusters"))

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p2 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainCB"))
p1 + p2


gila <- AddModuleScore( object = gila, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(gila, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainCB"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainCB\n source:CR"))
p1 + p2

FeaturePlot(gila, features = c("Slc17a7","Lingo1","Sv2b", "Fibcd1", "Ccdc88c","Lmo1","Tbc1d1","Mpped1"), ncol = 4, reduction = "umap", order=T, label=T) 
VlnPlot( gila, features = c("Slc17a7","Lingo1","Sv2b", "Fibcd1"), ncol = 1,pt.size = 0 )
VlnPlot( gila, features = c( "Ccdc88c","Lmo1","Tbc1d1","Mpped1"), ncol = 1,pt.size = 0 )




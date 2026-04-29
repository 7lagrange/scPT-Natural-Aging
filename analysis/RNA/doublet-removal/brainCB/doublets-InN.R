
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

setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/doublet/brainCB")

public <- readRDS("data/brainCB_public_10k.rds")
Idents(public) <- "cell_type"
DimPlot(public, group.by = "cell_type", reduction = "umap", label = TRUE,  label.size = 5 ) 

setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/doublet/brainCB")

rds_path <- "data/mouse_all.rds"
mouse <- readRDS(rds_path)

gila <- subset(
  mouse,
  subset = seurat_clusters %in% c(4, 5, 10, 16, 23, 26)
)

DefaultAssay(object = gila) <- "RNA"
gila = NormalizeData(gila, normalization.method = "LogNormalize", scale.factor = 10000)
gila = FindVariableFeatures(gila, selection.method = "vst", nfeatures = 3000)  
gila = ScaleData(gila, vars.to.regress = c("nFeature_RNA", "percent_mt"))
gila <- RunPCA(gila, features = VariableFeatures(object = gila))
reduction = "pca"
ndims = 1:20 # make it smaller when like like a star
gila <- FindNeighbors(gila, dims = ndims, reduction=reduction, k.param=15)
gila <- FindClusters(gila, resolution = 0.5)
gila <- RunUMAP(gila, dims = ndims, reduction=reduction)



p1 <- DimPlot( gila, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4 )
p2 <- DimPlot( gila, reduction = "umap", group.by = "cell_anno", label = TRUE, label.size = 4 ) + ggtitle("cell_anno_GY")
p3 <- DimPlot( gila, reduction = "umap", group.by = "main_anno", label = TRUE, label.size = 4 ) + ggtitle("main_anno_GY")

p1 + (p2/p3)
table(gila$seurat_clusters, gila$cell_anno)

saveRDS(gila, file = "data/mouse_InN.rds")

#########

p1 <- DimPlot(gila, reduction = "umap", group.by="seurat_clusters",  label = T, label.size = 5) + ggtitle("brainCB InN re-cluster")
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

all_markers <- c("Gabra6", # GC
                 "Eomes",	 # UBC
                 "Prkcd", # 
                 "Klhl1",
                 "Ppp1r17",
                 "Lgi2", # C1qb
                 "Mog",
                 "Pdgfra",
                 "Apoe",
                 "C1qb",
                 "Tie1",
                 "Dcn"
)

Idents(gila) <- "seurat_clusters"
FeaturePlot(gila, features = all_markers, reduction = "umap", ncol = 4, order=T, label=T)
DotPlot(gila, features = all_markers) + RotatedAxis()
FeaturePlot(gila, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_ribo"), ncol = 2, reduction = "umap", order=T, label=T) 
VlnPlot( gila, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_ribo"), ncol = 1,pt.size = 0 )



DimPlot( gila, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label=T,
         cells.highlight = WhichCells( mouse, idents = '26' ), cols = "grey85", cols.highlight = "red" 
) + ggtitle(paste0("all-cells-cluster 26 in InN subcluster"))



############# cluster 6 ##########  
cluster_id <- "6"
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
  ifelse(gila$seurat_clusters %in% c(6), "exclude", "others")
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

############# cluster 10 ##########  
cluster_id <- "10"
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
  ifelse(gila$seurat_clusters %in% c(10), "exclude", "others")
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

FeaturePlot(gila, features = c('Klhl1',"Aldh1a3", "Slc6a5","Htr2a", "Edil3"), 
            ncol = 3, reduction = "umap", order=T, label=T) 
VlnPlot( gila, features = c('Klhl1',"Aldh1a3", "Slc6a5","Htr2a", "Edil3"), ncol = 1,pt.size = 0 )

FeaturePlot(gila, features = c("Prkcd", "Nxph1", "Cdh22"), ncol = 3, reduction = "umap", order=T, label=T) 
VlnPlot( gila, features = c("Prkcd", "Nxph1", "Cdh22"), ncol = 1,pt.size = 0 )

############# cluster 11 ##########  
cluster_id <- "11"
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
  ifelse(gila$seurat_clusters %in% c(11), "exclude", "others")
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

############# cluster 12 ##########  
cluster_id <- "12"
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
  ifelse(gila$seurat_clusters %in% c(12), "exclude", "others")
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
############# cluster 13 ##########  
cluster_id <- "13"
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
  ifelse(gila$seurat_clusters %in% c(13), "exclude", "others")
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


############# cluster 15 ##########  
cluster_id <- "15"
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
  ifelse(gila$seurat_clusters %in% c(15), "exclude", "others")
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

############# cluster 2 ##########  
cluster_id <- "2"
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
  ifelse(gila$seurat_clusters %in% c(2), "exclude", "others")
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

############# cluster 7 ##########  
cluster_id <- "7"
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
  ifelse(gila$seurat_clusters %in% c(7), "exclude", "others")
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

gila$double_manual <- "singlet"
doublet_clusters <-c(6)
gila$double_manual[gila$seurat_clusters %in% doublet_clusters] <- "doublet"
DimPlot(gila, reduction = "umap", group.by="double_manual",  label = T, label.size = 5) + ggtitle("InN sub cluster doublets")



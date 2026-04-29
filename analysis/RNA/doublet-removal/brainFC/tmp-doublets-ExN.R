
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
  "/storage/zhangyanxiaoLab/niuyuxiao/projects/AD/data/zhang_CR/GSE187332_FC_final_seurat_object.RDS"
)
DimPlot(public, group.by = "cell_type", reduction = "umap", label = TRUE,  label.size = 5 ) 

setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/doublet/brainFC")

rds_path <- "data/mouse_all.rds"
mouse <- readRDS(rds_path)



gila <- subset(
  mouse,
  subset = seurat_clusters %in% c(0, 3, 4, 5, 7, 8, 11, 12, 14, 18, 25, 26, 27, 39)
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


#########

p1 <- DimPlot(gila, reduction = "umap", group.by="seurat_clusters",  label = T, label.size = 5) + ggtitle("brainFC ExN re-cluster")
p2 <- DimPlot(gila, reduction = "umap", group.by="scDblFinder.class",  label = T, label.size = 5) + ggtitle("brainFC scDblFinder.class")

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

p1 + p2 / p3

all_markers <- c("Slc17a7","Gad2",	
                 "Apoe", # Asc
                 "C1qb",
                 "Mog", # C1qb
                 "Pdgfra", "Tie1","Foxd1")

Idents(gila) <- "seurat_clusters"
FeaturePlot(gila, features = all_markers, reduction = "umap", ncol = 4, order=T, label=T)
DotPlot(gila, features = all_markers) + RotatedAxis()

# potential-doublets-clusters
DotPlot( subset( gila, idents = c(12, 13) ), features = all_markers, group.by = "seurat_clusters" ) + RotatedAxis()
Idents(mouse) <- "cell_anno"

############# cluster 11 ##########
cluster_id <- "12"
FeaturePlot(subset(gila, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( gila, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( gila, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( gila, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(gila) <- "seurat_clusters"

gila$one_vs_others <- ifelse(
  gila$seurat_clusters == cluster_id, cluster_id,
  ifelse(gila$seurat_clusters %in% c(12, 13), "exclude", "others")
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
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))
p1 + p2


gila <- AddModuleScore( object = gila, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(gila, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2

############# cluster 11 ##########
cluster_id <- "13"
FeaturePlot(subset(gila, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                      "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( gila, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( gila, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( gila, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(gila) <- "seurat_clusters"

gila$one_vs_others <- ifelse(
  gila$seurat_clusters == cluster_id, cluster_id,
  ifelse(gila$seurat_clusters %in% c(12, 13), "exclude", "others")
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
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))
p1 + p2


gila <- AddModuleScore( object = gila, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(gila, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2


saveRDS(gila, file = "data/mouse_ExN.rds")











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
  subset = seurat_clusters %in% c(0, 1, 2, 3, 7, 12, 13, 19, 20)
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

p1 <- DimPlot(gila, reduction = "umap", group.by="seurat_clusters",  label = T, label.size = 5) + ggtitle("brainCB ExN re-cluster")
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
         cells.highlight = WhichCells( mouse, idents = '2' ), cols = "grey85", cols.highlight = "red" 
) + ggtitle(paste0("all-cells-cluster 2 in ExN subcluster"))

saveRDS(gila, file = "data/mouse_ExN.rds")









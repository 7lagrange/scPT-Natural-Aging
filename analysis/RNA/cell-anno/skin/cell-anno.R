
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
library(harmony)


setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/cell-anno/skin")


load("/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/results/skin-merged/mouse_after_cellanno.RData")


######### harmony ######## https://mp.weixin.qq.com/s/9Y3DA_3xaSKRUKe9Q481aQ

mouse$histone_age_sample <- paste0(mouse$histone, '-', mouse$age, '-', mouse$mouse)
mouse <- RunHarmony( object = mouse, group.by.vars = "histone_age_sample", reduction = "pca", reduction.save ="harmony",
                     dims.use = 1:20, assay.use = "RNA", verbose = TRUE, plot_convergence = TRUE )
mouse <- RunUMAP(mouse, 
                 dims = 1:20,
                 n.neighbors = 30,
                 min.dist = 0.3,
                 reduction = "harmony",
                 reduction.name = "umap_harmony"
)

mouse <- FindNeighbors(mouse, dims = 1:20, reduction="harmony", k.param=15)
mouse <- FindClusters(mouse, resolution = c(0.5))

DimPlot(mouse, reduction = "umap_harmony", group.by="bc1",  label = T, label.size = 5)
DimPlot(mouse, reduction = "umap_harmony", group.by="cell_type",  label = T, label.size = 5)

p1 <- DimPlot(mouse, reduction = "umap_harmony", group.by="RNA_snn_res.0.5",  label = T, label.size = 5, raster=FALSE)
p2 <- DimPlot(mouse, reduction = "umap_harmony", group.by="age",  label = T, label.size = 5, raster=FALSE)
p3 <- DimPlot(mouse, reduction = "umap_harmony", group.by="mouse",  label = T, label.size = 5, raster=FALSE)
p4 <- DimPlot(mouse, reduction = "umap_harmony", group.by="histone",  label = T, label.size = 5, raster=FALSE)
(p1 + p2)  / (p3+p4)

p1

########### add cell type from text ############
Idents(mouse) <- "seurat_clusters"
new.cluster.ids <- read.table(file = "cellAnno.txt", header = TRUE, sep='\t')$label
new.cluster.ids <- gsub("/", "_", new.cluster.ids)
new.cluster.ids <- gsub("\\(", "", new.cluster.ids)
new.cluster.ids <- gsub("\\)", "", new.cluster.ids)

names(new.cluster.ids) <- levels(mouse)
mouse <- RenameIdents(mouse, new.cluster.ids)
mouse$cell_type <- mouse@active.ident
mouse$cell_type_cluster <- paste0(mouse@active.ident, mouse$seurat_clusters)


p1 <- DimPlot(mouse, reduction = "umap_harmony", group.by="cell_type", label = T, label.size = 6,raster=FALSE)
p1

mouse <- subset(mouse, subset = cell_type != "low_qual")


meta <- mouse@meta.data
meta$cell_id <- rownames(meta)
table(meta$cell_type)

write.csv(meta, 
          file = "~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/metadata/data/skin-final-cell-type.csv", 
          row.names = FALSE)

saveRDS(mouse, file = "data/skin_cell_anno.rds")





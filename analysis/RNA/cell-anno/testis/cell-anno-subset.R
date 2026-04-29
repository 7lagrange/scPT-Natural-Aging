
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


setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/cell-anno/testis/")



public <- readRDS("data/SPC_RS.combination.rds")
public$celltype <- public$cluster
Idents(public) <- "celltype"
DimPlot(public, reduction = "umap", group.by="cluster",  label = T, label.size = 5)

########## keep subset cells ###########
mouse_old <- readRDS("data/testis_cell_anno.rds")
mouse <- subset(mouse_old, subset = seurat_clusters %in% c(0,3,4,7,9,1,2,8,14,18,21,10,22,20,6,16,11,12))
DimPlot(mouse, reduction = "umap_harmony", group.by="seurat_clusters",  label = T, label.size = 5)
cells.keep <- rownames(Embeddings(mouse, "umap_harmony"))[
  Embeddings(mouse, "umap_harmony")[,2] < 10
]
mouse <- subset(mouse, cells = cells.keep)

DimPlot(mouse, reduction = "umap_harmony", group.by="seurat_clusters",  label = T, label.size = 5)
DimPlot(mouse, reduction = "umap_harmony", group.by="cell_type",  label = T, label.size = 5)

#saveRDS(mouse, file = "data/testis_subset_anno.rds")


mouse <- readRDS("data/testis_subset_anno.rds")
cat("Min count:", min(GetAssayData(mouse, slot = "counts")), "\n")
cat("Max count:", max(GetAssayData(mouse, slot = "counts")), "\n")

################### try ambinet RNA removal ############
library(celda)  # DecontX 在 celda 包里

# 转成 SingleCellExperiment
library(SingleCellExperiment)
sce <- SingleCellExperiment(
  assays = list(counts = GetAssayData(mouse, slot = "counts"))
)
sce$cluster <- mouse$seurat_clusters

# 运行 DecontX
sce <- decontX(sce, z = sce$cluster)
#saveRDS(sce, file = "data/deconX_sce.rds")

# 查看污染比例
summary(sce$decontX_contamination)
hist(sce$decontX_contamination, breaks = 50,
     main = "DecontX contamination per cell",
     xlab = "Contamination fraction")

mouse[["RNA"]] <- CreateAssayObject(
  counts = round(decontXcounts(sce))
)

# 验证替换成功
cat("Min count:", min(GetAssayData(mouse, slot="counts")), "\n")
cat("Max count:", max(GetAssayData(mouse, slot="counts")), "\n")
cat("Dim:", dim(GetAssayData(mouse, slot="counts")), "\n")


# 校正前后对比某个已知 ambient 基因的表达
# 比如 Prm1（精子特异）在 Sertoli 细胞里应该几乎不表达
# 校正前（原始）
before <- mean(GetAssayData(mouse_old, slot="counts")["Prm1",colnames(mouse)[mouse$cell_type == "Sertoli"]])
# 校正后
after <- mean(GetAssayData(mouse, slot="counts")["Prm1",colnames(mouse)[mouse$cell_type == "Sertoli"]])
cat(sprintf("Prm1 in Sertoli — Before: %.3f  After: %.3f\n", before, after))
# 期望：after 明显低于 before


# ── Step 2: 过滤零细胞 ────────────────────────────────────────────────────────
cat("\nBefore filtering:", ncol(mouse), "cells\n")

# 方法：保留 colSums > 0 的细胞
cells_keep <- colnames(counts_mat)[Matrix::colSums(counts_mat) > 0]
mouse <- subset(mouse, cells = cells_keep)
cat("After filtering:", ncol(mouse), "cells\n")

# ── Step 3: 更新 nCount_RNA 和 nFeature_RNA（因为 counts 变了）─────────────
mouse$nCount_RNA   <- Matrix::colSums(GetAssayData(mouse, slot = "counts"))
mouse$nFeature_RNA <- Matrix::colSums(GetAssayData(mouse, slot = "counts") > 0)

cat("\nnCount_RNA range:", range(mouse$nCount_RNA), "\n")
cat("nFeature_RNA range:", range(mouse$nFeature_RNA), "\n")


mouse <- NormalizeData(mouse)
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 3000)   
mouse <- ScaleData(mouse, vars.to.regress = c("nFeature_RNA", "percent_mt", "percent_ribo"))


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

saveRDS(mouse, file = "data/testis_subset_deconX_anno.rds")
mouse <- readRDS("data/testis_subset_deconX_anno.rds")


DimPlot(mouse, reduction = "umap_harmony", group.by="bc1",  label = T, label.size = 5)
DimPlot(mouse, reduction = "umap_harmony", group.by="cell_type",  label = T, label.size = 5)

p1 <- DimPlot(mouse, reduction = "umap_harmony", group.by="RNA_snn_res.0.5",  label = T, label.size = 5)
p2 <- DimPlot(mouse, reduction = "umap_harmony", group.by="age",  label = T, label.size = 5)
p3 <- DimPlot(mouse, reduction = "umap_harmony", group.by="mouse",  label = T, label.size = 5)
p4 <- DimPlot(mouse, reduction = "umap_harmony", group.by="histone",  label = T, label.size = 5)
(p1 + p2)  / (p3+p4)

p1

DimPlot(mouse, reduction = "umap_harmony", group.by="proj",  label = T, label.size = 5)






#

mouse$seurat_clusters <- factor(
  mouse$seurat_clusters,
  levels = as.character(c(0,1,3,7,16,17,4,6,14,2,12,5,10,11,8,9,13,15))
)

FeaturePlot(mouse, features = c("nFeature_RNA", "nCount_RNA","percent_mt","percent_ribo" ), ncol = 2, reduction = "umap_harmony", order=T, label=T) 
DotPlot(mouse, features = c("nFeature_RNA", "nCount_RNA","percent_mt","percent_ribo" ), group.by = "seurat_clusters", scale = T) 
VlnPlot(mouse, features  = c("nFeature_RNA", "nCount_RNA","percent_mt","percent_ribo" ),
        group.by  = "seurat_clusters",   # 换成 "cell_type" 如果想按注释分组
        ncol      = 2,
        pt.size   = 0,                   # 不画点，cluster多时更清晰
        cols      = scales::hue_pal()(length(unique(mouse$seurat_clusters)))
) 

DimPlot(mouse, reduction = "umap_harmony", group.by="cell_type",  label = T, label.size = 5)

DimPlot(mouse, reduction = "umap_harmony", group.by="seurat_clusters", split.by = "proj", ncol =2, label = T, label.size = 5)

plot_list <- lapply(unique(mouse$proj), function(proj) {
  cells_use <- colnames(mouse)[mouse$proj == proj]
  sub_obj   <- subset(mouse, cells = cells_use)
  
  FeaturePlot(sub_obj,
              reduction = "umap_harmony",
              features  = "nFeature_RNA",
              label     = TRUE, label.size = 3,
              pt.size   = 0.4,
              min.cutoff = 10,
              max.cutoff = 1500) +
    scale_color_gradientn(
      colors = c("#1a1a2e", "#16213e", "#0f3460", "#533483",
                 "#e94560", "#f5a623", "#f8e71c"),
      limits = c(0, 1500),
      name   = "nFeature"
    ) +
    ggtitle(proj) +
    theme_classic() +
    theme(
      plot.title       = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.key.height = unit(0.8, "cm"),
      legend.text      = element_text(size = 8),
      axis.title       = element_text(size = 8),
      axis.text        = element_text(size = 7)
    )
})

(plot_list[[1]] | plot_list[[2]]) /(plot_list[[3]] | plot_list[[4]])

table(mouse$proj)


mouse.sub <- subset(mouse, subset = proj=="20250524-scPT-Natural-Aging-testis-H3K27ac-H3K27me3")
cells.use <- sample(colnames(mouse.sub), 20000)
mouse.sub <- subset(mouse, cells = cells.use)
DimPlot(mouse.sub, reduction = "umap_harmony", group.by="seurat_clusters",  label = T, label.size = 5)
saveRDS(mouse.sub, file = "data/testis_sub.rds")



# https://mp.weixin.qq.com/s/eK0T3Fwcv9Gs4LK7j2esiQ
# https://biomni.phylo.bio/projects/proj_144d36958f/tasks/sess_8d91aebb2bbe

##################### endo ################
Idents(mouse) <- "seurat_clusters"
genes = c("Flt1", "Icam2", "Ly6c1", "Pecam1") # , , , 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(19),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2




##################### Macrophage ################
Idents(mouse) <- "seurat_clusters"
genes = c("Cd74", "Csf1r", "Cd68", "Mrc1") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(15),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2


##################### Leydig ################
Idents(mouse) <- "seurat_clusters"
genes = c("Cyp11a1", "Cyp17a1", "Hsd17b3", "Hsd3b1", "Hsd3b6") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(5),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )


p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2


##################### Telocyte ################
Idents(mouse) <- "seurat_clusters"
genes = c("Dcn", "Col1a2", "Col3a1", "Col4a4", "Lama2") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(13),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2



##################### Myoid ################
Idents(mouse) <- "seurat_clusters"
genes = c("Acta2", "Myh11") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(17),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2

##################### cluster 21 ################

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(21),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2

##################### cluster 18 ################

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(18),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2


##################### Sertoli ################
Idents(mouse) <- "seurat_clusters"
genes = c("Sox9", "Aard", "Amhr2", "Basp1") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(10,12,16,20),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2



##################### cluster 6 ################
Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(6),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2



##################### cluster 11 ################
Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(11),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2


##################### cluster 22 ################
Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(22),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2


##################### Leptotene ################
Idents(mouse) <- "seurat_clusters"
genes = c("Stra8","Ccnb1ip1","Hormad1","Hormad2","Sycp2") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(4),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2


##################### PrelepSPC ################
Idents(mouse) <- "seurat_clusters"
genes = c("Stra8",  "Ccnb1ip1", "Tex101",
          "Dazl", "Rec8", "Esx1", "Setdb2") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 

##################### Zygotene ################
Idents(mouse) <- "seurat_clusters"
genes = c("Sycp1","Dmc1","Rad51","Meiob","Spo11") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 

##################### Pachytene ################
Idents(mouse) <- "seurat_clusters"
genes = c("Prdm9","Mlh1","Mlh3","Cdc25a","Hspa2") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### Diplotene ################
Idents(mouse) <- "seurat_clusters"
genes = c("Ctnna3","Sh3gl3","Camk4","Fus","Aurka","Ccnb1",
          "Gm45408","Snhg11","Eef2") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### EarlyRS ################
Idents(mouse) <- "seurat_clusters"
genes = c("Catsper3", "Izumo1", "Spaca1", "Spaca3") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 

##################### MidRS ################
Idents(mouse) <- "seurat_clusters"
genes = c("Catsper2", "Spem1", "Tssk2","Spata3",  "Ube2b", "Ddx4") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### LateRS ################
Idents(mouse) <- "seurat_clusters"
genes = c("Catsper4", "Pex5l", "Sufu", "Tnp2", "Prm2",
          "Smcp", "Odf1", "Odf2", "Spag4") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### ES ################
Idents(mouse) <- "seurat_clusters"
genes = c("Cby3", "Izumo2", "Tssk6", "Tnp2", "Reep6","Vwf", "Cd93") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### SSC ################
Idents(mouse) <- "seurat_clusters"
genes = c("Rhox10", "Mageb4", "Zbtb16", "Nanos3", "Rarg", "Afp") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 

##################### SPG ################
Idents(mouse) <- "seurat_clusters"
genes = c("Dmrt1", "Crabp1", "Msh2", "Map7d2", "Scml2", "Rhox13", "Prdm9", "Meiob") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### SPC ################
Idents(mouse) <- "seurat_clusters"
genes = c("Cdc42ep3", "Ggnbp1", "Mllt10", "Ankrd54", "Spo11") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### cluster 3 ################
Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(3),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2




##################### cluster 0 ################
Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(0),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2

########### add cell type from text ############
table(public$celltype)

Idents(mouse) <- "seurat_clusters"
new.cluster.ids <- read.table(file = "cellAnno-sub.txt", header = TRUE, sep='\t')$label
new.cluster.ids <- gsub("/", "_", new.cluster.ids)
new.cluster.ids <- gsub("\\(", "", new.cluster.ids)
new.cluster.ids <- gsub("\\)", "", new.cluster.ids)

names(new.cluster.ids) <- levels(mouse)
mouse <- RenameIdents(mouse, new.cluster.ids)
mouse$cell_type <- mouse@active.ident
mouse$cell_type_cluster <- paste0(mouse@active.ident, mouse$seurat_clusters)


p1 <- DimPlot(mouse, reduction = "umap_harmony", group.by="cell_type", label = T, label.size = 6)
p1


meta <- mouse@meta.data
meta$cell_id <- rownames(meta)
table(meta$cell_type)

write.csv(meta, 
          file = "testis-sub-meta.csv", 
          row.names = FALSE)




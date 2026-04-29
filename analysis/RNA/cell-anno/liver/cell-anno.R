
rm(list=ls()) ; gc();

set.seed(42)

library(Seurat)
library(ggplot2)
library(readr)
library(grid)
library(gridExtra)
library(edgeR)
#library(FactoMineR)
library(dplyr) 
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(harmony)


setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/cell-anno/liver")

public <- readRDS(
  "/storage/zhangyanxiaoLab/qihongjian/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/caojunyue/data/liver_caoyuejun_50k.rds"
)
Idents(public) <- "celltype"

load("/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/results/liver-merged/mouse_after_cellanno.RData")

p1 <- DimPlot(mouse, reduction = "umap", group.by="RNA_snn_res.0.5",  label = T, label.size = 5)
p2 <- DimPlot(mouse, reduction = "umap", group.by="cell_type",  label = T, label.size = 5)
p1 + p2 


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

p1 <- DimPlot(mouse, reduction = "umap_harmony", group.by="RNA_snn_res.0.5",  label = T, label.size = 5)
p2 <- DimPlot(mouse, reduction = "umap_harmony", group.by="age",  label = T, label.size = 5)
p3 <- DimPlot(mouse, reduction = "umap_harmony", group.by="mouse",  label = T, label.size = 5)
p4 <- DimPlot(mouse, reduction = "umap_harmony", group.by="histone",  label = T, label.size = 5)
(p1 + p2)  / (p3+p4)
                                

saveRDS(mouse, file = "data/liver_cell_anno.rds")

FeaturePlot(mouse, features = c("nFeature_RNA", "nCount_RNA","percent_mt","percent_ribo" ), ncol = 2, reduction = "umap_harmony", order=T, label=T) 
VlnPlot(mouse, features  = c("nFeature_RNA", "nCount_RNA","percent_mt","percent_ribo" ),
        group.by  = "seurat_clusters",   # 换成 "cell_type" 如果想按注释分组
        ncol      = 2,
        pt.size   = 0,                   # 不画点，cluster多时更清晰
        cols      = scales::hue_pal()(length(unique(mouse$seurat_clusters)))
) 


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
FeaturePlot(mouse, features = unique(markers$gene)[1:9], ncol = 3, reduction = "umap_harmony", order=T, label=T) 


#####################  Hepatocytes ################
Hepatocytes = c('Alb',"Ttr", "Tat")
mouse <- AddModuleScore( object = mouse, features = Hepatocytes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c('Alb',"Ttr", "Tat", "cluster_marker1") , reduction = "umap_harmony", 
            label = TRUE, order=T)
DotPlot(mouse, features = c('Alb',"Ttr", "Tat", "cluster_marker1"),group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(0,1,4,6,13,17),
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
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in liver\n source:Caoyunjue-science"))
p1 + p2


#####################  Hepatic_stellate_cells ################
Idents(mouse) <- "seurat_clusters"
genes = c("Reln", "Rbp1") # Hepatic_stellate_cells
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(3, 12),
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
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in liver\n source:Caoyunjue-science"))
p1 + p2


#####################  endo ################
Idents(mouse) <- "seurat_clusters"
genes = c("Clec4g", "Kdr", "Aqp1") # endo
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(2, 7, 14),
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
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in liver\n source:Caoyunjue-science"))
p1 + p2

#####################  Kupffer_cells ################
Idents(mouse) <- "seurat_clusters"
genes = c("Clec4f", "C1qa", "Csf1r") # Kupffer_cells
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
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in liver\n source:Caoyunjue-science"))
p1 + p2




#####################  Cholangiocytes ################
Idents(mouse) <- "seurat_clusters"
genes = c("Krt19", "Epcam", "Sox9") # Cholangiocytes
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
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in liver\n source:Caoyunjue-science"))
p1 + p2


#####################  Mesothelial ################
Idents(mouse) <- "seurat_clusters"
genes = c("Msln", "Upk3b", "Wt1", "Krt19") # Mesothelial
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(16),
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
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = T, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in liver\n source:Caoyunjue-science"))
p1 + p2


#####################  Bcell ################
Idents(mouse) <- "seurat_clusters"
genes = c("Ighm", "Cd79b") # Bcell
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(9),
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
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in liver\n source:Caoyunjue-science"))
p1 + p2


#####################  Tcell ################
Idents(mouse) <- "seurat_clusters"
genes = c("Trac", "Trbc1", "Trbc2") # Tcell
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(8),
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
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in liver\n source:Caoyunjue-science"))
p1 + p2


#####################  cluster10 ################


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(10),
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
FeaturePlot(mouse, features = unique(markers$gene)[1:9], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in liver\n source:Caoyunjue-science"))
p1 + p2


############### add cell type from text ############### 

Idents(mouse) <- "RNA_snn_res.0.5"

new.cluster.ids <- read.table(file = "cellAnno.txt", header = TRUE, sep='\t')$label
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

write.csv(meta, 
          file = "~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/metadata/data/liver-final-cell-type.csv", 
          row.names = FALSE)

saveRDS(mouse, file = "data/liver_cell_anno.rds")


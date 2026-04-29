
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


setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/cell-anno/heart")

public <- readRDS(
  "/storage/zhangyanxiaoLab/qihongjian/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/caojunyue/data/heart_caoyuejun_50k.rds"
)
Idents(public) <- "celltype"
DimPlot(public, reduction = "umap", group.by="celltype",  label = T, label.size = 5)

load("/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/results/heart-merged/mouse_after_cellanno.RData")

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
                                
p1

saveRDS(mouse, file = "data/heart_cell_anno.rds")
#mouse <- readRDS("data/lung_cell_anno.rds")



FeaturePlot(mouse, features = c("nFeature_RNA", "nCount_RNA","percent_mt","percent_ribo" ), ncol = 2, reduction = "umap_harmony", order=T, label=T) 
VlnPlot(mouse, features  = c("nFeature_RNA", "nCount_RNA","percent_mt","percent_ribo" ),
        group.by  = "seurat_clusters",   # 换成 "cell_type" 如果想按注释分组
        ncol      = 2,
        pt.size   = 0,                   # 不画点，cluster多时更清晰
        cols      = scales::hue_pal()(length(unique(mouse$seurat_clusters)))
) 

# https://mp.weixin.qq.com/s/eK0T3Fwcv9Gs4LK7j2esiQ
# https://biomni.phylo.bio/projects/proj_144d36958f/tasks/sess_8d91aebb2bbe

##################### Cardiomyocyte ################
Idents(mouse) <- "seurat_clusters"
genes = c("Tnnt2", "Tnni3", "Myh6", "Myh7", "Ryr2") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(1,8),
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




##################### Tcell ################
Idents(mouse) <- "seurat_clusters"
genes = c("Cd3e", "Cd3d", "Cd8a", "Cd4") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


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
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:4]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2




##################### Bcell ################ Cd19, Ms4a1(Cd20), Cd79a, Cd79b
Idents(mouse) <- "seurat_clusters"
genes = c("Cd19", "Cd79a", "Cd79b") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(12),
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




##################### Epicardial ################ Wt1, Tbx18, Tcf21, Krt8, Krt18
Idents(mouse) <- "seurat_clusters"
genes = c("Wt1", "Tbx18", "Krt8", "Krt18") # 
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
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2




##################### Mural ################ Pdgfrb, Rgs5, Cspg4, Mcam
Idents(mouse) <- "seurat_clusters"
genes = c("Pdgfrb", "Rgs5", "Cspg4", "Mcam") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(7),
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




##################### Fibroblast ################ Pdgfra, Col1a1, Col1a2, Dcn, Lum, Tcf21
Idents(mouse) <- "seurat_clusters"
genes = c("Pdgfra", "Col1a1", "Col1a2", "Dcn", "Lum") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(0,6,14,18),
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


##################### Adipocyte ################ Plin1, Adipoq, Fabp4, Lpl
Idents(mouse) <- "seurat_clusters"
genes = c("Plin1", "Adipoq", "Fabp4", "Lpl") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


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




##################### neural ################ Plp1, Mpz, S100b, Sox10, Ncam1
Idents(mouse) <- "seurat_clusters"
genes = c("Plp1", "Mpz", "S100b", "Sox10", "Ncam1") # 
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




##################### Endocardial_endo ################ 
Idents(mouse) <- "seurat_clusters"
genes = c("Emcn", "Npr3", "Nfatc1", "Tbx20") # 
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

#####################  Capillary EC ################ 
Idents(mouse) <- "seurat_clusters"
genes = c("Cd36", "Aplnr", "Rgcc", "Slc16a1", "Gpihbp1") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(2,5),
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

#####################  Lymphatic Endothelial Cell ################ 
Idents(mouse) <- "seurat_clusters"
genes = c("Prox1", "Pdpn", "Flt4", "Lyve1", "Mmrn1") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


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

##################### Macrophage ################  Adgre1(F4/80), Cd68, Csf1r, Cx3cr1, Mrc1(Cd206)
Idents(mouse) <- "seurat_clusters" 
genes = c("Adgre1", "Cd68", "Csf1r", "Cx3cr1", "Mrc1") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


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


#####################  cluster13 ################ 

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

#####################  cluster19 ################ 

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


#####################  Neutrophil ################ 

Idents(mouse) <- "seurat_clusters" 
genes = c("Ly6g", "S100a8", "S100a9", "Csf3r", "Retnlg") # 
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


#####################  cluster 21 doublets ################ 
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


#####################  cluster 16 doublets ################ 
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


#####################  cluster 20 doublets ################ 

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(20),
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

#####################  cluster 15 doublets ################ 

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


########### add cell type from text ############
table(public$celltype)

Idents(mouse) <- "seurat_clusters"
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
table(meta$cell_type)

write.csv(meta, 
          file = "~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/metadata/data/heart-final-cell-type.csv", 
          row.names = FALSE)

saveRDS(mouse, file = "data/heart_cell_anno.rds")


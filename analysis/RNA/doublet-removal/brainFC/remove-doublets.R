
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
rds_path <- "/storage/zhangyanxiaoLab/yiguo/projects/Multi-batch_RNA_eachregion_cocluster_anno/data_FC/FC_three_batches_annotated.rds"
mouse <- readRDS(rds_path)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4 )
p2 <- DimPlot( mouse, reduction = "umap", group.by = "cell_anno", label = TRUE, label.size = 4 )
p3 <- DimPlot( mouse, reduction = "umap", group.by = "main_anno", label = TRUE, label.size = 4 )

p1 + (p2/p3)
table(mouse$seurat_clusters, mouse$cell_anno)

sce <- as.SingleCellExperiment(mouse)
sce <- scDblFinder(sce, samples = "SimpleID", clusters = "seurat_clusters", dbr=0.1, BPPARAM=MulticoreParam(16))

tmp <- as.Seurat(sce)
tmp_meta <- tmp@meta.data
DimPlot(tmp, group.by = "scDblFinder.cluster", reduction = "UMAP", label = TRUE,  label.size = 5 ) +ggtitle("seurat cluster")
clusterStickiness(sce)

table(sce$scDblFinder.score)#查看双细胞得分
table(sce$scDblFinder.class)


mouse_scDblFinder <- as.Seurat(sce)
metadata(sce)$scDblFinder.stats
scDblFinder_meta <- mouse_scDblFinder@meta.data
plotDoubletMap(sce)

score_vector  <- setNames(scDblFinder_meta$scDblFinder.score,  rownames(scDblFinder_meta))
class_vector  <- setNames(scDblFinder_meta$scDblFinder.class,  rownames(scDblFinder_meta))
cluster_vector <- setNames(scDblFinder_meta$scDblFinder.cluster, rownames(scDblFinder_meta))
origin_vector  <- setNames(scDblFinder_meta$scDblFinder.mostLikelyOrigin, rownames(scDblFinder_meta))

# 2. 将所有列加入 Seurat metadata
mouse@meta.data$scDblFinder.score            <- score_vector[colnames(mouse)]
mouse@meta.data$scDblFinder.class            <- class_vector[colnames(mouse)]
mouse@meta.data$scDblFinder.cluster          <- cluster_vector[colnames(mouse)]
mouse@meta.data$scDblFinder.mostLikelyOrigin <- origin_vector[colnames(mouse)]

# 3. 检查
head(mouse@meta.data[, c(
  "scDblFinder.score",
  "scDblFinder.class",
  "scDblFinder.cluster",
  "scDblFinder.mostLikelyOrigin"
)])


meta_df <- mouse@meta.data
meta_df$cell_id <- rownames(meta_df)
write.csv(meta_df, file = "data/mouse_metadata_first_round.csv", quote = FALSE)





#########

p1 <- DimPlot(mouse, reduction = "umap", group.by="seurat_clusters",  label = T, label.size = 5) + ggtitle("brainFC cluster")
p2 <- DimPlot(mouse, reduction = "umap", group.by="scDblFinder.class",  label = T, label.size = 5) + ggtitle("brainFC scDblFinder.class")
p2 <- FeaturePlot(mouse, features = "scDblFinder.score", reduction = "umap", order=T, label=T)


# 查看每个cluster中双细胞分类的分布
cluster_doublet_table <- table(mouse$seurat_clusters, mouse$scDblFinder.class)
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

Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_ribo"), ncol = 2, reduction = "umap", order=T, label=T) 
VlnPlot( mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_ribo"), ncol = 1,pt.size = 0 )

FeaturePlot(mouse, features = all_markers, reduction = "umap", ncol = 4, order=T, label=T)
DotPlot(mouse, features = all_markers) + RotatedAxis()

# potential-doublets-clusters
DotPlot( subset( mouse, idents = c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40) ), features = all_markers, group.by = "seurat_clusters" ) + RotatedAxis()


############# cluster 11 ##########
cluster_id <- "11"
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                      "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters == cluster_id, cluster_id,
  ifelse(mouse$seurat_clusters %in% c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40), "exclude", "others")
)
Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(mouse, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
            cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2

############# cluster 39 ##########
cluster_id <- "39"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters == cluster_id, cluster_id,
  ifelse(mouse$seurat_clusters %in% c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40), "exclude", "others")
)
Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(mouse, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2

############# cluster 36 ##########
cluster_id <- "36"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

############# cluster 29 ##########
cluster_id <- "29"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters == cluster_id, cluster_id,
  ifelse(mouse$seurat_clusters %in% c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40), "exclude", "others")
)
Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(mouse, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2


############# cluster 37 ##########
cluster_id <- "37"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

############# cluster 40 ##########
cluster_id <- "40"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters == cluster_id, cluster_id,
  ifelse(mouse$seurat_clusters %in% c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40), "exclude", "others")
)
Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(mouse, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2

############# cluster 40 ##########
cluster_id <- "23"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters == cluster_id, cluster_id,
  ifelse(mouse$seurat_clusters %in% c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40), "exclude", "others")
)
Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(mouse, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2

############# cluster 38 ##########
cluster_id <- "38"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

############# cluster 38 ##########
cluster_id <- "24"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters == cluster_id, cluster_id,
  ifelse(mouse$seurat_clusters %in% c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40), "exclude", "others")
)
Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(mouse, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2

############# cluster 22 ##########
cluster_id <- "22"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters == cluster_id, cluster_id,
  ifelse(mouse$seurat_clusters %in% c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40), "exclude", "others")
)
Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(mouse, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2

############# cluster 35 ##########
cluster_id <- "35"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters == cluster_id, cluster_id,
  ifelse(mouse$seurat_clusters %in% c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40), "exclude", "others")
)
Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(mouse, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2
############# cluster 31 ##########
cluster_id <- "31"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters == cluster_id, cluster_id,
  ifelse(mouse$seurat_clusters %in% c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40), "exclude", "others")
)
Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(mouse, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2

############# cluster 32 ##########
cluster_id <- "32"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters == cluster_id, cluster_id,
  ifelse(mouse$seurat_clusters %in% c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40), "exclude", "others")
)
Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(mouse, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2

############# cluster 30 ##########
cluster_id <- "30"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

############# cluster 34 ##########
cluster_id <- "34"
table(subset(mouse, seurat_clusters == cluster_id)$cell_anno)
FeaturePlot(subset(mouse, seurat_clusters == cluster_id), features = c("Slc17a7","Gad2",	"Apoe","C1qb",
                                                                       "Mog", "Pdgfra", "Tie1","Foxd1"
), reduction = "umap", ncol = 4, order=T)

p1 <- DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( mouse, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( mouse, idents = cluster_id ), features = all_markers, group.by = "seurat_clusters" ) + coord_flip()
p1 + p2

Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters == cluster_id, cluster_id,
  ifelse(mouse$seurat_clusters %in% c(11, 22, 23, 24, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40), "exclude", "others")
)
Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  subset(mouse, one_vs_others != "exclude"),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T,
                  cols = c("grey90", "orange", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score for cluster", cluster_id, " in brainFC\n source:CR"))
p1 + p2




mouse$double_manual <- "singlet"
doublet_clusters <-c(22, 32, 36, 37, 38)
mouse$double_manual[mouse$seurat_clusters %in% doublet_clusters] <- "doublet"
DimPlot(mouse, reduction = "umap", group.by="double_manual",  label = T, label.size = 5) + ggtitle("brainFC all doublets")



saveRDS(mouse, file = "data/mouse_all.rds")










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
library(scDblFinder)
library(BiocParallel)


setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/doublet/heart")





rds_path <- "~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/cell-anno/heart/data/heart_cell_anno.rds"
mouse <- readRDS(rds_path)

p1 <- DimPlot( mouse, reduction = "umap_harmony", group.by = "seurat_clusters", label = TRUE, label.size = 4 )
p2 <- DimPlot( mouse, reduction = "umap_harmony", group.by = "cell_type", label = TRUE, label.size = 4 )
p1 + p2

mouse <- subset(mouse, subset = cell_type != "low_qual")
mouse$SimpleID <- mouse$orig.ident

table(mouse$seurat_clusters, mouse$cell_type)


sce <- as.SingleCellExperiment(mouse)
sce <- scDblFinder(sce, samples = "SimpleID", clusters = "seurat_clusters", dbr=0.1, BPPARAM=MulticoreParam(16))

tmp <- as.Seurat(sce)
tmp_meta <- tmp@meta.data
DimPlot(tmp, group.by = "scDblFinder.cluster", reduction = "UMAP_HARMONY", label = TRUE,  label.size = 5 ) +ggtitle("seurat cluster")
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

p1 <- DimPlot(mouse, reduction = "umap_harmony", group.by="seurat_clusters",  label = T, label.size = 5) + ggtitle("heart cluster")
DimPlot(mouse, reduction = "umap_harmony", group.by="scDblFinder.class",  label = T, label.size = 5) + ggtitle("heart scDblFinder.class")
p2 <- FeaturePlot(mouse, features = "scDblFinder.score", reduction = "umap_harmony", order=T, label=T)

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

Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), ncol = 2, 
            reduction = "umap_harmony", order=T, label=T) 
VlnPlot( mouse, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), ncol = 1,pt.size = 0 )




saveRDS(mouse, file = "data/mouse_all.rds")







# remove doublets cluster and predicted

mouse$doublet <- ifelse(
  mouse$scDblFinder.class == "doublet" |
    mouse$seurat_clusters %in% c(14,16,18,19,20,21,22),
  "doublet",
  "singlet"
)

DimPlot( mouse, reduction = "umap_harmony", group.by = "doublet", label = TRUE, label.size = 4 )


##################### reclustering ##########
mouse
mouse <- subset(mouse, subset = doublet == "singlet")
mouse

DefaultAssay(object = mouse) <- "RNA"
mouse = NormalizeData(mouse, normalization.method = "LogNormalize", scale.factor = 10000)
mouse = FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 3000)  
mouse = ScaleData(mouse, vars.to.regress = c("nFeature_RNA", "percent_mt"))

mouse <- RunPCA(mouse, features = VariableFeatures(object = mouse))
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

DimPlot(mouse, reduction = "umap_harmony", group.by="cell_type",  label = T, label.size = 5)


table(mouse$seurat_clusters,mouse$cell_type)

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
          file = "~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/metadata/data/heart-final-cell-type.csv", 
          row.names = FALSE)

saveRDS(mouse, file = "data/heart_cell_remode_doublets.rds")




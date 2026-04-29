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
mouse <- readRDS("data/mouse_all.rds")
DimPlot( mouse, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4 )

meta <- mouse@meta.data[mouse$seurat_clusters %in% c(11, 27,28,29,30), ]
meta$barcode <- sub("^[^:]+:", "", meta$full_barcode_with_batch)
meta$cell_id <- paste0(meta$SimpleID, ":", meta$barcode)
all_meta <- meta[,c("SimpleID", "barcode", "cell_id")]
all_meta$double_manual <- "double_manual_all_cells"

meta <- mouse@meta.data[mouse$scDblFinder.class == "doublet", ]
meta$barcode <- sub("^[^:]+:", "", meta$full_barcode_with_batch)
meta$cell_id <- paste0(meta$SimpleID, ":", meta$barcode)
all_meta_pred <- meta[,c("SimpleID", "barcode", "cell_id")]
all_meta_pred$double_manual <- "double_pred_all_cells"



gila <- readRDS("data/mouse_InN.rds")
DimPlot( gila, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4 )
meta <- gila@meta.data[gila$seurat_clusters %in% c(6), ]
meta$barcode <- sub("^[^:]+:", "", meta$full_barcode_with_batch)
meta$cell_id <- paste0(meta$SimpleID, ":", meta$barcode)
ihn_meta <- meta[,c("SimpleID", "barcode", "cell_id")]
ihn_meta$double_manual <- "double_manual_InN_cells"


gila <- readRDS("data/mouse_Glia_high_resolution.rds")
DimPlot( gila, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4 )

meta <- gila@meta.data[gila$seurat_clusters %in% c(19,20,22,23,25,5,12,13,15,21,24), ] 
meta$barcode <- sub("^[^:]+:", "", meta$full_barcode_with_batch)
meta$cell_id <- paste0(meta$SimpleID, ":", meta$barcode)
glia_meta <- meta[,c("SimpleID", "barcode", "cell_id", "double_manual")]
glia_meta$double_manual <- "double_manual_Glia_cells"

# 合并
final_meta <- bind_rows(all_meta, all_meta_pred, ihn_meta, glia_meta)
final_meta <- distinct(final_meta, cell_id, .keep_all = TRUE)
write.csv(final_meta, file = "brainCB_doublet_cells_combined.csv", row.names = FALSE)

# 查看
dim(final_meta)
head(final_meta)

mouse$barcode <- sub("^[^:]+:", "", mouse$full_barcode_with_batch)
mouse$cell_id <- paste0(mouse$SimpleID, ":", mouse$barcode)
mouse$in_doublet_list <- ifelse(mouse$cell_id %in% final_meta$cell_id,
                                "doublet",
                                "keep")

table(mouse$in_doublet_list)
DimPlot(mouse, group.by = "in_doublet_list") + ggtitle("final doublets identified")




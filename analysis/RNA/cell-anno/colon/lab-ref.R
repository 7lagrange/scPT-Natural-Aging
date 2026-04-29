
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


setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/cell-anno/colon")


load("/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/results/colon-merged/mouse_after_cellanno.RData")

tissue <- "colon"

final_cell_type <- paste0(
  "/storage/zhangyanxiaoLab/qihongjian/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/metadata/data/",
  "xcolon-final-cell-type.csv"
)

if (file.exists(final_cell_type)) {
  print(tissue)
  mouse$cell_id <- rownames(mouse@meta.data)
  cell_anno <- read.csv(final_cell_type, stringsAsFactors = FALSE)
  mouse$cell_type <- NA
  idx <- match(mouse$cell_id, cell_anno$cell_id)
  mouse$cell_type <- cell_anno$cell_type[idx]
  sum(is.na(mouse$cell_type))
  mouse <- subset(mouse, subset = cell_type %in% unique(cell_anno$cell_type))
}


DimPlot(mouse, group.by = "cell_type",label = T,repel = T)

sce <- mouse
sce$celltype <- sce$cell_type

# value between 0 to 1137, so it looks like not log normalized 
av <-AverageExpression(sce,
                       group.by = "celltype",
                       assays = "RNA") 
av=av[[1]]
head(av)


cg=names(tail(sort(apply(av, 1, sd)),1000))
pheatmap::pheatmap(cor(av[cg,]))

Ref=av
ref_sce=SingleCellExperiment::SingleCellExperiment(assays=list(counts=Ref))
ref_sce=scater::logNormCounts(ref_sce)


logcounts(ref_sce)[1:4,1:4]
colData(ref_sce)$Type=colnames(Ref)
table(colnames(Ref))
ref_sce
save(ref_sce,file = '~/github/zhanglab-code/scRNA_scripts/data/mouse_colon.lab_data_ref.Rdata') 

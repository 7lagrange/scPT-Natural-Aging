
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

library(dplyr)


df1 <- readRDS("data/SPC_RS.combination.rds")
DimPlot(df1, reduction = "umap", group.by="cluster",  label = T, label.size = 5)

df2 <- readRDS("data/SPC.diff_PRO.rds")
df3 <- readRDS("data/RS.diff_PRO.rds")


DimPlot(df2, reduction = "umap", group.by="cluster",  label = T, label.size = 5)
DimPlot(df3, reduction = "umap", group.by="cluster",  label = T, label.size = 5)

public <- merge(
  df1,
  y = c(df2, df3),
  add.cell.ids = c("SPG", "SPC", "RS"),
  project = "combined"
)
DimPlot(public, reduction = "umap", group.by="cluster",  label = T, label.size = 5)





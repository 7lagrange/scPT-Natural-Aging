
rm(list=ls()) ; gc();

set.seed(42)

library(Seurat)
library(ggplot2)
library(readr)
library(grid)
library(gridExtra)
library(edgeR)

library(dplyr) 
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(tibble)

setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/analysis/H3K9me3/other-marks-zhuojie-regions")

H3K9me3_logfc <- read.csv("../H3K9me3-zhuojie-regions/data/H3K9me3_logfc_heatmap_matrix.csv",check.names = FALSE)
ct_annotation <- read_csv("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/analysis/H3K9me3/H3K9me3-zhuojie-regions/data/ct_annotation_annotated.csv")
ct_annotation <- as.data.frame(ct_annotation)
rownames(ct_annotation) <- ct_annotation$cell_type


antibody <- 'H3K9me3'
antibody <- 'H3K27me3'
antibody <- 'H3K36me3'
antibody <- 'H3K27ac'
antibody <- 'H3K4me1'
antibody <- 'H3K4me3'

tissue_colors <- c(
  "kidney"   = "#1f77b4",  # 蓝
  "liver"    = "#ff7f0e",  # 橙
  "heart"    = "#2ca02c",  # 绿
  "brainFC"  = "#d62728",  # 红
  "brainHip" = "#9467bd",  # 紫
  "brainCB"  = "#8c564b",  # 棕
  "BAT"      = "#e377c2",  # 粉
  "testis"   = "#7f7f7f",  # 灰
  "muscle"   = "#bcbd22",  # 黄绿
  "colon"    = "#17becf",  # 青
  "spleen"   = "#ff9896",  # 浅红
  "lung"     = "#98df8a"   # 浅绿
)


cluster_colors <- c(
  "1"      = "#d62728",  # 红
  "2"      = "#1f77b4",  # 蓝
  "3"      = "#ff7f0e",  # 橙
  "4"      = "#2ca02c",  # 绿
  "others" = "#7f7f7f"   # 灰
)

lineage_colors <- c(
  "ExN"         = "#17becf",  # 青
  "Epithelial"  = "#1f77b4",  # 蓝
  "Muscle"      = "#9467bd",  # 紫
  "Stromal"     = "#ff7f0e",  # 橙
  "Germ"    = "#7f7f7f",  # 灰
  "InN"         = "#bcbd22",  # 黄绿
  
  "Endothelial" = "#2ca02c",  # 绿
  "Immune"      = "#d62728",  # 红
  "Glia"        = "#8c564b",  # 棕
  "Neuronal"     = "#e377c2"  # 粉
  
)



########## read in depth info #########

DE_edgeR <- read.csv(paste0("../H3K9me3-zhuojie-regions/data/", antibody, "_edgeR_27M_over_3M_zhuojie.csv"))
DE_edgeR$cell_type <- paste0(DE_edgeR$tissue, "-", DE_edgeR$cell_type)
df_anno <- read.csv("/storage/zhangyanxiaoLab/qihongjian/github/zhanglab-code/projects/paired_seq_tag/code/zhuojie/kmeans_annotation.csv")
anno_split <- strsplit(df_anno$X, "[:-]")
df_anno$chrom <- sapply(anno_split, "[", 1)
df_anno$start <- as.numeric(sapply(anno_split, "[", 2))
df_anno$end   <- as.numeric(sapply(anno_split, "[", 3))

df_merged <- full_join(DE_edgeR, df_anno, by = c("chrom", "start", "end")) %>% mutate(cluster = ifelse(is.na(cluster), "others", cluster))

df_heat <- df_merged %>% dplyr::select(region = gene, cell_type, logFC = avg_logFC, cluster)

heatmap_matrix <- df_heat %>%
  dplyr::select(-cluster) %>%   # 丢掉 cluster 列
  pivot_wider(names_from = cell_type, values_from = logFC, values_fill = 0) %>%
  column_to_rownames(var = "region") %>%
  as.matrix()

heatmap_matrix <- heatmap_matrix[H3K9me3_logfc[[1]], colnames(H3K9me3_logfc)[-1]]
dim(heatmap_matrix)


row_annotation <- df_merged %>%
  dplyr::select(region = gene, cluster) %>%
  distinct() %>%
  column_to_rownames("region") %>%
  rownames_to_column("region") %>%
  mutate(
    chrom = str_extract(region, "^chr[^-]+")   # 提取 chr 开头到第一个 "-" 之前
  ) %>%
  column_to_rownames("region")

cluster_levels <- sort(unique(row_annotation$cluster))
cluster_colors <- setNames(c("red", "green", "blue", "purple", "black"), cluster_levels)

ct_annotation <- as.data.frame(ct_annotation[c("tissue", "lineage")])

pheatmap(heatmap_matrix,
         breaks = seq(-3, 3, length.out = 101),
         annotation_col = ct_annotation,  # 使用修改后的注释
         annotation_row = row_annotation,
         annotation_colors = list(
           tissue = tissue_colors,
           lineage = lineage_colors,  # 使用新的颜色
           cluster = cluster_colors
         ),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = FALSE,
         show_colnames = TRUE,
         scale = "none",
         filename = paste0("data/", antibody, "-zhuojie-regions-logfc.pdf"),
         width = 32,
         height = 20,
         angle_col = 90,
         fontsize = 20,        # 全局字体
         fontsize_col = 20,
         main = paste0(antibody, " aging logFC of large ", dim(heatmap_matrix)[1], " H3K9me3 peaks(>200kb) from zhuojie")
)



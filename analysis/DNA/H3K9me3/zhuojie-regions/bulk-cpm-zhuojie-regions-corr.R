
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

setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/analysis/H3K9me3/H3K9me3-cell-type-corr")


########### this script compare H3K9me3 3M mice samples correlation in zhuojie regions


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
  "Epithelial"  = "#1f77b4",  # 蓝
  "Stromal"     = "#ff7f0e",  # 橙
  "Endothelial" = "#2ca02c",  # 绿
  "Immune"      = "#d62728",  # 红
  "Muscle"      = "#9467bd",  # 紫
  "Adipose"     = "#e377c2",  # 粉
  "Glia"        = "#8c564b",  # 棕
  "ExN"         = "#17becf",  # 青
  "InN"         = "#bcbd22",  # 黄绿
  "Germline"    = "#7f7f7f",  # 灰
  "Stemcell"    = "black"  # 灰
)


########## read in depth info #########
antibody <- 'H3K9me3'

tissues <- c( "brainFC", "brainHip", "brainCB", "liver", "lung", "heart", 
             "kidney", "BAT", "muscle", "testis", "colon", "spleen" )

ct_annotation <- read_csv("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/analysis/H3K9me3/H3K9me3-zhuojie-regions/data/ct_annotation_annotated.csv")

DE_edgeR <- read.csv("data/H3K9me3_bulk_cpm_zhuojie.csv", check.names = FALSE)
rownames(DE_edgeR) <- DE_edgeR$gene

min_depth = 2

exclude_cts <- ct_annotation[
  ct_annotation$total_reads < min_depth | ct_annotation$lineage %in% c("Doublet", "Germ"),
]$cell_type
DE_edgeR <- DE_edgeR[, !colnames(DE_edgeR) %in% exclude_cts ]
keep_cts <- setdiff(colnames(DE_edgeR), "gene") 

DE_edgeR$gene <- rownames(DE_edgeR) 
anno_split <- strsplit(DE_edgeR$gene, "[:-]")
DE_edgeR$chrom <- sapply(anno_split, "[", 1)
DE_edgeR$start <- as.numeric(sapply(anno_split, "[", 2))
DE_edgeR$end   <- as.numeric(sapply(anno_split, "[", 3))
DE_edgeR$region_size <- DE_edgeR$end - DE_edgeR$start
DE_edgeR <- DE_edgeR[DE_edgeR$region_size > 200000,]
length(unique(DE_edgeR$gene))



head(DE_edgeR)

cor_mat <- cor(
  DE_edgeR[, keep_cts],
  use = "pairwise.complete.obs",
  method = "pearson"
)

rowSums(cor_mat)

head(cor_mat)


cell_anno <- data.frame(sample = colnames(cor_mat)) %>%
  tidyr::separate(
    sample,
    into = c("tissue", "cell_type"),
    sep = "-",
    remove = FALSE
  ) 


cell_anno <- cell_anno %>%
  left_join(
    ct_annotation %>% 
      dplyr::select(cell_type, lineage),
    by = c("sample" = "cell_type")
  )
rownames(cell_anno) <- cell_anno$sample

ph <- pheatmap(
  cor_mat,
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  annotation_col = cell_anno[c("tissue", "lineage")],
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  border_color = NA,
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 8,
  show_rownames = FALSE,  
  main = paste0("Cell-type pearson correlation of H3K9me3 cpm across 1704 peaks"),
  width = 36,   
  height = 36,   
  filename = "data/H3K9me3_zhuojiepeaks_cell_type_corr.pdf"
)

ph


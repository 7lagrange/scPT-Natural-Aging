
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
library(stringi)


setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/analysis/H3K9me3/H3K9me3-zhuojie-regions")

antibody <- 'H3K9me3'

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
ct_annotation <- read_csv("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/analysis/H3K9me3/H3K9me3-zhuojie-regions/data/ct_annotation_annotated.csv")
ct_annotation <- as.data.frame(ct_annotation)
rownames(ct_annotation) <- ct_annotation$cell_type

tissues <- c( "brainFC", "brainHip", "brainCB", "liver", "lung", "heart", "kidney", "BAT", "muscle", "testis", "colon", "spleen" )

DE_edgeR <- read.csv("data/H3K9me3_edgeR_27M_over_3M_zhuojie.csv")
DE_edgeR$cell_type <- paste0(DE_edgeR$tissue, "-", DE_edgeR$cell_type)

######## based the above plot, exclude samples that with few depth
min_depth = 2
exclude_cts <- ct_annotation[
  ct_annotation$total_reads < min_depth | ct_annotation$lineage %in% c("Doublet"),
]$cell_type
DE_edgeR <- DE_edgeR[!DE_edgeR$cell_type %in% exclude_cts, ]
DE_edgeR$logFC <- DE_edgeR$avg_logFC

###################### large region with zhuojie cluster ###########
df_anno <- read.csv("/storage/zhangyanxiaoLab/qihongjian/github/zhanglab-code/projects/paired_seq_tag/code/zhuojie/kmeans_annotation.csv")
anno_split <- strsplit(df_anno$X, "[:-]")
df_anno$chrom <- sapply(anno_split, "[", 1)
df_anno$start <- as.numeric(sapply(anno_split, "[", 2))
df_anno$end   <- as.numeric(sapply(anno_split, "[", 3))

df_merged <- full_join(DE_edgeR, df_anno, by = c("chrom", "start", "end")) %>% mutate(cluster = ifelse(is.na(cluster), "others", cluster))
df_merged$region_size <- df_merged$end - df_merged$start

# Assuming df_merged is your dataframe
library(dplyr)


df_c1_avglogfc <- df_merged %>%
  filter(cluster == 1) %>%  # Filter for cluster 1
  group_by(tissue, cell_type) %>%  # Group by tissue and cell_type
  dplyr::summarise(
    avg_logFC = mean(logFC, na.rm = TRUE),
    num_up_peaks   = sum(p_val_adj < 0.05 & logFC >  log2(1.1), na.rm = TRUE),
    num_down_peaks = sum(p_val_adj < 0.05 & logFC < -log2(1.1), na.rm = TRUE),
    .groups = "drop"
  ) %>%  # Calculate average logFC
  arrange(desc(-avg_logFC))  # Order by avg_logFC (large to small)

write.csv(df_c1_avglogfc, paste0("data/zhuojie_c1_avglogfc.csv"), quote = FALSE, row.names = FALSE)

km1_logfc_order <- df_c1_avglogfc$cell_type


df_merged <- df_merged[df_merged$region_size > 200000,]
length(unique(df_merged$gene))

df_region_count <- df_merged %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    n_sig_up   = sum(p_val_adj <= 0.05 & avg_logFC > 0, na.rm = TRUE),
    n_sig_down = sum(p_val_adj <= 0.05 & avg_logFC < 0, na.rm = TRUE),
    n_sig_total = n_sig_up + n_sig_down
  ) %>%
  arrange(desc(n_sig_total))

region_negative <- df_region_count %>% filter(n_sig_down >= 10) %>%  filter(n_sig_up < 10)  %>% pull(gene)
region_positive <- df_region_count %>% filter(n_sig_up >= 10) %>%  filter(n_sig_down < 10)  %>% pull(gene)
region_mix <- df_region_count %>% filter(n_sig_up < 10) %>%  filter(n_sig_down < 10)  %>% pull(gene)
regions <- unique(df_merged$gene)

selected_regions <- region_negative
selected_regions <- region_positive
selected_regions <- region_mix
selected_regions <- regions

df_heat <- df_merged[df_merged$gene %in% selected_regions,] %>% dplyr::select(region = gene, cell_type, logFC = avg_logFC, cluster)

heatmap_matrix <- df_heat %>%
  dplyr::select(-cluster) %>%   # 丢掉 cluster 列
  pivot_wider(names_from = cell_type, values_from = logFC, values_fill = 0) %>%
  column_to_rownames(var = "region") %>%
  as.matrix()

heatmap_matrix <- heatmap_matrix[, km1_logfc_order]
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
p <- pheatmap(heatmap_matrix,
              breaks = seq(-3, 3, length.out = 101),
              annotation_col = ct_annotation,
              annotation_row = row_annotation,
              annotation_colors = list(
                tissue = tissue_colors,
                lineage = lineage_colors,
                cluster = cluster_colors
              ),
              color = colorRampPalette(c("blue", "white", "red"))(100),
              cluster_rows = TRUE,
              cluster_cols = FALSE,
              show_rownames = FALSE,
              show_colnames = TRUE,
              scale = "none",
              filename = "data/combined_H3K9me3_heatmap2_logfc.pdf",
              width = 32,
              height = 24,
              angle_col = 90,
              fontsize = 20,
              fontsize_col = 22,
              fontsize_main = 30,
              main = paste0("H3K9me3 aging logFC of large ",
                            dim(heatmap_matrix)[1],
                            " H3K9me3 peaks(>200kb) from zhuojie")
)
p

row_order <- p$tree_row$order
heatmap_matrix_ordered <- heatmap_matrix[row_order, ]
write.csv(heatmap_matrix_ordered,
          file = "data/H3K9me3_logfc_heatmap_matrix.csv",
          quote = FALSE)



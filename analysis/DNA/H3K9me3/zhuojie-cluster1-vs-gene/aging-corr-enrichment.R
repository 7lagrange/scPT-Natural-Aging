
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

setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/analysis/H3K9me3/H3K9me3-gene-TE-corr")

df_corr <- read.csv("data/H3K9me3_gene_logFC_corr.csv")
genes <- unique(df_corr$gene)


########################## 2. 选择相关性显著且为负的基因 ########################## 
select_genes <- df_corr %>%
  filter(spearman_fdr < 0.05, spearman_r < 0) %>%
  pull(gene) %>%
  unique()

select_genes <- df_corr %>%
  filter(spearman_p < 0.05, spearman_r < 0) %>%
  pull(gene) %>%
  unique()


egobp_up <- clusterProfiler::enrichGO(
  keyType = "SYMBOL",
  gene = select_genes,
  universe = genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",  # biological process
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  readable = FALSE
)
egobp_up
barplot(egobp_up, showCategory = 50)



########################## 3. 选择相关性显著且为负的基因 ########################## 
select_genes <- df_corr %>%
  filter(spearman_fdr < 0.05, spearman_r > 0) %>%
  pull(gene) %>%
  unique()



egobp_up <- clusterProfiler::enrichGO(
  keyType = "SYMBOL",
  gene = select_genes,
  universe = genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",  # biological process
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  readable = FALSE
)
egobp_up

barplot(egobp_up, showCategory = 20)


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


# 构造映射表
lineage_map <- tribble(
  ~tissue,   ~cell_type,                              ~lineage,
  # BAT
  "BAT",     "BAT-Endothelial_cell",                  "Endothelial",
  "BAT",     "BAT-ASPC",                              "Stromal",
  "BAT",     "BAT-Adipocyte_C1",                      "Adipose",
  "BAT",     "BAT-Adipocyte_C2",                      "Adipose",
  "BAT",     "BAT-Epithelial_cell",                   "Epithelial",
  "BAT",     "BAT-Macrophage",                        "Immune",
  
  # brainCB
  "brainCB", "brainCB-GC_ExN",                        "ExN",
  "brainCB", "brainCB-Ogc_Glia",                      "Glia",
  "brainCB", "brainCB-Asc_Glia",                      "Glia",
  "brainCB", "brainCB-MLI1_InN",                      "InN",
  "brainCB", "brainCB-PLI_InN",                       "InN",
  "brainCB", "brainCB-cerebellar_Golgi_cell",        "InN",
  "brainCB", "brainCB-Bergmann_glial_cell",          "Glia",
  "brainCB", "brainCB-inhibitory_interneuron",        "InN",
  "brainCB", "brainCB-astrocyte_of_the_cerebellum", "Glia",
  "brainCB", "brainCB-cerebellar_granule_cell",      "ExN",
  
  # brainFC
  "brainFC", "brainFC-L6_ExN",                        "ExN",
  "brainFC", "brainFC-L4_ExN",                        "ExN",
  "brainFC", "brainFC-L2_3_ExN",                     "ExN",
  "brainFC", "brainFC-L5.Deptor_ExN",                "ExN",
  "brainFC", "brainFC-Claustrum_ExN",                 "ExN",
  "brainFC", "brainFC-InN.Six3_InN",                 "InN",
  "brainFC", "brainFC-L5.Parm1_ExN",                 "ExN",
  "brainFC", "brainFC-InN.Meis2_InN",                "InN",
  "brainFC", "brainFC-InN.Lrrc38.Plcxd3_InN",       "InN",
  "brainFC", "brainFC-InN.Sst.Pvalb",                "InN",
  "brainFC", "brainFC-InN.Vip.Npy",                   "InN",
  "brainFC", "brainFC-L5.Deptor",                     "ExN",
  "brainFC", "brainFC-InN.Six3",                      "InN",
  "brainFC", "brainFC-Claustrum",                     "ExN",
  "brainFC", "brainFC-InN.Meis2",                     "InN",
  "brainFC", "brainFC-Ogc",                           "Glia",
  "brainFC", "brainFC-Asc",                           "Glia",
  "brainFC", "brainFC-L4",                            "ExN",
  "brainFC", "brainFC-L6",                            "ExN",
  "brainFC", "brainFC-L2.3",                          "ExN",
  "brainFC", "brainFC-L5.Parm1",                      "ExN",
  
  # brainHip
  "brainHip","brainHip-Asc_Glia",                     "Glia",
  "brainHip","brainHip-Ogc_Glia",                     "Glia",
  "brainHip","brainHip-DG_ExN",                       "ExN",
  "brainHip","brainHip-SubEnt_ExN",                  "ExN",
  "brainHip","brainHip-CA1_ExN",                     "ExN",
  "brainHip","brainHip-CAsub2_ExN",                  "ExN",
  "brainHip","brainHip-CA3_ExN",                     "ExN",
  "brainHip","brainHip-InN_InN",                      "InN",
  "brainHip","brainHip-Inh",                          "InN",
  "brainHip","brainHip-CA2_3",                        "ExN",
  "brainHip","brainHip-Sub_Ent",                      "ExN",
  "brainHip","brainHip-CA1",                          "ExN",
  "brainHip","brainHip-DG",                           "ExN",
  "brainHip","brainHip-Asc",                          "Glia",
  "brainHip","brainHip-Peri",                         "Endothelial",
  
  # colon
  "colon",   "colon-Macrophages",                     "Immune",
  "colon",   "colon-Goblet",                          "Epithelial",
  "colon",   "colon-StemCells",                       "Stemcell",
  "colon",   "colon-Epithelial",                      "Epithelial",
  
  # heart
  "heart",   "heart-Endocardial_endothelial_cells",   "Endothelial",
  "heart",   "heart-Vascular_endothelial_cells",      "Endothelial",
  "heart",   "heart-Mural_cells",                     "Stromal",
  "heart",   "heart-Fibroblasts",                     "Stromal",
  "heart",   "heart-Lymphoid_cells_T_cells",          "Immune",
  "heart",   "heart-Lymphoid_cells_B_cells",          "Immune",
  "heart",   "heart-Myeloid_cells",                   "Immune",
  "heart",   "heart-Ventricular_cardiomyocytes",      "Muscle",
  
  # kidney
  "kidney",  "kidney-proximal_convoluted_tubule",     "Epithelial",
  "kidney",  "kidney-proximal_straight_tubule",       "Epithelial",
  "kidney",  "kidney-endo",                           "Endothelial",
  "kidney",  "kidney-Stroma",                         "Stromal",
  "kidney",  "kidney-Tcell",                         "Immune",
  "kidney",  "kidney-marco",                         "Immune",
  
  # liver
  "liver",   "liver-endo",                            "Endothelial",
  "liver",   "liver-HepaticStellateCells",            "Stromal",
  "liver",   "liver-Hepatocytes",                     "Epithelial",
  "liver",   "liver-KupfferCells",                    "Immune",
  
  # lung
  "lung",    "lung-vascular_endothelial_cells",       "Endothelial",
  "lung",    "lung-Type1_pneumocytes",                "Epithelial",
  "lung",    "lung-Type_2_pneumocytes",               "Epithelial",
  "lung",    "lung-Lipofibroblast",                   "Stromal",
  "lung",    "lung-Ciliated_cells",                   "Epithelial",
  "lung",    "lung-Club_cells",                       "Epithelial",
  "lung",    "lung-CD8+_T_cells",                       "Immune",
  "lung",    "lung-non-classical_monocyte_Ly6c2-",                       "Immune",
  "lung",    "lung-Eosinophils",                       "Immune",
  "lung",    "lung-B_cells",                       "Immune",
  "lung",    "lung-Interstitial_macrophages",                       "Immune",
  "lung",    "lung-Alveolar_macrophage",                       "Immune",
  "lung",    "lung-Mesothelial_cells",                       "Stromal",
  
  
  "spleen",    "spleen-B_cell",                       "Immune",
  
  
  
  # muscle
  "muscle",  "muscle-endothelial_cell",               "Endothelial",
  "muscle",  "muscle-cell_of_skeletal_muscle",        "Muscle",
  "muscle",  "muscle-macrophage",                     "Immune",
  "muscle",  "muscle-mesenchymal_stem_cell",          "Stromal",
  
  # testis
  "testis",  "testis-diplotene",                      "Germline",
  "testis",  "testis-late_round_spermatid",         "Germline",
  "testis",  "testis-leptotene_zygotene",           "Germline",
  "testis",  "testis-sertoli",                        "Stromal"
)


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

tissues <- c( "brainFC", "brainHip", "brainCB", "liver", "lung", "heart", "kidney", "BAT", "muscle", "testis", "colon", "spleen" )


all_depth <- list()
for (tissue in tissues) {
  file_path <- paste0(
    "~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/analysis/H3K9me3/H3K9me3-zhuojie-regions/data/", tissue, "_", antibody, "depth.csv"
  )
  df <- read.csv(file_path, quote = "", row.names = NULL)
  df$cell_type <- paste0(tissue, "-", df$cell_type)
  df$tissue <- tissue
  all_depth[[tissue]] <- df
}
dna_reads_per_cell_df<- do.call(rbind, all_depth)

ct_annotation <- dna_reads_per_cell_df %>%
  dplyr::select(tissue, cell_type, total_reads) %>%
  arrange(desc(total_reads)) %>%   # 按 depth 从高到低
  as.data.frame()

ct_annotation <- ct_annotation %>% left_join(lineage_map, by = c("tissue", "cell_type")) 

write.csv(ct_annotation, paste0("data/ct_annotation.csv"), quote = FALSE, row.names = FALSE)

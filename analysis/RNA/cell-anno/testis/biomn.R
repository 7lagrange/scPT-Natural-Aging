library(Seurat)
library(ggplot2)

# ============================================================
# 1. 读取数据
# ============================================================
seurat_obj <- readRDS("testis_sub.rds")

# ============================================================
# 2. 定义 Marker Gene List
# ============================================================
# SSC → SPG(undiff) → SPG(diff) → SPC → Spermatid(RS) → Sperm
# Preleptotene → Leptotene → Zygotene → Pachytene → Diplotene

marker_genes <- list(
  "SSC"            = c("Rhox10", "Mageb4", "Zbtb16", "Nanos3", "Rarg", "Afp"),
  "Early_diff_SPG" = c("Dmrt1", "Crabp1", "Msh2", "Map7d2", "Scml2"),
  "Late_diff_SPG"  = c("Rhox13", "Prdm9", "Meiob"),
  "EarlySPC"       = c("Cdc42ep3", "Ggnbp1", "Mllt10"),        # Pou5f2 不在数据中
  "LateSPC"        = c("Ankrd54", "Spo11"),
  
  "PrelepSPC" = c("Stra8", "Meiosin", "Ccnb1ip1", "Tex101",
                  "Dazl", "Rec8", "Esx1", "Setdb2"),
  
  
  "Leptotene"        = c("Stra8","Ccnb1ip1","Hormad1","Hormad2","Sycp2"),
  "Zygotene"         = c("Sycp1","Dmc1","Rad51","Meiob","Spo11"),
  "Pachytene"        = c("Prdm9","Mlh1","Mlh3","Cdc25a","Hspa2"),
  "Diplotene"        = c("Ctnna3","Sh3gl3","Camk4","Fus","Aurka","Ccnb1",
                         "Gm45408","Snhg11","Eef2"),
  
  "EarlyRS"        = c("Catsper3", "Izumo1", "Spaca1", "Spaca3"),
 # "EarlyRS"   = c("Acrv1", "Acr", "Dpy19l2", "Spata16", "Catsper1", "Tex36", "Tnp1", "Prm1"),
  "MidRS"     = c("Catsper2", "Spem1", "Tssk1", "Tssk2",
                  "Spata3", "Hmgb4", "Ube2b", "Ddx4"),
  "LateRS"    = c("Catsper4", "Pex5l", "Sufu", "Tnp2", "Prm2",
                  "Smcp", "Odf1", "Odf2", "Spag4"),

  "ES_s1"          = c("Cby3", "Izumo2", "Tssk6"),              # Tnp1 不在数据中
  "ES_s2"          = c("Tnp2", "Reep6"),
  "ES_s3"          = c("Vwf", "Cd93"),
  
  
  "Endothelial"    = c("Flt1", "Icam2", "Ly6c1", "Pecam1"),
  "Leydig"         = c("Cyp11a1", "Cyp17a1", "Hsd17b3", "Hsd3b1", "Hsd3b6"),
  "Sertoli"        = c("Sox9", "Aard", "Amhr2", "Basp1"),
  "Myoid"          = c("Acta2", "Myh11"),
  "Telocyte"       = c("Dcn", "Col1a2", "Col3a1", "Col4a4", "Lama2"),
  "Macrophage"     = c("Cd74", "Csf1r", "Cd68", "Mrc1"),
  "T_cells"        = c("Cd2", "Cd3e", "Cd3g", "Ms4a4b", "Trac", "Trbc2")
)

# 过滤掉数据集中不存在的基因
marker_genes_clean <- lapply(marker_genes, function(g) g[g %in% rownames(seurat_obj)])

# ============================================================
# 3. 计算每个 cluster 的 marker 平均表达量
# ============================================================
Idents(seurat_obj) <- "seurat_clusters"
avg_exp <- AverageExpression(seurat_obj, 
                             features = unlist(marker_genes_clean), 
                             assay = "RNA", 
                             layer = "data")$RNA

# 修正列名（Seurat 会自动加 'g' 前缀）
colnames(avg_exp) <- gsub("^g", "", colnames(avg_exp))

# 计算每个 cluster 对每种细胞类型的平均分
cluster_scores <- data.frame(cluster = colnames(avg_exp))
for (ct in names(marker_genes_clean)) {
  genes <- marker_genes_clean[[ct]]
  if (length(genes) > 0) {
    cluster_scores[[ct]] <- colMeans(avg_exp[genes, , drop = FALSE])
  }
}

# 查看每个 cluster 得分最高的细胞类型
ct_cols <- names(marker_genes_clean)
cluster_scores$top_celltype <- apply(cluster_scores[, ct_cols], 1, function(x) names(which.max(x)))
cluster_scores$cluster <- as.integer(cluster_scores$cluster)
print(cluster_scores[order(cluster_scores$cluster), c("cluster", "top_celltype")])

# ============================================================
# 4. 定义最终注释（结合分数 + cluster 19 已知为 Endothelial）
# ============================================================
annotations <- c(
  "0"  = "EarlySPC",
  "1"  = "LateRS",
  "2"  = "ES_s2",
  "3"  = "EarlySPC",
  "4"  = "EarlySPC",
  "5"  = "Leydig",
  "6"  = "EarlySPC",
  "7"  = "EarlySPC",
  "8"  = "EarlySPC",
  "9"  = "EarlySPC",
  "10" = "Sertoli",
  "11" = "EarlySPC",
  "12" = "Sertoli",
  "13" = "Telocyte",
  "14" = "Early_diff_SPG",
  "15" = "Macrophage",
  "16" = "Sertoli",
  "17" = "Myoid",
  "18" = "EarlySPC",
  "19" = "Endothelial",
  "20" = "Sertoli",
  "21" = "Telocyte",
  "22" = "Early_diff_SPG",
  "23" = "Leydig"
)

# 写入 metadata
seurat_obj$cell_type_new <- plyr::mapvalues(
  seurat_obj$seurat_clusters,
  from = names(annotations),
  to   = unname(annotations)
)

# ============================================================
# 5. UMAP 可视化
# ============================================================
cell_type_colors <- c(
  "EarlySPC"       = "#4E79A7",
  "LateRS"         = "#F28E2B",
  "ES_s2"          = "#E15759",
  "Leydig"         = "#76B7B2",
  "Sertoli"        = "#59A14F",
  "Telocyte"       = "#EDC948",
  "Early_diff_SPG" = "#B07AA1",
  "Macrophage"     = "#FF9DA7",
  "Myoid"          = "#9C755F",
  "Endothelial"    = "#BAB0AC"
)

Idents(seurat_obj) <- "cell_type_new"

p_umap <- DimPlot(seurat_obj,
                  reduction = "umap_harmony",
                  group.by  = "cell_type_new",
                  label     = TRUE,
                  label.size = 3.5,
                  repel     = TRUE,
                  cols      = cell_type_colors,
                  pt.size   = 0.3) +
  ggtitle("Mouse Testis - Cell Type Annotation") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("umap_annotation.png", p_umap, width = 9, height = 7, dpi = 200)

# ============================================================
# 6. Dot Plot 验证
# ============================================================
dot_markers <- c(
  "Zbtb16", "Nanos3",       # SSC
  "Dmrt1", "Crabp1",        # Early_diff_SPG
  "Cdc42ep3", "Mllt10",     # EarlySPC
  "Catsper3", "Izumo1",     # EarlyRS
  "Catsper4",               # LateRS
  "Reep6",                  # ES_s2
  "Cyp11a1", "Hsd3b1",      # Leydig
  "Sox9", "Aard",           # Sertoli
  "Acta2", "Myh11",         # Myoid
  "Dcn", "Col1a2",          # Telocyte
  "Pecam1", "Flt1",         # Endothelial
  "Cd68", "Csf1r"           # Macrophage
)
dot_markers <- dot_markers[dot_markers %in% rownames(seurat_obj)]

ct_order <- c("Early_diff_SPG", "EarlySPC", "LateRS", "ES_s2",
              "Leydig", "Sertoli", "Myoid", "Telocyte", "Endothelial", "Macrophage")
seurat_obj$cell_type_new <- factor(seurat_obj$cell_type_new, levels = ct_order)
Idents(seurat_obj) <- "cell_type_new"

p_dot <- DotPlot(seurat_obj,
                 features = dot_markers,
                 group.by = "cell_type_new",
                 dot.scale = 5) +
  scale_color_gradient2(low = "lightblue", mid = "white", high = "darkred", midpoint = 0) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marker Gene Validation") +
  xlab("") + ylab("")

ggsave("dotplot_validation.png", p_dot, width = 11, height = 5.5, dpi = 200)

# ============================================================
# 7. 保存结果
# ============================================================
write.csv(
  data.frame(Cluster = as.integer(names(annotations)), Cell_Type = unname(annotations)),
  "cluster_annotation.csv", row.names = FALSE
)
saveRDS(seurat_obj, "testis_annotated.rds")

cat("Done! Files saved: umap_annotation.png, dotplot_validation.png,",
    "cluster_annotation.csv, testis_annotated.rds\n")
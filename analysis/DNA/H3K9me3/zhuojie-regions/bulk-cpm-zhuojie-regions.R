
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
library(EnhancedVolcano)


antibody <- 'H3K9me3'

# colon need a seperated process
proj <- "20240713-scPT-Natural-Aging-colon-H3K9me3-H3K36me3"
setwd(paste0('~/projects/paired_seq_tag/results/', proj))
tissue <- sapply(strsplit(proj, "-"), `[`, 5)

mouse_meta_data <- read.csv(paste0("/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/results/",tissue, "-merged/rna-qc/mouse_metadata.csv"))

dna_path <- paste0("/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/data/merge/", proj, "/dna_bam/histone_dna_bam/H3K9me3_zhuojie")

dna <- CreateSeuratObject(
  counts = Read10X(dna_path, gene.column = 1),
  project = proj,
  min.cells = 1,
  min.features = 1
)
dna
head(dna@meta.data)

####### add info from meta data  #######
cell_ids <- mouse_meta_data$cell_id
ind <- match(rownames(dna@meta.data), cell_ids)
dna$cell_type <- mouse_meta_data[ind, ][["cell_type"]]
dna$histone <- mouse_meta_data[ind, ][["histone"]]
dna$mouse <- mouse_meta_data[ind, ][["mouse"]]
dna$age <- mouse_meta_data[ind, ][["age"]]
dna$mouse_age <- paste0(dna$mouse, '-', dna$age)
dna$orig.ident <- sub(":.*", "", colnames(dna))

dna
dna <- subset(x = dna, subset = cell_type %in% unique(mouse_meta_data$cell_type))
dna

dna1 <- subset(x = dna, subset = mouse != "mouse137")


proj <- "20241121-scPT-Natural-Aging-colon-H3K9me3-H3K36me3"
setwd(paste0('~/projects/paired_seq_tag/results/', proj))
tissue <- sapply(strsplit(proj, "-"), `[`, 5)


mouse_meta_data <- read.csv(paste0("/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/results/",
                                   tissue, "-merged/rna-qc/mouse_metadata.csv"))

dna_path <- paste0("/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/data/merge/", proj, "/dna_bam/histone_dna_bam/H3K9me3_zhuojie")

dna <- CreateSeuratObject(
  counts = Read10X(dna_path, gene.column = 1),
  project = proj,
  min.cells = 1,
  min.features = 1
)
dna
head(dna@meta.data)

####### add info from meta data  #######
cell_ids <- mouse_meta_data$cell_id
ind <- match(rownames(dna@meta.data), cell_ids)
dna$cell_type <- mouse_meta_data[ind, ][["cell_type"]]
dna$histone <- mouse_meta_data[ind, ][["histone"]]
dna$mouse <- mouse_meta_data[ind, ][["mouse"]]
dna$age <- mouse_meta_data[ind, ][["age"]]
dna$mouse_age <- paste0(dna$mouse, '-', dna$age)
dna$orig.ident <- sub(":.*", "", colnames(dna))

dna
dna <- subset(x = dna, subset = cell_type %in% unique(mouse_meta_data$cell_type))
dna

dna2 <- subset(x = dna, subset = mouse != "mouse128")


dna <- merge(
  x = dna1,
  y = dna2,
  add.cell.ids = c("dna1", "dna2"),  # 给细胞加前缀，防止重复
  project = "DNA_merged"
)


table(dna$cell_type)
table(dna$orig.ident)
table(dna$mouse)
table(dna$mouse, dna$age)
table(dna$mouse, dna$histone)
table(dna$cell_type, dna$age)

DefaultAssay(dna) <- "RNA"
counts <- GetAssayData(dna, slot = "counts")
dna$cell_type <- paste0(tissue, '-', dna$cell_type)
dna$cell_type <- as.factor(dna$cell_type)

meta <- dna@meta.data
meta$cell <- colnames(dna)
meta <- meta[match(colnames(counts), meta$cell), ]

cell_types <- unique(dna$cell_type )
ct_counts <- do.call(cbind, lapply(cell_types, function(ct) {
  idx <- which(meta$cell_type == ct)  # 找到哪些列属于该cell_type
  Matrix::rowSums(counts[, idx, drop = FALSE])
}))

colnames(ct_counts) <- cell_types

all_ct_counts <- list() 
all_ct_counts[[proj]] <- ct_counts

# colon need a seperated process
projs <- c(
  "20241027-scPT-Natural-Aging-brainHip-H3K9me3-H3K36me3",
  "20241117-scPT-Natural-Aging-brainCB-H3K9me3-H3K36me3",
  "20240507-scPT-Natural-Aging-brainFC-H3K9me3-H3K36me3",
  
  "20240704-scPT-Natural-Aging-liver-H3K9me3-H3K36me3",
  "20240819-scPT-Natural-Aging-lung-H3K9me3-H3K36me3",
  "20241012-scPT-Natural-Aging-heart-H3K9me3-H3K36me3",
  "20240715-scPT-Natural-Aging-kidney-H3K9me3-H3K36me3",
  
  "20241125-scPT-Natural-Aging-BAT-H3K9me3-H3K36me3",
  "20250107-scPT-Natural-Aging-muscle-H3K9me3-H3K36me3",
  "20250410-scPT-Natural-Aging-testis-H3K9me3-H3K36me3",
  "20251013-scPT-Natural-Aging-spleen-H3K9me3-H3K36me3"
  
)




for (proj in projs) {
    cat("Processing:", proj, "\n")

    
    setwd(paste0('~/projects/paired_seq_tag/results/', proj))
    tissue <- sapply(strsplit(proj, "-"), `[`, 5)
    

    dna_path <- paste0("/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/data/merge/", proj, "/dna_bam/histone_dna_bam/H3K9me3_zhuojie")
    
    dna <- CreateSeuratObject(
      counts = Read10X(dna_path, gene.column = 1),
      project = proj,
      min.cells = 1,
      min.features = 1
    )
    dna
    head(dna@meta.data)
    
    ####### add info from meta data  #######
    

    mouse_meta_data <- read.csv(paste0("/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/results/",
                                       tissue, "-merged/rna-qc/mouse_metadata.csv"))

    cell_ids <- mouse_meta_data$cell_id
    ind <- match(rownames(dna@meta.data), cell_ids)
    dna$cell_type <- mouse_meta_data[ind, ][["cell_type"]]
    dna$histone <- mouse_meta_data[ind, ][["histone"]]
    dna$mouse <- mouse_meta_data[ind, ][["mouse"]]
    dna$age <- mouse_meta_data[ind, ][["age"]]
    dna$mouse_age <- paste0(dna$mouse, '-', dna$age)
    dna$orig.ident <- sub(":.*", "", colnames(dna))
    
    dna
    dna <- subset(x = dna, subset = cell_type %in% unique(mouse_meta_data$cell_type))
    dna
    
    
    table(dna$cell_type)
    table(dna$orig.ident)
    table(dna$mouse)
    table(dna$mouse, dna$age)
    table(dna$mouse, dna$histone)
    table(dna$cell_type, dna$age)
    
    
    # should not affect differential peaks
    DefaultAssay(object = dna) <- "RNA"
    
    ###  differential test 
    dna <- subset(x = dna, subset = cell_type != "mix")
    dna <- subset(x = dna, subset = cell_type != "doublet")
    
    
    DefaultAssay(dna) <- "RNA"
    counts <- GetAssayData(dna, slot = "counts")
    dna$cell_type <- paste0(tissue, '-', dna$cell_type)
    dna$cell_type <- as.factor(dna$cell_type)
    
    meta <- dna@meta.data
    meta$cell <- colnames(dna)
    meta <- meta[match(colnames(counts), meta$cell), ]
    
    cell_types <- unique(dna$cell_type )
    ct_counts <- do.call(cbind, lapply(cell_types, function(ct) {
      idx <- which(meta$cell_type == ct)  # 找到哪些列属于该cell_type
      Matrix::rowSums(counts[, idx, drop = FALSE])
    }))
    colnames(ct_counts) <- cell_types
    all_ct_counts[[proj]] <- ct_counts
    

}

library(Matrix)

# 假设 all_ct_counts 已经是 list，存了每个proj的 ct_counts

# 1. 所有 gene 的并集
all_genes <- sort(unique(unlist(lapply(all_ct_counts, rownames))))

# 2. 对每个 matrix 对齐行，不存在的基因补 0
all_ct_counts_aligned <- lapply(all_ct_counts, function(mat) {
  missing_genes <- setdiff(all_genes, rownames(mat))
  
  # 补全缺失行，值为0
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix(0, nrow = length(missing_genes), ncol = ncol(mat),
                       dimnames = list(missing_genes, colnames(mat)))
    mat <- rbind(mat, zero_mat)
  }
  
  # 按 all_genes 排序行
  mat[all_genes, , drop = FALSE]
})

# 3. 按列合并
combined_ct_counts <- do.call(cbind, all_ct_counts_aligned)

dim(combined_ct_counts)
head(colnames(combined_ct_counts))


dge <- DGEList(counts = combined_ct_counts)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

cpm_mat <- as.data.frame(edgeR::cpm(dge, log = FALSE)) 
cpm_mat$gene <- rownames(cpm_mat)

setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/analysis/H3K9me3/H3K9me3-cell-type-corr")

write.csv(cpm_mat, paste0("data/H3K9me3_bulk_cpm_zhuojie.csv"), quote = FALSE, row.names = FALSE)

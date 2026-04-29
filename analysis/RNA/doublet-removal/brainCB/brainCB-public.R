rm(list = ls())

options(stringsAsFactors = F)

library(reticulate)


setwd("~/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/RNA/doublet/brainCB")

use_python("~/miniconda3/bin/python", required = TRUE)
py_config()

library("Seurat")
library("anndata")
data <- read_h5ad("data/CB_Kozareva.h5ad")
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
data

data@meta.data
table(data$cell_type)

sce <- data
sce$celltype <- sce$cell_type

max(sce@assays[["RNA"]]@counts)


sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce = FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000)  
sce = ScaleData(sce, vars.to.regress = c("nFeature_RNA"))

cells.use <- sample(colnames(sce), 100000)
sce_50k <- subset(sce, cells = cells.use)

table(sce_50k$celltype)

sce_50k <- RunPCA(sce_50k, features = VariableFeatures(object = sce_50k))
reduction = "pca"
ndims = 1:20 # make it smaller when like like a star
sce_50k <- FindNeighbors(sce_50k, dims = ndims, reduction=reduction, k.param=15)
sce_50k <- FindClusters(sce_50k, resolution = 0.5)
sce_50k <- RunUMAP(sce_50k, dims = ndims, reduction=reduction)
p1 <- DimPlot( sce_50k, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4 )
p2 <- DimPlot( sce_50k, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 4 ) 
p2
p1+p2


head(rownames(sce_50k))
tmp <- read.table(file = '/mnt/transposon1/zhangyanxiaoLab/qihongjian/projects/aging_setdb1/data/bulk_rna/rna/bulk_rna_counts.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
leaned_gene_ids <- sub("\\..*", "", tmp$Geneid)  # Remove version number from ENSEMBL IDs
ensembl_gene_names <- setNames(tmp$gene_name, leaned_gene_ids)

gene_ids <- sub("\\..*", "", rownames(sce_50k))
new_names <- ensembl_gene_names[gene_ids]
new_names[is.na(new_names)] <- gene_ids[is.na(new_names)]
new_names_unique <- make.unique(new_names)



rownames(sce_50k[["RNA"]]@counts) <- new_names_unique
rownames(sce_50k[["RNA"]]@data)   <- new_names_unique
#rownames(sce_50k[["RNA"]]@scale.data)  diff gene set

saveRDS(sce_50k, file = "data/brainCB_public_10k.rds")

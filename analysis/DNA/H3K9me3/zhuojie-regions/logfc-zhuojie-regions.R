


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



antibody <- 'H3K9me3'
antibody <- 'H3K27me3'
antibody <- 'H3K36me3'
antibody <- 'H3K27ac'
antibody <- 'H3K4me1'
antibody <- 'H3K4me3'


tissues <- c( "brainFC", "brainHip", "brainCB", "liver", "lung", "heart", "kidney", 
              "BAT", "muscle", "testis", "colon", "spleen" )

all_results <- list()

for (tissue in tissues) {
  
  dna_path <- paste0("/storage/zhangyanxiaoLab/qihongjian/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/DNA/counts/data/",
                     tissue, "-", antibody, "-H3K9me3_zhuojie")
  dna <- CreateSeuratObject(
    counts = Read10X(dna_path, gene.column = 1),
    project = tissue,
    min.cells = 1,
    min.features = 1
  )
  dna
  head(dna@meta.data)
  
  mouse_meta_data <- read.csv(paste0("/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag/results/",
                                     tissue, "-merged/rna-qc/mouse_metadata.csv"))
  
  cell_ids <- mouse_meta_data$cell_id
  ind <- match(rownames(dna@meta.data), cell_ids)
  dna$cell_type <- mouse_meta_data[ind, ][["cell_type"]]
  dna$histone <- mouse_meta_data[ind, ][["histone"]]
  dna$mouse <- mouse_meta_data[ind, ][["mouse"]]
  dna$age <- mouse_meta_data[ind, ][["age"]]
  dna$proj <- mouse_meta_data[ind, ][["proj"]]
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
  dna <- subset(x = dna, subset = cell_type != "Doublet")
  dna <- subset(x = dna, subset = cell_type != "Doublets")
  
  pseudobulks = Libra::to_pseudobulk(input = dna, 
                                     min_reps = 1,
                                     replicate_col = "mouse", cell_type_col = "cell_type", 
                                     label_col = "age",external = T)
  
  cell_types <- names(pseudobulks)
  results_list <- list()
  for (i in seq_along(pseudobulks)) {
    cell_type <- cell_types[i]
    
    x <- pseudobulks[[i]]
    targets <- data.frame(group_sample = colnames(x)) %>% mutate(group = gsub(".*\\:", "", group_sample))
    targets$group <- factor(targets$group, levels = c("age3M", "age27M"))
    design <- model.matrix(~group, data = targets)
    y <- DGEList(counts = x, group = targets$group)
    
    
    keep <- rowSums(edgeR::cpm(y) > 1) >= 2
    
    sum(keep)
    sum(keep) / length(keep)
    
    y <- y[keep, , keep.lib.sizes = FALSE] # lib.size is new total counts
    y <- calcNormFactors(y, method = "TMM")
    design <- model.matrix(~ y$samples$group)
    y <- tryCatch({
      estimateDisp(y, design, robust = TRUE)
    }, error = function(e) {
      message("Skipping iteration ", i, " due to error: ", conditionMessage(e))
      return(NULL)  # return something safe so the loop can continue
    })
    
    # Skip rest of the iteration if estimateDisp failed
    if (is.null(y)) next
    
    cat("Processing:", cell_type, y$samples$norm.factors, "\n")
    y <- estimateDisp(y, design, robust = TRUE)
    
    exp_cpms <- as.data.frame(edgeR::cpm(y))
    fit <- glmFit(y, design = design)
    test <- glmLRT(fit)
    res <- topTags(test, n = Inf) %>% as.data.frame() %>% rownames_to_column("gene")
    
    exp_cpms$gene <- rownames(exp_cpms)
    x$gene <- rownames(x)
    res <- merge(res, exp_cpms, by = "gene")
    res <- merge(res, x, by = "gene", suffixes = c("_cpm", "_counts"))
    res$cell_type <- cell_type
    results_list[[i]] <- res
    print(cell_type)
    
  }
  DE_edgeR <- bind_rows(results_list)
  
  
  # save results
  gene_split <- strsplit(DE_edgeR$gene, "-")
  DE_edgeR$chrom <- sapply(gene_split, "[", 1)
  DE_edgeR$start <- as.numeric(sapply(gene_split, "[", 2))
  DE_edgeR$end <- as.numeric(sapply(gene_split, "[", 3))
  DE_edgeR$length <- DE_edgeR$start - DE_edgeR$end 
  DE_edgeR$p_val <- DE_edgeR$PValue
  DE_edgeR$p_val_adj <- DE_edgeR$FDR
  DE_edgeR$avg_logFC <- DE_edgeR$logFC
  

  DE_edgeR <- DE_edgeR[, c(
    "gene", "logCPM", "cell_type",
    "chrom", "start", "end", "length",
    "p_val", "p_val_adj", "avg_logFC"
  )]
  
  DE_edgeR$tissue <- tissue
  
  all_results[[tissue]] <- DE_edgeR
  }

DE_edgeR <- bind_rows(all_results)

write.csv(DE_edgeR, paste0("/storage/zhangyanxiaoLab/qihongjian/github/zhanglab-code/projects/paired_seq_tag/code/machine-experiments/analysis/H3K9me3/H3K9me3-zhuojie-regions/data/",
                           antibody,  "_edgeR_27M_over_3M_zhuojie.csv"), quote = FALSE, row.names = FALSE)

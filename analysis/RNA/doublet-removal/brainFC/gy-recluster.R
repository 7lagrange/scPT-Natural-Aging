


gy <- readRDS("/storage/zhangyanxiaoLab/yiguo/projects/Multi-batch_RNA_eachregion_cocluster_anno/data_FC/FC_three_batches_clean_annotated.rds")
Idents(gy) <- "RNA_snn_res.0.9"
DimPlot(gy, group.by = "RNA_snn_res.0.9", reduction = "umap", label = TRUE,  label.size = 5 ) 

cluster_id <- '29'
gy$seurat_clusters_ct <- paste0(gy$RNA_snn_res.0.9,":",gy$cell_anno)
p1 <- DimPlot( gy, reduction = "umap", group.by = "seurat_clusters_ct", raster=FALSE, label=T,
               cells.highlight = WhichCells( gy, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster29 in GY recluster"))
p2 <- DotPlot( subset( gy, idents = cluster_id ), features = c("Slc17a7","Gad2", "Apoe","C1qb",
                                                                 "Mog", "Pdgfra", "Tie1","Foxd1", "Dcn")
               , group.by = "seurat_clusters" ) + coord_flip()
p1 + p2


cells_9 <- WhichCells(gy, idents = "29")
cells_9_mouse <- intersect(cells_9, colnames(gila))
table(gila$seurat_clusters[cells_9_mouse])
gila$seurat_clusters_ct <- paste0(gila$seurat_clusters,":",gila$cell_anno)
DimPlot( gila, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label=T,
         cells.highlight = cells_9_mouse, cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("GY cluster 29 highlighted in Gila sub clusters"))


cluster_id <- '7'
p1 <- DimPlot( gila, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, 
               cells.highlight = WhichCells( gila, idents = cluster_id ), cols = "grey85", cols.highlight = "red" ) + ggtitle(paste0("cluster: ", cluster_id))
p2 <- DotPlot( subset( gila, idents = cluster_id ), features = c("Slc17a7","Gad2", "Apoe","C1qb",
                                                                 "Mog", "Pdgfra", "Tie1","Foxd1", "Dcn")
               , group.by = "seurat_clusters" ) + coord_flip()
p1 + p2





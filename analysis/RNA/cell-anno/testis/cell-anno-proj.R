

mouse <- subset(mouse, subset = proj=="20250524-scPT-Natural-Aging-testis-H3K27ac-H3K27me3")
mouse
DimPlot(mouse, reduction = "umap_harmony", group.by="seurat_clusters",  label = T, label.size = 5)

# https://mp.weixin.qq.com/s/eK0T3Fwcv9Gs4LK7j2esiQ
# https://biomni.phylo.bio/projects/proj_144d36958f/tasks/sess_8d91aebb2bbe




##################### Leydig ################
Idents(mouse) <- "seurat_clusters"
genes = c("Cyp11a1", "Cyp17a1", "Hsd17b3", "Hsd3b1", "Hsd3b6") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 



##################### Telocyte ################
Idents(mouse) <- "seurat_clusters"
genes = c("Dcn", "Col1a2", "Col3a1", "Col4a4", "Lama2") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 




##################### Myoid ################
Idents(mouse) <- "seurat_clusters"
genes = c("Acta2", "Myh11") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### Sertoli ################
Idents(mouse) <- "seurat_clusters"
genes = c("Sox9", "Aard", "Amhr2", "Basp1") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 





##################### cluster 17 ################
FeaturePlot( subset(mouse, idents = 10), features = c("nFeature_RNA", "nCount_RNA","percent_mt","percent_ribo"),
             ncol = 2, reduction = "umap_harmony", order = TRUE, label = TRUE )


##################### EarlyRS ################
Idents(mouse) <- "seurat_clusters"
genes = c("Catsper3", "Izumo1", "Spaca1", "Spaca3") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 

##################### MidRS ################
Idents(mouse) <- "seurat_clusters"
genes = c("Catsper2", "Spem1", "Tssk2","Spata3",  "Ube2b", "Ddx4") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### LateRS ################
Idents(mouse) <- "seurat_clusters"
genes = c("Catsper4", "Pex5l", "Sufu", "Tnp2", "Prm2",
          "Smcp", "Odf1", "Odf2", "Spag4") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### ES ################
Idents(mouse) <- "seurat_clusters"
genes = c("Cby3", "Izumo2", "Tssk6", "Tnp2", "Reep6", "Vwf", "Cd93") # 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### SPC ################
Idents(mouse) <- "seurat_clusters"
genes = c("Cdc42ep3", "Ggnbp1", "Mllt10", # EarlySPC
          "Ankrd54", "Spo11") # lateSPC
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


##################### SSC ################
Idents(mouse) <- "seurat_clusters"
genes = c("Dmrt1", "Crabp1", "Msh2", "Map7d2", "Scml2",
          "Rhox13", "Prdm9", "Meiob"
          ) # lateSPC
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(4),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2


##################### cluster20 ################
Idents(mouse) <- "seurat_clusters"
genes = c("Rhox10", "Mageb4", "Zbtb16", "Nanos3", "Rarg", "Afp") #  , , , , , 
mouse <- AddModuleScore( object = mouse, features = genes, ctrl = 100,  name = "cluster_marker" )
FeaturePlot(mouse, features = c(genes, "cluster_marker1"), reduction = "umap_harmony", label = TRUE, order=T)
DotPlot(mouse, features = c(genes, "cluster_marker1"), group.by = "seurat_clusters", scale = T) 


Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(20),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2




##################### cluster 7 ################
Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(7),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2




##################### cluster 15 ################
Idents(mouse) <- "seurat_clusters"
mouse$one_vs_others <- ifelse(
  mouse$seurat_clusters %in% c(15),
  "marker", "others")

Idents(mouse) <- "one_vs_others"
markers <- FindAllMarkers(
  mouse,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)
markers <- markers %>% dplyr::filter(cluster == "marker") %>% dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
unique(markers$gene)[1:6]
Idents(mouse) <- "seurat_clusters"
FeaturePlot(mouse, features = unique(markers$gene)[1:6], ncol = 3, reduction = "umap_harmony", order=T, label=T) 

mouse <- AddModuleScore( object = mouse, features = list(unique(markers$gene)[1:6]), ctrl = 100,  name = "cluster_marker" )
p1 <- FeaturePlot(mouse, features = "cluster_marker1", reduction = "umap_harmony", label = TRUE, order=T,
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score"))

public <- AddModuleScore( object = public, features = list(unique(markers$gene)[1:6]), ctrl = 100, name = "cluster_marker" )
p2 <- FeaturePlot(public, features = "cluster_marker1", reduction = "umap", label = TRUE, order=T, 
                  cols = c("grey90", "blue", "red", "darkred")) + ggtitle(paste0("Top 6 Markers Module Score in heart\n source:Caoyunjue-science"))
p1 + p2

########### add cell type from text ############
table(public$celltype)

Idents(mouse) <- "seurat_clusters"
new.cluster.ids <- read.table(file = "cellAnno.txt", header = TRUE, sep='\t')$label
new.cluster.ids <- gsub("/", "_", new.cluster.ids)
new.cluster.ids <- gsub("\\(", "", new.cluster.ids)
new.cluster.ids <- gsub("\\)", "", new.cluster.ids)

names(new.cluster.ids) <- levels(mouse)
mouse <- RenameIdents(mouse, new.cluster.ids)
mouse$cell_type <- mouse@active.ident
mouse$cell_type_cluster <- paste0(mouse@active.ident, mouse$seurat_clusters)


p1 <- DimPlot(mouse, reduction = "umap_harmony", group.by="cell_type", label = T, label.size = 6)
p1


meta <- mouse@meta.data
meta$cell_id <- rownames(meta)
table(meta$cell_type)

write.csv(meta, 
          file = "testis-first-round-meta.csv", 
          row.names = FALSE)




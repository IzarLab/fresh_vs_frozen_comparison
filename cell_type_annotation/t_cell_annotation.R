library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(grid)
library(stringr)
library(dplyr)

#define folder locations for cutaneous melanoma and uveal primary data
foldersList = c("",
  "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata")
integrated_name_arr = c("BI5","cpoi-uvealprimarydata")

for (i in 1:length(foldersList)) {
  #load in rds objects
  if (i==1)
  {
    seu = readRDS("/data/fresh_vs_frozen_all_reannotate_BI5.rds")
  }
  else
  {
    system(paste0("aws s3 cp ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
    seu = readRDS(paste0("/data/",integrated_name_arr[i],"_integrated.rds"))
  }

  #subset to just T and NK cells, reintegrate these cells
  seu = subset(seu, manual_annotation_label %in% c("T cells/NK cells","NK cells"))

  DefaultAssay(seu) = "RNA"

  object.list = c()
  unique_idents = unique(seu$orig.ident)
  for (i1 in 1:length(unique_idents))
  {
    object.list = c(object.list, subset(seu, orig.ident==unique_idents[i1]))
  }

  dim_num = 25

  anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:dim_num)

  seu <- IntegrateData(anchorset = anchors, dims = 1:dim_num)

  seu <- ScaleData(object = seu)
  seu <- RunPCA(object = seu)
  seu <- FindNeighbors(seu, dims = 1:dim_num)
  seu <- FindClusters(seu)
  seu <- RunUMAP(object = seu, dims = 1:dim_num)

  #exclude some clusters that are likely not T/NK cells
  if (i==1)
  {
    seu$placeholder = TRUE
    seu$placeholder[seu$seurat_clusters %in% c(10)] = FALSE
    seu = subset(seu, placeholder)
  }
  if (i==2)
  {
    seu$placeholder = TRUE
    seu$placeholder[seu$seurat_clusters %in% c(10,12)] = FALSE
    seu = subset(seu, placeholder)
  }

  #reintegrate cells after excluding above clusters
  object.list = c()
  unique_idents = unique(seu$orig.ident)
  for (i1 in 1:length(unique_idents))
  {
    object.list = c(object.list, subset(seu, orig.ident==unique_idents[i1]))
  }

  dim_num = 10

  anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:dim_num)

  seu <- IntegrateData(anchorset = anchors, dims = 1:dim_num)

  seu <- ScaleData(object = seu)
  seu <- RunPCA(object = seu)
  seu <- FindNeighbors(seu, dims = 1:dim_num)
  seu <- FindClusters(seu)
  seu <- RunUMAP(object = seu, dims = 1:dim_num)

  #print umaps of seurat clusters and sample names
  pdf(paste0("fresh_vs_frozen_tcells_reannotate/fresh_vs_frozen_tcells_reannotate_",integrated_name_arr[i],".pdf"),width=12,height=12)
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "seurat_clusters", repel = T, label.size = 3) + ggtitle(paste0('seurat_clusters')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "orig.ident", repel = T, label.size = 3) + ggtitle(paste0('orig.ident')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))

  #reannotate t-cells based on marker gene expression below to finer classification
  seu$manual_annotation_label_tcell = "unknown"
  if (i==1)
  {
    seu$manual_annotation_label_tcell[seu$seurat_clusters %in% c(7)] = "Cycling T-cells"
    seu$manual_annotation_label_tcell[seu$seurat_clusters %in% c(0,1,4,6,8,9,11)] = "CD8+ T-cells"
    seu$manual_annotation_label_tcell[seu$seurat_clusters %in% c(3)] = "CD4+ T-cells"
    seu$manual_annotation_label_tcell[seu$seurat_clusters %in% c(2)] = "T-regs"
    seu$manual_annotation_label_tcell[seu$seurat_clusters %in% c(5)] = "TFH cells"
    seu$manual_annotation_label_tcell[seu$seurat_clusters %in% c(10)] = "NK cells"
  }
  if (i==2)
  {
    seu$manual_annotation_label_tcell[seu$seurat_clusters %in% c(7)] = "Cycling T-cells"
    seu$manual_annotation_label_tcell[seu$seurat_clusters %in% c(0,1,2,3,4,5,9,10,12)] = "CD8+ T-cells"
    seu$manual_annotation_label_tcell[seu$seurat_clusters %in% c(6)] = "CD4+ T-cells"
    seu$manual_annotation_label_tcell[seu$seurat_clusters %in% c(8)] = "T-regs"
    seu$manual_annotation_label_tcell[seu$seurat_clusters %in% c(11)] = "NK cells"
  }

  #print umaps of finer reannotation, QC metrics, and marker gene expression
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "manual_annotation_label_tcell", repel = T, label.size = 3) + ggtitle(paste0('manual_annotation_label_tcell')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))

  DefaultAssay(seu) = "RNA"

  print(FeaturePlot(seu, features = c("nFeature_RNA","nCount_RNA","percent.mt","ScrubDoublet_score"), min.cutoff = "q9", max.cutoff = "q90"))

  canonical_markers = c("CD4", "CD8A", "CD8B", "TCF7", "TOX", "FOXP3", "CXCL13", "NCAM1", "TOP2A", "MKI67", "IL7R","CD14","MRC1","XCL1","XCL2","TYROBP","KLRD1")
  for (j in seq(1,length(canonical_markers),4))
  {
    markers_to_plot = canonical_markers[j:min(j+3,length(canonical_markers))]
    print(FeaturePlot(seu, features = c(markers_to_plot), min.cutoff = "q9", max.cutoff = "q90"))
  }
  DefaultAssay(seu) = "integrated"
  
  #print heatmap of top 10 differential markers for each cluster
  allmarkers = FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  print(DoHeatmap(seu, features = top10$gene) + NoLegend() + ggtitle(paste0("DEG Heatmap")))
  write.table(allmarkers,paste0("fresh_vs_frozen_tcells_reannotate/fresh_vs_frozen_tcells_reannotate_",integrated_name_arr[i],"_markers.csv"),sep=",",quote=F,row.names=F,col.names=T)

  dev.off()

  saveRDS(seu,paste0("/data/fresh_vs_frozen_tcells_reannotate_",integrated_name_arr[i],".rds"))
}
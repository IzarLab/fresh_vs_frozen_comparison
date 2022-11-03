library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(grid)
library(stringr)

recluster_ribas_tcells = TRUE
recluster_ribas_myeloid = FALSE

#if options are set, select only t-cells or myeloid from ribas dataset to subcluster
if (recluster_ribas_tcells)
{
  outputcsv = "reannotate_ribas_melanoma_merged_tcells_markers.csv"
  outputpdf = "reannotate_ribas_melanoma_merged_tcells.pdf"
  outputrds = "/data/reannotate_ribas_melanoma_merged_tcells.rds"
  pdf(outputpdf,height=20,width=20)

  seu = readRDS("/data/reannotate_ribas_melanoma_merged.rds")
  seu = subset(seu, seurat_clusters %in% c(1,3,12,13,15))
} else if (recluster_ribas_myeloid) {
  outputcsv = "reannotate_ribas_melanoma_merged_myeloid_markers.csv"
  outputpdf = "reannotate_ribas_melanoma_merged_myeloid.pdf"
  outputrds = "/data/reannotate_ribas_melanoma_merged_myeloid.rds"
  pdf(outputpdf,height=20,width=20)

  seu = readRDS("/data/reannotate_ribas_melanoma_merged.rds")
  seu = subset(seu, seurat_clusters %in% c(4,15,19))
} else {
  #otherwise, download entire data for each ribas sample, score for cell cycle markers, and merge them into one large rds object
  outputcsv = "reannotate_ribas_melanoma_merged_markers.csv"
  outputpdf = "reannotate_ribas_melanoma_merged.pdf"
  outputrds = "/data/reannotate_ribas_melanoma_merged.rds"
  
  ribas_samples = c("ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh","ribas_310_on_later_previd_3_GEX_titrate_thresh","ribas_310_pre_GEX_5pv2_S26_L004_titrate_thresh")
  cluster_numbers = list(c(2,5,9,10),c(0,2,3,5,9,11),c(11))
  pdf(outputpdf,height=20,width=20)
  object.list = list()
  for (z in 1:length(ribas_samples)) {
    ribas_sample = ribas_samples[z]
    system(paste0("aws s3 cp s3://melanoma-ribas/ribas1/Seurat/",ribas_sample,"/",ribas_sample,"_cb.rds /data/",ribas_sample,"_cb.rds"))
    temp_rds = readRDS(paste0("/data/",ribas_sample,"_cb.rds"))

    temp_rds = CellCycleScoring(temp_rds, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)

    # print(DimPlot(temp_rds, reduction = "umap", label = T, group.by = "celltype_bped_main", repel = T, label.size = 3, shuffle = T) + ggtitle(paste0(ribas_sample,' celltype_bped_main')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)) + ggtitle(paste0(ribas_sample," celltype_bped_main")))
    # print(DimPlot(temp_rds, reduction = "umap", label = T, group.by = "seurat_clusters", repel = T, label.size = 3, shuffle = T) + ggtitle(paste0(ribas_sample,' seurat_cluster')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)) + ggtitle(paste0(ribas_sample," orig.ident")))
    # print(FeaturePlot(temp_rds, features = c("S.Score"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(ribas_sample," S.Score")))
    # print(FeaturePlot(temp_rds, features = c("G2M.Score"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(ribas_sample," G2M.Score")))

    # temp_rds = subset(temp_rds, seurat_clusters %in% cluster_numbers[[z]])

    if (length(colnames(temp_rds))>=100)
    {
      object.list = c(object.list, temp_rds)
    }
    system(paste0("rm /data/",ribas_sample,"_cb.rds"))
  }

  source("merge.SCTAssay.R")

  seu = merge.SCTAssay(x=object.list[[1]], y=object.list[[2]], merge.data = T, na.rm = T)
  seu = merge.SCTAssay(x=seu, y=object.list[[3]], merge.data = T, na.rm = T)
}

#recluster data using 25 dimensions in UMAP
dim_num = 25

seu <- ScaleData(object = seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:dim_num)
seu <- FindClusters(seu)
seu <- RunUMAP(object = seu, dims = 1:dim_num)

#if reclustering t cells, exclude a few clusters that are likely not t cells, and rerun umap
if (recluster_ribas_tcells)
{
  seu$placeholder = TRUE
  seu$placeholder[seu$seurat_clusters %in% c(6,13)] = FALSE
  seu = subset(seu, placeholder)

  seu <- ScaleData(object = seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- RunPCA(object = seu)
  seu <- FindNeighbors(seu, dims = 1:dim_num)
  seu <- FindClusters(seu)
  seu <- RunUMAP(object = seu, dims = 1:dim_num)

  seu$placeholder = TRUE
  seu$placeholder[seu$seurat_clusters %in% c(7)] = FALSE
  seu = subset(seu, placeholder)

  seu <- ScaleData(object = seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- RunPCA(object = seu)
  seu <- FindNeighbors(seu, dims = 1:dim_num)
  seu <- FindClusters(seu)
  seu <- RunUMAP(object = seu, dims = 1:dim_num)
}

seu$orig.ident[seu$orig.ident=="ribas_310_on_GEX_5pv2_S27_L004"] = "ribas_310_on"
seu$orig.ident[seu$orig.ident=="ribas_310_on_later_previd_3_GEX"] = "ribas_310_on_later"
seu$orig.ident[seu$orig.ident=="ribas_310_pre_GEX_5pv2_S26_L004"] = "ribas_310_pre"

#print umaps of seurat clusters and sample name
print(DimPlot(seu, reduction = "umap", label = T, group.by = "orig.ident", repel = T, label.size = 3, shuffle = T) + ggtitle('orig.ident') + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)) + ggtitle("merge object orig.ident"))
print(DimPlot(seu, reduction = "umap", label = T, group.by = "seurat_clusters", repel = T, label.size = 3, shuffle = T) + ggtitle('seurat_cluster') + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))

# anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:dim_num)

# seu <- IntegrateData(anchorset = anchors, dims = 1:dim_num)

# seu <- ScaleData(object = seu)
# seu <- RunPCA(object = seu)
# seu <- FindNeighbors(seu, dims = 1:dim_num)
# seu <- FindClusters(seu)
# seu <- RunUMAP(object = seu, dims = 1:dim_num)

#run cell cycle scoring, and print scatterplot of cell cycle scores
seu = CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
tempdf = data.frame(S=seu$S.Score, G2M=seu$G2M.Score)
print(ggplot(tempdf) + geom_point(aes(x=S, y=G2M)) + ggtitle("cell_cycle scoring"))

#print violin plots of QC metrics by seurat clusters
print(VlnPlot(seu, features = c("nCount_RNA"), group.by = "seurat_clusters", pt.size = 0))
print(FeaturePlot(seu, features = c("nCount_RNA"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0("nCount_RNA")))
print(VlnPlot(seu, features = c("nFeature_RNA"), group.by = "seurat_clusters", pt.size = 0))
print(FeaturePlot(seu, features = c("nFeature_RNA"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0("nFeature_RNA")))
print(VlnPlot(seu, features = c("percent.mt"), group.by = "seurat_clusters", pt.size = 0))
print(FeaturePlot(seu, features = c("percent.mt"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0("percent.mt")))
print(VlnPlot(seu, features = c("ScrubDoublet_score"), group.by = "seurat_clusters", pt.size = 0))
print(FeaturePlot(seu, features = c("ScrubDoublet_score"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0("ScrubDoublet_score")))

#if option is set, take only non-cycling cells based on cell cycle scoring, and recluster
recluster_by_cell_cycle = FALSE
if (recluster_by_cell_cycle)
{
  seu = subset(seu, G2M.Score<.25 & S.Score<.5)

  seu <- ScaleData(object = seu)
  seu <- RunPCA(object = seu)
  seu <- FindNeighbors(seu, dims = 1:dim_num)
  seu <- FindClusters(seu)
  seu <- RunUMAP(object = seu, dims = 1:dim_num)

  seu$orig.ident[seu$orig.ident=="ribas_310_on_GEX_5pv2_S27_L004"] = "ribas_310_on"
  seu$orig.ident[seu$orig.ident=="ribas_310_on_later_previd_3_GEX"] = "ribas_310_on_later"
  seu$orig.ident[seu$orig.ident=="ribas_310_pre_GEX_5pv2_S26_L004"] = "ribas_310_pre"

  print(DimPlot(seu, reduction = "umap", label = T, group.by = "seurat_clusters", repel = T, label.size = 3, shuffle = T) + ggtitle(paste0(ribas_sample,' seurat_cluster')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "orig.ident", repel = T, label.size = 3, shuffle = T) + ggtitle(paste0(ribas_sample,' orig.ident')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)) + ggtitle("integrated object orig.ident"))
}

#print umaps of marker genes
if (recluster_ribas_tcells) {
  canonical_markers = c("CD4", "CD8A", "CD8B", "TCF7", "TOX", "FOXP3", "CXCL13", "NCAM1", "TOP2A", "IL7R","CD14","MRC1")
} else if (recluster_ribas_myeloid) {
  canonical_markers = c("CD4", "CD8A", "CD8B", "TCF7", "TOX", "FOXP3", "CXCL13", "NCAM1", "TOP2A", "IL7R","CD14","MRC1")
} else {
  canonical_markers = c("CD4", "CD8A", "CD8B", "TCF7", "TOX", "FOXP3", "CXCL13", "NCAM1", "TOP2A", "IL7R","CD14","MRC1","PECAM1","DCUN1D3","FAP","ACTA2","MS4A1")
}
for (marker in canonical_markers)
{
  print(FeaturePlot(seu, features = c(marker), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(marker)))
}

#calculate markers for all clusters, and print heatmap of top 10 markers
allmarkers = FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
print(DoHeatmap(seu, features = top10$gene) + NoLegend() + ggtitle(paste0("DEGs")))
write.table(allmarkers,paste0(outputcsv),sep=",",quote=F,row.names=F,col.names=T)

dev.off()

saveRDS(seu,paste0(outputrds))

#if option is set, add a finer manual annotation of tcells, and save in rds object
add_manual_annotation = TRUE
pdf("reannotate_ribas_melanoma_merged_tcells_manual_annotation_label.pdf",width=7,height=7)
if (add_manual_annotation) {
  seu = readRDS("/data/reannotate_ribas_melanoma_merged_tcells.rds")
  seu$manual_annotation_label = ""
  seu$manual_annotation_label[seu$seurat_clusters %in% c(0,1,3,4,6,7,9,10,11,12,14)] = "CD8+ T-cells"
  seu$manual_annotation_label[seu$seurat_clusters %in% c(13)] = "NK cells"
  seu$manual_annotation_label[seu$seurat_clusters %in% c(8)] = "T-regs"
  seu$manual_annotation_label[seu$seurat_clusters %in% c(5)] = "TFH cells"
  seu$manual_annotation_label[seu$seurat_clusters %in% c(2)] = "Cycling T cells"
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "manual_annotation_label", repel = T, label.size = 5, shuffle = T) + ggtitle('manual_annotation_label') + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "orig.ident", repel = T, label.size = 3, shuffle = T) + ggtitle('orig.ident') + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))
}
dev.off()
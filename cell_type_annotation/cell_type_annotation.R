library(ggplot2)
library(Seurat)
library(stringr)
library(rlang)
library(rlist)
library(dplyr)

#load in BI5 scrna and snrna-seq samples, integrate, cluster, and store in rds file
integrate_BI5 = TRUE
if (integrate_BI5) {
  foldersList = c("s3://fresh-vs-frozen-comparison-ohio/BI5/scrna-seq",
    "s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq")
  integrated_name_arr = c("BI5_scrna-seq","BI5_snrna-seq")

  object.list = c()
  for (i in 1:2) {
    system(paste0("aws s3 cp ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
    integrated_rds = readRDS(paste0("/data/",integrated_name_arr[i],"_integrated.rds"))
    DefaultAssay(integrated_rds) = "RNA"
    object.list = c(object.list, integrated_rds)
  }

  anchors = FindIntegrationAnchors(object.list = object.list, dims = 1:20)
  seu = IntegrateData(anchorset = anchors, dims = 1:20)

  dimnum = 15
  seu = ScaleData(object = seu)
  seu = RunPCA(object = seu)
  seu = FindNeighbors(seu, dims = 1:15)
  seu = FindClusters(seu)
  seu = RunUMAP(object = seu, dims = 1:20)

  BI5_seu = seu
  saveRDS(BI5_seu,"/data/fresh_vs_frozen_all_reannotate_BI5.rds")
}

#load in ribas 310 samples, integrate, cluster, and store in rds file
integrate_ribas = TRUE
if (integrate_ribas) {
  ribas_samples = c("ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh","ribas_310_on_later_previd_3_GEX_titrate_thresh","ribas_310_pre_GEX_5pv2_S26_L004_titrate_thresh")
  object.list = list()
  for (z in 1:length(ribas_samples)) {
    ribas_sample = ribas_samples[z]
    system(paste0("aws s3 cp s3://melanoma-ribas/ribas1/Seurat/",ribas_sample,"/",ribas_sample,"_cb.rds /data/",ribas_sample,"_cb.rds"))
    temp_rds = readRDS(paste0("/data/",ribas_sample,"_cb.rds"))

    if (length(colnames(temp_rds))>=100)
    {
      object.list = c(object.list, temp_rds)
    }
    system(paste0("rm /data/",ribas_sample,"_cb.rds"))
  }

  anchors = FindIntegrationAnchors(object.list = object.list, dims = 1:20)
  seu = IntegrateData(anchorset = anchors, dims = 1:20)

  dimnum = 15
  seu = ScaleData(object = seu)
  seu = RunPCA(object = seu)
  seu = FindNeighbors(seu, dims = 1:15)
  seu = FindClusters(seu)
  seu = RunUMAP(object = seu, dims = 1:20)

  ribas_310_seu = seu
  saveRDS(ribas_310_seu,"/data/fresh_vs_frozen_all_reannotate_ribas_310.rds")
}

foldersList = c("",
  "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata",
  "s3://fresh-vs-frozen-comparison-ohio/nsclc",
  "",
  "s3://uveal-melanoma")
integrated_name_arr = c("BI5","cpoi-uvealprimarydata","nsclc","ribas_310","um_all")
calculateMarkers = FALSE

for (i in 1:length(foldersList)) {
  if (integrated_name_arr[i]=="BI5")
  {
    seu = readRDS("/data/fresh_vs_frozen_all_reannotate_BI5.rds")
  }
  else if (integrated_name_arr[i]=="ribas_310")
  {
    seu = readRDS("/data/fresh_vs_frozen_all_reannotate_ribas_310.rds")
  }
  else
  {
    system(paste0("aws s3 cp ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
    seu = readRDS(paste0("/data/",integrated_name_arr[i],"_integrated.rds"))
  }

  pdf(paste0("fresh_vs_frozen_all_reannotate/fresh_vs_frozen_all_reannotate_",integrated_name_arr[i],".pdf"))
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "seurat_clusters", repel = T, label.size = 3) + ggtitle(paste0('seurat_clusters')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "orig.ident", repel = T, label.size = 3) + ggtitle(paste0('orig.ident')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))

  if (integrated_name_arr[i]!="um_all")
  {
    seu$manual_annotation_label = "unknown"
  }
  if (integrated_name_arr[i]=="BI5")
  {
    seu$manual_annotation_label[seu$seurat_clusters %in% c(0,1,2,3,6,8,9,12,14)] = "Melanoma cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(7,11,13,18,20,23)] = "Myeloid"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(4,5,10,15,22,25)] = "T cells/NK cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(16)] = "B cells/Plasma cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(17)] = "Fibroblasts"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(19)] = "Neuronal cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(21)] = "Endothelial cells"
  }
  if (integrated_name_arr[i]=="cpoi-uvealprimarydata")
  {
    seu$manual_annotation_label[seu$seurat_clusters %in% c(0,1,2,4,5,6,8,10,16)] = "Melanoma cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(7,12,13,14,20)] = "Myeloid"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(3,9,11,14,15,17)] = "T cells/NK cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(18, 19)] = "Fibroblasts"
  }
  if (integrated_name_arr[i]=="nsclc")
  {
    seu$manual_annotation_label[seu$seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,20)] = "Lung cancer cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(10,15)] = "Myeloid"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(13)] = "T cells/NK cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(19)] = "NK cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(16,17)] = "B cells/Plasma cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(11,18)] = "Fibroblasts"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(12)] = "Endothelial cells"
  }
  if (integrated_name_arr[i]=="ribas"310")
  {
    seu$manual_annotation_label[seu$seurat_clusters %in% c(0,1,4,5,7,8,10,11,15)] = "Melanoma cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(6,12,16)] = "Myeloid" 
    seu$manual_annotation_label[seu$seurat_clusters %in% c(2,3,9,13)] = "T cells/NK cells"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(14)] = "Fibroblasts"
    seu$manual_annotation_label[seu$seurat_clusters %in% c(17,18)] = "Endothelial cells"
  }
  # old version for ribas_integrated_titrate_thresh
  # if (i==4)
  # {
  #   seu$manual_annotation_label[seu$seurat_clusters %in% c(0,1,2,4,5,6,9,11,12,13,16,17,21)] = "Melanoma cells"
  #   seu$manual_annotation_label[seu$seurat_clusters %in% c(3,20,24)] = "Myeloid"
  #   seu$manual_annotation_label[seu$seurat_clusters %in% c(7,8,10,15,19)] = "T cells/NK cells"
  #   seu$manual_annotation_label[seu$seurat_clusters %in% c(14,25)] = "Fibroblasts"
  #   seu$manual_annotation_label[seu$seurat_clusters %in% c(18)] = "Dendritic cells"
  #   seu$manual_annotation_label[seu$seurat_clusters %in% c(22)] = "Endothelial cells"
  # }
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "manual_annotation_label", repel = T, label.size = 3) + ggtitle(paste0('manual_annotation_label')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))

  DefaultAssay(seu) = "RNA"

  print(FeaturePlot(seu, features = c("nFeature_RNA","nCount_RNA","percent.mt","ScrubDoublet_score"), min.cutoff = "q9", max.cutoff = "q90"))

  #IG cells for B cells/plasma cells
  #Multiple NCAM1 XCL1 XCL2 TYROBP KLRD1 CD4 TCF4 COL4A1 KCNN3 DOCK2 PRDM1 SCT
  canonical_markers = c("CD3E",
  "NCAM1","XCL1","XCL2","TYROBP","KLRD1",
  "CD4",
  "CD8A", "CD8B",
  "FOXP3",
  "CD164", "VCAM1", "FCGR3A", "FCGR3B", "DOCK2",
  "PRDM1",
  "COL1A1",
  "PMEL", "MLANA", "MITF",
  "EPCAM",
  "TOP2A","MKI67",
  "GZMB",
  "VWF","ADAMTS9","PLVAP","FLT1","PECAM1","CALCRL","ADGRL4","CYYR1","TCF4","COL4A1","KCNN3",
  "NTM","DCLK1","SLC18A2","SAMD12","CUX2","KLHL13",
  "CD38","MZB1",
  "CLNK","CPVL","CLEC4C","SCT","C1QA","C1QB","C1QC","CD14","CX3CR1",
  "DCC","RGS7","DACH1","OOEP",
  "MACC1",
  "ELMO1")
  for (j in seq(1,length(canonical_markers),4))
  {
    markers_to_plot = canonical_markers[j:min(j+3,length(canonical_markers))]
    print(FeaturePlot(seu, features = c(markers_to_plot), min.cutoff = "q9", max.cutoff = "q90"))
  }
  DefaultAssay(seu) = "integrated"

  if (calculateMarkers)
  {
    allmarkers = FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    write.table(allmarkers,paste0("fresh_vs_frozen_all_reannotate/fresh_vs_frozen_all_reannotate_",integrated_name_arr[i],"_markers.csv"),sep=",",quote=F,row.names=F,col.names=T)
  }

  dev.off()

  seu$barebarcodes = unlist(lapply(strsplit(colnames(seu),"_"), function(x) x[1]))
  writeout_df = data.frame(barcode=seu$barebarcodes, orig.ident=seu$orig.ident, manual_annotation_label=seu$manual_annotation_label)
  if (integrated_name_arr[i]=="um_all")
  {
    writeout_df$manual_annotation_label[writeout_df$manual_annotation_label %in% c("Neuronal","Ribosomal","Ribosomal/Mitochondrial","Tumour","Tumour/Ribosomal")] = "Melanoma cells"
    writeout_df$manual_annotation_label[writeout_df$manual_annotation_label %in% c("Cell Cycle")] = "Cycling melanoma cells"
  }
  write.table(writeout_df,paste0(integrated_name_arr[i],"_manual_annotation_label.csv"),sep=",",col.names=T,row.names=F,quote=F)

  if (i==1)
  {
    saveRDS(seu,"/data/fresh_vs_frozen_all_reannotate_BI5.rds")
  }
  else if (i==4)
  {
    saveRDS(seu,"/data/fresh_vs_frozen_all_reannotate_ribas_310.rds")
  }
  else
  {
    saveRDS(seu, paste0("/data/",integrated_name_arr[i],"_integrated.rds"))
    system(paste0("aws s3 cp /data/",integrated_name_arr[i],"_integrated.rds ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds"))
  }
  system(paste0("rm /data/",integrated_name_arr[i],"_integrated.rds"))
}
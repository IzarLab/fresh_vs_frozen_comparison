library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(grid)
library(stringr)

#load uveal melanoma data, subset to t-cells, and save to rds file
seu = readRDS("/data/um_all_integrated.rds")
seu <- subset(seu, manual_annotation_label %in% c("T-cells"))

seu = ScaleData(object = seu)
seu = RunPCA(object = seu)
seu = FindNeighbors(seu, dims = 1:15)
seu = FindClusters(seu)
seu = RunUMAP(object = seu, dims = 1:20)

pdf("reannotate_uveal_melanoma_tcells.pdf")
print(DimPlot(seu, reduction = "umap", label = T, group.by = "integrated_snn_res.0.8", repel = T, label.size = 3) + ggtitle('celltype_bped_main') + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))
dev.off()

saveRDS(seu,"/data/reannotate_uveal_melanoma_tcells.rds")

#load uveal melanoma t-cells object, exclude cluster 7 which is unlikely to be t-cells
recluster_rds = readRDS("/data/reannotate_uveal_melanoma_tcells.rds")
recluster_df = data.frame(barebarcode = unlist(lapply(colnames(recluster_rds), function(x) {str_split(x,"_")[[1]][1]}))[recluster_rds$manual_annotation_label=="T-cells" & recluster_rds$RNA_snn_res.0.8!=7], orig.ident = recluster_rds$orig.ident[recluster_rds$manual_annotation_label=="T-cells" & recluster_rds$RNA_snn_res.0.8!=7])

pat_list = c("um_07_gk_on_S8_L001","um_07_gk_pre_S4_L001","uv003-uvme-snseq-3p-post","um_08_ar_on_S2_L001","um_08_ar_post_S3_L001","um_08_ar_pre_S1_L001","um_09_mw_on_S6_L001","um_09_mw_post_S7_L001","um_09_mw_pre_S5_L001","um_11_lc_on_S16_L002","um_11_lc_pre_S12_L002","um_12_ml_on_S10_L002","um_12_ml_post_S11_L002","um_12_ml_pre_S9_L002","um_15_lm_on_S14_L002","um_15_lm_post_S15_L002","um_15_lm_pre_S13_L002","um_16_rs_on_S18_L003","um_16_rs_post_S19_L003","um_16_rs_pre_S17_L003")

object.list = NULL
for(pat in pat_list){
  #load in rds files for each individual samples, subset to tcells, and add fields for mait and inkt evidence based on TR_frequency.csv file
  #store results in object.list array
  if (pat!="uv003-uvme-snseq-3p-post")
  {
    system(paste0("aws s3 cp s3://uveal-melanoma/Seurat/",pat,"/",pat,"_cb.rds /data/",pat,"_cb.rds"))
    temp_rds = readRDS(paste0('/data/',pat,'_cb.rds'))
    relevant_ids = recluster_df$barebarcode[recluster_df$orig.ident==pat]
    temp_rds$barcode = colnames(temp_rds)
    temp_rds = subset(temp_rds, barcode %in% relevant_ids)

    temp_rds$mait_evidence = ""
    temp_rds$inkt_evidence = ""
    csv_table = read.table(paste0(pat,"_TR_frequency.csv"),sep=",",header=T,quote=NULL)
    if (sum(names(csv_table)=="mait_evidence")!=0)
    {
      for (j in 1:length(csv_table$barcodes))
      {
	barcodes = unique(str_split(csv_table$barcodes[j],"\\|")[[1]])
	mait_evidence = csv_table$mait_evidence[j]
	inkt_evidence = csv_table$inkt_evidence[j]
	if (!is.na(mait_evidence))
	{
	  temp_rds$mait_evidence[colnames(temp_rds) %in% barcodes] = mait_evidence
	}
	if (!is.na(inkt_evidence))
	{
	  temp_rds$inkt_evidence[colnames(temp_rds) %in% barcodes] = inkt_evidence
	}
      }
    }

    if (length(colnames(temp_rds))>=100)
    {
      object.list = c(object.list, temp_rds)
    }
    system(paste0("rm /data/",pat,"_cb.rds"))
  }
}

#integrate samples in object.list array, using dimensions given in dim_nums array
dim_nums = c(25)
pdf("reannotate_uveal_melanoma_tcells_reintegrated.pdf",width=10,height=12)
for (dim_num in dim_nums)
{
  anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:dim_num)

  seu <- IntegrateData(anchorset = anchors, dims = 1:dim_num)

  seu <- ScaleData(object = seu)
  seu <- RunPCA(object = seu)
  seu <- FindNeighbors(seu, dims = 1:dim_num)
  seu <- FindClusters(seu)
  seu <- RunUMAP(object = seu, dims = 1:dim_num)

  #if number of dimensions = 25, exclude clusters 6 and 12 from the original clustering, then reintegrate with number of dimensions = 15
  if (dim_num==25)
  {
    seu$placeholder = TRUE
    seu$placeholder[seu$seurat_clusters %in% c(6,12)] = FALSE
    seu = subset(seu, placeholder)

    # seu <- ScaleData(object = seu)
    # seu <- RunPCA(object = seu)
    # seu <- FindNeighbors(seu, dims = 1:15)
    # seu <- FindClusters(seu)
    # seu <- RunUMAP(object = seu, dims = 1:15)

    recluster_df = data.frame(barebarcode = unlist(lapply(colnames(seu), function(x) {str_split(x,"_")[[1]][1]})), orig.ident = seu$orig.ident)

    object.list = NULL
    for(pat in unique(seu$orig.ident)){
      if (pat!="uv003-uvme-snseq-3p-post")
      {
    	system(paste0("aws s3 cp s3://uveal-melanoma/Seurat/",pat,"/",pat,"_cb.rds /data/",pat,"_cb.rds"))
    	temp_rds = readRDS(paste0('/data/',pat,'_cb.rds'))
    	relevant_ids = recluster_df$barebarcode[recluster_df$orig.ident==pat]
    	temp_rds$barcode = colnames(temp_rds)
    	temp_rds = subset(temp_rds, barcode %in% relevant_ids)

    	temp_rds$mait_evidence = ""
    	temp_rds$inkt_evidence = ""
    	csv_table = read.table(paste0(pat,"_TR_frequency.csv"),sep=",",header=T,quote=NULL)
        if (sum(names(csv_table)=="mait_evidence")!=0)
        {
	  for (j in 1:length(csv_table$barcodes))
	  {
	    barcodes = unique(str_split(csv_table$barcodes[j],"\\|")[[1]])
	    mait_evidence = csv_table$mait_evidence[j]
	    inkt_evidence = csv_table$inkt_evidence[j]
	    if (!is.na(mait_evidence))
	    {
	      temp_rds$mait_evidence[colnames(temp_rds) %in% barcodes] = mait_evidence
	    }
	    if (!is.na(inkt_evidence))
	    {
	      temp_rds$inkt_evidence[colnames(temp_rds) %in% barcodes] = inkt_evidence
	    }
	  }
	}

    	if (length(colnames(temp_rds))>=100)
    	{
    	  object.list = c(object.list, temp_rds)
    	}
    	system(paste0("rm /data/",pat,"_cb.rds"))
      }
    }
    
    anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:15)

    seu <- IntegrateData(anchorset = anchors, dims = 1:15)

    seu <- ScaleData(object = seu)
    seu <- RunPCA(object = seu)
    seu <- FindNeighbors(seu, dims = 1:15)
    seu <- FindClusters(seu)
    seu <- RunUMAP(object = seu, dims = 1:15)
  }

  dim_num_label = "25_then_15"

  #from here on, assume that samples were first reintegrated with 25 dimensions, then clusters 6 and 12 were removed, and then samples were reintegrated with 15 dimensions
  #print umaps of seurat clusters and sample name
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "integrated_snn_res.0.8", repel = T, label.size = 3) + ggtitle(paste0('dim_num ',dim_num_label,' seurat_cluster')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "orig.ident", repel = T, label.size = 3) + ggtitle(paste0('dim_num ',dim_num_label,' orig.ident')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))

  #print ummaps of marker genes
  DefaultAssay(seu) = "RNA"
  canonical_markers = c("CD4", "CD8A", "CD8B", "FOXP3", "CXCL13", "NCAM1", "TOP2A", "TOX")
  for (marker in canonical_markers)
  {
    print(FeaturePlot(seu, features = c(marker), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0("dim_num ",dim_num_label," ",marker)))
  }
  DefaultAssay(seu) = "integrated"

  #print umaps highlighting tcells with evidence of either MAIT or INKT status
  highlight_list = list(colnames(seu)[seu$mait_evidence!=""])
  names(highlight_list) = c("mait_evidence")
  print(DimPlot(seu, reduction = "umap", label = T, pt.size = 0.3, cells.highlight = highlight_list, cols.highlight = rev(rainbow(length(highlight_list))), sizes.highlight = 3) + ggtitle(paste0("dim_num ",dim_num_label," MAIT evidence")))
  highlight_list = list(colnames(seu)[seu$inkt_evidence!=""])
  names(highlight_list) = c("inkt_evidence")
  print(DimPlot(seu, reduction = "umap", label = T, pt.size = 0.3, cells.highlight = highlight_list, cols.highlight = rev(rainbow(length(highlight_list))), sizes.highlight = 3) + ggtitle(paste0("dim_num ",dim_num_label," INKT evidence")))

  #find marker genes for all markers, print heatmap of top 10 markers
  allmarkers = FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  print(DoHeatmap(seu, features = top10$gene) + NoLegend() + ggtitle(paste0("dim_num ",dim_num_label," DEG Heatmap")))
  write.table(allmarkers,paste0("reannotate_uveal_melanoma_tcells_reintegrated_dim_num_",dim_num,"_15_markers.csv"),sep=",",quote=F,row.names=F,col.names=T)

  saveRDS(seu,paste0("/data/reannotate_uveal_melanoma_tcells_reintegrated_dim_num_",dim_num,"_15.rds"))
}
dev.off()

#select only cluster 1 of reintegrated data, recluster with dimensions = 20
dim_num_label = "20"
pdf("reannotate_uveal_melanoma_tcells_reintegrated_cluster1.pdf",width=10,height=12)
seu_cluster1 = subset(seu, seurat_clusters==1)
seu_cluster1 <- ScaleData(object = seu_cluster1)
seu_cluster1 <- RunPCA(object = seu_cluster1)
seu_cluster1 <- FindNeighbors(seu_cluster1, dims = 1:20)
seu_cluster1 <- FindClusters(seu_cluster1, resolution = 1.2)
seu_cluster1 <- RunUMAP(object = seu_cluster1, dims = 1:20)

print(DimPlot(seu_cluster1, reduction = "umap", label = T, group.by = "seurat_clusters", repel = T, label.size = 3) + ggtitle(paste0('dim_num ',dim_num_label,' seurat_cluster1_cluster')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))
print(DimPlot(seu_cluster1, reduction = "umap", label = T, group.by = "orig.ident", repel = T, label.size = 3) + ggtitle(paste0('dim_num ',dim_num_label,' orig.ident')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))

DefaultAssay(seu_cluster1) = "RNA"
canonical_markers = c("CD4", "CD8A", "CD8B", "FOXP3", "CXCL13", "NCAM1", "TOP2A", "TOX")
for (marker in canonical_markers)
{
  print(FeaturePlot(seu_cluster1, features = c(marker), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0("dim_num ",dim_num_label," ",marker)))
}
DefaultAssay(seu_cluster1) = "integrated"

allmarkers = FindAllMarkers(seu_cluster1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
print(DoHeatmap(seu_cluster1, features = top10$gene) + NoLegend() + ggtitle(paste0("dim_num ",dim_num_label," DEG Heatmap")))
write.table(allmarkers,paste0("reannotate_uveal_melanoma_tcells_reintegrated_dim_num_",dim_num,"_15_markers_cluster1.csv"),sep=",",quote=F,row.names=F,col.names=T)
dev.off()

#relabel tcells in cluster 1 of original reintegrated data as either CD4 or CD8 tcells, based on subclustering of cluster 1 alone
seu$seurat_clusters_with_1_subclustered = as.character(seu$seurat_clusters)
for (i in 1:length(seu_cluster1$barcode)) {
  select_arr = (seu$barcode==seu_cluster1$barcode[i] & seu$orig.ident==seu_cluster1$orig.ident[i])
  if (seu_cluster1$seurat_clusters[i] %in% c(2,4,6,7,8,9))
  {
    seu$seurat_clusters_with_1_subclustered[select_arr] = "1_CD4"
  }
  else if (seu_cluster1$seurat_clusters[i] %in% c(0,1,3,5))
  {
    seu$seurat_clusters_with_1_subclustered[select_arr] = "1_CD8"
  }
}

dim_num_label = "25_then_15"
pdf("reannotate_uveal_melanoma_tcells_reintegrated_with_1_subclustered.pdf",width=7,height=7)

reloadprevRDS = TRUE
if (reloadprevRDS) {
  seu = readRDS(paste0("/data/reannotate_uveal_melanoma_tcells_reintegrated_with_1_subclustered_dim_num_",dim_num_label,".rds"))
}

#assign manual annotations of tcell status in final reintegrated data
seu$manual_annotation_label = ""
seu$manual_annotation_label[seu$seurat_clusters_with_1_subclustered %in% c("1_CD4","2","3","10")] = "CD4+ T-cells"
seu$manual_annotation_label[seu$seurat_clusters_with_1_subclustered %in% c("7")] = "T-regs"
seu$manual_annotation_label[seu$seurat_clusters_with_1_subclustered %in% c("0","1_CD8","4","6")] = "CD8+ T-cells"
seu$manual_annotation_label[seu$seurat_clusters_with_1_subclustered %in% c("8")] = "CD4+ T-cells TIGIT-high"
seu$manual_annotation_label[seu$seurat_clusters_with_1_subclustered %in% c("9")] = "NK cells NCAM-low"
seu$manual_annotation_label[seu$seurat_clusters_with_1_subclustered %in% c("5")] = "NK cells NCAM-high"

#print umaps of manual annotation and sample name
print(DimPlot(seu, reduction = "umap", label = T, group.by = "seurat_clusters_with_1_subclustered", repel = T, label.size = 3) + ggtitle(paste0('dim_num ',dim_num_label,' seurat_cluster1_cluster')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))
print(DimPlot(seu, reduction = "umap", label = T, group.by = "manual_annotation_label", repel = T, label.size = 5) + ggtitle(paste0('manual_annotation_label')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))
print(DimPlot(seu, reduction = "umap", label = T, group.by = "orig.ident", repel = T, label.size = 3) + ggtitle(paste0('orig.ident')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))

#print umaps of marker gene expression
DefaultAssay(seu) = "RNA"
canonical_markers = c("CD4", "CD8A", "CD8B", "FOXP3", "CXCL13", "NCAM1", "TOP2A", "TOX")
for (marker in canonical_markers)
{
  print(FeaturePlot(seu, features = c(marker), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0("dim_num ",dim_num_label," ",marker)))
}
DefaultAssay(seu) = "integrated"

#write table of marker genes, and save rds file containing reintegrated object
Idents(seu) = seu$seurat_clusters_with_1_subclustered
allmarkers = FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
print(DoHeatmap(seu, features = top10$gene) + NoLegend() + ggtitle(paste0("dim_num ",dim_num_label," DEG Heatmap")))
dev.off()

write.table(allmarkers,paste0("reannotate_uveal_melanoma_tcells_reintegrated_with_1_subclustered_dim_num_",dim_num_label,"_markers.csv"),sep=",",quote=F,row.names=F,col.names=T)

saveRDS(seu,paste0("/data/reannotate_uveal_melanoma_tcells_reintegrated_with_1_subclustered_dim_num_",dim_num_label,".rds"))
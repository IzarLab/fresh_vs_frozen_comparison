#!/usr/bin/env Rscript

### Seurat analysis for cellbender output provided as an argument
## needs to run in /home/ubuntu

print(paste('Start:',Sys.time()))

library(dplyr)
library(Seurat)
library(Matrix)
library(gplots)
library(ggplot2)
library(purrr)
library(cowplot)
#library(rscrublet)
library(DropletUtils)
library(SingleR)
library(SingleCellExperiment)
library(scater)
library(pheatmap)

pats = c("ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh")
for (z in 1:length(pats))
{
  pat = pats[z]
  patcellbender = "ribas_310_on_GEX_5pv2_S27_L004"
  doubs = .096
  maxCount = 60000
  maxFeature = 9000
  minCount = 100
  minFeature = 100
  maxMT = 20
  system(paste0("aws s3 cp s3://melanoma-ribas/ribas1/cellbender/",patcellbender,"/",patcellbender,"_filtered.h5 /data/",patcellbender,"_filtered.h5"))
  seu.data <- Read10X_h5(paste0("/data/",patcellbender,'_filtered.h5'), use.names = TRUE, unique.features = TRUE)

  system(paste0("aws s3 cp s3://melanoma-ribas/ribas1/Seurat/",pats[z],"/",pats[z],"_cb.rds /data/",pats[z],"_cb.rds"))
  seu <- readRDS(paste0("/data/",pat,'_cb.rds'))

  ### Initialize the Seurat object with the raw (non-normalized data)
  seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, min.cells = 1, min.features = 1)

  # Annotate MT genes
  seu_raw[["percent.mt"]] <- PercentageFeatureSet(seu_raw, pattern = "^MT-")

  # Annotate pat and facs info
  seu_raw[["patient"]]<-pat

  # identify doublets using scrublet
  matrixfilename = paste0(pat,"_MM.txt")
  scrubletoutputfile = paste0(pat,"_doublets.txt")
  writeMM(seu_raw@assays$RNA@counts,file=matrixfilename)
  system(paste0("python3 scrublet_code.py ",matrixfilename," ",doubs," ",scrubletoutputfile))
  scrubletoutput = read.table(scrubletoutputfile,header=T,quote=NULL)
  system(paste0("rm ",matrixfilename))
  system(paste0("rm ",scrubletoutputfile))
  #store doublet status, and score for strength of doublet prediction, in two new variables in seu_raw
  seu_raw[['ScrubDoublet']]<-scrubletoutput$predicted_doublets
  seu_raw[['ScrubDoublet_score']]<-scrubletoutput$doublet_scores

  stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 13))
  colnames(stats)<-c('sample','n_raw_features','n_raw_cells','n_ScrubDoublet','n_features','n_cells','median_features','median_counts','cutoff_features','cutoff_counts','cutoff_mt','virus_count_raw','virus_count')
  rownames(stats)<-pat
  stats$sample<-pat
  stats$n_raw_features<-dim(seu_raw@assays$RNA@counts)[1]
  stats$n_raw_cells<-dim(seu_raw@assays$RNA@counts)[2]
  stats$n_ScrubDoublet<-length(which(seu_raw@meta.data$ScrubDoublet ==T))
  stats$n_features<-dim(seu@assays$RNA@counts)[1]
  stats$n_cells<-dim(seu@assays$RNA@counts)[2]
  stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
  stats$median_counts<-round(median(seu@meta.data$nCount_RNA))
  stats$cutoff_features<-paste(minFeature,maxFeature)
  stats$cutoff_counts<-paste(minCount,maxCount)
  stats$cutoff_mt<-paste(maxMT)
  stats$virus_count_raw<-ifelse('SARS-COV-2' %in% rownames(seu_raw@assays$RNA@data),seu_raw@assays$RNA@counts['SARS-COV-2',] %>% sum,0)
  stats$virus_count<-ifelse('SARS-COV-2' %in% rownames(seu@assays$RNA@data),seu@assays$RNA@counts['SARS-COV-2',] %>% sum,0)

  ### Save objects
  system(paste0("mkdir /data/",pat))
  #saveRDS(seu, file = paste0("/data/",pat,"_titrate_thresh/",pat,'_titrate_thresh_cb.rds'))
  #save(seu, seu_raw,pred_bped_fine, file = paste0("data/",pat,"/data_",pat,"_full_cb.RData"))

  ### write pdf reports
  pdf(file = paste0("/data/",pat,"/plots_", pat,"_cb.pdf"))

  # stats
  textplot(t(stats),cex=1.2,halign='left')

  # plots raw data
  print(qplot(x=seu_raw@meta.data$nCount_RNA,y = seu_raw@meta.data$nFeature_RNA, col=seu_raw@meta.data$percent.mt, xlab = "nCount_RNA",
       ylab = "nFeature_RNA", main =paste0(seu_raw@meta.data$orig.ident[1], " raw data: nCount_RNA vs. nFeature_RNA")) + scale_colour_gradient(low="blue", high="green") + labs(color = "Percent MT") + theme_classic())
  print(VlnPlot(seu_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,split.by=NULL))

  # QC plot filtered
  print(qplot(x=seu@meta.data$nCount_RNA,y = seu@meta.data$nFeature_RNA, col=seu@meta.data$percent.mt, xlab = "nCount_RNA",
       ylab = "nFeature_RNA", main =paste0(seu@meta.data$orig.ident[1], " filtered: nCount_RNA vs. nFeature_RNA")) + scale_colour_gradient(low="blue", high="green") + labs(color = 'Percent MT') + theme_classic())
  print(VlnPlot(seu, features = 'nFeature_RNA'))
  print(VlnPlot(seu, features = "nCount_RNA"))
  print(VlnPlot(seu, features = "percent.mt"))

  # Determine metrics to plot present in seurat_control@meta.data
  metrics <-  c("nCount_RNA", "nFeature_RNA", "percent.mt")
  # Extract the UMAP coordinates for each cell and include information about the metrics to plot
  qc_data <- FetchData(seu, vars = c(metrics, "ident", "UMAP_1", "UMAP_2"))
  # Adding cluster label to center of cluster on UMAP
  umap_label <- FetchData(seu, vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
	  as.data.frame() %>% group_by(ident) %>% summarise(x=mean(UMAP_1), y=mean(UMAP_2))
  # Plot a UMAP plot for each metric
  map(metrics, function(qc){
	  ggplot(qc_data, aes(UMAP_1, UMAP_2)) +
			  geom_point(aes_string(color=qc), size=0.5) +
			  scale_color_gradientn(colours = rev(rainbow(4))) +theme_classic() +
			  geom_text(data=umap_label, aes(label=ident, x, y)) +
			  ggtitle(qc)+
			  theme(legend.text=element_text(size=6),legend.title=element_text(size=7),legend.key.width=unit(0.2,"cm"))
  }) %>% plot_grid(plotlist = .)

  # PCA
  print(DimPlot(seu, reduction = "pca"))
  print(DimHeatmap(seu, dims = 1:9, cells = 500, balanced = TRUE))
  print(JackStrawPlot(seu, dims = 1:20))
  print(ElbowPlot(seu))

  # UMAP
  print(DimPlot(seu, reduction = "umap",label = T,group.by = 'ident'))

  # features
  # if('SARS_COV_2' %in% rownames(seu@assays$RNA@data)){print(FeaturePlot(seu, features = 'SARS_COV_2', min.cutoff = "q9"))}
  # if('PTPRC' %in% rownames(seu@assays$RNA@data)){print(FeaturePlot(seu, features = 'PTPRC', min.cutoff = "q9"))}

  # if("CD4" %in% rownames(seu@assays$RNA@data) |
  #    "CD8A" %in% rownames(seu@assays$RNA@data) | 
  #    'IL7R'%in% rownames(seu@assays$RNA@data) | 
  #    "MKI67" %in% rownames(seu@assays$RNA@data)) {
  #   print(FeaturePlot(seu, features = c('CD4',"CD8A", 'IL7R', "MKI67"), min.cutoff = "q9"))}
  # if("CD68" %in% rownames(seu@assays$RNA@data) |
  #    "CD14" %in% rownames(seu@assays$RNA@data) | 
  #   'ITGAX'%in% rownames(seu@assays$RNA@data) | 
  #   "NKG7" %in% rownames(seu@assays$RNA@data)) {
  #     print(FeaturePlot(seu, features = c('CD68', # macrophage
  #                                             "CD14", # monocyte
  #                                             'ITGAX', # dendritic
  #                                             "NKG7"  # NK cell
  #                                             ), min.cutoff = "q9"))}

  ## bped


  print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine', repel = T) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=4)) + ggtitle("celltype_bped_fine"))

  dev.off()

  system(paste0("aws s3 sync /data/",pat," s3://melanoma-ribas/ribas1/Seurat/",pat))
  system(paste0("rm /data/",pat,"/",pat,"_cb.rds"))
  system(paste0("rm /data/",pat,"_filtered.h5"))

  print(paste('End:',Sys.time()))
}
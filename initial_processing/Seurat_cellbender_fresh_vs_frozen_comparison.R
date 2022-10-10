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
library(grid)

#setwd("/home/ubuntu/")

# foldersList = c("s3://fresh-vs-frozen-comparison/BI5/scrna-seq",
#   "s3://fresh-vs-frozen-comparison/BI5/snrna-seq",
#   "s3://fresh-vs-frozen-comparison/cpoi-uvealprimarydata",
#   "s3://fresh-vs-frozen-comparison/nsclc",
#   "s3://fresh-vs-frozen-comparison/sarcoma-sn")

# patsList = list(c("CD45negGEXBI5_S1_L001","CD45posGEXBI5_S1_L001"),
#   c("bi005-skcm-5snseq","bi005-skcm","skcm-bi005-5pv2-snseq"),
#   c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-E12","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-F1","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-GEX-F12"),
#   c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX","NSCL-NR001-5pv2-snseq"),
#   c("Sarcoma167GEX","Sarcoma322GEX","Sarcoma559GEX","Sarcoma708GEX"))

# maxCountsList = list(c(45000,25000),
#   c(30000,30000,16000),
#   c(50000,20000,25000),
#   c(50000,10000,7500,30000,40000,15000),
#   c(50000,25000,25000,50000))
# maxFeaturesList = list(c(6500,4500),
#   c(7500,7000,5000),
#   c(7500,4500,6000),
#   c(7000,4500,4000,7000,9000,6000),
#   c(8500,9000,6000,8000))

# downsampled parameters
foldersList = c("s3://fresh-vs-frozen-comparison-ohio/BI5/scrna-seq",
  "s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq",
  "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata",
  "s3://fresh-vs-frozen-comparison-ohio/nsclc",
  "s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")

patsList = list(c("CD45negGEXBI5_S1_L001","CD45posGEXBI5_S1_L001"),
  c("bi005-skcm-5snseq","bi005-skcm","skcm-bi005-5pv2-snseq"),
  c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-E12","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-F1","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-GEX-F12"),
  c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX","NSCL-NR001-5pv2-snseq"),
  c("BI5CST","BI5TST","NR1CST","NR1TST"))

maxCountsList = list(c(40000,18000),c(18000,10000,15000),c(45000,10000,15000),c(15000,13000,9000,12000,12000,13000),c(10000,10000,1000,800))
maxFeaturesList = list(c(6500,4000),c(6000,4500,6000),c(7000,4000,6000),c(4000,5000,4000,5000,5000,5000),c(4500,4500,700,600))

loadPrevRDS = FALSE
use_downsampled = TRUE

# for (z1 in 1:length(foldersList))
# {
#   pats = patsList[[z1]]
#   for (z in 1:length(pats))
#   {
#     if (use_downsampled)
#     {
#       system(paste0("aws s3 cp ",foldersList[z1],"/cellbender_downsampled/",pats[z],"/",pats[z],"_filtered.h5 /data/",pats[z],"_filtered.h5"))
#     }
#     else
#     {
#       system(paste0("aws s3 cp ",foldersList[z1],"/cellbender/",pats[z],"/",pats[z],"_filtered.h5 /data/",pats[z],"_filtered.h5"))
#     }
#     #pat<-commandArgs()[6]
#     #doubs<-as.numeric(commandArgs()[7])
#     pat = pats[z]#"uv001_uvmel_snseq_3p_pre"
#     doubs = .096
#     if (foldersList[z1]=="s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
#     {
#       doubs = .0047
#     }
#     print(pat)
#     print(doubs) #0.096 or 0.04

#     ##### Loading, merging, QC, dimension reduction #####
#     ### Load dataset
#     seu.data <- Read10X_h5(paste0('/data/',pat,'_filtered.h5'), use.names = TRUE, unique.features = TRUE)

#     ### Initialize the Seurat object with the raw (non-normalized data)
#     seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, min.cells = 1, min.features = 1)

#     # Annotate MT genes
#     seu_raw[["percent.mt"]] <- PercentageFeatureSet(seu_raw, pattern = "^MT-")

#     # Annotate pat and facs info
#     seu_raw[["patient"]]<-pat
#     #seu_raw[["organ"]]<-strsplit(pat,'_')[[1]][1]
#     #seu_raw[["group"]]<-strsplit(pat,'_')[[1]][3]
#     #seu_raw[["ID"]]<-strsplit(pat,'_')[[1]][2]
#     #seu_raw[["sequencing"]]<-'sn'

#     # identify doublets using scrublet
#     matrixfilename = paste0("/data/",pat,"_MM.txt")
#     scrubletoutputfile = paste0("/data/",pat,"_doublets.txt")
#     writeMM(seu_raw@assays$RNA@counts,file=matrixfilename)
#     system(paste0("python3 /mnt/vdb/home/ubuntu2/scrublet_code.py ",matrixfilename," ",doubs," ",scrubletoutputfile))
#     scrubletoutput = read.table(scrubletoutputfile,header=T,quote=NULL)
#     system(paste0("rm ",matrixfilename))
#     system(paste0("rm ",scrubletoutputfile))
#     #store doublet status, and score for strength of doublet prediction, in two new variables in seu_raw
#     seu_raw[['ScrubDoublet']]<-scrubletoutput$predicted_doublets
#     seu_raw[['ScrubDoublet_score']]<-scrubletoutput$doublet_scores

#     if (use_downsampled)
#     {
#       pdf(file = paste0("/data/plots_", pat,"_downsampled_raw_check.pdf"))
#     }
#     else
#     {
#       pdf(file = paste0("/data/plots_", pat,"_raw_check.pdf"))
#     }

#     print(qplot(x=seu_raw@meta.data$nCount_RNA,y = seu_raw@meta.data$nFeature_RNA, col=seu_raw@meta.data$percent.mt, xlab = "nCount_RNA",
# 	 ylab = "nFeature_RNA", main =paste0(seu_raw@meta.data$orig.ident[1], " raw data: nCount_RNA vs. nFeature_RNA")) + scale_colour_gradient(low="blue", high="green") + labs(color = "Percent MT") + theme_classic())
#     print(VlnPlot(seu_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,split.by=NULL))

#     dev.off()

#     system(paste0("rm /data/",pats[z],"_filtered.h5"))
#   }
# }

# nonsense = nonsense+1

justQualityClusters = FALSE

if (justQualityClusters) {
  pdf("fresh_vs_frozen_comparison_quality_clusters.pdf",width=14,height=126)
  pushViewport(viewport(layout = grid.layout(18,2)))
  vplayout = function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  rowcounter = 1
}

for (z1 in 1:length(foldersList))
{
  pats = patsList[[z1]]
  maxCounts = maxCountsList[[z1]]
  maxFeatures = maxFeaturesList[[z1]]
  for (z in 1:length(pats))
  {
    if (use_downsampled)
    {
      system(paste0("aws s3 cp ",foldersList[z1],"/cellbender_downsampled/",pats[z],"/",pats[z],"_filtered.h5 /data/",pats[z],"_filtered.h5"))
    }
    else
    {
      system(paste0("aws s3 cp ",foldersList[z1],"/cellbender/",pats[z],"/",pats[z],"_filtered.h5 /data/",pats[z],"_filtered.h5"))
    }
    #pat<-commandArgs()[6]
    #doubs<-as.numeric(commandArgs()[7])
    pat = pats[z]#"uv001_uvmel_snseq_3p_pre"
    doubs = .096
    if (foldersList[z1]=="s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
    {
      doubs = .0047
    }
    print(pat)
    print(doubs) #0.096 or 0.04

    ##### Loading, merging, QC, dimension reduction #####
    ### Load dataset
    seu.data <- Read10X_h5(paste0('/data/',pat,'_filtered.h5'), use.names = TRUE, unique.features = TRUE)

    ### Initialize the Seurat object with the raw (non-normalized data)
    seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, min.cells = 1, min.features = 1)

    # Annotate MT genes
    seu_raw[["percent.mt"]] <- PercentageFeatureSet(seu_raw, pattern = "^MT-")

    # Annotate pat and facs info
    seu_raw[["patient"]]<-pat
    #seu_raw[["organ"]]<-strsplit(pat,'_')[[1]][1]
    #seu_raw[["group"]]<-strsplit(pat,'_')[[1]][3]
    #seu_raw[["ID"]]<-strsplit(pat,'_')[[1]][2]
    #seu_raw[["sequencing"]]<-'sn'

    # identify doublets using scrublet
    matrixfilename = paste0("/data/",pat,"_MM.txt")
    scrubletoutputfile = paste0("/data/",pat,"_doublets.txt")
    writeMM(seu_raw@assays$RNA@counts,file=matrixfilename)
    system(paste0("python3 /mnt/vdb/home/ubuntu2/scrublet_code.py ",matrixfilename," ",doubs," ",scrubletoutputfile))
    scrubletoutput = read.table(scrubletoutputfile,header=T,quote=NULL)
    system(paste0("rm ",matrixfilename))
    system(paste0("rm ",scrubletoutputfile))
    #store doublet status, and score for strength of doublet prediction, in two new variables in seu_raw
    seu_raw[['ScrubDoublet']]<-scrubletoutput$predicted_doublets
    seu_raw[['ScrubDoublet_score']]<-scrubletoutput$doublet_scores

    if (loadPrevRDS)
    {
      system(paste0("aws s3 cp ",foldersList[z1],"/Seurat/",pats[z],"/",pats[z],"_cb.rds /data/",pats[z],"_cb.rds"))
      seu = readRDS(paste0("/data/",pats[z],"_cb.rds"))
      system(paste0("rm /data/",pats[z],"_cb.rds"))
    }
    else
    {
      ### subset 
      #Mesothelioma_12
      #minFeature<-200
      #maxFeature<- 2200#7500
      #minCount<- 100
      #maxCount<- 2800#40000
      #maxMT<-10
      #Mesothelioma_73
      minFeature = 300
      maxFeature = 1000000#maxFeatures[z]#2000
      minCount = 0#400
      maxCount = 1000000#maxCounts[z]#3000
      maxMT = 20
      seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & nCount_RNA > minCount & nCount_RNA < maxCount & percent.mt < maxMT & ScrubDoublet ==F)

      ### Workflow
      seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
      seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
      seu <- ScaleData(seu, features = rownames(seu))
      seu <- RunPCA(seu)
      seu <- JackStraw(seu, num.replicate = 100)
      seu <- ScoreJackStraw(seu, reduction = "pca",dims = 1:20)
      seu <- FindNeighbors(seu, dims = 1:15)
      seu <- FindClusters(seu)
      seu <- RunUMAP(seu, dims = 1:20)

      if (!justQualityClusters)
      {
	# ## Finding differentially expressed features (cluster biomarkers)
	# seu.markers <- FindAllMarkers(seu, only.pos = F, min.pct = 0.2, logfc.threshold = 0.25)
	# marker_table <- seu.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
	# write.csv(marker_table,paste0('markers_',pat,'_cellbender.csv'),row.names = F)

	# ### cellt type identification
	# seu_sce <- as.SingleCellExperiment(seu)

	# bped<-BlueprintEncodeData()
	# pred_bped_main <- SingleR(test = seu_sce, ref = bped, labels = bped$label.main)
	# pruneScores(pred_bped_main)
	# seu[['celltype_bped_main']]<-pred_bped_main$pruned.labels
	# pred_bped_fine <- SingleR(test = seu_sce, ref = bped, labels = bped$label.fine)
	# pruneScores(pred_bped_fine)
	# seu[['celltype_bped_fine']]<-pred_bped_fine$pruned.labels

	# iced<-DatabaseImmuneCellExpressionData()
	# pred_iced_main <- SingleR(test = seu_sce, ref = iced, labels = iced$label.main)
	# pruneScores(pred_iced_main)
	# seu[['celltype_iced_main']]<-pred_iced_main$pruned.labels
	# pred_iced_fine <- SingleR(test = seu_sce, ref = iced, labels = iced$label.fine)
	# pruneScores(pred_iced_fine)
	# seu[['celltype_iced_fine']]<-pred_iced_fine$pruned.labels

	# hpca<-HumanPrimaryCellAtlasData()
	# pred_hpca_main <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.main)
	# pruneScores(pred_hpca_main)
	# seu[['celltype_hpca_main']]<-pred_hpca_main$pruned.labels
	# print("HEREb")
	# pred_hpca_fine <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.fine)
	# print("HERE")
	# pruneScores(pred_hpca_fine)
	# print("HERE2")
	# seu[['celltype_hpca_fine']]<-pred_hpca_fine$pruned.labels
	# print("HERE3")

	# mid<-MonacoImmuneData()
	# pred_mid_main <- SingleR(test = seu_sce, ref = mid, labels = mid$label.main)
	# pruneScores(pred_mid_main)
	# seu[['celltype_mid_main']]<-pred_mid_main$pruned.labels
	# pred_mid_fine <- SingleR(test = seu_sce, ref = mid, labels = mid$label.fine)
	# pruneScores(pred_mid_fine)
	# seu[['celltype_mid_fine']]<-pred_mid_fine$pruned.labels
      }
    }

    ### stats
    stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 13))
    colnames(stats)<-c('sample','n_raw_features','n_raw_cells','n_ScrubDoublet','n_features','n_cells','median_features','median_counts','cutoff_features','cutoff_counts','cutoff_mt','virus_count_raw','virus_count')
    if (substr(pat,1,4)=="UMEL")
    {
      rownames(stats)<-substr(pat,12,length(pat))
      stats$sample<-substr(pat,12,length(pat))
    }
    else
    {
      rownames(stats)<-pat
      stats$sample<-pat
    }
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

    if (justQualityClusters)
    {
      seu_raw <- NormalizeData(seu_raw, normalization.method = "LogNormalize", scale.factor = 10000)
      seu_raw <- FindVariableFeatures(seu_raw, selection.method = "vst", nfeatures = 2000)
      seu_raw <- ScaleData(seu_raw, features = rownames(seu_raw))
      seu_raw <- RunPCA(seu_raw)
      seu_raw <- JackStraw(seu_raw, num.replicate = 100)
      seu_raw <- ScoreJackStraw(seu_raw, reduction = "pca",dims = 1:20)
      seu_raw <- FindNeighbors(seu_raw, dims = 1:15)
      seu_raw <- FindClusters(seu_raw)
      seu_raw <- RunUMAP(seu_raw, dims = 1:20)

      print(VlnPlot(seu_raw, features = 'nFeature_RNA', log = T) + ggtitle(paste0(pat," raw nFeature_RNA")),vp=vplayout(rowcounter,1))
      #print(VlnPlot(seu_raw, features = "nCount_RNA", log = T) + ggtitle(paste0(pat," raw nCount_RNA")),vp=vplayout(rowcounter,2))
      #print(VlnPlot(seu_raw, features = "percent.mt", log = T) + ggtitle(paste0(pat," raw percent.mt")),vp=vplayout(rowcounter,3))
      print(VlnPlot(seu, features = 'nFeature_RNA', log = T) + ggtitle(paste0(pat," with cutoffs nFeature_RNA")),vp=vplayout(rowcounter,2))
      #print(VlnPlot(seu, features = "nCount_RNA", log = T) + ggtitle(paste0(pat," with cutoffs nCount_RNA")),vp=vplayout(rowcounter,5))
      #print(VlnPlot(seu, features = "percent.mt", log = T) + ggtitle(paste0(pat," with cutoffs percent.mt")),vp=vplayout(rowcounter,6))

      rowcounter = rowcounter+1
    }
    else
    {
      ### Save objects
      system(paste0("mkdir /data/",pat,"_final_thresh"))
      saveRDS(seu, file = paste0("/data/",pat,"_final_thresh/",pat,'_final_thresh_cb.rds'))
      #save(seu, seu_raw,pred_bped_fine, file = paste0("data/",pat,"/data_",pat,"_full_cb.RData"))

      ### write pdf reports
      pdf(file = paste0("/data/",pat,"_final_thresh/plots_", pat,"_final_thresh_cb.pdf"))

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

      print("scores")

      # ## bped
      # plotScoreHeatmap(pred_bped_fine, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_bped_fine')
      # DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,label.size = 2.5) + 
      # 	ggtitle('celltype_bped_fine') +
      # 	guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
      # 	theme(legend.text=element_text(size=6))

      # ## bped
      # plotScoreHeatmap(pred_bped_main, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_bped_main')
      # DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_main',repel = T,label.size = 2.5) + 
      # 	ggtitle('celltype_bped_main') +
      # 	guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
      # 	theme(legend.text=element_text(size=6))

      # ## hpca
      # plotScoreHeatmap(pred_hpca_fine, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_hpca_fine')
      # DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_fine',repel = T,label.size = 2.5) + 
      # 	ggtitle('celltype_hpca_fine') +
      # 	guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
      # 	theme(legend.text=element_text(size=6))

      # ## hpca
      # plotScoreHeatmap(pred_hpca_main, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_hpca_main')
      # DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,label.size = 2.5) + 
      # 	ggtitle('celltype_hpca_main') +
      # 	guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
      # 	theme(legend.text=element_text(size=6))

      # print("dims")

      # print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine', repel = T) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=4)) + ggtitle("celltype_bped_fine"))
      # print("1")
      # print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_main', repel = T) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=4)) + ggtitle("celltype_bped_main"))
      # print("2")
      # #print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_fine', repel = T) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=4)) + ggtitle("celltype_hpca_fine"))
      # print("3")
      # print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main', repel = T) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=4)) + ggtitle("celltype_hpca_main"))
      # print("4")

      # print("off")

      dev.off()

      if (use_downsampled)
      {
        system(paste0("aws s3 sync /data/",pat,"_final_thresh ",foldersList[z1],"/Seurat_downsampled/",pat,"_final_thresh"))
      }
      else
      {
        system(paste0("aws s3 sync /data/",pat,"_final_thresh ",foldersList[z1],"/Seurat/",pat,"_final_thresh"))
      }
      system(paste0("rm /data/",pat,"_final_thresh/",pat,"_final_thresh_cb.rds"))
      system(paste0("rm /data/",pat,"_filtered.h5"))

      print(paste('End:',Sys.time()))
    }
  }
}

if (justQualityClusters) {
  dev.off()
}
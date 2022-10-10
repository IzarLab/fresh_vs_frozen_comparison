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

#setwd("/home/ubuntu/")

foldersList = c("s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")

patsList = list(c("BI5CST","BI5TST","NR1CST","NR1TST"))
#patsList = list(c("NR1CST","NR1TST"))

maxCountsList = list(c(70000,40000,8000,2000))
maxFeaturesList = list(c(11000,8000,3000,3500))
#maxCountsList = list(c(8000,2000))
#maxFeaturesList = list(c(3000,3500))

# for (z1 in 1:length(foldersList))
# {
#   pats = patsList[[z1]]
#   for (z in 1:length(pats))
#   {
#     system(paste0("aws s3 cp ",foldersList[z1],"/cellbender/",pats[z],"/",pats[z],"_filtered.h5 /data/",pats[z],"_filtered.h5"))
#     #pat<-commandArgs()[6]
#     #doubs<-as.numeric(commandArgs()[7])
#     pat = pats[z]#"uv001_uvmel_snseq_3p_pre"
#     doubs = .0047
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

#     pdf(file = paste0("/data/plots_", pat,"_raw_check.pdf"))

#     print(qplot(x=seu_raw@meta.data$nCount_RNA,y = seu_raw@meta.data$nFeature_RNA, col=seu_raw@meta.data$percent.mt, xlab = "nCount_RNA",
# 	 ylab = "nFeature_RNA", main =paste0(seu_raw@meta.data$orig.ident[1], " raw data: nCount_RNA vs. nFeature_RNA")) + scale_colour_gradient(low="blue", high="green") + labs(color = "Percent MT") + theme_classic())
#     print(VlnPlot(seu_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,split.by=NULL))

#     dev.off()

#     system(paste0("rm /data/",pats[z],"_filtered.h5"))
#   }
# }

# nonsense = nonsense+1

for (z1 in 1:length(foldersList))
{
  pats = patsList[[z1]]
  maxCounts = maxCountsList[[z1]]
  maxFeatures = maxFeaturesList[[z1]]
  for (z in 1:length(pats))
  {
    system(paste0("aws s3 cp ",foldersList[z1],"/cellbender/",pats[z],"/",pats[z],"_filtered.h5 /data/",pats[z],"_filtered.h5"))
    #pat<-commandArgs()[6]
    #doubs<-as.numeric(commandArgs()[7])
    pat = pats[z]#"uv001_uvmel_snseq_3p_pre"
    doubs = .0047
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

    ### subset 
    #Mesothelioma_12
    #minFeature<-200
    #maxFeature<- 2200#7500
    #minCount<- 100
    #maxCount<- 2800#40000
    #maxMT<-10
    #Mesothelioma_73
    minFeature = 100
    maxFeature = maxFeatures[z]#2000
    minCount = 200
    maxCount = maxCounts[z]#3000
    maxMT = 10
    seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & nCount_RNA > minCount & nCount_RNA < maxCount & percent.mt < maxMT & ScrubDoublet ==F)

    ### Workflow
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    seu <- ScaleData(seu, features = rownames(seu))
    seu <- RunPCA(seu, npcs=min(50,dim(seu)[2]-1))
    seu <- JackStraw(seu, num.replicate = 100)
    seu <- ScoreJackStraw(seu, reduction = "pca",dims = 1:min(20,dim(seu)[2]-1))
    seu <- FindNeighbors(seu, dims = 1:min(15,dim(seu)[2]-1))
    seu <- FindClusters(seu)
    seu <- RunUMAP(seu, dims = 1:min(20,dim(seu)[2]-1))

    # ## Finding differentially expressed features (cluster biomarkers)
    # seu.markers <- FindAllMarkers(seu, only.pos = F, min.pct = 0.2, logfc.threshold = 0.25)
    # marker_table <- seu.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
    # write.csv(marker_table,paste0('markers_',pat,'_cellbender.csv'),row.names = F)

    ### cellt type identification
    seu_sce <- as.SingleCellExperiment(seu)

    bped<-BlueprintEncodeData()
    pred_bped_main <- SingleR(test = seu_sce, ref = bped, labels = bped$label.main)
    pruneScores(pred_bped_main)
    seu[['celltype_bped_main']]<-pred_bped_main$pruned.labels
    pred_bped_fine <- SingleR(test = seu_sce, ref = bped, labels = bped$label.fine)
    pruneScores(pred_bped_fine)
    seu[['celltype_bped_fine']]<-pred_bped_fine$pruned.labels

    iced<-DatabaseImmuneCellExpressionData()
    pred_iced_main <- SingleR(test = seu_sce, ref = iced, labels = iced$label.main)
    pruneScores(pred_iced_main)
    seu[['celltype_iced_main']]<-pred_iced_main$pruned.labels
    pred_iced_fine <- SingleR(test = seu_sce, ref = iced, labels = iced$label.fine)
    pruneScores(pred_iced_fine)
    seu[['celltype_iced_fine']]<-pred_iced_fine$pruned.labels

    hpca<-HumanPrimaryCellAtlasData()
    pred_hpca_main <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.main)
    pruneScores(pred_hpca_main)
    seu[['celltype_hpca_main']]<-pred_hpca_main$pruned.labels
    print("HEREb")
    pred_hpca_fine <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.fine)
    print("HERE")
    pruneScores(pred_hpca_fine)
    print("HERE2")
    seu[['celltype_hpca_fine']]<-pred_hpca_fine$pruned.labels
    print("HERE3")

    mid<-MonacoImmuneData()
    pred_mid_main <- SingleR(test = seu_sce, ref = mid, labels = mid$label.main)
    pruneScores(pred_mid_main)
    seu[['celltype_mid_main']]<-pred_mid_main$pruned.labels
    pred_mid_fine <- SingleR(test = seu_sce, ref = mid, labels = mid$label.fine)
    pruneScores(pred_mid_fine)
    seu[['celltype_mid_fine']]<-pred_mid_fine$pruned.labels

    ### stats
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
    saveRDS(seu, file = paste0("/data/",pat,"/",pat,'_cb.rds'))
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
    for(metric in metrics)
    {
      print(FeaturePlot(seu, features = c(metric), min.cutoff = "q9", max.cutoff = "q90"))
    }

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

    ## bped
    plotScoreHeatmap(pred_bped_fine, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_bped_fine')
    DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,label.size = 2.5) + 
      ggtitle('celltype_bped_fine') +
      guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
      theme(legend.text=element_text(size=6))

    ## bped
    plotScoreHeatmap(pred_bped_main, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_bped_main')
    DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_main',repel = T,label.size = 2.5) + 
      ggtitle('celltype_bped_main') +
      guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
      theme(legend.text=element_text(size=6))

    ## hpca
    plotScoreHeatmap(pred_hpca_fine, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_hpca_fine')
    DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_fine',repel = T,label.size = 2.5) + 
      ggtitle('celltype_hpca_fine') +
      guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
      theme(legend.text=element_text(size=6))

    ## hpca
    plotScoreHeatmap(pred_hpca_main, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_hpca_main')
    DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,label.size = 2.5) + 
      ggtitle('celltype_hpca_main') +
      guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
      theme(legend.text=element_text(size=6))

    print("dims")

    print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine', repel = T) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=4)) + ggtitle("celltype_bped_fine"))
    print("1")
    print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_main', repel = T) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=4)) + ggtitle("celltype_bped_main"))
    print("2")
    #print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_fine', repel = T) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=4)) + ggtitle("celltype_hpca_fine"))
    print("3")
    print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main', repel = T) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=4)) + ggtitle("celltype_hpca_main"))
    print("4")

    print("off")

    dev.off()

    system(paste0("aws s3 sync /data/",pat," ",foldersList[z1],"/Seurat/",pat))
    system(paste0("rm /data/",pat,"/",pat,"_cb.rds"))
    system(paste0("rm /data/",pat,"_filtered.h5"))

    print(paste('End:',Sys.time()))
  }
}
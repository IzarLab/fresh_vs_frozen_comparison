#!/usr/bin/env Rscript
print(Sys.time())

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(purrr)
library(cowplot)
library(stringr)


outputnamearr = c("um_all_merged")#,"um_on_integrated","um_post_integrated","um_pre_integrated")
pat_list_arr = list(c("um_07_gk_on_S8_L001","um_07_gk_pre_S4_L001","um_08_ar_on_S2_L001","um_08_ar_post_S3_L001","um_08_ar_pre_S1_L001","um_09_mw_on_S6_L001","um_09_mw_post_S7_L001","um_09_mw_pre_S5_L001","um_11_lc_on_S16_L002","um_11_lc_pre_S12_L002","um_12_ml_on_S10_L002","um_12_ml_post_S11_L002","um_12_ml_pre_S9_L002","um_15_lm_on_S14_L002","um_15_lm_post_S15_L002","um_15_lm_pre_S13_L002","um_16_rs_on_S18_L003","um_16_rs_post_S19_L003","um_16_rs_pre_S17_L003"))#,c("um_07_gk_on_S8_L001","um_08_ar_on_S2_L001","um_09_mw_on_S6_L001","um_11_lc_on_S16_L002","um_12_ml_on_S10_L002","um_15_lm_on_S14_L002","um_16_rs_on_S18_L003"),c("uv003_uvme_snseq_3p_post","um_08_ar_post_S3_L001","um_09_mw_post_S7_L001","um_12_ml_post_S11_L002","um_15_lm_post_S15_L002","um_16_rs_post_S19_L003"),c("um_07_gk_pre_S4_L001","um_08_ar_pre_S1_L001","um_09_mw_pre_S5_L001","um_11_lc_pre_S12_L002","um_12_ml_pre_S9_L002","um_15_lm_pre_S13_L002","um_16_rs_pre_S17_L003"))
#outputnamearr = c("um_07_gk","um_08_ar","um_09_mw","um_11_lc")#c("um_12_ml","um_15_lm","um_16_rs")
#pat_list_arr = list(c("um_07_gk_on_S8_L001","um_07_gk_pre_S4_L001","uv003_uvme_snseq_3p_post"),c("um_08_ar_on_S2_L001","um_08_ar_post_S3_L001","um_08_ar_pre_S1_L001"),c("um_09_mw_on_S6_L001","um_09_mw_post_S7_L001","um_09_mw_pre_S5_L001"),c("um_11_lc_on_S16_L002","um_11_lc_pre_S12_L002"))#list(c("um_12_ml_on_S10_L002","um_12_ml_post_S11_L002","um_12_ml_pre_S9_L002"),c("um_15_lm_on_S14_L002","um_15_lm_post_S15_L002","um_15_lm_pre_S13_L002"),c("um_16_rs_on_S18_L003","um_16_rs_post_S19_L003","um_16_rs_pre_S17_L003"))

for (z in 1:length(pat_list_arr))
{
  pat_list = pat_list_arr[[z]]
  outputname = outputnamearr[z]

  # load objects
  object.list<-NULL
  for(pat in pat_list){
    system(paste0("aws s3 cp s3://uveal-melanoma/Seurat/",pat,"/",pat,"_cb.rds /data/",pat,"_cb.rds"))
    temp_obj = readRDS(paste0('/data/',pat,'_cb.rds'))
    object.list<-c(object.list,temp_obj)
    system(paste0("rm /data/",pat,"_cb.rds"))
  }

  if (length(grep("merged",outputname))==1)
  {
    source("merge.SCTAssay.R")
    seu = merge.SCTAssay(x=object.list[[1]], y=object.list[[2]], merge.data = T, na.rm = T)
    object.list[[1]] = NULL
    for (z1 in 3:length(object.list))
    {
      seu = merge.SCTAssay(x=seu, y=object.list[[z1]], merge.data = T, na.rm = T)
    }
  }
  else
  {
    # find anchors
    anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:20)

    # integrate data sets
    seu <- IntegrateData(anchorset = anchors, dims = 1:20)
  }

  object.list = list()

  # normal workflow
  seu <- ScaleData(object = seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- RunPCA(object = seu)
  seu <- FindNeighbors(seu, dims = 1:15)
  seu <- FindClusters(seu)
  seu <- RunUMAP(object = seu, dims = 1:20)

  DefaultAssay(seu) <- "RNA"
  seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
  DefaultAssay(seu) <- "integrated"

  # markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) %>% write.csv(paste0('/data/markers_',outputname,'.csv'),row.names = F)

  # save object
  # saveRDS(seu, file = paste0('/data/',outputname,'.rds'))

  # ### stats
  stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
  colnames(stats)<-c('sample','n_features','n_cells','median_features','median_counts')
  rownames(stats)<-'Lungs_integrated'
  stats$sample<-'lungs_all'
  stats$n_features<-dim(seu@assays$integrated@data)[1]
  stats$n_cells<-dim(seu@assays$integrated@data)[2]
  stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
  stats$median_counts<-round(median(seu@meta.data$nCount_RNA))

  pdf(file = paste0("/data/plots_",outputname,".pdf"))
  textplot(t(stats),cex=1.2,halign='left')
  print(DimPlot(seu, reduction = "pca",group.by = 'patient'))
  print(DimPlot(seu, reduction = "pca",group.by = 'Phase'))
  #print(DimPlot(seu, reduction = "pca",group.by = 'group'))
  print(DimHeatmap(seu, dims = 1:9, cells = 500, balanced = TRUE))
  print(DimHeatmap(seu, dims = 10:18, cells = 500, balanced = TRUE))
  print(DimPlot(seu, reduction = "umap",label = T,group.by = 'ident'))
  print(DimPlot(seu, reduction = "umap",label = F,group.by = 'patient'))
  print(DimPlot(seu, reduction = "umap",label = F,group.by = 'Phase'))
  #print(DimPlot(seu, reduction = "umap",label = F,group.by = 'group'))

  DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_main',repel = T,label.size = 3) + 
    ggtitle('celltype_bped_main') +
    guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
    theme(legend.text=element_text(size=10))

  DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,label.size = 3) + 
    ggtitle('celltype_bped_fine') +
    guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
    theme(legend.text=element_text(size=7))

  DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,label.size = 2.5) + 
    ggtitle('celltype_hpca_main') +
    guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
    theme(legend.text=element_text(size=7))

  # distribution
  cti<-data.frame(seu@meta.data$patient,seu@meta.data$celltype_bped_main)
  ggplot(cti, aes(x=seu.meta.data.patient, fill=seu.meta.data.celltype_bped_main))+
    geom_bar(colour="black",position="fill",size=0.25)+ ggtitle('Fractions of cell types in patients (bped)') + 
    theme_classic() +xlab("Patient") + ylab("Fraction (%)")+ labs(fill = "Cell type")+
    guides(fill = guide_legend(ncol = 1)) +
    theme(text = element_text(size=8),axis.text.x=element_text(angle=45,hjust=1),legend.key.size = unit(0.3, "cm"))

  cti<-data.frame(seu@meta.data$patient,seu@meta.data$celltype_hpca_main)
  ggplot(cti, aes(x=seu.meta.data.patient, fill=seu.meta.data.celltype_hpca_main))+
    geom_bar(colour="black",position="fill",size=0.25)+ ggtitle('Fractions of cell types in patients (hpca)') + 
    theme_classic() +xlab("Patient") + ylab("Fraction (%)")+ labs(fill = "Cell type")+
    guides(fill = guide_legend(ncol = 1)) +
    theme(text = element_text(size=8),axis.text.x=element_text(angle=45,hjust=1),legend.key.size = unit(0.3, "cm"))

  print(FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt', "ScrubDoublet_score"), min.cutoff = "q9",max.cutoff = "q90"))

  dev.off()
}
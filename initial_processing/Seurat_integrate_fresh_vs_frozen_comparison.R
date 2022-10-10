#!/usr/bin/env Rscript
print(Sys.time())

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(purrr)
library(cowplot)
library(stringr)
#library(celldex)

# foldersList = c("s3://fresh-vs-frozen-comparison-ohio/BI5/scrna-seq",
#   "s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq",
#   "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata",
#   "s3://fresh-vs-frozen-comparison-ohio/nsclc",
#   "s3://fresh-vs-frozen-comparison-ohio/sarcoma-sn")
# pat_list_arr = list(c("CD45negGEXBI5_S1_L001","CD45posGEXBI5_S1_L001"),
#   c("bi005-skcm-5snseq","bi005-skcm","skcm-bi005-5pv2-snseq"),
#   c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-E12","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-F1","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-GEX-F12"),
#   c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX","NSCL-NR001-5pv2-snseq"),
#   c("Sarcoma167GEX","Sarcoma322GEX","Sarcoma559GEX","Sarcoma708GEX"))
# outputnamearr = c("BI5_scrna-seq","BI5_snrna-seq","cpoi-uvealprimarydata","nsclc","sarcoma-sn")
# short_pat_list_arr = list(c("CD45neg","CD45pos"),
#   c("5snseq","3snseq","5pv2-snseq"),
#   c("SCRNA-5P-NA-E12","SCRNA-5P-NA-F1","SNRNA-5P-WI-F12"),
#   c("SCRNA_5P_NA","SNSEQ_3P_NI","SNSEQ_3P_WI","SNSEQ_5P_NI","SNSEQ_5P_WI","5pv2-snseq"),
#   c("Sarcoma167","Sarcoma322","Sarcoma559","Sarcoma708"))
# fresh_frozen_list_arr = list(c("fresh","fresh"),
#   c("frozen","frozen","frozen"),
#   c("fresh","fresh","frozen"),
#   c("fresh","frozen","frozen","frozen","frozen","5pv2-snseq"),
#   c("frozen","frozen","frozen","frozen"))
# outputnamearr = c("BI5_scrna-seq","BI5_snrna-seq","cpoi-uvealprimarydata","nsclc","sarcoma-sn")

# foldersList = c("s3://fresh-vs-frozen-comparison-ohio/nsclc")
# pat_list_arr = list(c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX","NSCL-NR001-5pv2-snseq"))
# outputnamearr = c("nsclc")
# short_pat_list_arr = list(c("SCRNA_5P_NA","SNSEQ_3P_NI","SNSEQ_3P_WI","SNSEQ_5P_NI","SNSEQ_5P_WI","5pv2-snseq"))
# fresh_frozen_list_arr = list(c("fresh","frozen","frozen","frozen","frozen","5pv2-snseq"))
# outputnamearr = c("nsclc")

# foldersList = c("s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
# pat_list_arr = list(c("BI5CST","BI5TST","NR1CST","NR1TST"))
# short_pat_list_arr = list(c("BI5CST","BI5TST","NR1CST","NR1TST"))
# fresh_frozen_list_arr = list(c("frozen","frozen","frozen","frozen"))
# outputnamearr = c("slyper_pipeline")

# foldersList = c("s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline","s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
# pat_list_arr = list(c("BI5CST","BI5TST"),c("NR1CST","NR1TST"))
# short_pat_list_arr = list(c("BI5CST","BI5TST"),c("NR1CST","NR1TST"))
# fresh_frozen_list_arr = list(c("frozen","frozen"),c("frozen","frozen"))
# outputnamearr = c("BI5","NR1")

use_downsampled = TRUE

# parameters for downsampled integration, excludes sarcoma samples
foldersList = c("s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline","s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
pat_list_arr = list(c("BI5CST","BI5TST"),c("NR1CST","NR1TST"))
short_pat_list_arr = list(c("BI5CST","BI5TST"),c("NR1CST","NR1TST"))
fresh_frozen_list_arr = list(c("frozen","frozen"),c("frozen","frozen"))
outputnamearr = c("BI5","NR1")

# foldersList = c("s3://fresh-vs-frozen-comparison-ohio/BI5/scrna-seq",
#   "s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq",
#   "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata",
#   "s3://fresh-vs-frozen-comparison-ohio/nsclc")
# pat_list_arr = list(c("CD45negGEXBI5_S1_L001","CD45posGEXBI5_S1_L001"),
#   c("bi005-skcm-5snseq","bi005-skcm","skcm-bi005-5pv2-snseq"),
#   c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-E12","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-F1","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-GEX-F12"),
#   c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX","NSCL-NR001-5pv2-snseq"))
# outputnamearr = c("BI5_scrna-seq","BI5_snrna-seq","cpoi-uvealprimarydata","nsclc")
# short_pat_list_arr = list(c("CD45neg","CD45pos"),
#   c("5snseq","3snseq","5pv2-snseq"),
#   c("SCRNA-5P-NA-E12","SCRNA-5P-NA-F1","SNRNA-5P-WI-F12"),
#   c("SCRNA_5P_NA","SNSEQ_3P_NI","SNSEQ_3P_WI","SNSEQ_5P_NI","SNSEQ_5P_WI","5pv2-snseq"))
# fresh_frozen_list_arr = list(c("fresh","fresh"),
#   c("frozen","frozen","frozen"),
#   c("fresh","fresh","frozen"),
#   c("fresh","frozen","frozen","frozen","frozen","5pv2-snseq"))
# outputnamearr = c("BI5_scrna-seq","BI5_snrna-seq","cpoi-uvealprimarydata","nsclc")

for (z in 1:length(pat_list_arr))
{
  pat_list = pat_list_arr[[z]]
  short_pat_list = short_pat_list_arr[[z]]
  fresh_frozen_list = fresh_frozen_list_arr[[z]]
  outputname = outputnamearr[z]

  # load objects
  for(i in 1:length(pat_list)){
    pat = pat_list[i]
    cbrds_name = paste0(pat,"_final_thresh_cb.rds")
    cbfolder_name = paste0(pat,"_final_thresh")
    if ((outputname=="slyper" || outputname=="BI5" || outputname=="NR1") && !(use_downsampled))
    {
      cbrds_name = paste0(pat,"_cb.rds")
      cbfolder_name = pat
    }
    if (use_downsampled)
    {
      system(paste0("aws s3 cp ",foldersList[z],"/Seurat_downsampled/",cbfolder_name,"/",cbrds_name," /data/",cbrds_name))
    }
    else
    {
      system(paste0("aws s3 cp ",foldersList[z],"/Seurat/",cbfolder_name,"/",cbrds_name," /data/",cbrds_name))
    }
    patunderscore = str_replace_all(pat,'-','_')
    assign(patunderscore,readRDS(paste0('/data/',cbrds_name)))
    eval(parse(text = paste0(patunderscore,"$orig.ident = \"",short_pat_list[i],"\"")))
    eval(parse(text = paste0(patunderscore,"$patient = \"",short_pat_list[i],"\"")))
    eval(parse(text = paste0(patunderscore,"$fresh_frozen = \"",fresh_frozen_list[i],"\"")))
  }

  # create object.list
  object.list<-NULL
  for(obj in pat_list){
    objunderscore = str_replace_all(obj,'-','_')
    object.list<-c(object.list,eval(parse(text = objunderscore)))
  }

  if (outputname=="BI5" || outputname=="NR1")
  {
    for (objidx in 1:length(object.list))
    {
      anobj = object.list[[objidx]]
      anobj$manual_annotation_label = "unknown"
      object.list[[objidx]] = anobj
    }
  }

  if (outputname=="BI5")
  {
    #system("aws s3 cp s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq/Seurat/bi005-skcm_final_thresh/bi005-skcm_final_thresh_cb.rds /data/bi005-skcm_final_thresh_cb.rds")
    if (use_downsampled)
    {
      system("aws s3 cp s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq/Seurat_downsampled/bi005-skcm_final_thresh/bi005-skcm_final_thresh_cb.rds /data/bi005-skcm_final_thresh_cb.rds")
      bi5obj = readRDS("/data/bi005-skcm_final_thresh_cb.rds")
      system("rm /data/bi005-skcm_final_thresh_cb.rds")
    }
    else
    {
      bi5obj = readRDS("/data/fresh_vs_frozen_all_reannotate_BI5.rds")
      bi5obj = subset(bi5obj, orig.ident=="3snseq")
    }
    object.list = c(object.list,bi5obj)
  }

  if (outputname=="NR1")
  {
    #system("aws s3 cp s3://fresh-vs-frozen-comparison-ohio/nsclc/Seurat/NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX_final_thresh/NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX_final_thresh_cb.rds /data/NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX_final_thresh_cb.rds")
    if (use_downsampled)
    {
      system("aws s3 cp s3://fresh-vs-frozen-comparison-ohio/nsclc/Seurat_downsampled/integrated/nsclc_integrated.rds /data/nsclc_integrated.rds")
    }
    else
    {
      system("aws s3 cp s3://fresh-vs-frozen-comparison-ohio/nsclc/Seurat/integrated/nsclc_integrated.rds /data/nsclc_integrated.rds")
    }
    nsclcobj = readRDS("/data/nsclc_integrated.rds")
    system("rm /data/nsclc_integrated.rds")
    object.list = c(object.list,subset(nsclcobj, orig.ident=="SNSEQ_3P_NI"))
    object.list = c(object.list,subset(nsclcobj, orig.ident=="SNSEQ_3P_WI"))
  }

  # find anchors
  anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:20)

  # integrate data sets
  seu <- IntegrateData(anchorset = anchors, dims = 1:20)

  # normal workflow
  seu <- ScaleData(object = seu)
  seu <- RunPCA(object = seu)
  seu <- FindNeighbors(seu, dims = 1:15)
  seu <- FindClusters(seu)
  seu <- RunUMAP(object = seu, dims = 1:20)

  DefaultAssay(seu) <- "RNA"
  #seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
  DefaultAssay(seu) <- "integrated"

  #markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) %>% write.csv(paste0('/data/markers_',outputname,'.csv'),row.names = F)

  # save object
  saveRDS(seu, file = paste0('/data/',outputname,'_integrated.rds'))

  ### stats
  stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
  colnames(stats)<-c('sample','n_features','n_cells','median_features','median_counts')
  rownames(stats)<-outputname
  stats$sample<-outputname
  stats$n_features<-dim(seu@assays$integrated@data)[1]
  stats$n_cells<-dim(seu@assays$integrated@data)[2]
  stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
  stats$median_counts<-round(median(seu@meta.data$nCount_RNA))

  pdf(file = paste0("/data/plots_",outputname,"_integrated.pdf"))
  textplot(t(stats),cex=1.2,halign='left')
  print(DimPlot(seu, reduction = "pca",group.by = 'patient'))
  print(DimPlot(seu, reduction = "pca",group.by = 'fresh_frozen'))
  print(DimPlot(seu, reduction = "umap",group.by = 'patient'))
  print(DimPlot(seu, reduction = "umap",group.by = 'fresh_frozen'))
  #print(DimPlot(seu, reduction = "pca",group.by = 'Phase'))
  #print(DimPlot(seu, reduction = "pca",group.by = 'group'))
  print(DimHeatmap(seu, dims = 1:9, cells = 500, balanced = TRUE))
  print(DimHeatmap(seu, dims = 10:18, cells = 500, balanced = TRUE))
  print(DimPlot(seu, reduction = "umap",label = T,group.by = 'ident'))
  print(DimPlot(seu, reduction = "umap",label = F,group.by = 'patient'))
  #print(DimPlot(seu, reduction = "umap",label = F,group.by = 'Phase'))
  #print(DimPlot(seu, reduction = "umap",label = F,group.by = 'group'))

  if (!use_downsampled)
  {
    print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_main',repel = T,label.size = 3) + 
      ggtitle('celltype_bped_main') +
      guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
      theme(legend.text=element_text(size=10)))

    print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,label.size = 3) + 
      ggtitle('celltype_bped_fine') +
      guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
      theme(legend.text=element_text(size=7)))

    print(DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,label.size = 2.5) + 
      ggtitle('celltype_hpca_main') +
      guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
      theme(legend.text=element_text(size=7)))

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
  }

  print(FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt', "ScrubDoublet_score"), min.cutoff = "q9",max.cutoff = "q90"))

  dev.off()

  if (use_downsampled)
  {
    system(paste0("aws s3 cp /data/plots_",outputname,"_integrated.pdf ",foldersList[z],"/Seurat_downsampled/integrated/plots_",outputname,"_integrated.pdf"))
    system(paste0("aws s3 cp /data/",outputname,"_integrated.rds ",foldersList[z],"/Seurat_downsampled/integrated/",outputname,"_integrated.rds"))
  }
  else
  {
    system(paste0("aws s3 cp /data/plots_",outputname,"_integrated.pdf ",foldersList[z],"/Seurat/integrated/plots_",outputname,"_integrated.pdf"))
    system(paste0("aws s3 cp /data/",outputname,"_integrated.rds ",foldersList[z],"/Seurat/integrated/",outputname,"_integrated.rds"))
  }
  for (pat in pat_list) {
    system(paste0("rm /data/",cbrds_name))
  }
  system(paste0("rm /data/",outputname,"_integrated.rds"))
}
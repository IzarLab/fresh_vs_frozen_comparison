library(ggplot2)
library(Seurat)
library(stringr)
library(rlang)
library(rlist)
library(dplyr)

foldersList = c("",
  "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata",
  "s3://fresh-vs-frozen-comparison-ohio/nsclc",
  "s3://melanoma-ribas/ribas1",
  "s3://uveal-melanoma")
integrated_name_arr = c("BI5","cpoi-uvealprimarydata","nsclc","ribas_integrated_titrate_thresh","um_all")
pdf("fresh_vs_frozen_all_reannotate/fresh_vs_frozen_all_reannotate_heatmap.pdf",height=25,width=25)
for (i in 1:length(integrated_name_arr)) {
  if (i==1)
  {
    seu = readRDS("/data/fresh_vs_frozen_all_reannotate_BI5.rds")
  }
  else
  {
    system(paste0("aws s3 cp ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
    seu = readRDS(paste0("/data/",integrated_name_arr[i],"_integrated.rds"))
  }
  
  DefaultAssay(seu)<-'RNA'
  seu<-NormalizeData(seu) %>% ScaleData()
  DefaultAssay(seu)<-'integrated'
  
  allmarkers = read.table(paste0("fresh_vs_frozen_all_reannotate/fresh_vs_frozen_all_reannotate_",integrated_name_arr[i],"_markers.csv"), header=T, quote=NULL, sep=",")
  if (i==3)
  {
    allmarkers = allmarkers[allmarkers$cluster!=0,]
    seu = subset(seu, seurat_clusters!=0)
  }
  else (i==4 || i==5)
  {
    allmarkers = allmarkers[!(allmarkers$cluster %in% c(0,1,2)),]
    seu$placeholder = T
    seu$placeholder[seu$seurat_clusters %in% c(0,1,2)] = F
    seu = subset(seu, placeholder)
  }
  top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  print(DoHeatmap(seu, features = top10$gene, raster = T, assay = "RNA") + NoLegend() + ggtitle(paste0(integrated_name_arr[i]," DEG Heatmap")))

  system(paste0("rm /data/",integrated_name_arr[i],"_integrated.rds"))
}
dev.off()
library(ggplot2)
library(Seurat)
library(stringr)
library(rlang)
library(rlist)

### title: Plotting violin plots of quality metrics for cutaneous melanoma and NSCLC datasets, using data downsampled to obtain comparable UMI counts across samples
### author: Yiping Wang date: 11/08/2022

integrated_name_arr_underscore = c("BI5","NR1")
integrated_name_arr = c("BI5","NR1")
fresh_idents = c("CD45negGEXBI5_S1_L001","CD45posGEXBI5_S1_L001","NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX")
identsarr = list(c("CD45negGEXBI5_S1_L001","CD45posGEXBI5_S1_L001","bi005-skcm","bi005-skcm-5snseq","skcm-bi005-5pv2-snseq","BI5CST","BI5TST"),c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX","NSCL-NR001-5pv2-snseq","NR1CST","NR1TST"))

for (i in 1:length(integrated_name_arr_underscore)) {

  source("merge.SCTAssay.R")
  integrated_rds = double_filtered_seurat_arr[[identsarr[[i]][1]]]
  for (j in 2:length(identsarr[[i]]))
  {
    integrated_rds = merge.SCTAssay(x=integrated_rds, y=double_filtered_seurat_arr[[identsarr[[i]][j]]], na.rm = T)
  }

  integrated_rds$fresh_frozen = "frozen"
  integrated_rds$fresh_frozen[integrated_rds$orig.ident %in% fresh_idents] = "fresh"

  source("fresh_vs_frozen_comparison/QC_metrics/rename_sample_IDs.R")
  integrated_rds = rename_IDs(integrated_rds,integrated_name_arr_underscore[i])

  pdf(paste0("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_QC/",integrated_name_arr_underscore[i],"_fresh_vs_frozen_QC_downsampled_automatic.pdf"),height=7,width=5.5/3*length(identsarr[[i]]))

  integrated_rds$orig.ident = paste0(integrated_rds$orig.ident,"\n(",integrated_rds$fresh_frozen,")")

  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_nFeature_RNA = integrated_rds$nFeature_RNA")))

  if (sum(integrated_rds$fresh_frozen=="fresh")!=0)
  {
    integrated_rds$orig.ident[integrated_rds$fresh_frozen=="fresh"] = paste0("1_",integrated_rds$orig.ident[integrated_rds$fresh_frozen=="fresh"])
  }
  integrated_rds$orig.ident[integrated_rds$fresh_frozen=="frozen"] = paste0("2_",integrated_rds$orig.ident[integrated_rds$fresh_frozen=="frozen"])

  relabel_list = unique(integrated_rds$orig.ident)
  names(relabel_list) = relabel_list
  for (i1 in 1:length(relabel_list))
  {
    relabel_list[i1] = substring(relabel_list[i1],3)
  }

  uniqueidents = unique(integrated_rds$orig.ident)

  aplot = VlnPlot(integrated_rds, features = c(paste0(integrated_name_arr_underscore[i],"_nFeature_RNA")), group.by = "orig.ident",pt.size = 0)
  uniqueidents = sort(unique(integrated_rds$orig.ident))
  ymax_arr = c(13000)
  ylabs_arr = c("# detected genes/cell")
  for (z2 in 1:1)
  {
    colorsarr = rep("",length(uniqueidents))
    colorsarr[grep("fresh",uniqueidents)] = "blue"
    colorsarr[grep("frozen",uniqueidents)] = "red"
    if (integrated_name_arr[i]=="BI5")
    {
      colorsarr[grep("BI5CST",uniqueidents)] = "green"
      colorsarr[grep("BI5TST",uniqueidents)] = "green"
    }
    if (integrated_name_arr[i]=="NR1")
    {
      colorsarr[grep("NR1CST",uniqueidents)] = "green"
      colorsarr[grep("NR1TST",uniqueidents)] = "green"
    }
    #reorder nsclc and ribas samples, so that fresh samples are plotted first
    if (integrated_name_arr[i]=="BI5")
    {
      uniqueidents_reorder = uniqueidents
      uniqueidents_reorder = uniqueidents_reorder[c(1,2,5,6,7,3,4)]
    }
    if (integrated_name_arr[i]=="NR1")
    {
      uniqueidents_reorder = uniqueidents
      uniqueidents_reorder = uniqueidents_reorder[c(1,8,4,5,6,7,2,3)]
    }
    aplot[[z2]] = aplot[[z2]] + scale_fill_manual(values = colorsarr, limits = uniqueidents) + scale_x_discrete(limits = uniqueidents_reorder, labels=relabel_list) + ylim(0,ymax_arr[z2]) + ylab(ylabs_arr[z2]) + theme(axis.text.x = element_text(angle=0, hjust=0.5), legend.position = "none") + geom_boxplot(fill="white",width=0.1)
  }
  AugmentPlot(aplot, dpi = 300)
  print(aplot)
  dev.off()
}

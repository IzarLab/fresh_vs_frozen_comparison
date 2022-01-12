library(Seurat)
library(infercnv)
library(ggplot2)
library(pheatmap)
library(grid)
library(rlist)
library(stringr)

setwd("/data/")

selectstep = 20
use_one_cell_type = TRUE
one_cell_type = "Melanocytes"
sort_by_cell_type = TRUE
use_subcluster_cell_types = TRUE
if (use_subcluster_cell_types) {
  sort_by_cell_type = TRUE
}
useRollingAverage = FALSE
allindividualsamples = TRUE

useUvealAll = FALSE
if (useUvealAll) {
  system("aws s3 cp s3://uveal-melanoma/Seurat/integrated/um_all_integrated.rds /data/um_all_integrated.rds")
  um_all_integrated = readRDS("/data/um_all_integrated.rds")
  system("rm /data/um_all_integrated.rds")

  sub_um_all_integrated = subset(um_all_integrated, manual_annotation_label %in% c("Neuronal","Ribosomal","Ribosomal/Mitochondrial","Tumour","Tumour/Ribosomal"))

  origpatslist = list(c("um_07_gk_pre_S4_L001","um_07_gk_on_S8_L001","uv003-uvm3-snseq-3p-post","um_08_ar_pre_S1_L001","um_08_ar_on_S2_L001","um_08_ar_post_S3_L001","um_09_mw_pre_S5_L001","um_09_mw_on_S6_L001","um_09_mw_post_S7_L001","um_11_lc_pre_S12_L002","um_11_lc_on_S16_L002","um_12_ml_post_S11_L002","um_12_ml_on_S10_L002","um_12_ml_pre_S9_L002","um_15_lm_pre_S13_L002","um_15_lm_on_S14_L002","um_15_lm_post_S15_L002","um_16_rs_pre_S17_L003","um_16_rs_on_S18_L003","um_16_rs_post_S19_L003"))
  patslist = list(paste0(origpatslist[[1]],"_manual_annotation_label"))
  cbpatslist = list(c("um_07_gk_pre_S4_L001","um_07_gk_on_S8_L001","uv003-uvme-snseq-3p-post","um_08_ar_pre_S1_L001","um_08_ar_on_S2_L001","um_08_ar_post_S3_L001","um_09_mw_pre_S5_L001","um_09_mw_on_S6_L001","um_09_mw_post_S7_L001","um_11_lc_pre_S12_L002","um_11_lc_on_S16_L002","um_12_ml_post_S11_L002","um_12_ml_on_S10_L002","um_12_ml_pre_S9_L002","um_15_lm_pre_S13_L002","um_15_lm_on_S14_L002","um_15_lm_post_S15_L002","um_16_rs_pre_S17_L003","um_16_rs_on_S18_L003","um_16_rs_post_S19_L003"))
  s3folderlist = rep("uveal-melanoma",length(cbpatslist[[1]]))
  cancertissuelist = list(rep("Melanocytes",length(cbpatslist[[1]])))
  outputnames = c("um_all")

  for (i in 1:length(origpatslist[[1]]))
  {
    pat = origpatslist[[1]][i]
    manual_pat = patslist[[1]][i]
    s3folder = s3folderlist[i]
    system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/without_plots/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
    system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
    infercnv_obj_temp = readRDS(paste0("/data/",pat,"_infercnv.rds"))

    sub_um_all_integrated$barebarcode = unlist(lapply(strsplit(colnames(sub_um_all_integrated),"_"), function(x) x[[1]][1]))
    names(sub_um_all_integrated$barebarcode) = sub_um_all_integrated$barebarcode
    if (pat=="uv003-uvm3-snseq-3p-post")
    {
      match_idxs = na.omit(match(sub_um_all_integrated$barebarcode[sub_um_all_integrated$orig.ident=="uv003-uvme-snseq-3p-post"],colnames(infercnv_obj_temp@expr.data)))
    }
    else
    {
      match_idxs = na.omit(match(sub_um_all_integrated$barebarcode[sub_um_all_integrated$orig.ident==pat],colnames(infercnv_obj_temp@expr.data)))
    }
    infercnv_obj_temp@tumor_subclusters$subclusters$Melanocytes[["manual"]] = match_idxs
    names(infercnv_obj_temp@tumor_subclusters$subclusters$Melanocytes[["manual"]]) = colnames(infercnv_obj_temp@expr.data)[match_idxs]
    infercnv_obj_temp@observation_grouped_cell_indices$Melanocytes = match_idxs

    orig_subclusters = names(infercnv_obj_temp@tumor_subclusters$subclusters$Melanocytes)
    orig_subclusters = orig_subclusters[grep("Melanocytes",orig_subclusters)]
    for (z in 1:length(orig_subclusters))
    {
      infercnv_obj_temp@tumor_subclusters$subclusters$Melanocytes[[orig_subclusters[z]]] = NULL
    }

    saveRDS(infercnv_obj_temp,paste0("/data/",manual_pat,"_infercnv.rds"))
    system(paste0("aws s3 cp /data/",manual_pat,"_infercnv.rds s3://",s3folder,"/inferCNV/without_plots/inferCNV_cellbender_subcluster_",manual_pat,"/run.final.infercnv_obj"))
  }

  selectstep = 20
  use_one_cell_type = TRUE
  one_cell_type = "Melanocytes"
  sort_by_cell_type = FALSE
  use_subcluster_cell_types = FALSE
  if (use_subcluster_cell_types) {
    sort_by_cell_type = TRUE
  }
  useRollingAverage = FALSE
  allindividualsamples = FALSE

  rm(um_all_integrated)
  rm(sub_um_all_integrated)
}

useNSCLCAll = FALSE
if (useNSCLCAll) {
  patslist = list(c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX_final_thresh","NSCL-NR001-5pv2-snseq_final_thresh"))
  cbpatslist = list(c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX_final_thresh","NSCL-NR001-5pv2-snseq_final_thresh"))
  s3folderlist = c("fresh-vs-frozen-comparison/nsclc")#"uveal-melanoma"
  outputnames = c("nsclc_melanocytes")
}

useUvealIndividual = FALSE
useFreshVsFrozenIndividual = FALSE
useRibasIndividual = FALSE
useRibasFromLinyueIndividual = FALSE
useUvealKmeansIndividual = TRUE
if (useUvealIndividual || useFreshVsFrozenIndividual || useRibasIndividual || useRibasFromLinyueIndividual || useUvealKmeansIndividual) {
  allindividualsamples = TRUE
}
if (allindividualsamples) {
  if (useUvealIndividual)
  {
    patslist = list(c("um_07_gk_pre_S4_L001"),c("um_07_gk_on_S8_L001"),c("uv003-uvm3-snseq-3p-post"),c("um_08_ar_pre_S1_L001"),c("um_08_ar_on_S2_L001"),c("um_08_ar_post_S3_L001"),c("um_09_mw_pre_S5_L001"),c("um_09_mw_on_S6_L001"),c("um_09_mw_post_S7_L001"),c("um_11_lc_pre_S12_L002"),c("um_11_lc_on_S16_L002"),c("um_12_ml_post_S11_L002"),c("um_12_ml_on_S10_L002"),c("um_12_ml_pre_S9_L002"),c("um_15_lm_pre_S13_L002"),c("um_15_lm_on_S14_L002"),c("um_15_lm_post_S15_L002"),c("um_16_rs_pre_S17_L003"),c("um_16_rs_on_S18_L003"),c("um_16_rs_post_S19_L003"))
     cbpatslist = list(c("um_07_gk_pre_S4_L001"),c("um_07_gk_on_S8_L001"),c("uv003-uvme-snseq-3p-post"),c("um_08_ar_pre_S1_L001"),c("um_08_ar_on_S2_L001"),c("um_08_ar_post_S3_L001"),c("um_09_mw_pre_S5_L001"),c("um_09_mw_on_S6_L001"),c("um_09_mw_post_S7_L001"),c("um_11_lc_pre_S12_L002"),c("um_11_lc_on_S16_L002"),c("um_12_ml_post_S11_L002"),c("um_12_ml_on_S10_L002"),c("um_12_ml_pre_S9_L002"),c("um_15_lm_pre_S13_L002"),c("um_15_lm_on_S14_L002"),c("um_15_lm_post_S15_L002"),c("um_16_rs_pre_S17_L003"),c("um_16_rs_on_S18_L003"),c("um_16_rs_post_S19_L003"))
  }

  if (useFreshVsFrozenIndividual)
  {
    patslist = list(c("CD45negGEXBI5_S1_L001_final_thresh"),c("CD45posGEXBI5_S1_L001_final_thresh"),c("bi005-skcm_final_thresh"),c("bi005-skcm-5snseq_final_thresh"),c("skcm-bi005-5pv2-snseq_final_thresh"),c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-E12_final_thresh"),c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-F1_final_thresh"),c("UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-GEX-F12_final_thresh"),c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX_final_thresh"),c("NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX_final_thresh"),c("NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX_final_thresh"),c("NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX_final_thresh"),c("NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX_final_thresh"),c("NSCL-NR001-5pv2-snseq_final_thresh"),c("Sarcoma167GEX_final_thresh"),c("Sarcoma322GEX_final_thresh"),c("Sarcoma559GEX_final_thresh"),c("Sarcoma708GEX_final_thresh"))

    cbpatslist = list(c("CD45negGEXBI5_S1_L001_final_thresh"),c("CD45posGEXBI5_S1_L001_final_thresh"),c("bi005-skcm_final_thresh"),c("bi005-skcm-5snseq_final_thresh"),c("skcm-bi005-5pv2-snseq_final_thresh"),c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-E12_final_thresh"),c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-F1_final_thresh"),c("UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-GEX-F12_final_thresh"),c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX_final_thresh"),c("NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX_final_thresh"),c("NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX_final_thresh"),c("NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX_final_thresh"),c("NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX_final_thresh"),c("NSCL-NR001-5pv2-snseq_final_thresh"),c("Sarcoma167GEX_final_thresh"),c("Sarcoma322GEX_final_thresh"),c("Sarcoma559GEX_final_thresh"),c("Sarcoma708GEX_final_thresh"))

    s3folderlist = c("fresh-vs-frozen-comparison/BI5/scrna-seq","fresh-vs-frozen-comparison/BI5/scrna-seq","fresh-vs-frozen-comparison/BI5/snrna-seq","fresh-vs-frozen-comparison/BI5/snrna-seq","fresh-vs-frozen-comparison/BI5/snrna-seq","fresh-vs-frozen-comparison/cpoi-uvealprimarydata","fresh-vs-frozen-comparison/cpoi-uvealprimarydata","fresh-vs-frozen-comparison/cpoi-uvealprimarydata","fresh-vs-frozen-comparison/nsclc","fresh-vs-frozen-comparison/nsclc","fresh-vs-frozen-comparison/nsclc","fresh-vs-frozen-comparison/nsclc","fresh-vs-frozen-comparison/nsclc","fresh-vs-frozen-comparison/nsclc","fresh-vs-frozen-comparison/sarcoma-sn","fresh-vs-frozen-comparison/sarcoma-sn","fresh-vs-frozen-comparison/sarcoma-sn","fresh-vs-frozen-comparison/sarcoma-sn")

    outputnames = c("CD45negGEXBI5_S1_L001_final_thresh","CD45posGEXBI5_S1_L001_final_thresh","bi005-skcm_final_thresh","bi005-skcm-5snseq_final_thresh","skcm-bi005-5pv2-snseq_final_thresh","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-E12_final_thresh","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-F1_final_thresh","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-GEX-F12_final_thresh","NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX_final_thresh","NSCL-NR001-5pv2-snseq_final_thresh","Sarcoma167GEX_final_thresh","Sarcoma322GEX_final_thresh","Sarcoma559GEX_final_thresh","Sarcoma708GEX_final_thresh")

    cancertissuelist = list(c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),
    c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),
    c("Epithelial cells","Neurons"),c("Epithelial cells","Neurons"),c("Epithelial cells","Neurons"),c("Epithelial cells","Neurons"),c("Epithelial cells","Neurons"),c("Epithelial cells","Neurons"),
    c("Chondrocytes","Fibroblasts","Smooth muscle"),c("Chondrocytes","Fibroblasts","Smooth muscle"),c("Chondrocytes","Fibroblasts","Smooth muscle"),c("Chondrocytes","Fibroblasts","Smooth muscle"))
  }

  if (useRibasIndividual)
  {
    patslist = c("ribas_204_pre_GEX_titrate_thresh","ribas_294_on_GEX_titrate_thresh","ribas_308_pre_GEX_titrate_thresh","ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh","ribas_310_on_later_previd_3_GEX_titrate_thresh","ribas_310_pre_GEX_5pv2_S26_L004_titrate_thresh","ribas_319_on_previd_2_GEX_titrate_thresh","ribas_319_pre_previd_1_GEX_titrate_thresh","ribas_328_on_GEX_titrate_thresh","ribas_329_on_GEX_titrate_thresh","ribas_334_pre_GEX_titrate_thresh","ribas_354_pre_GEX_titrate_thresh")
    cbpatslist = c("ribas_204_pre_GEX_titrate_thresh","ribas_294_on_GEX_titrate_thresh","ribas_308_pre_GEX_titrate_thresh","ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh","ribas_310_on_later_previd_3_GEX_titrate_thresh","ribas_310_pre_GEX_5pv2_S26_L004_titrate_thresh","ribas_319_on_previd_2_GEX_titrate_thresh","ribas_319_pre_previd_1_GEX_titrate_thresh","ribas_328_on_GEX_titrate_thresh","ribas_329_on_GEX_titrate_thresh","ribas_334_pre_GEX_titrate_thresh","ribas_354_pre_GEX_titrate_thresh")
    outputnames = c("ribas_204_pre_GEX_titrate_thresh","ribas_294_on_GEX_titrate_thresh","ribas_308_pre_GEX_titrate_thresh","ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh","ribas_310_on_later_previd_3_GEX_titrate_thresh","ribas_310_pre_GEX_5pv2_S26_L004_titrate_thresh","ribas_319_on_previd_2_GEX_titrate_thresh","ribas_319_pre_previd_1_GEX_titrate_thresh","ribas_328_on_GEX_titrate_thresh","ribas_329_on_GEX_titrate_thresh","ribas_334_pre_GEX_titrate_thresh","ribas_354_pre_GEX_titrate_thresh")
    s3folderlist = c("melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1")
    cancertissuelist = list(c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),c("Melanocytes"),c("Melanocytes"))
  }

  if (useRibasFromLinyueIndividual)
  {
    patslist = c("ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh_kmeans","ribas_310_on_later_previd_3_GEX_titrate_thresh_kmeans","ribas_310_pre_GEX_5pv2_S26_L004_titrate_thresh_kmeans")
    origpatslist = c("ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh","ribas_310_on_later_previd_3_GEX_titrate_thresh","ribas_310_pre_GEX_5pv2_S26_L004_titrate_thresh")
    cbpatslist = c("ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh","ribas_310_on_later_previd_3_GEX_titrate_thresh","ribas_310_pre_GEX_5pv2_S26_L004_titrate_thresh")
    outputnames = c("ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh_kmeans","ribas_310_on_later_previd_3_GEX_titrate_thresh_kmeans","ribas_310_pre_GEX_5pv2_S26_L004_titrate_thresh_kmeans")
    s3folderlist = c("melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1")
    cancertissuelist = list(c("Melanocytes"),c("Melanocytes"),c("Melanocytes"))
    numberlist = c("3","4","5")
    kmeans_table = read.table("/mnt/vdb/home/ubuntu2/r310_kmeans_clones.csv",sep=",",header=T,quote=NULL)
  }
  if (useUvealKmeansIndividual)
  {
    # patslist = c("um_12_ml_on_S10_L002_kmeans","um_12_ml_pre_S9_L002_kmeans")
    # origpatslist = c("um_12_ml_on_S10_L002","um_12_ml_pre_S9_L002")
    # cbpatslist = c("um_12_ml_on_S10_L002","um_12_ml_pre_S9_L002")
    # outputnames = c("um_12_ml_on_S10_L002_kmeans","um_12_ml_pre_S9_L002_kmeans")
    # s3folderlist = c("uveal-melanoma","uveal-melanoma")
    # cancertissuelist = list(c("Melanocytes"),c("Melanocytes"))
    # kmeans_table = read.table("/mnt/vdb/home/ubuntu2/um_12_kmeans_clones.csv",sep="\t",header=T,quote=NULL)

    patslist = c("um_11_lc_on_S16_L002_kmeans","um_11_lc_pre_S12_L002_kmeans")
    origpatslist = c("um_11_lc_on_S16_L002","um_11_lc_pre_S12_L002")
    cbpatslist = c("um_11_lc_on_S16_L002","um_11_lc_pre_S12_L002")
    outputnames = c("um_11_lc_on_S16_L002_kmeans","um_11_lc_pre_S12_L002_kmeans")
    s3folderlist = c("uveal-melanoma","uveal-melanoma")
    cancertissuelist = list(c("Melanocytes"),c("Melanocytes"))
    kmeans_table = read.table("/mnt/vdb/home/ubuntu2/um_11_kmeans_clones.csv",sep="\t",header=T,quote=NULL)
  }

  selectstep = 1
}

remove_ref_group = T

if (useRibasFromLinyueIndividual || useUvealKmeansIndividual) {
  for (largeindex in 1:length(patslist))
  {
    pats = patslist[[largeindex]]
    origpats = origpatslist[[largeindex]]
    cbpats = cbpatslist[[largeindex]]
    s3folder = s3folderlist[largeindex]
    if (useRibasFromLinyueIndividual)
    {
      number_id = numberlist[largeindex]
    }
    for (i in 1:length(origpats))
    {
      pat = origpats[i]
      kmeans_pat = pats[i]
      system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/without_plots/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
      system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
      infercnv_obj_temp = readRDS(paste0("/data/",pat,"_infercnv.rds"))

      if (useUvealKmeansIndividual)
      {
        kmeans_samples = kmeans_table$sample
	kmeans_barcodes = kmeans_table$barcode

        unique_kmeans = unique(kmeans_table$kmean)
	all_kmeans_idxs = c()
	for (unique_kmean in unique_kmeans)
	{
	  kmeans_idxs = match(kmeans_barcodes[kmeans_samples==pat & kmeans_table$kmean==unique_kmean],colnames(infercnv_obj_temp@expr.data))
	  kmeans_idxs_names = kmeans_barcodes[kmeans_samples==pat & kmeans_table$kmean==unique_kmean]
	  infercnv_obj_temp@tumor_subclusters$subclusters$Melanocytes[[paste0("kmeans_",as.character(unique_kmean))]] = kmeans_idxs
	  names(infercnv_obj_temp@tumor_subclusters$subclusters$Melanocytes[[paste0("kmeans_",as.character(unique_kmean))]]) = kmeans_idxs_names
	  all_kmeans_idxs = c(all_kmeans_idxs, kmeans_idxs)
	}
	infercnv_obj_temp@observation_grouped_cell_indices$Melanocytes = all_kmeans_idxs
      }
      if (useRibasFromLinyueIndividual)
      {
	linyue_sample_ids = unlist(lapply(strsplit(kmeans_table$barcode,"-"), function(x) {x[[3]]}))
	linyue_barcodes = unlist(lapply(strsplit(kmeans_table$barcode,"-"), function(x) {paste0(x[[1]],"-1")}))

	unique_kmeans = unique(kmeans_table$kmeans.inferCNV.clones)
	for (unique_kmean in unique_kmeans)
	{
	  linyue_idxs = match(linyue_barcodes[linyue_sample_ids==number_id & kmeans_table$kmeans.inferCNV.clones==unique_kmean],colnames(infercnv_obj_temp@expr.data))
	  linyue_idxs_names = linyue_barcodes[linyue_sample_ids==number_id & kmeans_table$kmeans.inferCNV.clones==unique_kmean]
	  infercnv_obj_temp@tumor_subclusters$subclusters$Melanocytes[[paste0("kmeans_",as.character(unique_kmean))]] = linyue_idxs
	  names(infercnv_obj_temp@tumor_subclusters$subclusters$Melanocytes[[paste0("kmeans_",as.character(unique_kmean))]]) = linyue_idxs_names
	}
      }
      
      orig_subclusters = names(infercnv_obj_temp@tumor_subclusters$subclusters$Melanocytes)
      orig_subclusters = orig_subclusters[grep("Melanocytes",orig_subclusters)]#c("Melanocytes.1.1.1.1","Melanocytes.1.1.1.2","Melanocytes.1.1.2.1","Melanocytes.1.1.2.2","Melanocytes.1.2.1.1","Melanocytes.1.2.1.2","Melanocytes.1.2.2.1","Melanocytes.1.2.2.2")
      for (z in 1:length(orig_subclusters))
      {
        infercnv_obj_temp@tumor_subclusters$subclusters$Melanocytes[[orig_subclusters[z]]] = NULL
      }

      saveRDS(infercnv_obj_temp,paste0("/data/",kmeans_pat,"_infercnv.rds"))
      system(paste0("aws s3 cp /data/",kmeans_pat,"_infercnv.rds s3://",s3folder,"/inferCNV/without_plots/inferCNV_cellbender_subcluster_",kmeans_pat,"/run.final.infercnv_obj"))
      #nonsense = nonsense+1
    }
  }
}

#nonsense = nonsense+1

for (largeindex in 1:length(patslist))
{

  pats = patslist[[largeindex]]
  cbpats = cbpatslist[[largeindex]]
  s3folder = s3folderlist[largeindex]
  if (use_one_cell_type)
  {
    cancertissues = cancertissuelist[[largeindex]]
  }
  else
  {
    cancertissues = c("not one cell type")
  }

  for (cancertissue in cancertissues)
  {
    one_cell_type = cancertissue
    if (allindividualsamples)
    {
      #pdf(paste0(pats[1],"_infercnv.pdf"),width=15,height=30)
    }
    else
    {
      #pdf(paste0("integrated_cnv_expr_data.pdf"),width=15,height=30)
    }
    all_gene_names = c()
    all_cell_names = c()
    for (patidx in 1:length(pats))
    {
      pat = pats[patidx]
      cbpat = cbpats[patidx]
      system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/without_plots/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
      system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
      if (!(useRibasFromLinyueIndividual || useUvealKmeansIndividual || useUvealAll))
      {
        cbpat = pat
      }
      if (pat=="uv003-uvm3-snseq-3p-post")
      {
	cbpat = "uv003-uvme-snseq-3p-post"
      }
      system(paste0("aws s3 cp s3://",s3folder,"/Seurat/",cbpat,"/",cbpat,"_cb.rds /data/",pat,"_cb.rds"))

      # system(paste0("aws s3 cp s3://melanoma-ribas/ribas1/inferCNV/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
      # cbpat = pat
      # system(paste0("aws s3 cp s3://melanoma-ribas/ribas1/Seurat/",cbpat,"/",cbpat,"_cb.rds /data/",pat,"_cb.rds"))

      infercnv_obj_temp = readRDS(paste0("/data/",pat,"_infercnv.rds"))
      colnames(infercnv_obj_temp@expr.data) = paste0(colnames(infercnv_obj_temp@expr.data),"_",pat)
      if (use_one_cell_type)
      {
	eval(parse(text=paste0("all_adipocyte_idxs = infercnv_obj_temp@observation_grouped_cell_indices$\"",one_cell_type,"\"")))
	infercnv_obj_temp@expr.data = infercnv_obj_temp@expr.data[,all_adipocyte_idxs]
	for (aname in names(infercnv_obj_temp@observation_grouped_cell_indices))
	{
	  if (aname!=one_cell_type)
	  {
	    infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	    infercnv_obj_temp@tumor_subclusters$subclusters[[aname]] = NULL
	    infercnv_obj_temp@tumor_subclusters$hc[[aname]] = NULL
	  }
	  else
	  {
	    infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	    infercnv_obj_temp@observation_grouped_cell_indices[[pat]] = 1:dim(infercnv_obj_temp@expr.data)[2]
	    for (subcluster in names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]]))
	    {
	      infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]] = match(paste0(names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]]),"_",pat),colnames(infercnv_obj_temp@expr.data))
	    }
	  }
	}
	if (length(infercnv_obj_temp@reference_grouped_cell_indices)==1)
	{
	  infercnv_obj_temp@reference_grouped_cell_indices[[1]] = c()
	}
      }

      all_gene_names = union(all_gene_names, rownames(infercnv_obj_temp@expr.data))
      all_cell_names = c(all_cell_names, colnames(infercnv_obj_temp@expr.data))

      system(paste0("rm /data/",pat,"_cb.rds"))
    }

    all_expr_data = matrix(0, length(all_gene_names), length(all_cell_names))
    all_gene_order = data.frame(chr = rep("",length(all_gene_names)))
    all_gene_order$start = 1
    all_gene_order$stop = 1
    all_reference_grouped_cell_indices = c()
    all_observation_grouped_cell_indices = list()
    all_cell_type_bped_main = data.frame(cell_type_bped_main = rep("",length(all_cell_names)))
    cluster_cnv_profiles = matrix(0, length(all_gene_names), 0)

    rownames(all_expr_data) = all_gene_names
    colnames(all_expr_data) = all_cell_names
    rownames(all_gene_order) = all_gene_names
    rownames(all_cell_type_bped_main) = all_cell_names
    rownames(cluster_cnv_profiles) = all_gene_names

    if (sort_by_cell_type)
    {
      for (patidx in 1:length(pats))
      {
        pat = pats[patidx]
	cbpat = cbpats[patidx]
	system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/without_plots/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
	system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
	if (!(useRibasFromLinyueIndividual || useUvealKmeansIndividual || useUvealAll))
	{
	  cbpat = pat
	}
	if (pat=="uv003-uvm3-snseq-3p-post")
	{
	  cbpat = "uv003-uvme-snseq-3p-post"
	}
	system(paste0("aws s3 cp s3://",s3folder,"/Seurat/",cbpat,"/",cbpat,"_cb.rds /data/",pat,"_cb.rds"))

	# system(paste0("aws s3 cp s3://melanoma-ribas/ribas1/inferCNV/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
	# cbpat = pat
	# system(paste0("aws s3 cp s3://melanoma-ribas/ribas1/Seurat/",cbpat,"/",cbpat,"_cb.rds /data/",pat,"_cb.rds"))

	infercnv_obj_temp = readRDS(paste0("/data/",pat,"_infercnv.rds"))
	if (use_one_cell_type)
	{
	  eval(parse(text=paste0("all_adipocyte_idxs = infercnv_obj_temp@observation_grouped_cell_indices$\"",one_cell_type,"\"")))
	  infercnv_obj_temp@expr.data = infercnv_obj_temp@expr.data[,all_adipocyte_idxs]
	  colnames(infercnv_obj_temp@expr.data) = paste0(colnames(infercnv_obj_temp@expr.data),"_",pat)
	  for (aname in names(infercnv_obj_temp@observation_grouped_cell_indices))
	  {
	    if (aname!=one_cell_type)
	    {
	      infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	      infercnv_obj_temp@tumor_subclusters$subclusters[[aname]] = NULL
	      infercnv_obj_temp@tumor_subclusters$hc[[aname]] = NULL
	    }
	    else
	    {
	      infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	      infercnv_obj_temp@observation_grouped_cell_indices[[pat]] = 1:dim(infercnv_obj_temp@expr.data)[2]
	      for (subcluster in names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]]))
	      {
		infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]] = match(paste0(names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]]),"_",pat),colnames(infercnv_obj_temp@expr.data))
	      }
	    }
	  }
	  if (length(infercnv_obj_temp@reference_grouped_cell_indices)==1)
	  {
	    infercnv_obj_temp@reference_grouped_cell_indices[[1]] = c()
	  }
	}
	orig_obj_temp = readRDS(paste0("/data/",pat,"_cb.rds"))
	orig_obj_temp_colnames = paste0(colnames(orig_obj_temp),"_",pat)

	if (use_subcluster_cell_types)
	{
	  for (subcluster in names(infercnv_obj_temp@tumor_subclusters$subclusters))
	  {
	    for (subcluster1 in names(infercnv_obj_temp@tumor_subclusters$subclusters[[subcluster]]))
	    {
	      print(subcluster1)
	      if (subcluster!="tcell" && subcluster!="tcell_and_bcell")
	      {
		all_cell_type_bped_main$cell_type_bped_main[infercnv_obj_temp@tumor_subclusters$subclusters[[subcluster]][[subcluster1]]] = subcluster1
	      }
	    }
	  }
	}
	else
	{
	  matchidxs = match(orig_obj_temp_colnames,all_cell_names)
	  matchidxs = matchidxs[!is.na(matchidxs)]
	  all_cell_type_bped_main$cell_type_bped_main[matchidxs] = orig_obj_temp$celltype_bped_main
	}
	#browser()

	system(paste0("rm /data/",pat,"_cb.rds"))
      }

      all_cell_names = all_cell_names[order(all_cell_type_bped_main$cell_type_bped_main)]
      all_cell_type_bped_main = data.frame(cell_type_bped_main = all_cell_type_bped_main$cell_type_bped_main[order(all_cell_type_bped_main$cell_type_bped_main)])
      colnames(all_expr_data) = all_cell_names
      rownames(all_cell_type_bped_main) = all_cell_names
    }

    for (pat in pats)
    {
      infercnv_obj_temp = readRDS(paste0("/data/",pat,"_infercnv.rds"))
      colnames(infercnv_obj_temp@expr.data) = paste0(colnames(infercnv_obj_temp@expr.data),"_",pat)

      #browser()

      if (use_one_cell_type)
      {
	eval(parse(text=paste0("all_adipocyte_idxs = infercnv_obj_temp@observation_grouped_cell_indices$\"",one_cell_type,"\"")))
	infercnv_obj_temp@expr.data = infercnv_obj_temp@expr.data[,all_adipocyte_idxs]
	for (aname in names(infercnv_obj_temp@observation_grouped_cell_indices))
	{
	  if (aname!=one_cell_type)
	  {
	    infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	    infercnv_obj_temp@tumor_subclusters$subclusters[[aname]] = NULL
	    infercnv_obj_temp@tumor_subclusters$hc[[aname]] = NULL
	  }
	  else
	  {
	    infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	    infercnv_obj_temp@observation_grouped_cell_indices[[pat]] = 1:dim(infercnv_obj_temp@expr.data)[2]
	    for (subcluster in names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]]))
	    {
	      infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]] = match(paste0(names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]]),"_",pat),colnames(infercnv_obj_temp@expr.data))
	    }
	  }
	}
	if (length(infercnv_obj_temp@reference_grouped_cell_indices)==1)
	{
	  infercnv_obj_temp@reference_grouped_cell_indices = list(c())
	}
      }

      #browser()

      if (!(use_one_cell_type) || use_subcluster_cell_types)
      {
	for (aname in names(infercnv_obj_temp@observation_grouped_cell_indices))
	{
	  infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	}
	for (a_cell_type in unique(all_cell_type_bped_main$cell_type_bped_main))
	{
	  infercnv_obj_temp@observation_grouped_cell_indices[[a_cell_type]] = match(rownames(all_cell_type_bped_main)[all_cell_type_bped_main$cell_type_bped_main==a_cell_type],colnames(infercnv_obj_temp@expr.data))
	}
      }

      #browser()

      all_expr_data[match(rownames(infercnv_obj_temp@expr.data),all_gene_names),match(colnames(infercnv_obj_temp@expr.data),all_cell_names)] = infercnv_obj_temp@expr.data

      all_reference_grouped_cell_indices = c(all_reference_grouped_cell_indices,match(colnames(infercnv_obj_temp@expr.data)[infercnv_obj_temp@reference_grouped_cell_indices[[1]]],all_cell_names))

      for (celltype in names(infercnv_obj_temp@observation_grouped_cell_indices))
      {
	if (!(celltype %in% names(all_observation_grouped_cell_indices)))
	{
	  eval(parse(text=paste0("all_observation_grouped_cell_indices = list.append(all_observation_grouped_cell_indices, \"",celltype,"\" = c())")))
	  print(length(all_observation_grouped_cell_indices))
	}
	eval(parse(text=paste0("an_observation_grouped_cell_indices = all_observation_grouped_cell_indices$`",celltype,"`")))
	eval(parse(text=paste0("an_observation_grouped_cell_indices_temp = infercnv_obj_temp@observation_grouped_cell_indices$`",celltype,"`")))
	an_observation_grouped_cell_indices = c(an_observation_grouped_cell_indices,match(colnames(infercnv_obj_temp@expr.data)[an_observation_grouped_cell_indices_temp],all_cell_names))
	eval(parse(text=paste0("all_observation_grouped_cell_indices[[\"",celltype,"\"]] = an_observation_grouped_cell_indices")))
      }

      all_gene_order$chr[match(rownames(infercnv_obj_temp@gene_order),all_gene_names)] = infercnv_obj_temp@gene_order$chr
      all_gene_order$start[match(rownames(infercnv_obj_temp@gene_order),all_gene_names)] = infercnv_obj_temp@gene_order$start
      all_gene_order$stop[match(rownames(infercnv_obj_temp@gene_order),all_gene_names)] = infercnv_obj_temp@gene_order$stop

      # main_cell_types = names(infercnv_obj_temp@tumor_subclusters$subclusters)
      # for (i in 1:length(main_cell_types))
      # {
      # 	cluster_names = eval(parse(text=paste0("names(infercnv_obj_temp@tumor_subclusters$subclusters$`",main_cell_types[i],"`)")))
      # 	for (j in 1:length(cluster_names))
      # 	{
      # 	  cellidxs = eval(parse(text=paste0("infercnv_obj_temp@tumor_subclusters$subclusters$`",main_cell_types[i],"`$`",cluster_names[j],"`")))
      # 	  if (length(cellidxs)>1)
      # 	  {
      # 	    cnv_profile = rowMeans(infercnv_obj_temp@expr.data[,cellidxs])
      # 	  }
      # 	  else
      # 	  {
      # 	    cnv_profile = infercnv_obj_temp@expr.data[,cellidxs]
      # 	  }

      # 	  cnv_profile_expanded = matrix(0,length(all_gene_names),1)
      # 	  rownames(cnv_profile_expanded) = all_gene_names
      # 	  cnv_profile_expanded[match(names(cnv_profile),rownames(cnv_profile_expanded))] = cnv_profile
      # 	  cluster_cnv_profiles = cbind(cluster_cnv_profiles,cnv_profile_expanded)
      # 	  colnames(cluster_cnv_profiles)[dim(cluster_cnv_profiles)[2]] = paste0(pat,"_",cluster_names[j])
      # 	  if (j==2)
      # 	  {
      # 	    #nonsense = nonsense+1
      # 	  }
      # 	  print(dim(cluster_cnv_profiles))
      # 	}
      # }
    }

    all_expr_data = all_expr_data[order(all_gene_order$chr,all_gene_order$start),]
    all_gene_order = all_gene_order[order(as.numeric(all_gene_order$chr),as.numeric(all_gene_order$start)),]

    infercnv_obj = infercnv_obj_temp
    infercnv_obj@expr.data = all_expr_data
    infercnv_obj@gene_order = all_gene_order
    infercnv_obj@reference_grouped_cell_indices = list(all_reference_grouped_cell_indices)
    infercnv_obj@observation_grouped_cell_indices = all_observation_grouped_cell_indices
    infercnv_obj@tumor_subclusters = NULL

    common_gene_names = all_gene_names
    for (pat2 in pats)
    {
      infercnv_obj_temp2 = readRDS(paste0("/data/",pat2,"_infercnv.rds"))
      common_gene_names = intersect(common_gene_names,rownames(infercnv_obj_temp2@gene_order))
    }

    infercnv_obj@expr.data = infercnv_obj@expr.data[match(common_gene_names,rownames(infercnv_obj@expr.data)),]
    infercnv_obj@gene_order = infercnv_obj@gene_order[match(common_gene_names,rownames(infercnv_obj@gene_order)),]

    selectidxs = seq(1,dim(infercnv_obj@expr.data)[2],selectstep)
    dfselect = data.frame(origidxs = selectidxs, newidxs = 1:length(selectidxs))
    infercnv_obj@expr.data = infercnv_obj@expr.data[,selectidxs]
    infercnv_obj@reference_grouped_cell_indices[[1]] = dfselect$newidxs[match(infercnv_obj@reference_grouped_cell_indices[[1]][infercnv_obj@reference_grouped_cell_indices[[1]] %in% selectidxs],dfselect$origidxs)]
    observation_grouped_cell_indices_temp = infercnv_obj@observation_grouped_cell_indices
    for (aname in names(infercnv_obj@observation_grouped_cell_indices))
    {
      testsub = dfselect$newidxs[match(infercnv_obj@observation_grouped_cell_indices[[aname]][infercnv_obj@observation_grouped_cell_indices[[aname]] %in% selectidxs],dfselect$origidxs)]
      #testsub = infercnv_obj@observation_grouped_cell_indices[[aname]][dfselect$newidxs[match(infercnv_obj@observation_grouped_cell_indices[[aname]] %in% selectidxs,dfselect$origidxs)]]
      print(length(testsub))
      print(aname)
    }
    for (aname in names(infercnv_obj@observation_grouped_cell_indices))
    {
      testsub = dfselect$newidxs[match(infercnv_obj@observation_grouped_cell_indices[[aname]][infercnv_obj@observation_grouped_cell_indices[[aname]] %in% selectidxs],dfselect$origidxs)]
      if (length(testsub)<2)
      {
	observation_grouped_cell_indices_temp[[aname]] = NULL
      }
      else
      {
	observation_grouped_cell_indices_temp[[aname]] = testsub
      }
    }
    infercnv_obj@observation_grouped_cell_indices = observation_grouped_cell_indices_temp

    obs_annotations_groups = rep(-1, length(colnames(infercnv_obj@expr.data)))
    names(obs_annotations_groups) = colnames(infercnv_obj@expr.data)
    obs_index_groupings = infercnv_obj@observation_grouped_cell_indices
    counter <- 1
    for (obs_index_group in obs_index_groupings) {
	obs_annotations_groups[ obs_index_group ] <- counter
	counter <- counter + 1
    }
    obs_annotations_groups[infercnv_obj@reference_grouped_cell_indices[[1]]] = counter+1
    infercnv_obj@observation_grouped_cell_indices$negidxs = which(obs_annotations_groups==-1)

    #fudges for single cell type plotting only
    if (use_one_cell_type)
    {
      #print(length(infercnv_obj@reference_grouped_cell_indices[[1]])==0)
      #print(sum(is.na(infercnv_obj@reference_grouped_cell_indices[[1]]))>0)
      #print(length(infercnv_obj@reference_grouped_cell_indices[[1]])==0 || sum(is.na(infercnv_obj@reference_grouped_cell_indices[[1]]))>0)
      if (length(infercnv_obj@reference_grouped_cell_indices[[1]])==0 || sum(is.na(infercnv_obj@reference_grouped_cell_indices[[1]]))>0)
      {
	infercnv_obj@reference_grouped_cell_indices[[1]] = c(1)
      }
    }
    if (length(infercnv_obj@observation_grouped_cell_indices$negidxs)==0)
    {
      infercnv_obj@observation_grouped_cell_indices$negidxs = NULL
    }

    chrs = unique(infercnv_obj@gene_order$chr)
    expr.data.temp = matrix(0,dim(infercnv_obj@expr.data)[1],dim(infercnv_obj@expr.data)[2])
    rownames(expr.data.temp) = rownames(infercnv_obj@expr.data)
    colnames(expr.data.temp) = colnames(infercnv_obj@expr.data)
    #gene_order_temp = data.frame(chr=character(), start = integer(), end = integer())
    #rownames(gene_order_temp) = rownames(infercnv_obj@gene_order)

    if (useRollingAverage)
    {
      for (z1 in 1:dim(infercnv_obj@expr.data)[2])
      {
	for (chr in chrs)
	{
	  expr.data.temp[infercnv_obj@gene_order$chr==chr,z1] = rollmean(infercnv_obj@expr.data[infercnv_obj@gene_order$chr==chr,z1],k=11,fill=NA)
	}
      }
      infercnv_obj@gene_order = infercnv_obj@gene_order[!is.na(expr.data.temp[,1]),]
      #expr.data.temp = na.omit(expr.data.temp)
      #infercnv_obj@expr.data = expr.data.temp
      infercnv_obj@expr.data = expr.data.temp[!is.na(expr.data.temp[,1]),]
      print("roll")
    }

    #lineardata = infercnv_obj@expr.data[1:(dim(infercnv_obj@expr.data)[1]*dim(infercnv_obj@expr.data)[2])]

    oldmedian = median(median(infercnv_obj@expr.data))
    infercnv_obj@expr.data[infercnv_obj@expr.data>.88 & infercnv_obj@expr.data<1.12] = oldmedian

    if (remove_ref_group)
    {
      if (length(infercnv_obj@reference_grouped_cell_indices[[1]])>=10)
      {
        infercnv_obj@reference_grouped_cell_indices[[1]] = infercnv_obj@reference_grouped_cell_indices[[1]][1:10]
      }
    }

    #nonsense = nonsense+1
    source("/mnt/vdb/home/ubuntu2/uveal-melanoma-code/plot_cnv_fullcode.R")
    plot_cnv(infercnv_obj)
    if (use_one_cell_type)
    {
      system(paste0("mv infercnv.pdf ",outputnames[largeindex],"_",str_replace_all(one_cell_type," ","_"),"_infercnv.pdf"))
    }
    else
    {
      system(paste0("mv infercnv.pdf ",outputnames[largeindex],"_infercnv.pdf"))
    }
  }
}
library(ggplot2)
library(Seurat)
library(stringr)
library(rlang)
library(rlist)
library(dplyr)

### title: Print out csv files of count matrices for single-nucleus and single-cell sequencing data across all datasets
### author: Yiping Wang date: 11/08/2022

# rename_tcr_arr = c("Mel_sc_5_CD45+","Mel_sn_5","UM_sc_5","UM_sc_5_CD45+","UM_sn_5_inhib","NSCLC_sc_5","NSCLC_sn_5","NSCLC_sn_5_inhib","ribas_pre","ribas_on","ribas_on_later","UMEL_1_1","UMEL_1_2","UMEL_2_1","UMEL_2_2","UMEL_2_3","UMEL_3_1","UMEL_3_2","UMEL_3_3","UMEL_4_1","UMEL_4_2","UMEL_5_1","UMEL_5_2","UMEL_5_3","UMEL_6_1","UMEL_6_2","UMEL_6_3","UMEL_7_1","UMEL_7_2","UMEL_7_3")

# names(rename_tcr_arr) = c("TCRBI5_S1_L001","bi005-skcm-5snseq-TCR","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F2","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F3","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-TCR-F9","NSCL_NR001_SCRNA_5P_NA_BRAIN_TCR","NSCL_NR001_SNSEQ_5P_NI_BRAIN_TCR","NSCL_NR001_SNSEQ_5P_WI_BRAIN_TCR","ribas1_pre_tcr_S35_L004","ribas1_on_tcr_S36_L004","ribas_310_on_later_previd_3_TCR","um_07_gk_pre_S4_L001","um_07_gk_on_S8_L001","um_08_ar_pre_S1_L001","um_08_ar_on_S2_L001","um_08_ar_post_S3_L001","um_09_mw_pre_S5_L001","um_09_mw_on_S6_L001","um_09_mw_post_S7_L001","um_11_lc_pre_S12_L002","um_11_lc_on_S16_L002","um_12_ml_pre_S9_L002","um_12_ml_on_S10_L002","um_12_ml_post_S11_L002","um_15_lm_pre_S13_L002","um_15_lm_on_S14_L002","um_15_lm_post_S15_L002","um_16_rs_pre_S17_L003","um_16_rs_on_S18_L003","um_16_rs_post_S19_L003")

# for (i in 1:length(rename_tcr_arr)) {
#   aname = names(rename_tcr_arr)[i]
#   system(paste0("cp /mnt/vdb/home/ubuntu2/",aname,"/all_contig_annotations.csv /data/GEO/processed/",rename_tcr_arr[i],"_all_contig_annotations.csv"))
#   system(paste0("cp /mnt/vdb/home/ubuntu2/",aname,"/filtered_contig_annotations.csv /data/GEO/processed/",rename_tcr_arr[i],"_filtered_contig_annotations.csv"))
#   system(paste0("cp /mnt/vdb/home/ubuntu2/",aname,"/clonotypes.csv /data/GEO/processed/",rename_tcr_arr[i],"_clonotypes.csv"))
#   system(paste0("cp /mnt/vdb/home/ubuntu2/",aname,"/consensus_annotations.csv /data/GEO/processed/",rename_tcr_arr[i],"_consensus_annotations.csv"))
# }

# nonsense = nonsense+1

# foldersList = c("",
#   "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata",
#   "s3://fresh-vs-frozen-comparison-ohio/nsclc",
#   "",
#   "s3://uveal-melanoma")
# integrated_name_arr = c("BI5","cpoi-uvealprimarydata","nsclc","ribas_310","um_all")
# integrated_name_arr_underscore = c("Mel","UM","NSCLC","ribas","UMEL")

foldersList = c("s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline",
  "s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
integrated_name_arr = c("BI5","NR1")
integrated_name_arr_underscore = c("BI5","NR1")

counts_dim_arr = c()
for (i in 1:length(foldersList)) {
  if (integrated_name_arr[i]=="BI5" && foldersList[i]!="s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
  {
    seu = readRDS("/data/fresh_vs_frozen_all_reannotate_BI5.rds")
  }
  else if (integrated_name_arr[i]=="ribas_310")
  {
    seu = readRDS("/data/fresh_vs_frozen_all_reannotate_ribas_310.rds")
  }
  else if (foldersList[i]=="s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
  {
    system(paste0("aws s3 cp ",foldersList[i],"/Seurat_downsampled/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
    seu = readRDS(paste0("/data/",integrated_name_arr[i],"_integrated.rds"))
  }
  else
  {
    system(paste0("aws s3 cp ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
    seu = readRDS(paste0("/data/",integrated_name_arr[i],"_integrated.rds"))
  }

  print(integrated_name_arr[i])
  print(dim(seu))
  print(length(unique(seu$orig.ident)))

  if (integrated_name_arr[i]=="BI5" && foldersList[i]!="s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
  {
    seu$orig.ident[seu$orig.ident=="CD45neg"] = "Mel_sc_5_CD45-"
    seu$orig.ident[seu$orig.ident=="CD45pos"] = "Mel_sc_5_CD45+"
    seu$orig.ident[seu$orig.ident=="3snseq"] = "Mel_sn_3"
    seu$orig.ident[seu$orig.ident=="5pv2-snseq"] = "Mel_sn_5v2"
    seu$orig.ident[seu$orig.ident=="5snseq"] = "Mel_sn_5"
  }
  else if (integrated_name_arr[i]=="cpoi-uvealprimarydata")
  {
    seu$orig.ident[seu$orig.ident=="SCRNA-5P-NA-E12"] = "UM_sc_5"
    seu$orig.ident[seu$orig.ident=="SCRNA-5P-NA-F1"] = "UM_sc_5_CD45+"
    seu$orig.ident[seu$orig.ident=="SNRNA-5P-WI-F12"] = "UM_sn_5_inhib"
  }
  else if (integrated_name_arr[i]=="nsclc")
  {
    seu$orig.ident[seu$orig.ident=="SCRNA_5P_NA"] = "NSCLC_sc_5"
    seu$orig.ident[seu$orig.ident=="5pv2-snseq"] = "NSCLC_sn_5v2"
    seu$orig.ident[seu$orig.ident=="SNSEQ_3P_NI"] = "NSCLC_sn_3"
    seu$orig.ident[seu$orig.ident=="SNSEQ_3P_WI"] = "NSCLC_sn_3_inhib"
    seu$orig.ident[seu$orig.ident=="SNSEQ_5P_NI"] = "NSCLC_sn_5"
    seu$orig.ident[seu$orig.ident=="SNSEQ_5P_WI"] = "NSCLC_sn_5_inhib"
  }
  else if (integrated_name_arr[i]=="ribas_310")
  {
    seu$orig.ident[seu$orig.ident=="ribas_310_pre_GEX_5pv2_S26_L004"] = "ribas_pre"
    seu$orig.ident[seu$orig.ident=="ribas_310_on_GEX_5pv2_S27_L004"] = "ribas_on"
    seu$orig.ident[seu$orig.ident=="ribas_310_on_later_previd_3_GEX"] = "ribas_on_later"
  }
  else if (integrated_name_arr[i]=="um_all")
  {
    seu$orig.ident[seu$orig.ident=="um_07_gk_pre_S4_L001"] = "UMEL_1_1"
    seu$orig.ident[seu$orig.ident=="um_07_gk_on_S8_L001"] = "UMEL_1_2"
    seu$orig.ident[seu$orig.ident=="uv003-uvme-snseq-3p-post"] = "UMEL_1_3"
    seu$orig.ident[seu$orig.ident=="um_08_ar_pre_S1_L001"] = "UMEL_2_1"
    seu$orig.ident[seu$orig.ident=="um_08_ar_on_S2_L001"] = "UMEL_2_2"
    seu$orig.ident[seu$orig.ident=="um_08_ar_post_S3_L001"] = "UMEL_2_3"
    seu$orig.ident[seu$orig.ident=="um_09_mw_pre_S5_L001"] = "UMEL_3_1"
    seu$orig.ident[seu$orig.ident=="um_09_mw_on_S6_L001"] = "UMEL_3_2"
    seu$orig.ident[seu$orig.ident=="um_09_mw_post_S7_L001"] = "UMEL_3_3"
    seu$orig.ident[seu$orig.ident=="um_11_lc_pre_S12_L002"] = "UMEL_4_1"
    seu$orig.ident[seu$orig.ident=="um_11_lc_on_S16_L002"] = "UMEL_4_2"
    seu$orig.ident[seu$orig.ident=="um_12_ml_pre_S9_L002"] = "UMEL_5_1"
    seu$orig.ident[seu$orig.ident=="um_12_ml_on_S10_L002"] = "UMEL_5_2"
    seu$orig.ident[seu$orig.ident=="um_12_ml_post_S11_L002"] = "UMEL_5_3"
    seu$orig.ident[seu$orig.ident=="um_15_lm_pre_S13_L002"] = "UMEL_6_1"
    seu$orig.ident[seu$orig.ident=="um_15_lm_on_S14_L002"] = "UMEL_6_2"
    seu$orig.ident[seu$orig.ident=="um_15_lm_post_S15_L002"] = "UMEL_6_3"
    seu$orig.ident[seu$orig.ident=="um_16_rs_pre_S17_L003"] = "UMEL_7_1"
    seu$orig.ident[seu$orig.ident=="um_16_rs_on_S18_L003"] = "UMEL_7_2"
    seu$orig.ident[seu$orig.ident=="um_16_rs_post_S19_L003"] = "UMEL_7_3"
  }
  else if (integrated_name_arr[i]=="BI5" && foldersList[i]=="s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
  {
    seu = subset(seu, orig.ident %in% c("BI5CST","BI5TST"))
  }
  else if (integrated_name_arr[i]=="NR1" && foldersList[i]=="s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
  {
    seu = subset(seu, orig.ident %in% c("NR1CST","NR1TST"))
  }

  if ((integrated_name_arr[i]=="BI5" || integrated_name_arr[i]=="NR1") && foldersList[i]=="s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
  {
    write.csv(as.data.frame(seu@assays$integrated@data),paste0('/data/GEO/processed/',integrated_name_arr_underscore[i],'_slyper_processed_data.csv'))
  }
  else
  {
    write.csv(as.data.frame(seu@assays$integrated@data),paste0('/data/GEO/processed/',integrated_name_arr_underscore[i],'_processed_data.csv'))
  }

  for (anident in unique(seu$orig.ident))
  {
    su = subset(seu, orig.ident==anident)
    system(paste0('mkdir /data/GEO/raw/',anident))
    write.csv(as.data.frame(su@assays$RNA@counts),paste0('/data/GEO/raw/',anident,'/',anident,'_raw_counts.csv'))
    counts_dim_arr = c(counts_dim_arr, dim(su@assays$RNA@counts)[2])
  }

  seu$barcode = colnames(seu)
  metadata_df<- seu@meta.data %>% select(c('barcode','orig.ident','ScrubDoublet_score'))

  if ((integrated_name_arr[i]=="BI5" || integrated_name_arr[i]=="NR1") && foldersList[i]=="s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
  {
    write.csv(metadata_df,paste0('/data/GEO/processed/',integrated_name_arr_underscore[i],'_slyper_meta_data.csv'),row.names = F)
  }
  else
  {
    write.csv(metadata_df,paste0('/data/GEO/processed/',integrated_name_arr_underscore[i],'_meta_data.csv'),row.names = F)
  }
}
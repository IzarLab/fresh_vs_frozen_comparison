library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(purrr)
library(cowplot)
library(stringr)

system("aws s3 cp s3://uveal-melanoma/Seurat/integrated/um_all_integrated.rds /data/um_all_integrated.rds")
seu = readRDS("/data/um_all_integrated.rds")
system("rm /data/um_all_integrated.rds")

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

outdf = data.frame(barcode = seu$barebarcodes, sample = seu$orig.ident, manual_annotation_label = seu$manual_annotation_label)
outdf = outdf[order(outdf$sample),]
write.table(outdf, "/data/uveal_melanoma_metastasis_cell_type_annotation.csv", row.names = F, col.names = T, sep = ",", quote=FALSE)

nonsense = nonsense+1

pat_list_arr = c("UMEL_1_1","UMEL_1_2","UMEL_1_3","UMEL_2_1","UMEL_2_2","UMEL_2_3","UMEL_3_1","UMEL_3_2","UMEL_3_3","UMEL_4_1","UMEL_4_2","UMEL_5_1","UMEL_5_2","UMEL_5_3","UMEL_6_1","UMEL_6_2","UMEL_6_3","UMEL_7_1","UMEL_7_2","UMEL_7_3")

system("rm -r /data/uveal_melanoma_metastasis_matrices")
system("mkdir /data/uveal_melanoma_metastasis_matrices")

for(pat in pat_list_arr){
  system(paste0("aws s3 cp s3://fresh-vs-frozen-comparison-ohio/GEO/raw/",pat,"/",pat,"_raw_counts.csv /data/uveal_melanoma_metastasis_matrices/",pat,"_raw_counts.csv"))
}

system("zip -r /data/uveal_melanoma_metastasis_matrices.zip /data/uveal_melanoma_metastasis_matrices")

system("rm -r /data/uveal_melanoma_metastasis_matrices")

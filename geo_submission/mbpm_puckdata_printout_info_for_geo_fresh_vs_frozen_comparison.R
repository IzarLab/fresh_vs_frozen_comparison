library(Seurat)
library(ggplot2)
library(rlist)
library(grid)

### title: Print out csv files of count matrices for spatial sequencing data for sequential cutaneous melanoma samples
### author: Yiping Wang date: 11/08/2022

#pats = c("MPM08_pre_slide","MPM08_on_slide","MPM08_on_later_slide")
pats = c("MPM08_on_slide","MPM08_on_later_slide")

for (pat in pats) {
  system(paste0("aws s3 cp s3://melanoma-brain-mets/SlideSeqV2/mbpm_puckdata_add_rctd_and_sigs/",pat,"_with_rctd_sigs.rds /data/mbpm_puckdata_add_rctd_and_sigs/",pat,"_with_rctd_sigs.rds"))
  puck = readRDS(paste0("/data/mbpm_puckdata_add_rctd_and_sigs/",pat,"_with_rctd_sigs.rds"))
  system(paste0("rm /data/mbpm_puckdata_add_rctd_and_sigs/",pat,"_with_rctd_sigs.rds"))
  output_pat = pat
  if (pat=="puck5")
  {
    output_pat="MPM01_rep1_slide"
  }
  if (pat=="puck6final")
  {
    output_pat="MPM01_rep2_slide"
  }
  if (pat=="puck7_20_feature_threshold")
  {
    output_pat="MBM05_rep3_slide"
  }
  if (pat=="puck8_20_feature_threshold")
  {
    output_pat="MBM11_rep3_slide"
  }
  output_pat = str_replace(output_pat, "MPM08", "ribas")

  spatialdf = data.frame(barcode = colnames(puck), xcoord = puck$image@coordinates$x, ycoord = puck$image@coordinates$y)

  write.csv(spatialdf, paste0("/data/mbpm_puckdata_printout_info_for_geo_fresh_vs_frozen_comparison/",output_pat,"_spatial_coordinates.csv"))
  system(paste0("aws s3 cp /data/mbpm_puckdata_printout_info_for_geo_fresh_vs_frozen_comparison/",output_pat,"_spatial_coordinates.csv s3://fresh-vs-frozen-comparison-ohio/GEO/raw/",output_pat,"/",output_pat,"_spatial_coordinates.csv"))
  system(paste0("rm /data/mbpm_puckdata_printout_info_for_geo_fresh_vs_frozen_comparison/",output_pat,"_spatial_coordinates.csv"))
  
  write.csv(as.data.frame(puck@assays$Spatial@counts),paste0("/data/mbpm_puckdata_printout_info_for_geo_fresh_vs_frozen_comparison/",output_pat,"_raw_counts.csv"))
  system(paste0("aws s3 cp /data/mbpm_puckdata_printout_info_for_geo_fresh_vs_frozen_comparison/",output_pat,"_raw_counts.csv s3://fresh-vs-frozen-comparison-ohio/GEO/raw/",output_pat,"/",output_pat,"_raw_counts.csv"))
  system(paste0("rm /data/mbpm_puckdata_printout_info_for_geo_fresh_vs_frozen_comparison/",output_pat,"_raw_counts.csv"))
}
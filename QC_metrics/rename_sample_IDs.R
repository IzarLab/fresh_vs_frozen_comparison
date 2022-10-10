rename_IDs <- function(integrated_rds, dataset_name) {
  if (dataset_name=="Mel" || dataset_name=="BI5")
  {
    integrated_rds$orig.ident[integrated_rds$orig.ident=="CD45neg"] = "Mel_sc_5_CD45-"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="CD45negGEXBI5_S1_L001"] = "Mel_sc_5_CD45-"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="CD45pos"] = "Mel_sc_5_CD45+"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="CD45posGEXBI5_S1_L001"] = "Mel_sc_5_CD45+"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="3snseq"] = "Mel_sn_3"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="bi005-skcm"] = "Mel_sn_3"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="5pv2-snseq"] = "Mel_sn_5v2"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="skcm-bi005-5pv2-snseq"] = "Mel_sn_5v2"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="5snseq"] = "Mel_sn_5"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="bi005-skcm-5snseq"] = "Mel_sn_5"
  }
  if (dataset_name=="UM")
  {
    integrated_rds$orig.ident[integrated_rds$orig.ident=="SCRNA-5P-NA-E12"] = "UM_sc_5"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="SCRNA-5P-NA-F1"] = "UM_sc_5_CD45+"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="SNRNA-5P-WI-F12"] = "UM_sn_5_inhib"
  }
  if (dataset_name=="NSCLC" || dataset_name=="NR1")
  {
    integrated_rds$orig.ident[integrated_rds$orig.ident=="SCRNA_5P_NA"] = "NSCLC_sc_5"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX"] = "NSCLC_sc_5"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="5pv2-snseq"] = "NSCLC_sn_5v2"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="NSCL-NR001-5pv2-snseq"] = "NSCLC_sn_5v2"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="SNSEQ_3P_NI"] = "NSCLC_sn_3"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX"] = "NSCLC_sn_3"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="SNSEQ_3P_WI"] = "NSCLC_sn_3_inhib"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX"] = "NSCLC_sn_3_inhib"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="SNSEQ_5P_NI"] = "NSCLC_sn_5"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX"] = "NSCLC_sn_5"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="SNSEQ_5P_WI"] = "NSCLC_sn_5_inhib"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX"] = "NSCLC_sn_5_inhib"
  }
  if (dataset_name=="ribas")
  {
    integrated_rds$orig.ident[integrated_rds$orig.ident=="ribas_310_pre_GEX_5pv2_S26_L004"] = "pre"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="ribas_310_on_GEX_5pv2_S27_L004"] = "on"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="ribas_310_on_later_previd_3_GEX"] = "on_later"
  }
  if (dataset_name=="UMEL")
  {
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_07_gk_pre_S4_L001"] = "UMEL_1_1"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_07_gk_on_S8_L001"] = "UMEL_1_2"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="uv003-uvme-snseq-3p-post"] = "UMEL_1_3"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_08_ar_pre_S1_L001"] = "UMEL_2_1"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_08_ar_on_S2_L001"] = "UMEL_2_2"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_08_ar_post_S3_L001"] = "UMEL_2_3"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_09_mw_pre_S5_L001"] = "UMEL_3_1"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_09_mw_on_S6_L001"] = "UMEL_3_2"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_09_mw_post_S7_L001"] = "UMEL_3_3"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_11_lc_pre_S12_L002"] = "UMEL_4_1"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_11_lc_on_S16_L002"] = "UMEL_4_2"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_12_ml_pre_S9_L002"] = "UMEL_5_1"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_12_ml_on_S10_L002"] = "UMEL_5_2"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_12_ml_post_S11_L002"] = "UMEL_5_3"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_15_lm_pre_S13_L002"] = "UMEL_6_1"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_15_lm_on_S14_L002"] = "UMEL_6_2"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_15_lm_post_S15_L002"] = "UMEL_6_3"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_16_rs_pre_S17_L003"] = "UMEL_7_1"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_16_rs_on_S18_L003"] = "UMEL_7_2"
    integrated_rds$orig.ident[integrated_rds$orig.ident=="um_16_rs_post_S19_L003"] = "UMEL_7_3"
  }
  return(integrated_rds)
}
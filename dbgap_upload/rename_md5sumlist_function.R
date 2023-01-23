library(stringr)

rename_md5sumlist_function <- function(afilename) {
  longnameslist = list(c("CD45negGEXBI5","CD45posGEXBI5","BI5","skcm-bi005-5pv2-snseq","bi005-skcm-5snseq"),
    c("TCRBI5","bi005-skcm-5snseq-TCR"),
    c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-E12","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-F1","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-GEX-F12"),
    c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F2","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F3","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-TCR-F9"),
    c("NSCL-NR001-SCRNA-5P-NA-BRAIN-GEX","NSCL-NR001-SNSEQ-3P-NI-BRAIN-GEX","NSCL-NR001-SNSEQ-3P-WI-BRAIN-GEX","NSCL-NR001-SNSEQ-5P-NI-BRAIN-GEX","NSCL-NR001-SNSEQ-5P-WI-BRAIN-GEX","NSCL-NR001-5pv2-snseq"),
    c("NSCL-NR001-SCRNA-5P-NA-BRAIN-TCR","NSCL-NR001-SNSEQ-5P-NI-BRAIN-TCR","NSCL-NR001-SNSEQ-5P-WI-BRAIN-TCR"),
    c("ribas1_pre_5pv2","ribas1_on_5pv2","Ribas3103GEX"),
    c("Ribas-1Pre_S102","Ribas-1ON_S103","Ribas310dot3_S105"),
    c("um_07_gk_pre_S4_L001","um_07_gk_on_S8_L001","uv003-uvme-snseq-3p-post","um_08_ar_pre_S1_L001","um_08_ar_on_S2_L001","um_08_ar_post_S3_L001","um_09_mw_pre_S5_L001","um_09_mw_on_S6_L001","um_09_mw_post_S7_L001","um_11_lc_pre_S12_L002","um_11_lc_on_S16_L002","um_12_ml_pre_S9_L002","um_12_ml_on_S10_L002","um_12_ml_post_S11_L002","um_15_lm_pre_S13_L002","um_15_lm_on_S14_L002","um_15_lm_post_S15_L002","um_16_rs_pre_S17_L003","um_16_rs_on_S18_L003","um_16_rs_post_S19_L003"),
    c("UM_07_GK_pre_TCR","UM_07_GK_on_TCR","UM_08_AR_pre_TCR","UM_08_AR_on_TCR","UM_08_AR_post_TCR","UM_09_MW_pre_TCR","UM_09_MW_on_TCR","UM_09_MW_post_TCR","UM_11_LC_pre_TCR","UM_11_LC_on_TCR","UM_12_ML_pre_TCR","UM_12_ML_on_TCR","UM_12_ML_post_TCR","UM_15_LM_pre_TCR","UM_15_LM_on_TCR","UM_15_LM_post_TCR","UM_16_RS_pre_TCR","UM_16_RS_on_TCR","UM_16_RS_post_TCR"),
    c("uM_07_GK_Pre","uM_07_GK_ON","uM_08_AR_Pre","uM_08_AR_ON","uM_08_AR_POST","uM_09_MW_Pre","uM_09_MW_ON","uM_09_MW_POST","uM_11_LC_Pre","uM_11_LC_ON","uM_12_nL_Pre","uM_12_ML_ON","uM_12_ML_POST","uM_15_LM_Pre","uM_15_LM_ON","uM_15_LM_POST","uM_16_RS_Pre","uM_16_RS_ON","uM_16_RS_POST"),
    c("BI5CST","BI5TST","NR1CST","NR1TST"),
    c("JB_BI_demultiplexing_ECM08_pre_S11_L004","JB_BI_demultiplexing_ECM08_on_S12_L004","JB_BI_demultiplexing_ECM08_on_later_S13_L004"),
    c("HM-baseline_TGACCA"))

  shortnameslist = list(c("Mel_sc_5p_CD45_minus","Mel_sc_5p_CD45_plus","Mel_sn_3p","Mel_sn_5p_v2","Mel_sn_5p"),
    c("Mel_sc_5p_CD45_plus_TCR","Mel_sn_5p_TCR"),
    c("UM_sc_5p","UM_sc_5p_CD45_plus","UM_sn_5p_inhib"),
    c("UM_sc_5p_TCR","UM_sc_5p_CD45_plus_TCR","UM_sn_5p_inhib_TCR"),
    c("NSCLC_sc_5p","NSCLC_sn_3p","NSCLC_sn_3p_inhib","NSCLC_sn_5p","NSCLC_sn_5p_inhib","NSCLC_sn_5p_v2"),
    c("NSCLC_sc_5p_TCR","NSCLC_sn_5p_TCR","NSCLC_sn_5p_inhib_TCR"),
    c("ribas_310_pre","ribas_310_on","ribas_310_on_later"),
    c("ribas_310_pre_usWGS","ribas_310_on_usWGS","ribas_310_on_later_usWGS"),
    c("UMEL_1_1","UMEL_1_2","UMEL_1_3","UMEL_2_1","UMEL_2_2","UMEL_2_3","UMEL_3_1","UMEL_3_2","UMEL_3_3","UMEL_4_1","UMEL_4_2","UMEL_5_1","UMEL_5_2","UMEL_5_3","UMEL_6_1","UMEL_6_2","UMEL_6_3","UMEL_7_1","UMEL_7_2","UMEL_7_3"),
    c("UMEL_1_1_TCR","UMEL_1_2_TCR","UMEL_2_1_TCR","UMEL_2_2_TCR","UMEL_2_3_TCR","UMEL_3_1_TCR","UMEL_3_2_TCR","UMEL_3_3_TCR","UMEL_4_1_TCR","UMEL_4_2_TCR","UMEL_5_1_TCR","UMEL_5_2_TCR","UMEL_5_3_TCR","UMEL_6_1_TCR","UMEL_6_2_TCR","UMEL_6_3_TCR","UMEL_7_1_TCR","UMEL_7_2_TCR","UMEL_7_3_TCR"),
    c("UMEL_1_1_usWGS","UMEL_1_2_usWGS","UMEL_2_1_usWGS","UMEL_2_2_usWGS","UMEL_2_3_usWGS","UMEL_3_1_usWGS","UMEL_3_2_usWGS","UMEL_3_3_usWGS","UMEL_4_1_usWGS","UMEL_4_2_usWGS","UMEL_5_1_usWGS","UMEL_5_2_usWGS","UMEL_5_3_usWGS","UMEL_6_1_usWGS","UMEL_6_2_usWGS","UMEL_6_3_usWGS","UMEL_7_1_usWGS","UMEL_7_2_usWGS","UMEL_7_3_usWGS"),
    c("BI5CST","BI5TST","NR1CST","NR1TST"),
    c("ribas_310_pre_slide","ribas_310_on_slide","ribas_310_on_later_slide"),
    c("ribas_310_pre_WES"))

  ashortname = ""
  for (i in 1:length(longnameslist)) {
    for (j in 1:length(longnameslist[[i]]))
    {
      if (length(grep(paste0("^",longnameslist[[i]][j]),afilename))!=0)
      {
	ashortname = shortnameslist[[i]][j]
      }
    }
  }

  return(ashortname)
}
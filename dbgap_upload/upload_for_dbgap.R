library(stringr)

foldersList = c("s3://fresh-vs-frozen-comparison-ohio/BI5/scrna-seq/raw_data/raw_BI5",
  "s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq/raw_data",
  "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata/raw_data",
  "s3://fresh-vs-frozen-comparison-ohio/nsclc/raw_data",
  "s3://melanoma-ribas/raw_data/snRNA-seq",
  "s3://melanoma-ribas/raw_data/TCR-seq",
  "s3://melanoma-ribas/raw_data/WGS/fastq",
  "s3://uveal-melanoma/raw_data",
  "s3://uveal-melanoma/raw_data/new_seq_run",
  "s3://uveal-melanoma/raw_data/WGS/fastqs",
  "s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline/raw_data",
  "s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline/raw_data",
  "s3://fresh-vs-frozen-comparison-ohio/melanoma-ribas/raw_data",
  "s3://fresh-vs-frozen-comparison-ohio/melanoma-ribas/raw_data")

sampleFolderList = list(c("CD45negGEXBI5_S1_L001","CD45posGEXBI5_S1_L001","TCRBI5_S1_L001"),
  c("bi005-skcm-5snseq","bi005-skcm","skcm-bi005-5pv2-snseq","bi005-skcm-5snseq-TCR"),
  c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-E12/fastq/201208_A00521_0191_AHLGFCDSXY","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-F1/fastq/201208_A00521_0191_AHLGFCDSXY","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-GEX-F12/fastq/201208_A00521_0191_AHLGFCDSXY","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-TCR-F9/fastq/201208_A00521_0191_AHLGFCDSXY","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F2/fastq/201208_A00521_0191_AHLGFCDSXY","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F3/fastq/201208_A00521_0191_AHLGFCDSXY"),
  c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX","NSCL-NR001-5pv2-snseq","NSCL_NR001_SCRNA_5P_NA_BRAIN_TCR","NSCL_NR001_SNSEQ_5P_NI_BRAIN_TCR","NSCL_NR001_SNSEQ_5P_WI_BRAIN_TCR"),
  c("ribas1_on_5pv2","ribas1_pre_5pv2","Ribas3103GEX"),
  c("ribas1_on_tcr","ribas1_pre_tcr","Ribas3103TCRLib"),
  c("Ribas-1ON_S103","Ribas-1Pre_S102","Ribas310dot3_S105"),
  c("um_07_gk_on_S8_L001","um_07_gk_pre_S4_L001","um_08_ar_on_S2_L001","um_08_ar_post_S3_L001","um_08_ar_pre_S1_L001","um_09_mw_on_S6_L001","um_09_mw_post_S7_L001","um_09_mw_pre_S5_L001","um_11_lc_on_S16_L002","um_11_lc_pre_S12_L002","um_12_ml_on_S10_L002","um_12_ml_post_S11_L002","um_12_ml_pre_S9_L002","um_15_lm_on_S14_L002","um_15_lm_post_S15_L002","um_15_lm_pre_S13_L002","um_16_rs_on_S18_L003","um_16_rs_post_S19_L003","um_16_rs_pre_S17_L003","uv003-uvme-snseq-3p-post"),
  c("UM_07_GK_on_TCR","UM_07_GK_pre_TCR","UM_08_AR_on_TCR","UM_08_AR_post_TCR","UM_08_AR_pre_TCR","UM_09_MW_on_TCR","UM_09_MW_post_TCR","UM_09_MW_pre_TCR","UM_11_LC_on_TCR","UM_11_LC_pre_TCR","UM_12_ML_on_TCR","UM_12_ML_post_TCR","UM_12_ML_pre_TCR","UM_15_LM_on_TCR","UM_15_LM_post_TCR","UM_15_LM_pre_TCR","UM_16_RS_on_TCR","UM_16_RS_post_TCR","UM_16_RS_pre_TCR"),
  c("uM_07_GK_ON_S146_L001","uM_07_GK_Pre_S152_L001","uM_08_AR_ON_S158_L001","uM_08_AR_POST_S145_L001","uM_08_AR_Pre_S155_L001","uM_09_MW_ON_S148_L001","uM_09_MW_POST_S159_L001","uM_09_MW_Pre_S154_L001","uM_11_LC_ON_S149_L001","uM_11_LC_Pre_S150_L001","uM_12_ML_ON_S151_L001","uM_12_ML_POST_S153_L001","uM_12_nL_Pre_S156_L001","uM_15_LM_ON_S144_L001","uM_15_LM_POST_S157_L001","uM_15_LM_Pre_S161_L001","uM_16_RS_ON_S147_L001","uM_16_RS_POST_S162_L001","uM_16_RS_Pre_S160_L001"),
  c("BI5CST","BI5TST"),
  c("NR1CST","NR1TST"),
  c("JB_BI_demultiplexing_ECM08_pre_S11_L004","JB_BI_demultiplexing_ECM08_on_S12_L004","JB_BI_demultiplexing_ECM08_on_later_S13_L004"),
  c("PT0310_Baseline"))

# glacierrestorelist = read.table("glacier_restore_list_edited.txt",sep=" ",header=F,quote=NULL)
# restorelist = glacierrestorelist$V4
# for (i in 1:length(restorelist)) {
#   restoreparams = "'{\"Days\":25,\"GlacierJobParameters\":{\"Tier\":\"Standard\"}}'"
#   if (str_detect(restorelist[i],"ibas"))
#   {
#     system(paste0("aws s3api restore-object --bucket melanoma-ribas --key ",restorelist[i]," --restore-request ",restoreparams))
#     system(paste0("aws s3api head-object --bucket melanoma-ribas --key ",restorelist[i]))
#   }
#   else
#   {
#     system(paste0("aws s3api restore-object --bucket uveal-melanoma --key ",restorelist[i]," --restore-request ",restoreparams))
#     system(paste0("aws s3api head-object --bucket uveal-melanoma --key ",restorelist[i]))
#   }
# }

# nonsense = nonsense+1

uploaded_files = c()
#nonsense = nonsense+1
for (i in 7:7) {#1:length(sampleFolderList)) {
  for (j in 1:length(sampleFolderList[[i]])) {
    print(sampleFolderList[[i]][j])
    #system(paste0("aws s3 sync ",foldersList[i],"/",sampleFolderList[[i]][j]," ",sampleFolderList[[i]][j]))
    system(paste0("aws s3 sync ",foldersList[i],"/",sampleFolderList[[i]][j]," ",sampleFolderList[[i]][j]," --force-glacier-transfer"))
    #aws_list_of_files = system(paste0("aws s3 ls ",foldersList[i],"/",sampleFolderList[[i]][j]), intern = TRUE)
    list_of_files = system(paste0("ls ",sampleFolderList[[i]][j]), intern = TRUE)
    for (afile in list_of_files)
    {
      if (str_detect(afile,"_R1_") || str_detect(afile,"_R2_") || str_detect(afile,"R1.fastq") || str_detect(afile,"R2.fastq"))
      {
        uploaded_files = c(uploaded_files, afile)
        system(paste0("ascp -i private.openssh -Q -l 1000m -k 1 ",sampleFolderList[[i]][j],"/",afile," asp-sra@gap-submit.ncbi.nlm.nih.gov:protected"))
      }
    }
    system(paste0("rm -r ",sampleFolderList[[i]][j]))
  }
}

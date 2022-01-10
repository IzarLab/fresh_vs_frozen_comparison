library(ggplot2)
library(Seurat)
library(stringr)
library(rlang)
library(rlist)

foldersList = c("s3://fresh-vs-frozen-comparison-ohio/BI5/scrna-seq",
  "s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq",
  "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata",
  "s3://fresh-vs-frozen-comparison-ohio/nsclc",
  "s3://melanoma-ribas/ribas1",
  "s3://uveal-melanoma")
integrated_name_arr = c("BI5_scrna-seq","BI5_snrna-seq","cpoi-uvealprimarydata","nsclc","ribas_integrated_titrate_thresh","um_all")
integrated_name_arr_underscore = c("Mel_scrna_seq","Mel_snrna_seq","UM","NSCLC","ribas","UMEL")

source("/mnt/vdb/home/ubuntu2/uveal-melanoma-code/fresh_vs_frozen_QC_stress_sigs.R")

ribas_ifng_sig = read.table("ribas_ifng_sig.txt", header = F, sep = ",")
stress_sig_list = list.append(stress_sig_list, ribas_ifng_sig)
names(stress_sig_list) = c(names(stress_sig_list), "ribas_ifng_sig")

object.list = c()
for (i in 1:2) {
  system(paste0("aws s3 cp ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
  integrated_rds = readRDS(paste0("/data/",integrated_name_arr[i],"_integrated.rds"))
  DefaultAssay(integrated_rds) = "RNA"
  for (stress_sig in names(stress_sig_list))
  {
    integrated_rds = AddModuleScore(integrated_rds, features = list(na.omit(stress_sig_list[[stress_sig]])), name = stress_sig, assay = "RNA", search = T)
  }
  object.list = c(object.list, integrated_rds)
}

makePubFigures = TRUE

source("merge.SCTAssay.R")

pdf(paste0("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_QC/Mel_fresh_vs_frozen_QC.pdf"),height=7,width=5.5*length(unique(seu$orig.ident)))

uniqueidents = unique(seu$orig.ident)
seu$orig.ident[seu$orig.ident=="CD45neg"] = "Mel_sc_5_CD45-"
seu$orig.ident[seu$orig.ident=="CD45pos"] = "Mel_sc_5_CD45+"
seu$orig.ident[seu$orig.ident=="3snseq"] = "Mel_sn_3"
seu$orig.ident[seu$orig.ident=="5pv2-snseq"] = "Mel_sn_5v2"
seu$orig.ident[seu$orig.ident=="5snseq"] = "Mel_sn_5"

seu$orig.ident = paste0(seu$orig.ident,"\n(",seu$fresh_frozen,")")

seu$Mel_nCount_RNA = seu$nCount_RNA
seu$Mel_nFeature_RNA = seu$nFeature_RNA
seu$Mel_ScrubDoublet_score = seu$ScrubDoublet_score
seu$Mel_percent.mt = seu$percent.mt
seu$Mel_stress_sig_dysfunctional_cd81 = seu$stress_sig_dysfunctional_cd81
seu$Mel_stress_sig_nmeth_celseq1 = seu$stress_sig_nmeth_celseq1
seu$Mel_stress_sig_nmeth_sortseq_cluster11 = seu$stress_sig_nmeth_sortseq_cluster11
seu$Mel_stress_sig_nmeth_sortseq_cluster41 = seu$stress_sig_nmeth_sortseq_cluster41
seu$Mel_stress_sig_genomebiol_cryopreserve1 = seu$stress_sig_genomebiol_cryopreserve1
seu$Mel_stress_sig_brain_met1 = seu$stress_sig_brain_met1
seu$Mel_ribas_ifng_sig1 = seu$ribas_ifng_sig1
seu$Mel_num_zero_count_genes = colSums(seu@assays$RNA@counts==0)

seu$orig.ident[seu$fresh_frozen=="fresh"] = paste0("1_",seu$orig.ident[seu$fresh_frozen=="fresh"])
seu$orig.ident[seu$fresh_frozen=="frozen"] = paste0("2_",seu$orig.ident[seu$fresh_frozen=="frozen"])

relabel_list = unique(seu$orig.ident)
names(relabel_list) = relabel_list
for (i in 1:length(relabel_list))
{
  relabel_list[i] = substring(relabel_list[i],3)
}

aplot = VlnPlot(seu, features = c("Mel_nFeature_RNA","Mel_percent.mt","Mel_stress_sig_nmeth_celseq1"), group.by = "orig.ident",pt.size = 0)
uniqueidents = unique(seu$orig.ident)
ymax_arr = c(13000,25,3)
ylabs_arr = c("# detected genes/cell","percent mitochondrial reads","expression of signature")
for (z2 in 1:3) {
  colorsarr = rep("",length(uniqueidents))
  colorsarr[grep("fresh",uniqueidents)] = "blue"
  colorsarr[grep("frozen",uniqueidents)] = "red"
  aplot[[z2]] = aplot[[z2]] + scale_fill_manual(values = colorsarr, limits = uniqueidents) + scale_x_discrete(labels = relabel_list) + ylim(0,ymax_arr[z2]) + ylab(ylabs_arr[z2]) + theme(axis.text.x = element_text(angle=0, hjust=0.5))
}
AugmentPlot(aplot, dpi = 300)
print(aplot)

atable = table(seu$orig.ident)
cellnumbersdf = data.frame(sample = names(atable), count = atable)

testarrlong = list()
testarr = list()
fresh_ids = unique(seu$orig.ident[seu$fresh_frozen=="fresh"])
frozen_ids = unique(seu$orig.ident[seu$fresh_frozen=="frozen"])
testarr_labels = c()
for (i in 1:length(fresh_ids))
{
  for (j in 1:length(frozen_ids))
  {
    test1 = wilcox.test(seu$nFeature_RNA[seu$orig.ident==fresh_ids[i]],seu$nFeature_RNA[seu$orig.ident==frozen_ids[j]])
    test2 = wilcox.test(seu$percent.mt[seu$orig.ident==fresh_ids[i]],seu$percent.mt[seu$orig.ident==frozen_ids[j]])
    test3 = wilcox.test(seu$stress_sig_nmeth_celseq1[seu$orig.ident==fresh_ids[i]],seu$stress_sig_nmeth_celseq1[seu$orig.ident==frozen_ids[j]])
    testarr1 = c(test1$p.value, test2$p.value, test3$p.value)
    names(testarr1) = c("nfeature test","percentmt test","stress sig test")
    testarr = list.append(testarr, testarr1)
    testarr_labels = c(testarr_labels, paste0(fresh_ids[i]," ",frozen_ids[i]))
  }
}
names(testarr) = testarr_labels
print(testarr)
testarrlong = list.append(testarrlong, testarr)

dev.off()

qualitydf = data.frame(ident=character(), stat=character(), value=double())
uniqueidents = unique(seu$orig.ident)
for (uniqueident in uniqueidents)
{
  print(uniqueident)
  tempdf = data.frame(ident=uniqueident, stat="median_nFeature_RNA", value=median(seu$nFeature_RNA[seu$orig.ident==uniqueident]))
  qualitydf = rbind(qualitydf, tempdf)
  tempdf = data.frame(ident=uniqueident, stat="median_percent.mt", value=median(seu$percent.mt[seu$orig.ident==uniqueident]))
  qualitydf = rbind(qualitydf, tempdf)
  tempdf = data.frame(ident=uniqueident, stat="median_stress_sig", value=median(seu$stress_sig_nmeth_celseq1[seu$orig.ident==uniqueident]))
  qualitydf = rbind(qualitydf, tempdf)
}

for (i in 3:length(foldersList)) {
  system(paste0("aws s3 cp ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
  integrated_rds = readRDS(paste0("/data/",integrated_name_arr[i],"_integrated.rds"))
  DefaultAssay(integrated_rds) = "RNA"
  if (sum(names(integrated_rds@meta.data)=="fresh_frozen")==0)
  {
    integrated_rds$fresh_frozen = "frozen"
  }
  if (integrated_name_arr[i]=="nsclc")
  {
    integrated_rds$fresh_frozen[integrated_rds$fresh_frozen=="5pv2-snseq"] = "frozen"
  }

  if (integrated_name_arr_underscore[i]=="ribas")
  {
    integrated_rds$placeholder = FALSE
    integrated_rds$placeholder[grep("310",integrated_rds$orig.ident)] = TRUE
    integrated_rds = subset(integrated_rds, placeholder)
  }

  pdf(paste0("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_QC/",integrated_name_arr_underscore[i],"_fresh_vs_frozen_QC.pdf"),height=7,width=5.5*length(unique(integrated_rds$orig.ident)))

  uniqueidents = unique(integrated_rds$orig.ident)
  integrated_rds$orig.ident[integrated_rds$orig.ident=="SCRNA-5P-NA-E12"] = "UM_sc_5"
  integrated_rds$orig.ident[integrated_rds$orig.ident=="SCRNA-5P-NA-F1"] = "UM_sc_5_CD45+"
  integrated_rds$orig.ident[integrated_rds$orig.ident=="SNRNA-5P-WI-F12"] = "UM_sn_5_inhib"
  integrated_rds$orig.ident[integrated_rds$orig.ident=="SCRNA_5P_NA"] = "NSCLC_sc_5"
  integrated_rds$orig.ident[integrated_rds$orig.ident=="5pv2-snseq"] = "NSCLC_sn_5v2"
  integrated_rds$orig.ident[integrated_rds$orig.ident=="SNSEQ_3P_NI"] = "NSCLC_sn_3"
  integrated_rds$orig.ident[integrated_rds$orig.ident=="SNSEQ_3P_WI"] = "NSCLC_sn_3_inhib"
  integrated_rds$orig.ident[integrated_rds$orig.ident=="SNSEQ_5P_NI"] = "NSCLC_sn_5"
  integrated_rds$orig.ident[integrated_rds$orig.ident=="SNSEQ_5P_WI"] = "NSCLC_sn_5_inhib"
  integrated_rds$orig.ident[integrated_rds$orig.ident=="ribas_310_pre_GEX_5pv2_S26_L004"] = "pre"
  integrated_rds$orig.ident[integrated_rds$orig.ident=="ribas_310_on_GEX_5pv2_S27_L004"] = "on"
  integrated_rds$orig.ident[integrated_rds$orig.ident=="ribas_310_on_later_previd_3_GEX"] = "on_later"
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
  # for (i1 in 1:length(uniqueidents))
  # {
  #   newident = uniqueidents[i1]
  #   newident = str_replace(newident,"BI5_","")
  #   newident = str_replace(newident,"nsclc_","")
  #   newident = str_replace(newident,"cpoi-uvealprimarydata_","")
  #   newident = str_replace(newident,"Sarcoma","")
  #   newident = str_replace(newident,"GEX","")
  #   newident = str_replace(newident,"_GEX","")
  #   newident = str_replace(newident,"ribas_","")
  #   newident = str_replace(newident,"um_","")
  #   newident = str_replace(newident,"_S.*_L.*","")
  #   newident = str_replace(newident,"uv003-uvme-snseq-3p-post","07_gk_post")
  #   integrated_rds$orig.ident[integrated_rds$orig.ident==uniqueidents[i1]] = newident
  # }

  integrated_rds$orig.ident = paste0(integrated_rds$orig.ident,"\n(",integrated_rds$fresh_frozen,")")
  integrated_rds = AddModuleScore(integrated_rds, features = list(na.omit(stress_sig_dysfunctional_cd8)), name = "stress_sig_dysfunctional_cd8", assay = "RNA", search = T)
  integrated_rds = AddModuleScore(integrated_rds, features = list(na.omit(stress_sig_nmeth_celseq)), name = "stress_sig_nmeth_celseq", assay = "RNA", search = T)
  integrated_rds = AddModuleScore(integrated_rds, features = list(na.omit(stress_sig_nmeth_sortseq_cluster1)), name = "stress_sig_nmeth_sortseq_cluster1", assay = "RNA", search = T)
  integrated_rds = AddModuleScore(integrated_rds, features = list(na.omit(stress_sig_nmeth_sortseq_cluster4)), name = "stress_sig_nmeth_sortseq_cluster4", assay = "RNA", search = T)
  integrated_rds = AddModuleScore(integrated_rds, features = list(na.omit(stress_sig_genomebiol_cryopreserve)), name = "stress_sig_genomebiol_cryopreserve", assay = "RNA", search = T)
  integrated_rds = AddModuleScore(integrated_rds, features = list(na.omit(stress_sig_brain_met)), name = "stress_sig_brain_met", assay = "RNA", search = T)
  integrated_rds = AddModuleScore(integrated_rds, features = list(na.omit(ribas_ifng_sig)), name = "ribas_ifng_sig", assay = "RNA", search = T)

  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_nCount_RNA = integrated_rds$nCount_RNA")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_nFeature_RNA = integrated_rds$nFeature_RNA")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_ScrubDoublet_score = integrated_rds$ScrubDoublet_score")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_percent.mt = integrated_rds$percent.mt")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_stress_sig_dysfunctional_cd81 = integrated_rds$stress_sig_dysfunctional_cd81")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_stress_sig_nmeth_celseq1 = integrated_rds$stress_sig_nmeth_celseq1")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_stress_sig_nmeth_sortseq_cluster11 = integrated_rds$stress_sig_nmeth_sortseq_cluster11")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_stress_sig_nmeth_sortseq_cluster41 = integrated_rds$stress_sig_nmeth_sortseq_cluster41")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_stress_sig_genomebiol_cryopreserve1 = integrated_rds$stress_sig_genomebiol_cryopreserve1")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_stress_sig_brain_met1 = integrated_rds$stress_sig_brain_met1")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_ribas_ifng_sig1 = integrated_rds$ribas_ifng_sig1")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_num_zero_count_genes = colSums(integrated_rds@assays$RNA@counts==0)")))

  integrated_rds$orig.ident[integrated_rds$fresh_frozen=="fresh"] = paste0("1_",integrated_rds$orig.ident[integrated_rds$fresh_frozen=="fresh"])
  integrated_rds$orig.ident[integrated_rds$fresh_frozen=="frozen"] = paste0("2_",integrated_rds$orig.ident[integrated_rds$fresh_frozen=="frozen"])

  relabel_list = unique(integrated_rds$orig.ident)
  names(relabel_list) = relabel_list
  for (i1 in 1:length(relabel_list))
  {
    relabel_list[i1] = substring(relabel_list[i1],3)
  }

  aplot = VlnPlot(integrated_rds, features = c(paste0(integrated_name_arr_underscore[i],"_nFeature_RNA"),paste0(integrated_name_arr_underscore[i],"_percent.mt"),paste0(integrated_name_arr_underscore[i],"_stress_sig_nmeth_celseq1")), group.by = "orig.ident",pt.size = 0)
  uniqueidents = sort(unique(integrated_rds$orig.ident))
  ymax_arr = c(13000,25,3)
  ylabs_arr = c("# detected genes/cell","percent mitochondrial reads","expression of signature")
  for (z2 in 1:3)#length(aplot))
  {
    colorsarr = rep("",length(uniqueidents))
    colorsarr[grep("fresh",uniqueidents)] = "blue"
    colorsarr[grep("frozen",uniqueidents)] = "red"
    if (integrated_name_arr[i]=="um_all")
    {
      colorsarr[grep("UMEL_1_3",uniqueidents)] = "yellow"
    }
    if (integrated_name_arr[i]=="nsclc" || integrated_name_arr[i]=="ribas_integrated_titrate_thresh")
    {
      if (integrated_name_arr[i]=="nsclc")
      {
	uniqueidents_reorder = uniqueidents
	uniqueidents_reorder = uniqueidents_reorder[c(1,6,2,3,4,5)]
      }
      if (integrated_name_arr[i]=="ribas_integrated_titrate_thresh")
      {
	uniqueidents_reorder = uniqueidents
	uniqueidents_reorder = uniqueidents_reorder[c(3,1,2)]
      }
      aplot[[z2]] = aplot[[z2]] + scale_fill_manual(values = colorsarr, limits = uniqueidents) + scale_x_discrete(limits = uniqueidents_reorder, labels=relabel_list) + ylim(0,ymax_arr[z2]) + ylab(ylabs_arr[z2]) + theme(axis.text.x = element_text(angle=0, hjust=0.5))
    }
    else
    {
      aplot[[z2]] = aplot[[z2]] + scale_fill_manual(values = colorsarr, limits = uniqueidents) + ylim(0,ymax_arr[z2]) + ylab(ylabs_arr[z2]) + scale_x_discrete(labels = relabel_list) + theme(axis.text.x = element_text(angle=0, hjust=0.5))
    }
  }
  AugmentPlot(aplot, dpi = 300)
  print(aplot)

  system(paste0("rm /data/",integrated_name_arr[i],"_integrated.rds"))

  atable = table(integrated_rds$orig.ident)
  adf = data.frame(sample = names(atable), count = atable)
  cellnumbersdf = rbind(cellnumbersdf, adf)

  if (length(unique(integrated_rds$fresh_frozen))!=1)
  {
    testarr = list()
    fresh_ids = unique(integrated_rds$orig.ident[integrated_rds$fresh_frozen=="fresh"])
    frozen_ids = unique(integrated_rds$orig.ident[integrated_rds$fresh_frozen=="frozen"])
    testarr_labels = c()
    for (i1 in 1:length(fresh_ids))
    {
      for (j1 in 1:length(frozen_ids))
      {
	test1 = wilcox.test(integrated_rds$nFeature_RNA[integrated_rds$orig.ident==fresh_ids[i1]],integrated_rds$nFeature_RNA[integrated_rds$orig.ident==frozen_ids[j1]])
	test2 = wilcox.test(integrated_rds$percent.mt[integrated_rds$orig.ident==fresh_ids[i1]],integrated_rds$percent.mt[integrated_rds$orig.ident==frozen_ids[j1]])
	test3 = wilcox.test(integrated_rds$stress_sig_nmeth_celseq1[integrated_rds$orig.ident==fresh_ids[i1]],integrated_rds$stress_sig_nmeth_celseq1[integrated_rds$orig.ident==frozen_ids[j1]])
	testarr1 = c(test1$p.value, test2$p.value, test3$p.value)
	names(testarr1) = c("nfeature test","percentmt test","stress sig test")
	testarr = list.append(testarr, testarr1)
	testarr_labels = c(testarr_labels, paste0(fresh_ids[i1]," ",frozen_ids[j1]))
      }
    }
    names(testarr) = testarr_labels
    print(testarr)
    testarrlong = list.append(testarrlong, testarr)
  }
  dev.off()

  print(integrated_name_arr[i])
  uniqueidents = unique(integrated_rds$orig.ident)
  for (uniqueident in uniqueidents)
  {
    print(uniqueident)
    tempdf = data.frame(ident=uniqueident, stat="median_nFeature_RNA", value=median(integrated_rds$nFeature_RNA[integrated_rds$orig.ident==uniqueident]))
    qualitydf = rbind(qualitydf, tempdf)
    tempdf = data.frame(ident=uniqueident, stat="median_percent.mt", value=median(integrated_rds$percent.mt[integrated_rds$orig.ident==uniqueident]))
    qualitydf = rbind(qualitydf, tempdf)
    tempdf = data.frame(ident=uniqueident, stat="median_stress_sig", value=median(integrated_rds$stress_sig_nmeth_celseq1[integrated_rds$orig.ident==uniqueident]))
    qualitydf = rbind(qualitydf, tempdf)
  }
}

names(stress_medians) = stress_medians_names

pdf("fresh_vs_frozen_cell_count.pdf")
theme_set(theme_bw())
print(ggplot(cellnumbersdf, aes(y=count.Freq, x=sample)) + geom_bar(stat="identity") + ylab("Count") + theme(axis.text.x = element_text(angle=45, hjust=1)))

dev.off()
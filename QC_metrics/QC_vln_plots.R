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

#load in stress signature genes
source("/mnt/vdb/home/ubuntu2/uveal-melanoma-code/fresh_vs_frozen_QC_stress_sigs.R")

ribas_ifng_sig = read.table("ribas_ifng_sig.txt", header = F, sep = ",")
stress_sig_list = list.append(stress_sig_list, ribas_ifng_sig)
names(stress_sig_list) = c(names(stress_sig_list), "ribas_ifng_sig")

#load data for BI5 scrna and snrna-seq samples, calculate stress signatures
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

#merge together BI5 scrna and snrna-seq samples
source("merge.SCTAssay.R")
seu = merge.SCTAssay(x=object.list[1], y=object.list[2], merge.data = T, na.rm = T)

pdf(paste0("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_QC/Mel_fresh_vs_frozen_QC.pdf"),height=7,width=5.5*length(unique(seu$orig.ident)))

#rename BI5 sample names to shorter names for display in figure, add description of whether they are fresh or frozen
uniqueidents = unique(seu$orig.ident)
source("fresh_vs_frozen_rename_IDs.R")
seu = rename_IDs(seu)

seu$orig.ident = paste0(seu$orig.ident,"\n(",seu$fresh_frozen,")")

#add new QC metric metadata fields, with name of dataset attached
seu$Mel_nCount_RNA = seu$nCount_RNA
seu$Mel_nFeature_RNA = seu$nFeature_RNA
seu$Mel_ScrubDoublet_score = seu$ScrubDoublet_score
seu$Mel_percent.mt = seu$percent.mt
for (stress_sig in names(stress_sig_list))
{
  eval(parse(text=paste0("seu$Mel_",stress_sig,"1 = seu$",stress_sig,"1")))
}

#prepend either 1 or 2 to sample names, depending on fresh/frozen status, to ensure that fresh samples are displayed first
seu$orig.ident[seu$fresh_frozen=="fresh"] = paste0("1_",seu$orig.ident[seu$fresh_frozen=="fresh"])
seu$orig.ident[seu$fresh_frozen=="frozen"] = paste0("2_",seu$orig.ident[seu$fresh_frozen=="frozen"])

#clip off 1 or 2 from sample names, to create array of labels that will actually be displayed in figure
relabel_list = unique(seu$orig.ident)
names(relabel_list) = relabel_list
for (i in 1:length(relabel_list))
{
  relabel_list[i] = substring(relabel_list[i],3)
}

#violin plot figure, using manually encoded upper limits for violins, and fresh samples colored blue, frozen red
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

#count number of cells that are present in each sample
atable = table(seu$orig.ident)
cellnumbersdf = data.frame(sample = names(atable), count = atable)

#use wilcoxon test to test for differences in QC metrics between each pair of fresh and frozen samples
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

#store median QC metrics for each sample in a dataframe
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

#repeat above steps for remaining datasets
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

  #for ribas samples, only look at sample 310
  if (integrated_name_arr_underscore[i]=="ribas")
  {
    integrated_rds$placeholder = FALSE
    integrated_rds$placeholder[grep("310",integrated_rds$orig.ident)] = TRUE
    integrated_rds = subset(integrated_rds, placeholder)
  }

  pdf(paste0("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_QC/",integrated_name_arr_underscore[i],"_fresh_vs_frozen_QC.pdf"),height=7,width=5.5*length(unique(integrated_rds$orig.ident)))

  uniqueidents = unique(integrated_rds$orig.ident)
  source("fresh_vs_frozen_rename_IDs.R")
  integrated_rds = rename_IDs(integrated_rds)

  integrated_rds$orig.ident = paste0(integrated_rds$orig.ident,"\n(",integrated_rds$fresh_frozen,")")
  for (stress_sig in names(stress_sig_list))
  {
    integrated_rds = AddModuleScore(integrated_rds, features = list(na.omit(stress_sig_list[[stress_sig]])), name = stress_sig, assay = "RNA", search = T)
  }

  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_nCount_RNA = integrated_rds$nCount_RNA")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_nFeature_RNA = integrated_rds$nFeature_RNA")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_ScrubDoublet_score = integrated_rds$ScrubDoublet_score")))
  eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_percent.mt = integrated_rds$percent.mt")))
  for (stress_sig in names(stress_sig_list))
  {
    eval(parse(text=paste0("integrated_rds$",integrated_name_arr_underscore[i],"_",stress_sig,"1 = integrated_rds$",stress_sig,"1")))
  }
  
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
    #label 3p sample in UMEL dataset as yellow
    if (integrated_name_arr[i]=="um_all")
    {
      colorsarr[grep("UMEL_1_3",uniqueidents)] = "yellow"
    }
    #reorder nsclc and ribas samples, so that fresh samples are plotted first
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

#plot total number of cells in each sample
pdf("fresh_vs_frozen_cell_count.pdf")
theme_set(theme_bw())
print(ggplot(cellnumbersdf, aes(y=count.Freq, x=sample)) + geom_bar(stat="identity") + ylab("Count") + theme(axis.text.x = element_text(angle=45, hjust=1)))

dev.off()
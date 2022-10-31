library(ggplot2)
library(Seurat)
library(stringr)
library(rlang)
library(rlist)

#define list of s3 folders to download data from, names of each dataset, and alternate names that replace dashes with underscores
foldersList = c("s3://fresh-vs-frozen-comparison-ohio/BI5/scrna-seq",
  "s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq",
  "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata",
  "s3://fresh-vs-frozen-comparison-ohio/nsclc",
  "s3://melanoma-ribas/ribas1",
  "s3://uveal-melanoma",
  "s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline",
  "s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
integrated_name_arr = c("BI5_scrna-seq","BI5_snrna-seq","cpoi-uvealprimarydata","nsclc","ribas_integrated_titrate_thresh","um_all","BI5","NR1")
integrated_name_arr_underscore = c("Mel_scrna_seq","Mel_snrna_seq","UM","NSCLC","ribas","UMEL","BI5","NR1")

#set options for using downsampled data from s3, and showing all stress signatures versus just one
use_downsampled = FALSE
possible_downsampled_suffix = ""
if (use_downsampled) {
  possible_downsampled_suffix = "_downsampled"
}
heightParam = 7
show_all_stress_sigs = TRUE
if (show_all_stress_sigs) {
  heightParam = 21
}

#read in table of ensembl gene ids mapped to hgnc ids, remove ensembl ids that map to more than one hgnc id
ensembl_gene_to_hgnc = read.table("mart_export_ensembl_gene_to_hgnc.txt",quote=NULL,header=T,sep="\t")
dup_ids_temp = table(ensembl_gene_to_hgnc$HGNC_symbol)
dup_ids = names(dup_ids_temp)[dup_ids_temp>1]
ensembl_gene_to_hgnc = subset(ensembl_gene_to_hgnc, !(HGNC_symbol %in% dup_ids))
rownames(ensembl_gene_to_hgnc) = ensembl_gene_to_hgnc[["HGNC_symbol"]]
ensembl_gene_to_hgnc[["HGNC_symbol"]] = NULL

#load in stress signature genes
source("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_comparison/QC_metrics/load_stress_sigs.R")

#either show all stress signature results, or just results from nature methods paper for celseq
if (show_all_stress_sigs) {
  stress_sigs_to_plot = stress_sig_names
} else {
  stress_sigs_to_plot = c("stress_sig_nmeth_celseq")
}

#read in interferon-gamma signature from ribas data and add to stress_signature list
ribas_ifng_sig = read.table("ribas_ifng_sig.txt", header = F, sep = ",")
stress_sig_list = list.append(stress_sig_list, ribas_ifng_sig)
names(stress_sig_list)[length(stress_sig_list)] = "ribas_ifng_sig"
stress_sig_list_orig = stress_sig_list
#if using downsampled data, and therefore measuring ribas stress signature, change ensembl gene ids to hgnc
if (use_downsampled)
{
  for (stress_sig in names(stress_sig_list))
  {
    if (stress_sig=="ribas_ifng_sig")
    {
      stress_sig_list[[stress_sig]] = ensembl_gene_to_hgnc[stress_sig_list[[stress_sig]]$V1,]
    }
    else
    {
      stress_sig_list[[stress_sig]] = ensembl_gene_to_hgnc[stress_sig_list[[stress_sig]],]
    }
  }
}

cell_type_median_features = c()
cell_type_median_counts = c()
cell_type_median_percentmt = c()
cell_type_median_stress_sig = c()
cell_type_total_features = c()

#load data for BI5 scrna and snrna-seq samples, calculate stress signatures
object.list = c()
for (i in 1:2) {
  if (use_downsampled)
  {
    system(paste0("aws s3 cp ",foldersList[i],"/Seurat_downsampled/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
  }
  else
  {
    system(paste0("aws s3 cp ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
  }
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

pdf(paste0("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_QC/Mel_fresh_vs_frozen_QC",possible_downsampled_suffix,".pdf"),height=heightParam,width=5.5*length(unique(seu$orig.ident)))

#rename BI5 sample names to shorter names for display in figure, add description of whether they are fresh or frozen
uniqueidents = unique(seu$orig.ident)
source("fresh_vs_frozen_comparison/QC_metrics/rename_sample_IDs.R")
seu = rename_IDs(seu,"Mel")

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
aplot = VlnPlot(seu, features = c("Mel_nFeature_RNA","Mel_percent.mt",paste0("Mel_",stress_sigs_to_plot,"1")), group.by = "orig.ident",pt.size = 0)# + geom_boxplot()
uniqueidents = unique(seu$orig.ident)
ymax_arr = c(13000,25,3)
ylabs_arr = c("# detected genes/cell","percent mitochondrial reads","expression of signature")
for (z2 in 1:(2+length(stress_sigs_to_plot))) {
  colorsarr = rep("",length(uniqueidents))
  colorsarr[grep("fresh",uniqueidents)] = "blue"
  colorsarr[grep("frozen",uniqueidents)] = "red"
  aplot[[z2]] = aplot[[z2]] + scale_fill_manual(values = colorsarr, limits = uniqueidents) + scale_x_discrete(labels = relabel_list) + ylim(0,ymax_arr[z2]) + ylab(ylabs_arr[z2]) + theme(axis.text.x = element_text(angle=0, hjust=0.5)) + geom_boxplot(fill="white",width=0.1)
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

seu = readRDS("/data/fresh_vs_frozen_all_reannotate_BI5.rds")
for (stress_sig in names(stress_sig_list))
{
  seu = AddModuleScore(seu, features = list(na.omit(stress_sig_list[[stress_sig]])), name = stress_sig, assay = "RNA", search = T)
}

#calculate QC median metrics for each cell type in rds object, as well as for all non-tumor cells generally
for (an_orig_ident in unique(seu$orig.ident)) {
  seu_sub = subset(seu, orig.ident==an_orig_ident)
  cell_types_arr = unique(seu_sub$manual_annotation_label)
  cell_types_arr = c(cell_types_arr, "non-tumor")
  for (a_cell_type in cell_types_arr)
  {
    if (a_cell_type=="non-tumor")
    {
      seu_sub$temparr = !(seu_sub$manual_annotation_label %in% c("Melanoma cells","Lung cancer cells"))
      seu_sub2 = subset(seu_sub, temparr)
    }
    else
    {
      seu_sub2 = subset(seu_sub, manual_annotation_label==a_cell_type)
    }
    cell_type_median_features = c(cell_type_median_features, median(seu_sub2$nFeature_RNA))
    names(cell_type_median_features)[length(cell_type_median_features)] = paste0(an_orig_ident,"_",a_cell_type)

    cell_type_median_counts = c(cell_type_median_counts, median(seu_sub2$nCount_RNA))
    names(cell_type_median_counts)[length(cell_type_median_counts)] = paste0(an_orig_ident,"_",a_cell_type)

    cell_type_median_percentmt = c(cell_type_median_percentmt, median(seu_sub2$percent.mt))
    names(cell_type_median_percentmt)[length(cell_type_median_percentmt)] = paste0(an_orig_ident,"_",a_cell_type)

    cell_type_median_stress_sig = c(cell_type_median_stress_sig, median(seu_sub2$stress_sig_nmeth_celseq1))
    names(cell_type_median_stress_sig)[length(cell_type_median_stress_sig)] = paste0(an_orig_ident,"_",a_cell_type)

    cell_type_total_features = c(cell_type_total_features, sum(rowSums(seu_sub2@assays$RNA@counts)!=0))
    names(cell_type_total_features)[length(cell_type_total_features)] = paste0(an_orig_ident,"_",a_cell_type)
  }
}

#make violin plot of cell-type specific quality metrics
seu$orig.ident.bare = unlist(lapply(strsplit(seu$orig.ident,"\n"), function(x) {x[[1]][1]}))
seu$cell_type_and_sample = paste0(seu$orig.ident.bare,"\n",seu$manual_annotation_label)
pdf(paste0("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_QC/Mel_cell_type_QC.pdf"),height=3*heightParam,width=3*length(unique(seu$cell_type_and_sample)))
aplot = VlnPlot(seu, features = c("nFeature_RNA", "percent.mt", paste0(stress_sigs_to_plot,"1")), group.by = "cell_type_and_sample", pt.size = 0, ncol=1)
for (z2 in 1:3)
{
  aplot[[z2]] = aplot[[z2]] + ylim(0,ymax_arr[z2]) + ylab(ylabs_arr[z2]) + theme(axis.text.x = element_text(angle=0, hjust=0.5))
}
print(aplot)
dev.off()

source("merge.SCTAssay.R")
#repeat above steps for remaining datasets
for (i in 3:(length(foldersList))) {
  if (use_downsampled)
  {
    system(paste0("aws s3 cp ",foldersList[i],"/Seurat_downsampled/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
  }
  else
  {
    system(paste0("aws s3 cp ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
  }
  integrated_rds = readRDS(paste0("/data/",integrated_name_arr[i],"_integrated.rds"))

  if (integrated_name_arr_underscore[i]=="NR1" || integrated_name_arr_underscore[i]=="BI5")
  {
    integrated_rds$fresh_frozen = "frozen"
  }

  if (integrated_name_arr[i]=="BI5")
  {
    integrated_rds = subset(integrated_rds, (orig.ident=="BI5CST" | orig.ident=="BI5TST"))
    #system("aws s3 cp s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq/Seurat/bi005-skcm_final_thresh/bi005-skcm_final_thresh_cb.rds /data/bi005-skcm_final_thresh_cb.rds")
    if (use_downsampled)
    {
      system("aws s3 cp s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq/Seurat_downsampled/bi005-skcm_final_thresh/bi005-skcm_final_thresh_cb.rds /data/bi005-skcm_final_thresh_cb.rds")
      bi5obj = readRDS("/data/bi005-skcm_final_thresh_cb.rds")
      system("rm /data/bi005-skcm_final_thresh_cb.rds")
    }
    else
    {
      bi5obj = readRDS("/data/fresh_vs_frozen_all_reannotate_BI5.rds")
    }
    integrated_rds = merge.SCTAssay(x=integrated_rds, y=bi5obj, merge.data = T, na.rm = T)
  }

  if (integrated_name_arr[i]=="NR1")
  {
    integrated_rds = subset(integrated_rds, (orig.ident=="NR1CST" | orig.ident=="NR1TST"))
    #system("aws s3 cp s3://fresh-vs-frozen-comparison-ohio/nsclc/Seurat/NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX_final_thresh/NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX_final_thresh_cb.rds /data/NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX_final_thresh_cb.rds")
    if (use_downsampled)
    {
      system("aws s3 cp s3://fresh-vs-frozen-comparison-ohio/nsclc/Seurat_downsampled/integrated/nsclc_integrated.rds /data/nsclc_integrated.rds")
    }
    else
    {
      system("aws s3 cp s3://fresh-vs-frozen-comparison-ohio/nsclc/Seurat/integrated/nsclc_integrated.rds /data/nsclc_integrated.rds")
    }
    nsclcobj = readRDS("/data/nsclc_integrated.rds")
    system("rm /data/nsclc_integrated.rds")
    integrated_rds = merge.SCTAssay(x=integrated_rds, y=nsclcobj, merge.data = T, na.rm = T)
  }

  DefaultAssay(integrated_rds) = "RNA"
  if (sum(names(integrated_rds@meta.data)=="fresh_frozen")==0)
  {
    integrated_rds$fresh_frozen = "frozen"
  }
  if (integrated_name_arr[i]=="nsclc" || integrated_name_arr[i]=="NR1")
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

  pdf(paste0("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_QC/",integrated_name_arr_underscore[i],"_fresh_vs_frozen_QC",possible_downsampled_suffix,".pdf"),height=heightParam,width=5.5*length(unique(integrated_rds$orig.ident)))

  uniqueidents = unique(integrated_rds$orig.ident)
  source("fresh_vs_frozen_comparison/QC_metrics/rename_sample_IDs.R")
  integrated_rds = rename_IDs(integrated_rds,integrated_name_arr_underscore[i])

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

  aplot = VlnPlot(integrated_rds, features = c(paste0(integrated_name_arr_underscore[i],"_nFeature_RNA"),paste0(integrated_name_arr_underscore[i],"_percent.mt"),paste0(integrated_name_arr_underscore[i],"_",stress_sigs_to_plot,"1")), group.by = "orig.ident",pt.size = 0)
  uniqueidents = sort(unique(integrated_rds$orig.ident))
  ymax_arr = c(13000,25,3)
  ylabs_arr = c("# detected genes/cell","percent mitochondrial reads","expression of signature")
  for (z2 in 1:(2+length(stress_sigs_to_plot)))#length(aplot))
  {
    colorsarr = rep("",length(uniqueidents))
    colorsarr[grep("fresh",uniqueidents)] = "blue"
    colorsarr[grep("frozen",uniqueidents)] = "red"
    #label 3p sample in UMEL dataset as yellow
    if (integrated_name_arr[i]=="um_all")
    {
      colorsarr[grep("UMEL_1_3",uniqueidents)] = "yellow"
    }
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
    if (integrated_name_arr[i]=="nsclc" || integrated_name_arr[i]=="ribas_integrated_titrate_thresh" || integrated_name_arr[i]=="BI5" || integrated_name_arr[i]=="NR1")
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
      aplot[[z2]] = aplot[[z2]] + scale_fill_manual(values = colorsarr, limits = uniqueidents) + scale_x_discrete(limits = uniqueidents_reorder, labels=relabel_list) + ylim(0,ymax_arr[z2]) + ylab(ylabs_arr[z2]) + theme(axis.text.x = element_text(angle=0, hjust=0.5)) + geom_boxplot(fill="white",width=0.1)
    }
    else
    {
      aplot[[z2]] = aplot[[z2]] + scale_fill_manual(values = colorsarr, limits = uniqueidents) + ylim(0,ymax_arr[z2]) + ylab(ylabs_arr[z2]) + scale_x_discrete(labels = relabel_list) + theme(axis.text.x = element_text(angle=0, hjust=0.5)) + geom_boxplot(fill="white",width=0.1)
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

  for (an_orig_ident in unique(integrated_rds$orig.ident)) {
    integrated_rds_sub = subset(integrated_rds, orig.ident==an_orig_ident)
    cell_types_arr = unique(integrated_rds_sub$manual_annotation_label)
    cell_types_arr = c(cell_types_arr, "non-tumor")
    for (a_cell_type in cell_types_arr)
    {
      if (a_cell_type=="non-tumor")
      {
        integrated_rds_sub$temparr = !(integrated_rds_sub$manual_annotation_label %in% c("Melanoma cells","Lung cancer cells"))
        integrated_rds_sub2 = subset(integrated_rds_sub, temparr)
      }
      else
      {
        integrated_rds_sub2 = subset(integrated_rds_sub, manual_annotation_label==a_cell_type)
      }
      cell_type_median_features = c(cell_type_median_features, median(integrated_rds_sub2$nFeature_RNA))
      names(cell_type_median_features)[length(cell_type_median_features)] = paste0(an_orig_ident,"_",a_cell_type)

      cell_type_median_counts = c(cell_type_median_counts, median(integrated_rds_sub2$nCount_RNA))
      names(cell_type_median_counts)[length(cell_type_median_counts)] = paste0(an_orig_ident,"_",a_cell_type)

      cell_type_median_percentmt = c(cell_type_median_percentmt, median(integrated_rds_sub2$percent.mt))
      names(cell_type_median_percentmt)[length(cell_type_median_percentmt)] = paste0(an_orig_ident,"_",a_cell_type)

      cell_type_median_stress_sig = c(cell_type_median_stress_sig, median(integrated_rds_sub2$stress_sig_nmeth_celseq1))
      names(cell_type_median_stress_sig)[length(cell_type_median_stress_sig)] = paste0(an_orig_ident,"_",a_cell_type)

      cell_type_total_features = c(cell_type_total_features, sum(rowSums(integrated_rds_sub2@assays$RNA@counts)!=0))
      names(cell_type_total_features)[length(cell_type_total_features)] = paste0(an_orig_ident,"_",a_cell_type)
    }
  }

  integrated_rds$orig.ident.bare = unlist(lapply(strsplit(str_sub(integrated_rds$orig.ident,start=3),"\n"), function(x) {x[[1]][1]}))
  if (integrated_name_arr_underscore[i]=="BI5")
  {
    integrated_rds$temparr = (integrated_rds$orig.ident.bare=="BI5CST" | integrated_rds$orig.ident.bare=="BI5TST")
    integrated_rds = subset(integrated_rds, temparr)
  }
  if (integrated_name_arr_underscore[i]=="NR1")
  {
    integrated_rds$temparr = (integrated_rds$orig.ident.bare=="NR1CST" | integrated_rds$orig.ident.bare=="NR1TST")
    integrated_rds = subset(integrated_rds, temparr)
  }
  integrated_rds$cell_type_and_sample = paste0(integrated_rds$orig.ident.bare,"\n",integrated_rds$manual_annotation_label)
  pdf(paste0("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_QC/",integrated_name_arr_underscore[i],"_cell_type_QC.pdf"),height=3*heightParam,width=3*length(unique(integrated_rds$cell_type_and_sample)))
  aplot = VlnPlot(integrated_rds, features = c("nFeature_RNA", "percent.mt", paste0(stress_sigs_to_plot,"1")), group.by = "cell_type_and_sample", pt.size = 0, ncol=1)
  for (z2 in 1:3)
  {
    aplot[[z2]] = aplot[[z2]] + ylim(0,ymax_arr[z2]) + ylab(ylabs_arr[z2]) + theme(axis.text.x = element_text(angle=0, hjust=0.5))
  }
  print(aplot)
  dev.off()
}

cell_type_table = data.frame(cell_type_and_sample = names(cell_type_median_features), median_features = cell_type_median_features, median_counts = cell_type_median_counts, median_percentmt = cell_type_median_percentmt, median_stress_sig = cell_type_median_stress_sig, total_features = cell_type_total_features)
write.table(cell_type_table, "cell_type_QC_metrics.csv", sep=",", row.names=F, col.names=T)

names(stress_medians) = stress_medians_names

#plot total number of cells in each sample
pdf(paste0("fresh_vs_frozen_cell_count",possible_downsampled_suffix,".pdf"))
theme_set(theme_bw())
print(ggplot(cellnumbersdf, aes(y=count.Freq, x=sample)) + geom_bar(stat="identity") + ylab("Count") + theme(axis.text.x = element_text(angle=45, hjust=1)))

dev.off()
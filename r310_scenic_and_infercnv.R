library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(infercnv)
library(stringr)
library(rlist)
library(fgsea)
library(stringr)

pats = c("ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh_kmeans","ribas_310_on_later_previd_3_GEX_titrate_thresh_kmeans","ribas_310_pre_GEX_5pv2_S26_L004_titrate_thresh_kmeans")
cbpats = c("ribas_310_on_GEX_5pv2_S27_L004_titrate_thresh","ribas_310_on_later_previd_3_GEX_titrate_thresh","ribas_310_pre_GEX_5pv2_S26_L004_titrate_thresh")
shortpats = c("on","on_later","pre")
s3folderlist = c("melanoma-ribas/ribas1","melanoma-ribas/ribas1","melanoma-ribas/ribas1")
cancertissuelist = list(c("Melanocytes"),c("Melanocytes"),c("Melanocytes"))

obj_arr = c()
for (i in 1:length(pats)) {
  #download Seurat rds objects and infercnv objects for each sample, and subset to cancer cells
  pat = pats[i]
  cbpat = cbpats[i]
  system(paste0("aws s3 cp s3://",s3folderlist[i],"/Seurat/",cbpat,"/",cbpat,"_cb.rds /data/",cbpat,"_cb.rds"))
  system(paste0("aws s3 cp s3://",s3folderlist[i],"/inferCNV/without_plots/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))

  orig_obj = readRDS(paste0("/data/",cbpat,"_cb.rds"))
  infercnv_obj = readRDS(paste0("/data/",pat,"_infercnv.rds"))
  orig_obj$infercnv = "not inferred"
  system(paste0("rm /data/",cbpat,"_cb.rds"))
  system(paste0("rm /data/",pat,"_infercnv.rds"))

  orig_obj = subset(orig_obj, (celltype_bped_main %in% cancertissuelist[[i]]))

  #add infercnv subcluster information, containing kmeans results, to rds object
  for (subcluster in names(infercnv_obj@tumor_subclusters$subclusters))
  {
    if (subcluster %in% cancertissuelist[[i]])
    {
      for (subcluster1 in names(infercnv_obj@tumor_subclusters$subclusters[[subcluster]]))
      {
	matchidxs = match(colnames(infercnv_obj@expr.data)[infercnv_obj@tumor_subclusters$subclusters[[subcluster]][[subcluster1]]],colnames(orig_obj))
	matchidxs = na.omit(matchidxs)
	orig_obj$infercnv[matchidxs] = subcluster1
      }

    }
  }

  #add field distinguishing between kmeans_2 and non_kmeans 2 cells
  orig_obj$infercnv_broad = "not inferred"
  orig_obj$infercnv_broad[orig_obj$infercnv %in% c("kmeans_2")] = "kmeans_2"
  orig_obj$infercnv_broad[orig_obj$infercnv %in% c("kmeans_0","kmeans_1","kmeans_3")] = "not_kmeans_2"

  obj_arr = c(obj_arr, orig_obj)
}

names(obj_arr) = shortpats

#load scenic results for each sample from csv file
scenic_arr = list()
for (shortpat in shortpats) {
  scenic_table = read.table(paste0("/data/r310_",shortpat,"_auc_mtx.csv"), sep=",", quote="", header=T)
  rownames(scenic_table) = scenic_table[,1]
  scenic_table = scenic_table[,-1]
  rownames(scenic_table) = unlist(lapply(str_split(rownames(scenic_table),"_"), function(x) {x[[1]]}))
  rownames(scenic_table) = str_replace_all(rownames(scenic_table), "\\.", "-")
  scenic_arr = list.append(scenic_arr, scenic_table)
}
names(scenic_arr) = shortpats

#perform wilcoxon test for AUCell score differences between clone 2 and non-clone 2 cells
wilc_pval_arr = c()
wilc_fc_arr = c()
wilc_nonzero_kmeans_2_arr = c()
wilc_nonzero_not_kmeans_2_arr = c()
for (shortpat in shortpats) {
  scenic_table = scenic_arr[[shortpat]]
  infercnv_broad_arr = obj_arr[[shortpat]]$infercnv_broad
  orig_obj_colnames = colnames(obj_arr[[shortpat]])
  for (atf in colnames(scenic_table))
  {
    tf_col = scenic_table[[atf]]
    names(tf_col) = rownames(scenic_table)
    tf_col_kmeans_2 = tf_col[orig_obj_colnames[infercnv_broad_arr=="kmeans_2"]]
    tf_col_not_kmeans_2 = tf_col[orig_obj_colnames[infercnv_broad_arr=="not_kmeans_2"]]

    tf_col_kmeans_2 = tf_col_kmeans_2[!is.na(names(tf_col_kmeans_2))]
    tf_col_not_kmeans_2 = tf_col_not_kmeans_2[!is.na(names(tf_col_not_kmeans_2))]

    test_res = wilcox.test(tf_col_kmeans_2, tf_col_not_kmeans_2)
    wilc_pval_arr = c(wilc_pval_arr, test_res$p.value)
    names(wilc_pval_arr)[length(wilc_pval_arr)] = paste0(shortpat,"_",atf)
    wilc_fc_arr = c(wilc_fc_arr, mean(tf_col_kmeans_2)/mean(tf_col_not_kmeans_2))
    names(wilc_fc_arr)[length(wilc_fc_arr)] = paste0(shortpat,"_",atf)
    wilc_nonzero_kmeans_2_arr = c(wilc_nonzero_kmeans_2_arr, sum(tf_col_kmeans_2!=0)/length(tf_col_kmeans_2))
    names(wilc_nonzero_kmeans_2_arr)[length(wilc_nonzero_kmeans_2_arr)] = paste0(shortpat,"_",atf)
    wilc_nonzero_not_kmeans_2_arr = c(wilc_nonzero_not_kmeans_2_arr, sum(tf_col_not_kmeans_2!=0)/length(tf_col_not_kmeans_2))
    names(wilc_nonzero_not_kmeans_2_arr)[length(wilc_nonzero_not_kmeans_2_arr)] = paste0(shortpat,"_",atf)
    # nonsense = nonsense+1
  }
}

#filter for tfs that are significant according to Bonferroni-corrected p-val, and have logfc > .25, across all three timepoints
all_wilc_tfs = names(wilc_pval_arr)
pre_tfs = c()
on_tfs = c()
on_later_tfs = c()
for (i in 1:length(all_wilc_tfs)) {
  if (wilc_pval_arr[[i]]<=.05/length(wilc_pval_arr) && !is.na(wilc_fc_arr[[i]]) && abs(log(wilc_fc_arr[[i]]))>=0.25)
  {
    words = str_split(all_wilc_tfs[i],"_")[[1]]
    if (words[1]=="pre")
    {
      pre_tfs = c(pre_tfs, paste0(words[length(words)-1],"_",words[length(words)]))
    }
    else if (words[1]=="on" && words[2]!="later")
    {
      on_tfs = c(on_tfs, paste0(words[length(words)-1],"_",words[length(words)]))
    }
    else if (words[1]=="on" && words[2]=="later")
    {
      on_later_tfs = c(on_later_tfs, paste0(words[length(words)-1],"_",words[length(words)]))
    }
  }
}

intersect_tfs = intersect(intersect(pre_tfs,on_tfs),on_later_tfs)
intersect_tfs_temp = c()
intersect_fcs = c()
for (i in 1:length(intersect_tfs))
{
  astring = intersect_tfs[i]

  intersect_fcs = c(intersect_fcs, log((wilc_fc_arr[[paste0("pre_",astring)]]+wilc_fc_arr[[paste0("on_",astring)]]+wilc_fc_arr[[paste0("on_later_",astring)]])/3))

  astring = str_sub(astring,end=length(astring)-6)
  intersect_tfs_temp = c(intersect_tfs_temp,astring)
}
intersect_tfs = intersect_tfs_temp
fgsea_input = intersect_fcs#rep(1,length(intersect_tfs))
names(fgsea_input) = intersect_tfs

#perform gsea on tfs significant in all three timepoints, using hallmarks pathways
Hfgsea = gmtPathways("h.all.v7.4.symbols.gmt")

result_fgsea = fgsea(Hfgsea,fgsea_input,minSize=-1,maxSize=Inf,scoreType="pos",nperm=1000)

dotplotdf = data.frame(timepoint = character(), tf = character(), logfc = double(), pval = double())

#construct dataframe containing all wilcoxon
for (i in 1:length(intersect_tfs))
{
  astring = intersect_tfs[i]

  times = c("pre","on","on_later")
  for (atime in times)
  {
    tempdf = data.frame(timepoint=atime, tf=astring, logfc=log(wilc_fc_arr[[paste0(atime,"_",astring,"_pos")]]), pval=wilc_pval_arr[[paste0(atime,"_",astring,"_pos")]])
    dotplotdf = rbind(dotplotdf, tempdf)
  }
}

half_min_nonzero = min(dotplotdf$pval[dotplotdf$pval!=0])/2
dotplotdf$pval[dotplotdf$pval==0] = half_min_nonzero

theme_set(theme_bw())

dotplot_avglogfc = c()
dotplot_tfs = unique(dotplotdf$tf)
for (atf in dotplot_tfs)
{
  dotplot_avglogfc = c(dotplot_avglogfc, mean(dotplotdf$logfc[dotplotdf$tf==atf]))
}
names(dotplot_avglogfc) = dotplot_tfs
dotplot_avglogfc = sort(dotplot_avglogfc, decreasing=F)

dotplotdftemp = data.frame(timepoint = character(), tf = character(), logfc = double(), pval = double())
for (atf in names(dotplot_avglogfc))
{
  times = c("pre","on","on_later")
  for (atime in times)
  {
    tempdf = subset(dotplotdf, timepoint==atime & tf==atf)
    dotplotdftemp = rbind(dotplotdftemp, tempdf)
  }
}
dotplotdf = dotplotdftemp

dotplotdf$timepoint = factor(dotplotdf$timepoint, levels = c("pre","on","on_later"))
dotplotdf$tf = factor(dotplotdf$tf, levels = unique(dotplotdf$tf))

dotplotdf$minus_log_pval = -log(dotplotdf$pval)

dotplotdf[["Log FC"]] = dotplotdf$logfc
dotplotdf[["-Log 2 P-val"]] = dotplotdf$minus_log_pval
ggplot(dotplotdf, aes(x = timepoint, y = tf)) + scale_colour_gradientn(colours = c("blue","blue","gray","red","red"), breaks = c(min(dotplotdf$logfc), quantile(dotplotdf$logfc,0.05)[[1]], 0, quantile(dotplotdf$logfc,0.95))[[1]], max(dotplotdf$logfc)) + geom_point(aes(color = `Log FC`, size = `-Log 2 P-val`)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("TF") + xlab("Time Point") + labs(color = "Log FC") + labs(minus_log_pval = "-Log 2 P-val")
ggsave(paste0("r310_scenic_and_infercnv.pdf"), width = 5, height = 7)

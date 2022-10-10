library(hypeR)
library(rlist)
library(fgsea)
library(zoo)
library(pheatmap)
library(stringr)
library(ggplot2)

useRollingAverage = TRUE

wgspats = c("uM-07-GK-ON_S3","uM-07-GK-Pre_S9",
"uM-08-AR-ON_S15","uM-08-AR-POST_S2","uM-08-AR-Pre_S12",
"uM-09-MW-ON_S5","uM-09-MW-POST_S16","uM-09-MW-Pre_S11",
"uM-11-LC-ON_S6","uM-11-LC-Pre_S7",
"uM-12-ML-ON_S8","uM-12-ML-POST_S10","uM-12-nL-Pre_S13",
"uM-15-LM-ON_S1","uM-15-LM-POST_S14",
"uM-16-RS-ON_S4")

tumorfractions = c(.2528,.8325,
.4263,.8426,.3277,
.4621,.06156,.6152,
.8629,.8366,
.6126,.03006,.02325,
.7217,.7818,
.9085)

infercnvpats = c("um_07_gk_on_S8_L001","um_07_gk_pre_S4_L001","um_08_ar_on_S2_L001","um_08_ar_post_S3_L001","um_08_ar_pre_S1_L001","um_09_mw_on_S6_L001","um_09_mw_post_S7_L001","um_09_mw_pre_S5_L001","um_11_lc_on_S16_L002","um_11_lc_pre_S12_L002","um_12_ml_on_S10_L002","um_12_ml_post_S11_L002","um_12_ml_pre_S9_L002","um_15_lm_on_S14_L002","um_15_lm_post_S15_L002","um_16_rs_on_S18_L003")

allsingleinfercnvarr_across_samples = c()
allsinglewgsarr_across_samples = c()
printoutinfercnvdf = data.frame(sample=character(), chr = character(), start = integer(), end = integer(), infercnv_mean_score = double(), wgs_score = double())

for (infercnvpat in infercnvpats)
{
  pat = infercnvpat
  if (pat=="uv003_uvme_snseq_3p_post")
  {
    pattemp = "uv003-uvm3-snseq-3p-post"
    #system(paste0("aws s3 cp s3://uveal-melanoma/Seurat/",pattemp,"/",pattemp,"_cb.rds ",pat,"_cb.rds"))
    system(paste0("aws s3 cp s3://uveal-melanoma/inferCNV/inferCNV_cellbender_subcluster_",pattemp,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
    system(paste0("aws s3 cp s3://uveal-melanoma/inferCNV/without_plots/inferCNV_cellbender_subcluster_",pattemp,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
  }
  else
  {
    #system(paste0("aws s3 cp s3://uveal-melanoma/Seurat/",pat,"/",pat,"_cb.rds ",pat,"_cb.rds"))
    system(paste0("aws s3 cp s3://uveal-melanoma/inferCNV/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
    system(paste0("aws s3 cp s3://uveal-melanoma/inferCNV/without_plots/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
  }
}

allgenes = c()
for (z in 1:length(wgspats))
{
  infercnvpat = infercnvpats[z]
  infercnv_obj = readRDS(paste0("/data/",infercnvpat,"_infercnv.rds"))
  allgenes = sort(unique(c(allgenes,rownames(infercnv_obj@gene_order))))
}

# gene_table_df_all = data.frame(genename=allgenes)
# for (z in 1:length(infercnvpats))
# {
#   tempdf = data.frame(x=rep("NA",length(allgenes)))
#   names(tempdf) = c(infercnvpats[z])
#   gene_table_df_all = cbind(gene_table_df_all,tempdf)
# }

# #nonsense = nonsense+1

# for (z in 1:length(wgspats))
# {
#   print(infercnvpats[z])
#   wgspat = wgspats[z]
#   infercnvpat = infercnvpats[z]

#   infercnv_obj = readRDS(paste0("/data/",infercnvpat,"_infercnv.rds"))

#   chrs = unique(infercnv_obj@gene_order$chr)
#   expr.data.temp = matrix(0,dim(infercnv_obj@expr.data)[1],dim(infercnv_obj@expr.data)[2])
#   rownames(expr.data.temp) = rownames(infercnv_obj@expr.data)
#   colnames(expr.data.temp) = colnames(infercnv_obj@expr.data)
#   #gene_order_temp = data.frame(chr=character(), start = integer(), end = integer())
#   #rownames(gene_order_temp) = rownames(infercnv_obj@gene_order)

#   if (useRollingAverage)
#   {
#     for (z1 in 1:dim(infercnv_obj@expr.data)[2])
#     {
#       for (chr in chrs)
#       {
# 	expr.data.temp[infercnv_obj@gene_order$chr==chr,z1] = rollmean(infercnv_obj@expr.data[infercnv_obj@gene_order$chr==chr,z1],k=11,fill=NA)
#       }
#     }
#     infercnv_obj@gene_order = infercnv_obj@gene_order[!is.na(expr.data.temp[,1]),]
#     #expr.data.temp = na.omit(expr.data.temp)
#     #infercnv_obj@expr.data = expr.data.temp
#     infercnv_obj@expr.data = expr.data.temp[!is.na(expr.data.temp[,1]),]
#     print("roll")
#   }

#   lineardata = infercnv_obj@expr.data[1:(dim(infercnv_obj@expr.data)[1]*dim(infercnv_obj@expr.data)[2])]

#   oldmedian = median(lineardata)
#   infercnv_obj@expr.data[infercnv_obj@expr.data>.88 & infercnv_obj@expr.data<1.12] = oldmedian

#   wgsdata = read.table(paste0("/mnt/vdb/home/ubuntu2/",wgspat,".seg"),header=T,sep="\t",quote=NULL)
#   infercnv_obj_melanocyte = infercnv_obj#readRDS(paste0("/mnt/vdb/home/ubuntu2/",infercnvpat,"_infercnv.rds"))
#   infercnv_obj_melanocyte@expr.data = infercnv_obj_melanocyte@expr.data[,infercnv_obj_melanocyte@observation_grouped_cell_indices[["Melanocytes"]]]
#   infercnv_obj_adipocyte = infercnv_obj#readRDS(paste0("/mnt/vdb/home/ubuntu2/",infercnvpat,"_infercnv.rds"))
#   infercnv_obj_adipocyte@expr.data = infercnv_obj_adipocyte@expr.data[,infercnv_obj_adipocyte@observation_grouped_cell_indices[["Adipocytes"]]]
#   infercnv_obj_rest = infercnv_obj#readRDS(paste0("/mnt/vdb/home/ubuntu2/",infercnvpat,"_infercnv.rds"))
#   restidxs = c()
#   for (celltype in names(infercnv_obj_rest@observation_grouped_cell_indices))
#   {
#     if (celltype!="Melanocytes")
#     {
#       restidxs = c(restidxs,infercnv_obj_rest@observation_grouped_cell_indices[[celltype]])
#     }
#   }
#   infercnv_obj_rest@expr.data = infercnv_obj_rest@expr.data[,restidxs]

#   #nonsense = nonsense+1

#   gene_table_df = data.frame(genename = rownames(infercnv_obj@gene_order),segment="NA")

#   allinfercnvarr = c()
#   allwgsarr = c()
#   allsingleinfercnvarr = c()
#   allsinglewgsarr = c()

#   for (i in 1:length(wgsdata$sample))
#   {
#     wgsarr = rep(wgsdata$median[i],wgsdata$bins[i])
#     infercnvarr_melanocyte = rowMeans(infercnv_obj_melanocyte@expr.data)[infercnv_obj_melanocyte@gene_order$chr==paste0("chr",wgsdata$chr[i]) & infercnv_obj_melanocyte@gene_order$start>wgsdata$start[i] & infercnv_obj_melanocyte@gene_order$stop<wgsdata$end[i]]
#     infercnvarr_adipocyte = rowMeans(infercnv_obj_adipocyte@expr.data)[infercnv_obj_adipocyte@gene_order$chr==paste0("chr",wgsdata$chr[i]) & infercnv_obj_adipocyte@gene_order$start>wgsdata$start[i] & infercnv_obj_adipocyte@gene_order$stop<wgsdata$end[i]]
#     infercnvarr_rest = rowMeans(infercnv_obj_rest@expr.data)[infercnv_obj_rest@gene_order$chr==paste0("chr",wgsdata$chr[i]) & infercnv_obj_rest@gene_order$start>wgsdata$start[i] & infercnv_obj_rest@gene_order$stop<wgsdata$end[i]]

#     gene_table_df$segment[infercnv_obj@gene_order$chr==paste0("chr",wgsdata$chr[i]) & infercnv_obj@gene_order$start>wgsdata$start[i] & infercnv_obj@gene_order$stop<wgsdata$end[i]] = paste0(wgsdata$chr[i],"_",wgsdata$start[i],"_",wgsdata$end[i])

#     #infercnvarr = infercnvarr_melanocyte*tumorfractions[z] + infercnvarr_rest*(1-tumorfractions[z])
#     infercnvarr = infercnvarr_melanocyte

#     if (length(wgsarr)!=0 && length(infercnvarr)!=0)
#     {
#       if (length(wgsarr)<length(infercnvarr))
#       {
# 	infercnvarr = infercnvarr[sample(1:length(infercnvarr),length(wgsarr),replace=F)]
#       }
#       if (length(infercnvarr)<length(wgsarr))
#       {
# 	wgsarr = wgsarr[sample(1:length(wgsarr),length(infercnvarr),replace=F)]
#       }
#       allwgsarr = c(allwgsarr,wgsarr)
#       allinfercnvarr = c(allinfercnvarr,infercnvarr)
#       allsinglewgsarr = c(allsinglewgsarr,median(wgsarr))
#       allsingleinfercnvarr = c(allsingleinfercnvarr,median(infercnvarr))
#       #nonsense = nonsense+1

#       tempdf = data.frame(sample=wgsdata$sample[i], chr = wgsdata$chr[i], start = wgsdata$start[i], end = wgsdata$end[i], infercnv_mean_score = median(infercnvarr), wgs_score = median(wgsarr))
#       printoutinfercnvdf = rbind(printoutinfercnvdf,tempdf)
#     }
#   }

#   gene_table_df_all[match(gene_table_df$genename,gene_table_df_all$genename),names(gene_table_df_all)==infercnvpats[z]] = gene_table_df$segment

#   #nonsense = nonsense+1

#   print(cor.test(allinfercnvarr,allwgsarr,method="spearman")$estimate)
#   print(cor.test(allsingleinfercnvarr,allsinglewgsarr,method="spearman")$estimate)
#   rho = cor.test(allsingleinfercnvarr,allsinglewgsarr,method="spearman")$estimate
#   tempdf = data.frame(x=allsingleinfercnvarr,y=allsinglewgsarr)
#   pdf(paste0("/mnt/vdb/home/ubuntu2/correlate_wgs_and_infercnv/",infercnvpat,"_wgs_vs_infercnv.pdf"),width=9,height=7)
#   print(ggplot(tempdf, aes(x, y)) + geom_point(size=3) + xlab("Mean InferCNV Score") + ylab("WGS Median Log Relative Copy Number") + ggtitle(paste0(infercnvpats[z]," Spearman's R^2 = ",rho^2)) + theme_set(theme_bw()) + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 20)))
#   dev.off()

#   allsingleinfercnvarr_across_samples = c(allsingleinfercnvarr_across_samples,allsingleinfercnvarr)
#   allsinglewgsarr_across_samples = c(allsinglewgsarr_across_samples,allsinglewgsarr)
# }

# nonsense = nonsense+1

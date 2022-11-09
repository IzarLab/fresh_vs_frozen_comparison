library(hypeR)
library(rlist)
library(fgsea)
library(zoo)
library(pheatmap)
library(stringr)
library(ggplot2)
library(infercnv)

### title: Correlate average arm-level infercnv and whole-genome sequencing copy number alterations in uveal melanoma metastatic dataset
### author: Yiping Wang date: 11/08/2022

useRollingAverage = FALSE
flattenToMedian = FALSE
skipIndividualCalc = FALSE

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

chrlengths = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)
centromere_locations = c(123400000,93900000,90900000,50000000,48800000,59800000,60100000,45200000,43000000,39800000,53400000,35500000,17700000,17200000,19000000,36800000,25100000,18500000,26200000,28100000,12000000,15000000)
names(chrlengths) = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
names(centromere_locations) = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")

if (!skipIndividualCalc) {
  infercnvarr_across_segments_median_across_samples = c()
  wgsarr_across_segments_median_across_samples = c()
  print_infercnvarr_wgsarr_across_segments_median_df = data.frame(sample=character(), chr = character(), arm = character(), start = integer(), end = integer(), infercnv_mean_score = double(), wgs_score = double())

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

  gene_to_segment_df_across_samples = data.frame(genename=allgenes)
  for (z in 1:length(infercnvpats))
  {
    tempdf = data.frame(x=rep("NA",length(allgenes)))
    names(tempdf) = c(infercnvpats[z])
    gene_to_segment_df_across_samples = cbind(gene_to_segment_df_across_samples,tempdf)
  }

  for (z in 1:length(wgspats))
  {
    print(infercnvpats[z])
    wgspat = wgspats[z]
    infercnvpat = infercnvpats[z]

    infercnv_obj = readRDS(paste0("/data/",infercnvpat,"_infercnv.rds"))

    chrs = unique(infercnv_obj@gene_order$chr)
    expr.data.temp = matrix(0,dim(infercnv_obj@expr.data)[1],dim(infercnv_obj@expr.data)[2])
    rownames(expr.data.temp) = rownames(infercnv_obj@expr.data)
    colnames(expr.data.temp) = colnames(infercnv_obj@expr.data)

    if (useRollingAverage)
    {
      for (z1 in 1:dim(infercnv_obj@expr.data)[2])
      {
	for (chr in chrs)
	{
	  expr.data.temp[infercnv_obj@gene_order$chr==chr,z1] = rollmean(infercnv_obj@expr.data[infercnv_obj@gene_order$chr==chr,z1],k=11,fill=NA)
	}
      }
      infercnv_obj@gene_order = infercnv_obj@gene_order[!is.na(expr.data.temp[,1]),]
      infercnv_obj@expr.data = expr.data.temp[!is.na(expr.data.temp[,1]),]
      print("roll")
    }

    lineardata = infercnv_obj@expr.data[1:(dim(infercnv_obj@expr.data)[1]*dim(infercnv_obj@expr.data)[2])]

    oldmedian = median(lineardata)
    if (flattenToMedian)
    {
      infercnv_obj@expr.data[infercnv_obj@expr.data>.88 & infercnv_obj@expr.data<1.12] = oldmedian
    }

    wgsdata = read.table(paste0("/mnt/vdb/home/ubuntu2/",wgspat,".seg"),header=T,sep="\t",quote=NULL)
    infercnv_obj_melanocyte = infercnv_obj
    infercnv_obj_melanocyte@expr.data = infercnv_obj_melanocyte@expr.data[,infercnv_obj_melanocyte@observation_grouped_cell_indices[["Melanocytes"]]]
    infercnv_obj_adipocyte = infercnv_obj
    infercnv_obj_adipocyte@expr.data = infercnv_obj_adipocyte@expr.data[,infercnv_obj_adipocyte@observation_grouped_cell_indices[["Adipocytes"]]]
    infercnv_obj_rest = infercnv_obj
    restidxs = c()
    for (celltype in names(infercnv_obj_rest@observation_grouped_cell_indices))
    {
      if (celltype!="Melanocytes")
      {
	restidxs = c(restidxs,infercnv_obj_rest@observation_grouped_cell_indices[[celltype]])
      }
    }
    infercnv_obj_rest@expr.data = infercnv_obj_rest@expr.data[,restidxs]

    gene_to_segment_df = data.frame(genename = rownames(infercnv_obj@gene_order),segment="NA")

    infercnvarr_across_segments = c()
    wgsarr_across_segments = c()
    infercnvarr_across_segments_median = c()
    wgsarr_across_segments_median = c()

    for (i in 1:length(wgsdata$sample))
    {
      wgsarr = rep(wgsdata$median[i],wgsdata$bins[i])
      infercnvarr_melanocyte = rowMeans(infercnv_obj_melanocyte@expr.data)[infercnv_obj_melanocyte@gene_order$chr==paste0("chr",wgsdata$chr[i]) & infercnv_obj_melanocyte@gene_order$start>wgsdata$start[i] & infercnv_obj_melanocyte@gene_order$stop<wgsdata$end[i]]
      infercnvarr_adipocyte = rowMeans(infercnv_obj_adipocyte@expr.data)[infercnv_obj_adipocyte@gene_order$chr==paste0("chr",wgsdata$chr[i]) & infercnv_obj_adipocyte@gene_order$start>wgsdata$start[i] & infercnv_obj_adipocyte@gene_order$stop<wgsdata$end[i]]
      infercnvarr_rest = rowMeans(infercnv_obj_rest@expr.data)[infercnv_obj_rest@gene_order$chr==paste0("chr",wgsdata$chr[i]) & infercnv_obj_rest@gene_order$start>wgsdata$start[i] & infercnv_obj_rest@gene_order$stop<wgsdata$end[i]]

      gene_to_segment_df$segment[infercnv_obj@gene_order$chr==paste0("chr",wgsdata$chr[i]) & infercnv_obj@gene_order$start>wgsdata$start[i] & infercnv_obj@gene_order$stop<wgsdata$end[i]] = paste0(wgsdata$chr[i],"_",wgsdata$start[i],"_",wgsdata$end[i])

      #infercnvarr = infercnvarr_melanocyte*tumorfractions[z] + infercnvarr_rest*(1-tumorfractions[z])
      infercnvarr = infercnvarr_melanocyte

      if (length(wgsarr)!=0 && length(infercnvarr)!=0)
      {
	if (length(wgsarr)<length(infercnvarr))
	{
	  infercnvarr = infercnvarr[sample(1:length(infercnvarr),length(wgsarr),replace=F)]
	}
	if (length(infercnvarr)<length(wgsarr))
	{
	  wgsarr = wgsarr[sample(1:length(wgsarr),length(infercnvarr),replace=F)]
	}
	wgsarr_across_segments = c(wgsarr_across_segments,wgsarr)
	infercnvarr_across_segments = c(infercnvarr_across_segments,infercnvarr)
	wgsarr_across_segments_median = c(wgsarr_across_segments_median,median(wgsarr))
	infercnvarr_across_segments_median = c(infercnvarr_across_segments_median,median(infercnvarr))

	centromere_location = centromere_locations[[paste0("chr",wgsdata$chr[i])]]
	if (wgsdata$start[i]<=centromere_location && wgsdata$end[i]>=centromere_location)
	{
	  anarm = paste0("chr",wgsdata$chr[i],"p"," and ","chr",wgsdata$chr[i],"q")
	}
	else if (wgsdata$start[i]<=centromere_location)
	{
	  anarm = paste0("chr",wgsdata$chr[i],"p")
	}
	else
	{
	  anarm = paste0("chr",wgsdata$chr[i],"q")
	}
	tempdf = data.frame(sample=wgsdata$sample[i], chr = wgsdata$chr[i], arm = anarm, start = wgsdata$start[i], end = wgsdata$end[i], infercnv_mean_score = median(infercnvarr), wgs_score = median(wgsarr))
	print_infercnvarr_wgsarr_across_segments_median_df = rbind(print_infercnvarr_wgsarr_across_segments_median_df,tempdf)
      }
    }

    gene_to_segment_df_across_samples[match(gene_to_segment_df$genename,gene_to_segment_df_across_samples$genename),names(gene_to_segment_df_across_samples)==infercnvpats[z]] = gene_to_segment_df$segment

    print(cor.test(infercnvarr_across_segments,wgsarr_across_segments,method="spearman")$estimate)
    print(cor.test(infercnvarr_across_segments_median,wgsarr_across_segments_median,method="spearman")$estimate)
    rho = cor.test(infercnvarr_across_segments_median,wgsarr_across_segments_median,method="spearman")$estimate
    tempdf = data.frame(x=infercnvarr_across_segments_median,y=wgsarr_across_segments_median)
    pdf(paste0("/mnt/vdb/home/ubuntu2/correlate_wgs_and_infercnv/",infercnvpat,"_wgs_vs_infercnv.pdf"),width=9,height=7)
    print(ggplot(tempdf, aes(x, y)) + geom_point(size=3) + xlab("Mean InferCNV Score") + ylab("WGS Median Log Relative Copy Number") + ggtitle(paste0(infercnvpats[z]," Spearman's R^2 = ",rho^2)) + theme_set(theme_bw()) + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 20)))
    dev.off()

    infercnvarr_across_segments_median_across_samples = c(infercnvarr_across_segments_median_across_samples,infercnvarr_across_segments_median)
    wgsarr_across_segments_median_across_samples = c(wgsarr_across_segments_median_across_samples,wgsarr_across_segments_median)
  }
}

armdf = data.frame(sample=character(), shortchr=character(), chr=character(), start=integer(), end=integer(), infercnv_mean_score=double(), wgs_score=double())
uniquesamples = unique(print_infercnvarr_wgsarr_across_segments_median_df$sample)
for (i in 1:length(uniquesamples))
{
  for (j in 1:length(chrlengths))
  {
    tempdf = data.frame(sample=uniquesamples[i],shortchr=substring(names(chrlengths)[j],4),chr=paste0(names(chrlengths)[j],"p"),start=1,end=centromere_locations[j])
    armdf = rbind(armdf, tempdf)
    tempdf = data.frame(sample=uniquesamples[i],shortchr=substring(names(chrlengths)[j],4),chr=paste0(names(chrlengths)[j],"q"),start=centromere_locations[j]+1,end=chrlengths[j])
    armdf = rbind(armdf, tempdf)
  }
}

for (i in 1:length(armdf$sample))
{
  printdf_sub = subset(print_infercnvarr_wgsarr_across_segments_median_df, sample==armdf$sample[i] & chr==armdf$shortchr[i])
  overlapfracs = c()
  infercnv_mean_scores = c()
  wgs_scores = c()
  for (j in 1:length(printdf_sub$sample))
  {
    overlapstart = max(armdf$start[i], printdf_sub$start[j])
    overlapend = min(armdf$end[i], printdf_sub$end[j])
    overlap = overlapend-overlapstart
    if (overlap>0)
    {
      overlapfracs = c(overlapfracs, overlap/(armdf$end[i]-armdf$start[i]))
      infercnv_mean_scores = c(infercnv_mean_scores, printdf_sub$infercnv_mean_score[j])
      wgs_scores = c(wgs_scores, printdf_sub$wgs_score[j])
    }
  }
  if (length(overlapfracs)>0)
  {
    armdf$infercnv_mean_score[i] = sum(overlapfracs*infercnv_mean_scores)
    armdf$wgs_score[i] = sum(overlapfracs*wgs_scores)
  }
  else
  {
    armdf$infercnv_mean_score[i] = 1
    armdf$wgs_score[i] = 0
  }
  # if (armdf$sample[i]=="uM-07-GK-ON_S3" && armdf$chr[i]=="chr13p")
  # {
  #   nonsense = nonsense+1
  # }
}

#remove arms that had no overlap with wgs data
armdf = subset(armdf, infercnv_mean_score!=1 & wgs_score!=0)

print(cor.test(infercnvarr_across_segments_median_across_samples,wgsarr_across_segments_median_across_samples,method="spearman")$estimate)
nonzerowgsarr = (wgsarr_across_segments_median_across_samples>.1 | wgsarr_across_segments_median_across_samples < (-.1))

rho = cor.test(infercnvarr_across_segments_median_across_samples,wgsarr_across_segments_median_across_samples,method="spearman")$estimate
rhononzero = cor.test(infercnvarr_across_segments_median_across_samples[nonzerowgsarr],wgsarr_across_segments_median_across_samples[nonzerowgsarr],method="spearman")$estimate

permute_sample_rho_arr = c()
uniquesamples = unique(armdf$sample)
print("ATHERE")
for (i in 1:10000)
{
  selectsamples = sample(uniquesamples, length(uniquesamples)-2)
  selectarr = (nonzerowgsarr & (print_infercnvarr_wgsarr_across_segments_median_df$sample %in% selectsamples))
  rhononzero_permute = cor.test(print_infercnvarr_wgsarr_across_segments_median_df$infercnv_mean_score[selectarr],print_infercnvarr_wgsarr_across_segments_median_df$wgs_score[selectarr],method="spearman")$estimate
  permute_sample_rho_arr = c(permute_sample_rho_arr, rhononzero_permute^2)
}

permute_chr_rho_arr = c()
uniquechr = unique(print_infercnvarr_wgsarr_across_segments_median_df$chr)
print("ATHERE")
for (i in 1:length(uniquechr))
{
  selectarr = (nonzerowgsarr & (print_infercnvarr_wgsarr_across_segments_median_df$chr!=uniquechr[i]))
  rhononzero_permute = cor.test(print_infercnvarr_wgsarr_across_segments_median_df$infercnv_mean_score[selectarr],print_infercnvarr_wgsarr_across_segments_median_df$wgs_score[selectarr],method="spearman")$estimate
  permute_chr_rho_arr = c(permute_chr_rho_arr, rhononzero_permute^2)
}

permute_arm_rho_arr = c()
uniquearm = unique(print_infercnvarr_wgsarr_across_segments_median_df$arm)
print("ATHERE")
for (i in 1:length(uniquearm))
{
  selectarr = (nonzerowgsarr)
  selectarr[grep(uniquearm[i],print_infercnvarr_wgsarr_across_segments_median_df$arm)] = FALSE
  rhononzero_permute = cor.test(print_infercnvarr_wgsarr_across_segments_median_df$infercnv_mean_score[selectarr],print_infercnvarr_wgsarr_across_segments_median_df$wgs_score[selectarr],method="spearman")$estimate
  permute_arm_rho_arr = c(permute_arm_rho_arr, rhononzero_permute^2)
}

tempdf = data.frame(x=infercnvarr_across_segments_median_across_samples,y=wgsarr_across_segments_median_across_samples,class="nearzero")
tempdf$color[nonzerowgsarr] = "not_nearzero"

pdf(paste0("/mnt/vdb/home/ubuntu2/correlate_wgs_and_infercnv/all_samples_wgs_vs_infercnv.pdf"))
print(ggplot(tempdf, aes(x=x, y=y, color=color)) + geom_point() + scale_colour_manual(values = c("black","grey"), breaks = c("not_nearzero", "nearzero")) + xlab("Mean InferCNV Score") + ylab("WGS Median Log Relative Copy Number (MLRCN)") + geom_hline(yintercept = .1, color = "red") + geom_hline(yintercept = -.1, color = "red") + ggtitle(paste0("Spearman's R^2 = ",rho^2,"\nSpearman's R^2 Only With WGS MLRCN > .1 or < -.1 = ",rhononzero^2)) + theme_set(theme_bw()))
print(hist(permute_sample_rho_arr, main="correlation after dropping 2 samples in each permutation", xlab="Spearman's R^2"))
print(hist(permute_chr_rho_arr, main="correlation after dropping each possible chromosome", xlab="Spearman's R^2"))
print(hist(permute_arm_rho_arr, main="correlation after dropping each possible arm", xlab="Spearman's R^2"))
dev.off()

print(cor.test(armdf$infercnv_mean_score,armdf$wgs_score,method="spearman")$estimate)
nonzerowgsarr = (armdf$wgs_score>.1 | armdf$wgs_score < (-.1))

rho = cor.test(armdf$infercnv_mean_score,armdf$wgs_score,method="spearman")$estimate
rhononzero = cor.test(armdf$infercnv_mean_score[nonzerowgsarr],armdf$wgs_score[nonzerowgsarr],method="spearman")$estimate

permute_sample_rho_arr = c()
uniquesamples = unique(armdf$sample)
print("ATHERE")
for (i in 1:10000)
{
  selectsamples = sample(uniquesamples, length(uniquesamples)-2)
  selectarr = (nonzerowgsarr & (armdf$sample %in% selectsamples))
  rhononzero_permute = cor.test(armdf$infercnv_mean_score[selectarr],armdf$wgs_score[selectarr],method="spearman")$estimate
  permute_sample_rho_arr = c(permute_sample_rho_arr, rhononzero_permute^2)
}

permute_chr_rho_arr = c()
uniquechr = unique(armdf$shortchr)
print("ATHERE")
for (i in 1:length(uniquechr))
{
  selectarr = (nonzerowgsarr & (armdf$shortchr!=uniquechr[i]))
  rhononzero_permute = cor.test(armdf$infercnv_mean_score[selectarr],armdf$wgs_score[selectarr],method="spearman")$estimate
  permute_chr_rho_arr = c(permute_chr_rho_arr, rhononzero_permute^2)
}

permute_arm_rho_arr = c()
uniquearm = unique(armdf$chr)
print("ATHERE")
for (i in 1:length(uniquearm))
{
  selectarr = (nonzerowgsarr & (armdf$chr!=uniquearm[i]))
  rhononzero_permute = cor.test(armdf$infercnv_mean_score[selectarr],armdf$wgs_score[selectarr],method="spearman")$estimate
  permute_arm_rho_arr = c(permute_arm_rho_arr, rhononzero_permute^2)
}

tempdf = data.frame(x=armdf$infercnv_mean_score,y=armdf$wgs_score,class="nearzero")
tempdf$color[nonzerowgsarr] = "not_nearzero"

pdf(paste0("/mnt/vdb/home/ubuntu2/correlate_wgs_and_infercnv/all_samples_wgs_vs_infercnv_armlevel.pdf"))
print(ggplot(tempdf, aes(x=x, y=y, color=color)) + geom_point() + scale_colour_manual(values = c("black","grey"), breaks = c("not_nearzero", "nearzero")) + xlab("Mean InferCNV Score") + ylab("WGS Median Log Relative Copy Number (MLRCN)") + geom_hline(yintercept = .1, color = "red") + geom_hline(yintercept = -.1, color = "red") + ggtitle(paste0("Spearman's R^2 = ",rho^2,"\nSpearman's R^2 Only With WGS MLRCN > .1 or < -.1 = ",rhononzero^2)) + theme_set(theme_bw()))
print(hist(permute_sample_rho_arr, main="correlation after dropping 2 samples in each permutation", xlab="Spearman's R^2"))
print(hist(permute_chr_rho_arr, main="correlation after dropping each possible chromosome", xlab="Spearman's R^2"))
print(hist(permute_arm_rho_arr, main="correlation after dropping each possible arm", xlab="Spearman's R^2"))
dev.off()

write.table(print_infercnvarr_wgsarr_across_segments_median_df,"/mnt/vdb/home/ubuntu2/correlate_wgs_and_infercnv/infercnv_region_estimates.txt",sep="\t",col.names=T,row.names=F,quote=F)

gene_to_segment_df_across_samples2 = data.frame(genename=gene_to_segment_df_across_samples$genename, segment="NA")
for (i in 2:dim(gene_to_segment_df_across_samples)[2])
{
  gene_to_segment_df_across_samples2$segment[gene_to_segment_df_across_samples[,i]!="NA"] = gene_to_segment_df_across_samples[gene_to_segment_df_across_samples[,i]!="NA",i]
}
write.table(gene_to_segment_df_across_samples,"/mnt/vdb/home/ubuntu2/correlate_wgs_and_infercnv/infercnv_genes_by_wgs_region.txt",sep="\t",col.names=T,row.names=F,quote=F)
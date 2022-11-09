library(rlist)
library(ggplot2)
library(stringr)

### title: Correlate average arm-level infercnv copy numter alterations in fresh vs. Frozen samples in cutaneous melanoma, uveal primary, and NSCLC datasets
### author: Yiping Wang date: 11/08/2022

#define lists of sample names, folder locations, cancer cell type annotations, and fresh/frozen status
patslist = list(c("CD45negGEXBI5_S1_L001_final_thresh","bi005-skcm-5snseq_final_thresh","skcm-bi005-5pv2-snseq_final_thresh"),c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-E12_final_thresh","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-F1_final_thresh","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-GEX-F12_final_thresh"),c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX_final_thresh","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX_final_thresh","NSCL-NR001-5pv2-snseq_final_thresh"),c("Sarcoma167GEX_final_thresh","Sarcoma322GEX_final_thresh","Sarcoma559GEX_final_thresh","Sarcoma708GEX_final_thresh"))

shortpatslist = list(c("CD45neg","5snseq","5pv2-snseq"),c("SCRNA-5P-NA-E12","SCRNA-5P-NA-F1","SNRNA-5P-WI-F12"),c("SCRNA_5P_NA","SCRNA_5P_NI","SCRNA_5P_WI","5pv2-snseq"),c("167","322","559","708"))

s3folderlist = list(c("fresh-vs-frozen-comparison-ohio/BI5/scrna-seq","fresh-vs-frozen-comparison-ohio/BI5/snrna-seq","fresh-vs-frozen-comparison-ohio/BI5/snrna-seq"),c("fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata","fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata","fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata"),c("fresh-vs-frozen-comparison-ohio/nsclc","fresh-vs-frozen-comparison-ohio/nsclc","fresh-vs-frozen-comparison-ohio/nsclc","fresh-vs-frozen-comparison-ohio/nsclc"),c("fresh-vs-frozen-comparison-ohio/sarcoma-sn","fresh-vs-frozen-comparison-ohio/sarcoma-sn","fresh-vs-frozen-comparison-ohio/sarcoma-sn","fresh-vs-frozen-comparison-ohio/sarcoma-sn"))

cancertissuelist = list(list(c("Melanocytes"),c("Melanocytes"),c("Melanocytes")),
list(c("Melanocytes"),c("Melanocytes"),c("Melanocytes")),
list(c("Epithelial cells","Neurons"),c("Epithelial cells","Neurons"),c("Epithelial cells","Neurons"),c("Epithelial cells","Neurons")),
list(c("Chondrocytes","Fibroblasts","Smooth muscle"),c("Chondrocytes","Fibroblasts","Smooth muscle"),c("Chondrocytes","Fibroblasts","Smooth muscle"),c("Chondrocytes","Fibroblasts","Smooth muscle")))

fresh_frozen_list = list(c("fresh","frozen","frozen"),c("fresh","fresh","frozen"),c("fresh","frozen","frozen","frozen"),c("frozen","frozen","frozen","frozen"))

title_names = c("BI5 melanoma brain met","cpoi uveal primary","nsclc")

#define lengths of chromosomes and locations of centromeres
chrlengths = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)
centromere_locations = c(123400000,93900000,90900000,50000000,48800000,59800000,60100000,45200000,43000000,39800000,53400000,35500000,17700000,17200000,19000000,36800000,25100000,18500000,26200000,28100000,12000000,15000000)
names(chrlengths) = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
names(centromere_locations) = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")

pdf("fresh_vs_frozen_correlate_infercnv.pdf")
doubledf = data.frame(x=double(), y=double(), comparison=character())
titleline = ""
#values in topermutearr control whether to run just a single comparison of fresh vs frozen copy number for a dataset,
#or run 1000 permutation, with each permutation taking a random 50% subset of the data
topermutearr = c("permutecells","nopermute")
infercnv_arr = list()
corres_arr = list()
for (permuteopt in topermutearr)
{
  if (permuteopt=="nopermute")
  {
    numiters = 1
  }
  else
  {
    numiters = 1000
  }
  for (aniter in 1:numiters)
  {
    for (largeindex in 1:length(patslist))
    {
      pats = patslist[[largeindex]]
      shortpats = shortpatslist[[largeindex]]
      s3folders = s3folderlist[[largeindex]]
      cancertissues = cancertissuelist[[largeindex]]
      fresh_frozen_arr = fresh_frozen_list[[largeindex]]
      fresh_armlevel_cnv_arr = list()
      frozen_armlevel_cnv_arr = list()
      for (i in 1:length(pats))
      {
	pat = pats[i]
	s3folder = s3folders[i]

	#if (!file.exists(paste0("/data/",pat,"_infercnv.rds")))
	#{
	#  system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
	#}
	#load in infercnv object for a particular sample, if not already present in infercnv_arr
	if (sum(pat==names(infercnv_arr))==0)
	{
	  system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/inferCNV_cellbender_subcluster_",pat,"/run.final.infercnv_obj /data/",pat,"_infercnv.rds"))
	  infercnv_obj = readRDS(paste0("/data/",pat,"_infercnv.rds"))
	  infercnv_arr[[pat]] = infercnv_obj
	}
	infercnv_obj = infercnv_arr[[pat]]
	
	canceridxs = c()
	cancertissuearr = cancertissues[[i]]

	#subset infercnv data to cells that are only in cancer tissues
        for (j in 1:length(cancertissuearr))
	{
	  canceridxs = c(canceridxs, infercnv_obj@observation_grouped_cell_indices[[cancertissuearr[j]]])
	}

        #if performing multiple permutations, select half of all cancer cells in each permutation
	if (permuteopt=="permutecells")
	{
	  canceridxs = sample(canceridxs,floor(length(canceridxs)/2))
	}

        #for each chromosome arm, calculate average infercnv score across all genes in that arm
	average_cnv = rowMeans(infercnv_obj@expr.data[,canceridxs])

	average_cnv_temp = c()
	for (achr in sort(unique(infercnv_obj@gene_order$chr)))
	{
	  startsteps = c(1,centromere_locations[achr])
	  stopsteps = c(centromere_locations[achr]+1, chrlengths[achr])

	  print(achr)
	  for (z in 1:length(startsteps))
	  {
	    included_genes = rownames(infercnv_obj@gene_order)[infercnv_obj@gene_order$chr==achr & infercnv_obj@gene_order$start>startsteps[z] & infercnv_obj@gene_order$stop<stopsteps[z]]
	    if (length(included_genes)>0)
	    {
	      average_cnv_temp = c(average_cnv_temp, mean(average_cnv[included_genes]))
	    }
	    else
	    {
	      average_cnv_temp = c(average_cnv_temp, 1)
	    }
	    if (length(names(average_cnv_temp))==0)
	    {
	      names(average_cnv_temp) = c(paste0(achr," ",startsteps[z]," ",stopsteps[z]))
	    }
	    else
	    {
	      names(average_cnv_temp)[length(names(average_cnv_temp))] = paste0(achr," ",startsteps[z]," ",stopsteps[z])
	    }
	  }
	}
	average_cnv = average_cnv_temp

	if (fresh_frozen_arr[i]=="fresh")
	{
	  fresh_armlevel_cnv_arr = list.append(fresh_armlevel_cnv_arr, average_cnv)
	}
	if (fresh_frozen_arr[i]=="frozen")
	{
	  frozen_armlevel_cnv_arr = list.append(frozen_armlevel_cnv_arr, average_cnv)
	}
	system(paste0("rm /data/",pat,"_infercnv.rds"))
      }

      if (length(fresh_armlevel_cnv_arr)>0 && length(frozen_armlevel_cnv_arr)>0)
      {
        #determine chromosome arms that have average cnv levels calculated for all fresh samples in dataset
	common_arms = names(fresh_armlevel_cnv_arr[[1]])
	for (i in 1:length(fresh_armlevel_cnv_arr))
	{
	  common_arms = intersect(common_arms, names(fresh_armlevel_cnv_arr[[i]]))
	}

        #determine average cnv levels across all arms that have cnv levels in all fresh samples
	fresh_average_cnv = fresh_armlevel_cnv_arr[[1]][common_arms]
	if (length(fresh_armlevel_cnv_arr)>1)
	{
	  for (i in 2:length(fresh_armlevel_cnv_arr))
	  {
	    fresh_average_cnv = fresh_average_cnv + fresh_armlevel_cnv_arr[[i]][common_arms]
	  }
	}
	fresh_average_cnv = fresh_average_cnv/length(fresh_armlevel_cnv_arr)
	
        #determine chromosome arms that have average cnv levels calculated for all frozen samples in dataset
	common_arms = names(frozen_armlevel_cnv_arr[[1]])
	for (i in 1:length(frozen_armlevel_cnv_arr))
	{
	  common_arms = intersect(common_arms, names(frozen_armlevel_cnv_arr[[i]]))
	}

        #determine average cnv levels across all arms that have cnv levels in all frozen samples
	frozen_average_cnv = frozen_armlevel_cnv_arr[[1]][common_arms]
	if (length(frozen_armlevel_cnv_arr)>1)
	{
	  for (i in 2:length(frozen_armlevel_cnv_arr))
	  {
	    frozen_average_cnv = frozen_average_cnv + frozen_armlevel_cnv_arr[[i]][common_arms]
	  }
	}
	frozen_average_cnv = frozen_average_cnv/length(frozen_armlevel_cnv_arr)

        #determine whether we are comparing one fresh sample vs. multiple frozen samples, or vice versa
	if (length(fresh_armlevel_cnv_arr)==1)
	{
	  singlet_gene_cnv = fresh_armlevel_cnv_arr[[1]]
	  multiple_gene_cnv = frozen_armlevel_cnv_arr
	  singlet_pats = shortpats[fresh_frozen_arr=="fresh"]
	  multiple_pats = shortpats[fresh_frozen_arr=="frozen"]
	  singlet_type = "fresh"
	  multiple_type = "frozen"
	}
	else
	{
	  singlet_gene_cnv = frozen_armlevel_cnv_arr[[1]]
	  multiple_gene_cnv = fresh_armlevel_cnv_arr
	  singlet_pats = shortpats[fresh_frozen_arr=="frozen"]
	  multiple_pats = shortpats[fresh_frozen_arr=="fresh"]
	  singlet_type = "frozen"
	  multiple_type = "fresh"
	}

        #for each sample in multiple_type samples, correlate their average cnv vs. singlet_type sample 
	for (z in 1:length(multiple_gene_cnv))
	{
	  common_arms = intersect(names(singlet_gene_cnv), names(multiple_gene_cnv[[z]]))
	  corres = cor.test(singlet_gene_cnv[common_arms], multiple_gene_cnv[[z]][common_arms], method="spearman")

          #if running permutations, save all correlation coefficients in array
          if (permuteopt=="permutecells")
	  {
	    corres_label = paste0(singlet_pats[1]," vs. ",multiple_pats[z])
	    if (sum(names(corres_arr)==corres_label)==0)
	    {
	      corres_arr[[corres_label]] = c()
	    }
	    corres_arr[[corres_label]] = c(corres_arr[[corres_label]], (corres$estimate)^2)
	    #print(corres_arr)
	  }

          #if not running permutations, this implies that 1000 permutations have already been run
	  #for each sample in multiple_type samples, plot scatter plot of arm-level cnv vs. single_type arm-level cnv
          if (permuteopt=="nopermute")
	  {
	    corres_label = paste0(singlet_pats[1]," vs. ",multiple_pats[z])

            tempdf = data.frame(x=singlet_gene_cnv[common_arms], y=multiple_gene_cnv[[z]][common_arms])
	    tempdf = tempdf[abs(tempdf$x-1)>.01 & abs(tempdf$y-1)>.01,]
	    theme_set(theme_bw())
	    if (singlet_type=="fresh")
	    {
	      print(ggplot(tempdf, aes(x, y)) + geom_point(alpha = 1) + xlab(paste0(singlet_pats[1]," (Fresh) copy number")) + ylab(paste0(multiple_pats[z]," (Frozen) copy number")) + xlim(c(0.8,1.2)) + ylim(c(0.8,1.2)) + ggtitle(paste0(title_names[largeindex],"\nSpearman's R^2: ",as.character(round(corres$estimate^2,2)))))
	    }
	    else
	    {
	      print(ggplot(tempdf, aes(x, y)) + geom_point(alpha = 1) + xlab(paste0(multiple_pats[z]," (Fresh) copy number")) + ylab(paste0(singlet_pats[1]," (Frozen) copy number")) + xlim(c(0.8,1.2)) + ylim(c(0.8,1.2)) + ggtitle(paste0(title_names[largeindex],"\n","\nSpearman's R^2: ",as.character(round(corres$estimate^2,2)))))
	    }
	    #print histogram of correlation coefficients from permutation analysis
	    print(hist(corres_arr[[corres_label]], main=paste0(corres_label," random selection of half of cells"), xlab="Spearman's R^2"))

            #if comparing cutaneous melanoma fresh vs. frozen, save results in dataframe for later plotting
	    if (singlet_pats[1]=="CD45neg" && (multiple_pats[z]=="5snseq" || multiple_pats[z]=="5pv2-snseq"))
	    {
	      tempdf$comparison = paste0(singlet_pats[1]," vs. ",multiple_pats[z])
	      doubledf = rbind(doubledf, tempdf)
	      if (titleline=="")
	      {
		titleline = paste0(titleline,multiple_pats[z]," vs. ",singlet_pats[1]," Spearman's R^2: ",as.character(round(corres$estimate^2,2)),"\n")
	      }
	      else
	      {
		titleline = paste0(titleline,multiple_pats[z]," vs. ",singlet_pats[1]," Spearman's R^2: ",as.character(round(corres$estimate^2,2)))
	      }
	    }
	  }
	}
      }
    }
  }
}
dev.off()

#plot results of cutaneous melanoma fresh vs. frozen comparison
pdf("BI5_correlate_infercnv.pdf",height=7,width=10)
doubledf_sub1 = subset(doubledf, comparison=="CD45neg vs. 5pv2-snseq")
lm1 = lm(y ~ x, data = doubledf_sub1)
doubledf_sub2 = subset(doubledf, comparison=="CD45neg vs. 5snseq")
lm2 = lm(y ~ x, data = doubledf_sub2)
lmdf1 = data.frame(pred = predict(lm1, doubledf_sub1), x=doubledf_sub1$x, comparison=doubledf_sub1$comparison)
lmdf2 = data.frame(pred = predict(lm2, doubledf_sub2), x=doubledf_sub2$x, comparison=doubledf_sub2$comparison)
print(ggplot(doubledf, aes(x=x, y=y, color=comparison)) + geom_point(alpha = 1) + geom_line(data=lmdf1, aes(y=pred, x=x)) + geom_line(data=lmdf2, aes(y=pred, x=x)) + xlab(paste0("CD45neg (Fresh) copy number")) + ylab(paste0("5snseq or 5pv2-snseq (Frozen) copy number")) + ggtitle(titleline) + xlim(c(0.8,1.2)) + ylim(c(0.8,1.2)))
dev.off()
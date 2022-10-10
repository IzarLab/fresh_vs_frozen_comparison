library(rlist)
library(scRepertoire)
library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(cowplot)
library(scater)
library(stringr)
library(limma)
library(reshape2)
library(circlize)

useRibas = FALSE
if (useRibas) {
  patslist = list(c("ribas1_pre_tcr_S35_L004","ribas1_on_tcr_S36_L004","ribas_310_on_later_previd_3_TCR"),c("ribas_319_pre_previd_1_TCR","ribas_319_on_previd_2_TCR"))
  #shortpatslist = list(c("ribas_310_pre","ribas_310_on","ribas_319_on_later"),c("ribas_319_pre","ribas_319_on"))
  shortpatslist = list(c("pre","on","on_later"),c("pre","on"))
  colorslist = list(c("black","black","black"),c("black","black"))
} else {
  patslist = list(c("TCRBI5_S1_L001","bi005-skcm-5snseq-TCR"),
    c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F2","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F3", "UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-TCR-F9"),
    c("NSCL_NR001_SCRNA_5P_NA_BRAIN_TCR","NSCL_NR001_SNSEQ_5P_NI_BRAIN_TCR","NSCL_NR001_SNSEQ_5P_WI_BRAIN_TCR"))
  #shortpatslist = list(c("CD45pos","5snseq"),c("SCRNA-5P-NA-F2","SCRNA-5P-NA-F3","SNRNA-5P-WI-F9"),c("SCRNA-5P-NA","SNRNA-5P-NI","SNRNA-5P-WI"))
  shortpatslist = list(c("Mel_sc_5_CD45+","Mel_sn_5"),c("UM_sc_5","UM_sc_5_CD45+","UM_sn_5_inhib"),c("NSCLC_sc_5","NSCLC_sn_5","NSCLC_sn_5_inhib"))
  colorslist = list(c("blue","red"),c("blue","blue","red"),c("blue","red","red"))
}

# patslist = list(c("PA005","PA019","PA025","PA034","PA042","PA043","PA054","PA056","PA060","PA067","PA070","PA076","PA080","PA104","PA125","PA141"))
# shortpatslist = patslist

fullcolors = TRUE
if (useRibas) {
  pdf("ribas_melanoma_TCR_circos.pdf")
} else {
  if (fullcolors) {
    pdf("fresh_vs_frozen_comparison_TCR_circos_full.pdf")
  } else {
    pdf("fresh_vs_frozen_comparison_TCR_circos.pdf")
  }
}

for (z in 1:length(patslist))
{
  print(z)
  pats = patslist[[z]]
  shortpats = shortpatslist[[z]]
  colors = colorslist[[z]]
  contig_list<-c()
  for(i in 1:length(pats)){
    contig_list[[i]] <- read.csv(paste0(pats[i],'/filtered_contig_annotations.csv'), stringsAsFactors = F)
  }

  combined <- combineTCR(contig_list, samples = pats, ID = pats, cells ="T-AB",filterMulti = F)
  names(combined)<- pats

  removepats = c()
  for (i in 1:length(combined))
  {
    uniqueCTgene = unique(combined[[i]]$CTgene)
    uniqueCTgene = uniqueCTgene[!is.na(uniqueCTgene)]
    if (length(uniqueCTgene)==0)
    {
      removepats = c(removepats,pats[i])
    }
  }
  if (length(removepats)>0)
  {
    for (i in 1:length(removepats))
    {
      combined[[removepats[i]]] = NULL
      corresponding_shortpat = shortpats[pats==removepats[i]]
      pats = setdiff(pats, removepats[i])
      shortpats = setdiff(shortpats, corresponding_shortpat)
    }
  }
  #nonsense = nonsense+1

  contig_list<-c()
  for(i in 1:length(pats)){
    contig_list[[i]] <- read.csv(paste0(pats[i],'/filtered_contig_annotations.csv'), stringsAsFactors = F)
  }

  combined <- combineTCR(contig_list, samples = pats, ID = pats, cells ="T-AB",filterMulti = F)
  names(combined)<- pats

  tcrdf = data.frame(sectors = rep(as.character(1:length(pats)),each=2), x = rep(c(0,3),times=length(pats)), y = rep(c(0,1),times=length(pats)))

  circos.par("track.height" = 0.1)
  circos.initialize(tcrdf$sectors, x = tcrdf$x)

  rlelist = list()
  frequentrlelist = list()
  rlelist_noncum = list()
  for (i in 1:length(combined))
  {
    rleresult = rle(sort(combined[[i]]$CTgene))$lengths
    names(rleresult) = rle(sort(combined[[i]]$CTgene))$values

    #nonsense = nonsense+1

    rleresult_noncum = rleresult#sort(rleresult, decreasing=T)[1:min(10000,length(rleresult))]
    rlelist_noncum = list.append(rlelist_noncum, rleresult_noncum)

    rleresult = rleresult*3/sum(rleresult)
    frequentrlelist = list.append(frequentrlelist, sort(rleresult, decreasing=T)[1:min(10000,length(rleresult))])
    rleresult = cumsum(rleresult)
    rlelist = list.append(rlelist,rleresult)
    print(length(rleresult))
  }
  names(rlelist) = pats
  names(frequentrlelist) = pats
  names(rlelist_noncum) = pats

  for (i in 1:length(rlelist_noncum))
  {
    for (j in 1:length(rlelist_noncum))
    {
      if (i!=j)
      {
	mutual_clonotypes = intersect(names(rlelist_noncum[[i]]), names(rlelist_noncum[[j]]))
	all_clonotypes = union(names(rlelist_noncum[[i]]), names(rlelist_noncum[[j]]))
	commonrle_ith = rep(0, length(all_clonotypes))
	names(commonrle_ith) = all_clonotypes
	commonrle_ith[names(rlelist_noncum[[i]])] = rlelist_noncum[[i]]
	commonrle_jth = rep(0, length(all_clonotypes))
	names(commonrle_jth) = all_clonotypes
	commonrle_jth[names(rlelist_noncum[[j]])] = rlelist_noncum[[j]]
	#print("multinomial")
	#print(length(all_clonotypes))
	if (length(all_clonotypes)!=0)
	{
	  #print(dmultinom(rlelist_noncum[[j]][mutual_clonotypes], prob=rlelist_noncum[[i]][mutual_clonotypes]))
	  #print(dmultinom(commonrle_jth[1:10], prob=commonrle_ith[1:10]))
	  #nonsense = nonsense+1
	}
	# else
	# {
	#   print("no shared clonotypes")
	# }
	print("hypergeom")
	print(pats[i])
	print(pats[j])
	print(length(mutual_clonotypes))
	print(length(names(rlelist_noncum[[i]])))
	print(length(all_clonotypes)-length(names(rlelist_noncum[[i]])))
	print(length(names(rlelist_noncum[[j]])))
	print(dhyper(x=length(mutual_clonotypes), m=length(names(rlelist_noncum[[i]])), n=length(all_clonotypes)-length(names(rlelist_noncum[[i]])), k=length(names(rlelist_noncum[[j]]))))
      }
    }
  }

  for (i in 1:length(rlelist_noncum))
  {
    writedf = data.frame(clonotype=names(rlelist_noncum[[i]]), frequency=rlelist_noncum[[i]])
    write.table(writedf,"rlelist_noncum.txt",sep="\t",quote=F,row.names=F,col.names=T)
    system("python3 /mnt/vdb/home/ubuntu2/uveal-melanoma-code/giniscript.py")
    gini_val = read.table("gini_val.txt",sep="\t",header=T,quote=NULL)$gini
    print(pats[i])
    print(gini_val)
  }

  common_clonotypes_list = list()
  common_clonotypes_names_list = c()
  for (i in 1:length(rlelist))
  {
    for (j in 1:length(rlelist))
    {
      if (i!=j)
      {
	intersectarr = c()
	for (z1 in 1:length(names(rlelist[[i]])))
	{
	  for (z2 in 1:length(names(rlelist[[j]])))
	  {
	    if (names(rlelist[[i]])[z1]==names(rlelist[[j]])[z2] && (sum(names(rlelist[[i]])[z1]==names(frequentrlelist[[i]]))!=0 || sum(names(rlelist[[j]])[z2]==names(frequentrlelist[[j]]))!=0))
	    {
	      intersectarr = c(intersectarr, names(rlelist[[i]])[z1])
	    }
	  }
	}
	common_clonotypes_list = list.append(common_clonotypes_list, intersectarr)
	common_clonotypes_names_list = c(common_clonotypes_names_list, paste0(i," ",j))
      }
    }
  }
  names(common_clonotypes_list) = common_clonotypes_names_list
    
    common_clonotypes_all = unique(c(unlist(common_clonotypes_list)))
    if (fullcolors)
    {
      common_colors_all = rand_color(length(common_clonotypes_all))
    }
    else
    {
      common_colors_all = rand_color(length(common_clonotypes_all), transparency = 0.4)
    }
    common_colors_list = list()
    for (i in 1:length(common_clonotypes_list))
    {
      common_colors_list = list.append(common_colors_list, common_colors_all[match(common_clonotypes_list[[i]],common_clonotypes_all)])
    }
  #}


  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
      xlim = CELL_META$xlim
      ylim = CELL_META$ylim

      circos.text(CELL_META$xcenter,
        CELL_META$cell.ylim[2] + mm_y(5),
	shortpats[CELL_META$sector.numeric.index],
	facing = "bending",
	niceFacing = T,
	col = colors[CELL_META$sector.numeric.index])
      rleresult = rlelist[[CELL_META$sector.numeric.index]]
      colorvec = rep("gray",length(rleresult))
      if (FALSE)#(length(pats)==2)
      {
        colorvec[match(common_clonotypes,names(rleresult))] = common_colors
      }
      else
      {
	matchidxs = match(common_clonotypes_all,names(rleresult))
        colorvec[matchidxs[!is.na(matchidxs)]] = common_colors_all[which(!is.na(matchidxs))]
      }
      circos.rect(c(0,rleresult[1:length(rleresult)-1]), rep(0, length(rleresult)), rleresult, rep(1, length(rleresult)), col = colorvec, border = "black")
  })

  for (i1 in 1:length(common_clonotypes_list))
  {
    if (length(common_clonotypes_list[[i1]]>0))
    {
      for (i in 1:length(common_clonotypes_list[[i1]]))
      {
	rleidxsstring = strsplit(names(common_clonotypes_list)[i1]," ")[[1]]
	rleidxs = as.integer(rleidxsstring)
	index1 = which(names(rlelist[[rleidxs[1]]])==common_clonotypes_list[[i1]][i])
	index2 = which(names(rlelist[[rleidxs[2]]])==common_clonotypes_list[[i1]][i])

        xcoord1higher = rlelist[[rleidxs[1]]][index1]
        if (index1!=1)
	{
	  xcoord1lower = rlelist[[rleidxs[1]]][index1-1]
	}
	else
	{
	  xcoord1lower = 0
	}
        xcoord2higher = rlelist[[rleidxs[2]]][index2]
        if (index2!=1)
	{
	  xcoord2lower = rlelist[[rleidxs[2]]][index2-1]
	}
	else
	{
	  xcoord2lower = 0
	}

        truesectororder = as.integer(sort(as.character(1:length(pats))))

	circos.link(truesectororder[rleidxs[1]], c(xcoord1lower, xcoord1higher), truesectororder[rleidxs[2]], c(xcoord2lower, xcoord2higher), col = common_colors_list[[i1]][i], border = "black")
      }
    }
  }
}
dev.off()

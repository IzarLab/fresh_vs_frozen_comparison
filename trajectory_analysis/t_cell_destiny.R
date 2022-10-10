library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(grid)
library(stringr)
library(rlist)
library(scRepertoire)

colDC <- c('#DE8C00', '#F564E3', '#7CAE00', '#00B4F0', '#00C08B')

prefix_arr = c("cd8")#"treg_and_tfh","treg","tfh","cd4")#,"cd8")
useRibas = FALSE
for (i in 1:length(prefix_arr))
{
  #specify figure layout. if plotting ribas non-cd8 t-cells, plot multiple diffusion maps on grid of subplots. otherwise, plot each diffusion map on a separate individual plot.
  if (useRibas)
  {
    if (prefix_arr[i]=="cd8")
    {
      pdf(paste0("ribas_melanoma_",prefix_arr[i],"_destiny_trajectory.pdf"),height=8,width=8)
    }
    else
    {
      pdf(paste0("ribas_melanoma_",prefix_arr[i],"_destiny_trajectory.pdf"),height=60,width=40)
    }
  }
  else
  {
    pdf(paste0("uveal_melanoma_",prefix_arr[i],"_destiny_trajectory.pdf"),height=8,width=8)
  }
  vplayout = function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

  figurecount = 1

  #load in data, extract t-cell subgroups to be used for diffusion map
  if (useRibas)
  {
    seu = readRDS("/data/reannotate_ribas_melanoma_merged_tcells.rds")
  }
  else
  {
    seu = readRDS("/data/reannotate_uveal_melanoma_tcells_reintegrated_with_1_subclustered_dim_num_25_then_15.rds")
  }
  if (useRibas)
  {
    if (prefix_arr[i]=="cd4")
    {
      seu = subset(seu, seurat_clusters %in% c(12))
    }
    if (prefix_arr[i]=="treg")
    {
      seu = subset(seu, seurat_clusters %in% c(7))
    }
    if (prefix_arr[i]=="tfh")
    {
      seu = subset(seu, seurat_clusters %in% c(5))
    }
    if (prefix_arr[i]=="treg_and_tfh")
    {
      seu = subset(seu, seurat_clusters %in% c(5,7))
    }
    if (prefix_arr[i]=="cd8")
    {
      seu = subset(seu, seurat_clusters %in% c(0,1,3,4,6,7,9,10,11,12,14))
    }
  }
  else
  {
    if (prefix_arr[i]=="cd4")
    {
      seu = subset(seu, seurat_clusters_with_1_subclustered %in% c("1_CD4","2","3","10"))
    }
    if (prefix_arr[i]=="cd8")
    {
      seu = subset(seu, seurat_clusters_with_1_subclustered %in% c("0","1_CD8","4","6","8"))
    } 
  }
  seu$treatment_group = unlist(lapply(strsplit(seu$orig.ident,"_"), function(x) {x[4]}))
  seu$celltype_bped_main[!(seu$celltype_bped_main %in% c("CD4+ T-cells","CD8+ T-cells","Melanocytes","NK cells"))] = "Other cell types"

  #load in t-cell state and differentiation signatures, and calculate enrichment scores
  azizi_signatures = read.table("/mnt/vdb/home/ubuntu2/uveal-melanoma-code/azizi_signatures.csv",header=T,sep=",",quote=NULL)
  azizi_signatures_tcell_exhaustion = read.table("/mnt/vdb/home/ubuntu2/uveal-melanoma-code/azizi_signatures_tcell_exhaustion.txt",header=T,sep="\t",quote=NULL)
  while (dim(azizi_signatures_tcell_exhaustion)[1]<dim(azizi_signatures)[1])
  {
    tempdf = data.frame(Terminal.exhaustion="",Precursor.exhaustion="")
    azizi_signatures_tcell_exhaustion = rbind(azizi_signatures_tcell_exhaustion, tempdf)
  }

  azizi_signatures$G1.S = NULL
  azizi_signatures$G2.M = NULL
  azizi_signatures$M1.Macrophage.Polarization = NULL
  azizi_signatures$M2.Macrophage.Polarization = NULL
  azizi_signatures$Terminal.exhaustion = azizi_signatures_tcell_exhaustion$Terminal.exhaustion
  azizi_signatures$Precursor.exhaustion = azizi_signatures_tcell_exhaustion$Precursor.exhaustion
  if (prefix_arr[i]!="cd8")
  {
    azizi_signatures$CD8.T.Cell.Activation = NULL
  }
  azizi_signatures[["just_TOX"]] = c("TOX",rep("",dim(azizi_signatures)[1]-1))
  azizi_signatures[["just_TCF7"]] = c("TCF7",rep("",dim(azizi_signatures)[1]-1))

  for (j in 1:length(names(azizi_signatures)))
  {
    seu = AddModuleScore(seu, features = list(na.omit(azizi_signatures[[names(azizi_signatures)[j]]])), name = names(azizi_signatures)[j], assay = "RNA", search = T)
  }

  #load in t-cell differentiation signatures and calculate enrichment scores. for negative marker genes, first multiply expression of gene by -1, calculate enrichment score, then multiply expression again by -1.
  diff_sigs = read.table("/mnt/vdb/home/ubuntu2/uveal-melanoma-code/azizi_differentiation_state_markers.txt",sep="\t",header=T,quote=NULL)
  for (j in 1:length(names(diff_sigs)))
  {
    sig_genenames = diff_sigs[,j]
    sig_genenames = sig_genenames[sig_genenames!=""]
    rev_genenames = c()
    for (k in 1:length(sig_genenames))
    {
      strlen = str_length(sig_genenames[k])

      threelastchars = substr(sig_genenames[k],strlen-2,strlen)
      if (threelastchars=="+/-")
      {
  	sig_genenames[k] = substr(sig_genenames[k],1,strlen-4)
      }
      else
      {
  	lastchar = substr(sig_genenames[k],strlen,strlen)
  	if (lastchar=="-")
  	{
  	  rev_genenames = c(rev_genenames, substr(sig_genenames[k],1,strlen-2))
  	}
  	sig_genenames[k] = substr(sig_genenames[k],1,strlen-2)
      }
    }

    for (k in 1:length(rev_genenames))
    {
      seu@assays$RNA@counts[rev_genenames[k],] = -seu@assays$RNA@counts[rev_genenames[k],]
      seu@assays$RNA@data[rev_genenames[k],] = -seu@assays$RNA@data[rev_genenames[k],]
    }

    seu = AddModuleScore(seu, features = list(na.omit(sig_genenames)), name = names(diff_sigs)[j], assay = "RNA", search = T)

    for (k in 1:length(rev_genenames))
    {
      seu@assays$RNA@counts[rev_genenames[k],] = -seu@assays$RNA@counts[rev_genenames[k],]
      seu@assays$RNA@data[rev_genenames[k],] = -seu@assays$RNA@data[rev_genenames[k],]
    }

    azizi_signatures[[names(diff_sigs[j])]] = c(sig_genenames,rep("",dim(azizi_signatures)[1]-length(sig_genenames)))
  }

  #if plotting cd8 t-cells, only project a subset of signatures onto diffusion map
  cd8_only = c("just_TOX","just_TCF7","TCell.Terminal.Differentiation","T.N","T.SCM","T.CM","T.EM","T.TE","Precursor.exhaustion")
  if (prefix_arr[i]=="cd8")
  {
    for (sig in names(azizi_signatures))
    {
      if (!(sig %in% cd8_only))
      {
        azizi_signatures[[sig]] = NULL
      }
    }
  }

  #score cells by cell cycle signatures, project scatter plot of cycling scores
  seu = CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
  tempdf = data.frame(S=seu$S.Score, G2M=seu$G2M.Score)
  if (prefix_arr[i]=="cd8")
  {
    print(ggplot(tempdf) + geom_point(aes(x=S, y=G2M)) + ggtitle(paste0(prefix_arr[i],"_cell_cycle")))
  }
  else
  {
    print(ggplot(tempdf) + geom_point(aes(x=S, y=G2M)) + ggtitle(paste0(prefix_arr[i],"_cell_cycle")), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
  }

  seu$original_seurat_clusters = seu$seurat_clusters
  seu = ScaleData(object = seu)
  seu = RunPCA(object = seu)
  seu = FindNeighbors(seu, dims = 1:15)
  seu = FindClusters(seu)
  seu = RunUMAP(object = seu, dims = 1:20)

  seu$nCount_RNA = 1.06*seu$nCount_RNA
  seu$nFeature_RNA = 1.06*seu$nFeature_RNA

  #for cd8 t-cells, assign clonality expansion information
  if (prefix_arr[i]=="cd8")
  {
    source("assign_TCR_clonality.R")
    seu = assign_TCR_clonality(seu, dataset_names[i])

    highlight_list = list()
    clone_numbers = sort(unique(seu$clonality))
    if (length(clone_numbers)!=0)
    {
      seu$dummy_index = 1:length(seu$orig.ident)
      for (i1 in 1:length(clone_numbers)) {
  	highlight_list = list.append(highlight_list,colnames(seu)[seu$clonality==clone_numbers[i1]])
      }
      display_numbers = as.character(clone_numbers)
      for (i1 in 1:length(display_numbers))
      {
  	if (str_length(display_numbers[i1])==1)
  	{
  	  display_numbers[i1] = paste0("0",display_numbers[i1])
  	}
      }
      names(highlight_list) = paste0("Clonality: ", display_numbers)
    }

    unique_clonality_groups = unique(seu$clonality_group)
    highlight_list = list()
    for (agroup in unique_clonality_groups)
    {
      highlight_list = list.append(highlight_list, colnames(seu)[seu$clonality_group==agroup])
    }
    names(highlight_list) = unique_clonality_groups

    seu$orig.ident[seu$orig.ident=="ribas1_on_tcr_S36_L004"] = "ribas_310_on"
    seu$orig.ident[seu$orig.ident=="ribas_310_on_later_previd_3_TCR"] = "ribas_310_on_later"
    seu$orig.ident[seu$orig.ident=="ribas1_pre_tcr_S35_L004"] = "ribas_310_pre"
  }

  #select 2000 most variable features, use them to run diffusion map
  seu = FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  if (useRibas && prefix_arr[i]=="cd8")
  {
    datamat = seu@assays$RNA@data[seu@assays$RNA@var.features,]
  }
  else
  {
    datamat = seu@assays$integrated@data
  }

  es <- as.ExpressionSet(as.data.frame(t(as.matrix(datamat))))
  es@phenoData@data <- seu@meta.data
  dm <- DiffusionMap(es, verbose = T, n_pcs = 30)

  #read in frequencies of clonotypes for ribas samples
  if (useRibas)
  {
    ribas_tcr_pats = c("ribas1_pre_tcr_S35_L004","ribas1_on_tcr_S36_L004","ribas_310_on_later_previd_3_TCR")
    ribas_short_pats = c("ribas_310_pre","ribas_310_on","ribas_310_on_later")
    contig_list = c()
    for(ribas_tcr_idx in 1:length(ribas_tcr_pats))
    {
      contig_list[[ribas_tcr_idx]] = read.csv(paste0(ribas_tcr_pats[ribas_tcr_idx],'/filtered_contig_annotations.csv'), stringsAsFactors = F)
    }

    combined = combineTCR(contig_list, samples = ribas_tcr_pats, ID = ribas_tcr_pats, cells = "T-AB", filterMulti = F)
    names(combined) = ribas_short_pats
  }
  else
  {
    for(unique_ident_idx in 1:length(unique_idents))
    {
      contig_list[[unique_ident_idx]] = read.csv(paste0(unique_idents[unique_ident_idx],'/filtered_contig_annotations.csv'), stringsAsFactors = F)
    }

    combined = combineTCR(contig_list, samples = unique_idents, ID = unique_idents, cells = "T-AB", filterMulti = F)
    names(combined) = unique_idents
  }

  #determine proportions of clonotypes in each sample
  proptable1 = table(combined[[1]]$CTgene)/sum(table(combined[[1]]$CTgene))
  proptable2 = table(combined[[2]]$CTgene)/sum(table(combined[[2]]$CTgene))
  proptable3 = table(combined[[3]]$CTgene)/sum(table(combined[[3]]$CTgene))

  #determine clonotypes that are shared between pre and on, and on and on_later samples
  #calculate differences in clonotype frequency
  shared_1_and_2 = intersect(unique(combined[[1]]$CTgene),unique(combined[[2]]$CTgene))
  shared_1_and_2_diffs = (proptable2[shared_1_and_2] - proptable1[shared_1_and_2])
  shared_1_and_2_diffs = shared_1_and_2_diffs[names(sort(abs(shared_1_and_2_diffs),decreasing=T))]

  shared_2_and_3 = intersect(unique(combined[[2]]$CTgene),unique(combined[[3]]$CTgene))
  shared_2_and_3_diffs = (proptable3[shared_2_and_3] - proptable2[shared_2_and_3])
  shared_2_and_3_diffs = shared_2_and_3_diffs[names(sort(abs(shared_2_and_3_diffs),decreasing=T))]

  #determine top 20 most common clonotypes in each sample, and combine them
  #top_diff_names = union(names(shared_2_and_3_diffs[1:10]),names(shared_1_and_2_diffs[1:10]))
  top_diff_names = union(names(sort(table(combined[[1]]$CTgene),decreasing=T)),names(sort(table(combined[[2]]$CTgene),decreasing=T))[1:20])
  top_diff_names = union(top_diff_names,names(sort(table(combined[[3]]$CTgene),decreasing=T))[1:20])

  #plot diffusion maps with QC metrics projected on top
  qc_fields = c("nCount_RNA","nFeature_RNA","percent.mt","ScrubDoublet_score")
  qc_titles = paste0(qc_fields," before quality filtering")
  if (prefix_arr[i]=="cd8")
  {
    for(qc_idx in 1:length(qc_fields))
    {
      par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
      palette(viridis(100))
      maintitle = qc_titles[qc_idx]
      maintitle = paste0(prefix_arr[i],"_",maintitle)
      print(plot.DiffusionMap(dm, c(1,2,3), col = es@phenoData@data[[qc_fields[qc_idx]]], pch = 20, main = maintitle))
    }
  }
  else
  {
    for(qc_idx in 1:length(qc_fields))
    {
      par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
      palette(viridis(100))
      maintitle = qc_titles[qc_idx]
      maintitle = paste0(prefix_arr[i],"_",maintitle)
      print(plot.DiffusionMap(dm, c(1,2), col = es@phenoData@data[[qc_fields[qc_idx]]], pch = 20) + ggtitle(maintitle), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
      figurecount = figurecount+1
    }
  }

  #plot values of diffusion components 1-10 on gene expression umaps of t-cells
  for(dcidx in 1:10)
  {
    eval(parse(text=paste0("seu$DC",dcidx," = dm$DC",dcidx)))
  }

  if (prefix_arr[i]=="cd8")
  {
    for (dcidx in 1:10)
    {
      print(FeaturePlot(seu, features = c(paste0("DC",dcidx)), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC",dcidx)))
    }
  }
  else
  {
    for (dcidx in 1:10)
    {
      print(FeaturePlot(seu, features = c(paste0("DC",dcidx)), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC",dcidx)), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
      figurecount = figurecount+1
    }
  }

  #project strength of t-cell state and differentiation signatures onto diffusion maps
  #print 3 signatures per row for non-cd8 t-cells, 5 signatures for cd8
  rownum = floor((figurecount-1)/9)+1
  colnum = ((figurecount-1) %% 9)+1
  rownum = rownum+1
  colnum = 1
  for (j in 1:length(names(azizi_signatures)))
  {
    if (prefix_arr[i]=="cd8")
    {
      if (j!=1 && (j %% 5)==1)
      {
  	rownum = rownum+1
  	colnum = 1
      }
    }
    else
    {
      if (j!=1 && (j %% 3)==1)
      {
  	rownum = rownum+1
  	colnum = 1
      }
    }
    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(100))

    #calculate correlation values of signature with DCs 1-4
    testres1 = cor.test(seu$DC1, seu[[paste0(names(azizi_signatures)[j],"1")]][[1]],method="spearman")
    test_pval1 = round(testres1$p.val,8)
    testres2 = cor.test(seu$DC2, seu[[paste0(names(azizi_signatures)[j],"1")]][[1]],method="spearman")
    test_pval2 = round(testres2$p.val,8)
    testres3 = cor.test(seu$DC3, seu[[paste0(names(azizi_signatures)[j],"1")]][[1]],method="spearman")
    test_pval3 = round(testres3$p.val,8)
    testres4 = cor.test(seu$DC4, seu[[paste0(names(azizi_signatures)[j],"1")]][[1]],method="spearman")
    test_pval4 = round(testres4$p.val,8)

    #print correlation values for DCs 1 and 2, and DCs 3 and 4, in title
    #or alternatively, do not print correlation values for non-cd8 t-cells
    maintitle = paste0(names(azizi_signatures)[j],"1")
    maintitle = paste0(prefix_arr[i],"_",maintitle)
    maintitle_first = paste0(maintitle,"\nSpearman p-val with DC1: ",test_pval1)
    maintitle_first = paste0(maintitle_first,"\nSpearman p-val with DC2: ",test_pval2)
    maintitle_second = paste0(maintitle,"\nSpearman p-val with DC3: ",test_pval3)
    maintitle_second = paste0(maintitle_second,"\nSpearman p-val with DC4: ",test_pval4)

    maintitle_third = names(azizi_signatures)[j]
    if (maintitle_third=="just_TOX")
    {
      maintitle_third = "TOX"
    }
    if (maintitle_third=="just_TCF7")
    {
      maintitle_third = "TCF7"
    }

    if (prefix_arr[i]=="cd8")
    {
      print(plot.DiffusionMap(dm, c(1,2,3), col = es@phenoData@data[[paste0(names(azizi_signatures)[j],"1")]], pch = 20, main = maintitle_third))
      print(colorlegend(es@phenoData@data[[paste0(names(azizi_signatures)[j],"1")]], viridis(length(es@phenoData@data[[paste0(names(azizi_signatures)[j],"1")]])), posx = c(0.88, 0.9), posy = c(0, 0.65)))
      print(plot.DiffusionMap(dm, c(1,2), col = es@phenoData@data[[paste0(names(azizi_signatures)[j],"1")]], pch = 20, main = maintitle_third) + ggtitle(maintitle_third))
      print(plot.DiffusionMap(dm, c(3,4), col = es@phenoData@data[[paste0(names(azizi_signatures)[j],"1")]], pch = 20, main = maintitle_third) + ggtitle(maintitle_third))
      print(FeaturePlot(seu, features = c(paste0(names(azizi_signatures)[j],"1")), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(names(azizi_signatures)[j],"1")))
    }
    else
    {
      print(plot.DiffusionMap(dm, c(1,2), col = es@phenoData@data[[paste0(names(azizi_signatures)[j],"1")]], pch = 20, main = maintitle_first), vp=vplayout(rownum, colnum))
      colnum = colnum+1
      print(plot.DiffusionMap(dm, c(3,4), col = es@phenoData@data[[paste0(names(azizi_signatures)[j],"1")]], pch = 20, main = maintitle_second), vp=vplayout(rownum, colnum))
      colnum = colnum+1
      print(FeaturePlot(seu, features = c(paste0(names(azizi_signatures)[j],"1")), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(names(azizi_signatures)[j],"1")), vp=vplayout(rownum, colnum))
      colnum = colnum+1
    }
  }

  #calculate correlation of all genes expression with first four DCs
  #write out results to csv file, excluding ribosomal and mitochondrial genes
  corr_df = data.frame(gene=character(), dc = integer(), rho = double(), pval = double())
  for (z1 in 1:4)
  {
    corr_df_temp = data.frame(gene=character(), dc = integer(), rho = double(), pval = double())
    eval(parse(text=paste0("dc_arr = seu$DC",z1)))
    for (i1 in 1:length(rownames(datamat)))
    {
      corr_res = cor.test(datamat[i1,],dc_arr,method="spearman")
      tempdf = data.frame(gene=rownames(datamat)[i1],dc=z1,rho=corr_res$estimate,pval=corr_res$p.val)
      corr_df_temp = rbind(corr_df_temp, tempdf)
    }
    corr_df_temp = corr_df_temp[order(corr_df_temp$rho),]
    corr_df = rbind(corr_df, corr_df_temp)
  }
  corr_df$placeholder = (1:length(corr_df$gene))
  corr_df = subset(corr_df, !(placeholder %in% grep("RPS|RPL|MT-",corr_df$gene)))
  write.table(corr_df,paste0(prefix_arr[i],"_dc_corr_genes.csv"),sep=",",quote=FALSE,row.names=F,col.names=T)

  clonality_palette = c("gray","blue",colorRampPalette(c("yellow","red"))(4))
  names(clonality_palette) = names(highlight_list)

  if (prefix_arr[i]=="cd8")
  {
    #plot clonal expansion on diffusion map if plotting cd8
    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(length(unique(seu$clonality_group))))
    maintitle = "Clonal expansion"
    color_factor = as.factor(es@phenoData@data$clonality_group)
    print(plot.DiffusionMap(dm, c(1,2,3), col = color_factor, pch = 20, main = maintitle, pal = clonality_palette[levels(color_factor)]))

    #for ribas samples, plot 20 largest clonotypes from each sample on diffusion map
    #color samples by timepoint: pre, on, and on_later
    #write out data to csv file, for use by rgl_script.R later
    if (useRibas)
    {
      writeout_df = data.frame(DC1=seu$DC1, DC2=seu$DC2, DC3=seu$DC3)
      for (top_diff_name in top_diff_names)
      {
	es@phenoData@data$clonality_over_time = ""
	for (ribas_short_pat in ribas_short_pats)
	{
	  combined_barebarcodes = unlist(lapply(str_split(combined[[ribas_short_pat]]$barcode,"_"), function (x) {x[length(x)]}))	
	  associated_barcodes = combined_barebarcodes[combined[[ribas_short_pat]]$CTgene==top_diff_name]
	  es@phenoData@data$clonality_over_time[es@phenoData@data$orig.ident==ribas_short_pat & es@phenoData@data$barebarcodes %in% associated_barcodes] = ribas_short_pat
	}

	par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
	clonality_over_time_palette = c("gray","blue","green","red")
	names(clonality_over_time_palette) = c("","ribas_310_pre","ribas_310_on","ribas_310_on_later")
	maintitle = paste0(top_diff_name,'\nclonality_over_time')
	maintitle = paste0(prefix_arr[i],"\n",maintitle)
	color_factor = as.factor(es@phenoData@data$clonality_over_time)
	print(plot.DiffusionMap(dm, c(1,2,3), col = color_factor, pch = 20, main = maintitle, pal = clonality_over_time_palette[levels(color_factor)], draw.legend = T))
	eval(parse(text=paste0("writeout_df$`",top_diff_name,"` = es@phenoData@data$clonality_over_time")))
      }
      write.table(writeout_df, "ribas_melanoma_cd8_clonality_over_time.csv", sep=",", quote=FALSE, row.names=F, col.names=T)
    }
  }
  else
  {
    #plot diffusion maps with QC metrics projected on top
    misc_fields = c("treatment_group","celltype_bped_main","original_seurat_clusters")
    misc_titles = misc_fields
    if (length(unique(seu$original_seurat_clusters))>1)
    {
      misc_fields = c(misc_fields,"original_seurat_clusters")
      misc_titles = misc_fields
    }
    
    rownum = rownum+1
    for (miscidx in 1:length(misc_fields))
    {
      colnum = miscidx
      par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
      palette(colDC)
      maintitle = misc_titles[miscidx]
      maintitle = paste0(prefix_arr[i],"_",maintitle)
      print(plot.DiffusionMap(dm, c(1,2), col = as.factor(es@phenoData@data[[misc_fields[miscidx]]]), pch = 20, main = maintitle), vp=vplayout(rownum, colnum))
      figurecount = figurecount+1
    }
  }
  dev.off()
}
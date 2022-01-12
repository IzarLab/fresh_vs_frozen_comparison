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

# colBP <- c('#A80D11', '#008DB8')
# colSCSN <- c('#E1AC24', '#288F56')
colDC <- c('#DE8C00', '#F564E3', '#7CAE00', '#00B4F0', '#00C08B')

prefix_arr = c("cd8")#"treg_and_tfh","treg","tfh","cd4")#,"cd8")
useRibas = FALSE
for (i in 1:length(prefix_arr))
{
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
  if (prefix_arr[i]=="cd8")
  {
    #pushViewport(viewport(layout = grid.layout(11,10)))
  }
  else
  {
    #pushViewport(viewport(layout = grid.layout(13,9)))
  }
  vplayout = function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

  figurecount = 1

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
      #seu@assays$RNA@scale.data[rev_genenames[k],] = -seu@assays$RNA@scale.data[rev_genenames[k],]
    }

    seu = AddModuleScore(seu, features = list(na.omit(sig_genenames)), name = names(diff_sigs)[j], assay = "RNA", search = T)

    for (k in 1:length(rev_genenames))
    {
      seu@assays$RNA@counts[rev_genenames[k],] = -seu@assays$RNA@counts[rev_genenames[k],]
      seu@assays$RNA@data[rev_genenames[k],] = -seu@assays$RNA@data[rev_genenames[k],]
      #seu@assays$RNA@scale.data[rev_genenames[k],] = -seu@assays$RNA@scale.data[rev_genenames[k],]
    }

    azizi_signatures[[names(diff_sigs[j])]] = c(sig_genenames,rep("",dim(azizi_signatures)[1]-length(sig_genenames)))
  }

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

  # #seu = subset(seu, G2M.Score<.03 & S.Score<.3)
  # #seu = subset(seu, G2M.Score<.25 & S.Score<.5)

  seu$original_seurat_clusters = seu$seurat_clusters
  seu = ScaleData(object = seu)
  seu = RunPCA(object = seu)
  seu = FindNeighbors(seu, dims = 1:15)
  seu = FindClusters(seu)
  seu = RunUMAP(object = seu, dims = 1:20)

  seu$nCount_RNA = 1.06*seu$nCount_RNA
  seu$nFeature_RNA = 1.06*seu$nFeature_RNA

  if (prefix_arr[i]=="cd8")
  {
    seu$barebarcodes = unlist(lapply(strsplit(colnames(seu),"_"), function(x) {x[1]}))

    if (useRibas)
    {
      seu$orig.ident[seu$orig.ident=="ribas_310_on"] = "ribas1_on_tcr_S36_L004"
      seu$orig.ident[seu$orig.ident=="ribas_310_on_later"] = "ribas_310_on_later_previd_3_TCR"
      seu$orig.ident[seu$orig.ident=="ribas_310_pre"] = "ribas1_pre_tcr_S35_L004"
    }
    unique_idents = unique(seu$orig.ident)
    seu$clonality = 0
    for (i1 in 1:length(unique_idents))
    {
      csv_table = read.table(paste0(unique_idents[i1],"_TR_frequency.csv"),sep=",",header=T,quote=NULL)
      for (j in 1:length(csv_table$barcodes))
      {
  	barcodes = unique(str_split(csv_table$barcodes[j],"\\|")[[1]])
  	seu$clonality[(seu$barebarcodes %in% barcodes) & (seu$orig.ident==unique_idents[i1])] = length(barcodes)
      }
    }

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

    seu$clonality_group = "Unmatched with TCR sequencing"
    seu$clonality_group[seu$clonality==1] = "Unexpanded clones"
    seu$clonality_group[seu$clonality==2] = "Expanded clones with clonality 2"
    seu$clonality_group[seu$clonality>2 & seu$clonality<=5] = "Expanded clones with clonality > 2 and <= 5"
    seu$clonality_group[seu$clonality>5 & seu$clonality<=20] = "Expanded clones with clonality > 5 and <= 20"
    seu$clonality_group[seu$clonality>20] = "Expanded clones with clonality > 20"

    highlight_list = list(colnames(seu)[seu$clonality_group=="Unmatched with TCR sequencing"],colnames(seu)[seu$clonality_group=="Unexpanded clones"],colnames(seu)[seu$clonality_group=="Expanded clones with clonality 2"],colnames(seu)[seu$clonality_group=="Expanded clones with clonality > 2 and <= 5"],colnames(seu)[seu$clonality_group=="Expanded clones with clonality > 5 and <= 20"],colnames(seu)[seu$clonality_group=="Expanded clones with clonality > 20"])
    names(highlight_list) = c("Unmatched with TCR sequencing","Unexpanded clones","Expanded clones with clonality 2","Expanded clones with clonality > 2 and <= 5","Expanded clones with clonality > 5 and <= 20","Expanded clones with clonality > 20")

    seu$orig.ident[seu$orig.ident=="ribas1_on_tcr_S36_L004"] = "ribas_310_on"
    seu$orig.ident[seu$orig.ident=="ribas_310_on_later_previd_3_TCR"] = "ribas_310_on_later"
    seu$orig.ident[seu$orig.ident=="ribas1_pre_tcr_S35_L004"] = "ribas_310_pre"
  }

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

  proptable1 = table(combined[[1]]$CTgene)/sum(table(combined[[1]]$CTgene))
  proptable2 = table(combined[[2]]$CTgene)/sum(table(combined[[2]]$CTgene))
  proptable3 = table(combined[[3]]$CTgene)/sum(table(combined[[3]]$CTgene))

  shared_1_and_2 = intersect(unique(combined[[1]]$CTgene),unique(combined[[2]]$CTgene))
  shared_1_and_2_diffs = (proptable2[shared_1_and_2] - proptable1[shared_1_and_2])
  shared_1_and_2_diffs = shared_1_and_2_diffs[names(sort(abs(shared_1_and_2_diffs),decreasing=T))]

  shared_2_and_3 = intersect(unique(combined[[2]]$CTgene),unique(combined[[3]]$CTgene))
  shared_2_and_3_diffs = (proptable3[shared_2_and_3] - proptable2[shared_2_and_3])
  shared_2_and_3_diffs = shared_2_and_3_diffs[names(sort(abs(shared_2_and_3_diffs),decreasing=T))]

  #top_diff_names = union(names(shared_2_and_3_diffs[1:10]),names(shared_1_and_2_diffs[1:10]))
  top_diff_names = union(names(sort(table(combined[[1]]$CTgene),decreasing=T)),names(sort(table(combined[[2]]$CTgene),decreasing=T))[1:20])
  top_diff_names = union(top_diff_names,names(sort(table(combined[[3]]$CTgene),decreasing=T))[1:20])

  if (prefix_arr[i]=="cd8")
  {
    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(100))
    maintitle = 'nCount_RNA before quality filtering'
    maintitle = paste0(prefix_arr[i],"_",maintitle)
    print(plot.DiffusionMap(dm, c(1,2,3), col = es@phenoData@data$nCount_RNA, pch = 20, main = maintitle))

    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(100))
    maintitle = 'nFeature_RNA before quality filtering'
    maintitle = paste0(prefix_arr[i],"_",maintitle)
    print(plot.DiffusionMap(dm, c(1,2,3), col = es@phenoData@data$nFeature_RNA, pch = 20, main = maintitle))

    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(100))
    maintitle = 'percent.mt before quality filtering'
    maintitle = paste0(prefix_arr[i],"_",maintitle)
    print(plot.DiffusionMap(dm, c(1,2,3), col = es@phenoData@data$percent.mt, pch = 20, main = maintitle))

    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(100))
    maintitle = 'ScrubDoublet_score before quality filtering'
    maintitle = paste0(prefix_arr[i],"_",maintitle)
    print(plot.DiffusionMap(dm, c(1,2,3), col = es@phenoData@data$ScrubDoublet_score, pch = 20, main = maintitle))
  }
  else
  {
    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(100))
    maintitle = 'nCount_RNA before quality filtering'
    maintitle = paste0(prefix_arr[i],"_",maintitle)
    print(plot.DiffusionMap(dm, c(1,2), col = es@phenoData@data$nCount_RNA, pch = 20) + ggtitle(maintitle), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1

    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(100))
    maintitle = 'nFeature_RNA before quality filtering'
    maintitle = paste0(prefix_arr[i],"_",maintitle)
    print(plot.DiffusionMap(dm, c(1,2), col = es@phenoData@data$nFeature_RNA, pch = 20) + ggtitle(maintitle), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1

    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(100))
    maintitle = 'percent.mt before quality filtering'
    maintitle = paste0(prefix_arr[i],"_",maintitle)
    print(plot.DiffusionMap(dm, c(1,2), col = es@phenoData@data$percent.mt, pch = 20) + ggtitle(maintitle), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1

    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(100))
    maintitle = 'ScrubDoublet_score before quality filtering'
    maintitle = paste0(prefix_arr[i],"_",maintitle)
    print(plot.DiffusionMap(dm, c(1,2), col = es@phenoData@data$ScrubDoublet_score, pch = 20) + ggtitle(maintitle), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
  }

  nfeature_thresh = quantile(seu$nFeature_RNA,.95)
  ncount_thresh = quantile(seu$nCount_RNA,.95)
  percentmt_thresh = quantile(seu$percent.mt,.95)
  scrubdoublet_thresh = quantile(seu$ScrubDoublet_score,.95)
  # seu = subset(seu, nFeature_RNA<nfeature_thresh & nCount_RNA<ncount_thresh & percent.mt<percentmt_thresh & ScrubDoublet_score<scrubdoublet_thresh)

  # seu = ScaleData(object = seu)
  # seu = RunPCA(object = seu)
  # seu = FindNeighbors(seu, dims = 1:15)
  # seu = FindClusters(seu)
  # seu = RunUMAP(object = seu, dims = 1:20)

  # seu = FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  # if (prefix_arr[i]=="cd8")
  # {
  #   datamat = seu@assays$RNA@data[seu@assays$RNA@var.features,]
  # }
  # else
  # {
  #   datamat = seu@assays$integrated@data
  # }

  # es <- as.ExpressionSet(as.data.frame(t(as.matrix(datamat))))
  # es@phenoData@data <- seu@meta.data
  # dm <- DiffusionMap(es, verbose = T, n_pcs = 30)

  seu$DC1 = dm$DC1
  seu$DC2 = dm$DC2
  seu$DC3 = dm$DC3
  seu$DC4 = dm$DC4
  seu$DC5 = dm$DC5
  seu$DC6 = dm$DC6
  seu$DC7 = dm$DC7
  seu$DC8 = dm$DC8
  seu$DC9 = dm$DC9
  seu$DC10 = dm$DC10

  if (prefix_arr[i]=="cd8")
  {
    print(FeaturePlot(seu, features = c("DC1"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC1")))
    print(FeaturePlot(seu, features = c("DC2"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC2")))
    print(FeaturePlot(seu, features = c("DC3"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC3")))
    print(FeaturePlot(seu, features = c("DC4"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC4")))
    print(FeaturePlot(seu, features = c("DC5"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC5")))
    print(FeaturePlot(seu, features = c("DC6"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC6")))
    print(FeaturePlot(seu, features = c("DC7"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC7")))
    print(FeaturePlot(seu, features = c("DC8"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC8")))
    print(FeaturePlot(seu, features = c("DC9"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC9")))
    print(FeaturePlot(seu, features = c("DC10"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC10")))
  }
  else
  {
    print(FeaturePlot(seu, features = c("DC1"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC1")), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
    print(FeaturePlot(seu, features = c("DC2"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC2")), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
    print(FeaturePlot(seu, features = c("DC3"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC3")), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
    print(FeaturePlot(seu, features = c("DC4"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC4")), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
    print(FeaturePlot(seu, features = c("DC5"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC5")), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
    print(FeaturePlot(seu, features = c("DC6"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC6")), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
    print(FeaturePlot(seu, features = c("DC7"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC7")), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
    print(FeaturePlot(seu, features = c("DC8"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC8")), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
    print(FeaturePlot(seu, features = c("DC9"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC9")), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
    print(FeaturePlot(seu, features = c("DC10"), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(prefix_arr[i]," DC10")), vp=vplayout(floor((figurecount-1)/9)+1, ((figurecount-1) %% 9)+1))
    figurecount = figurecount+1
  }

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
    testres1 = cor.test(seu$DC1, seu[[paste0(names(azizi_signatures)[j],"1")]][[1]],method="spearman")
    test_pval1 = round(testres1$p.val,8)
    testres2 = cor.test(seu$DC2, seu[[paste0(names(azizi_signatures)[j],"1")]][[1]],method="spearman")
    test_pval2 = round(testres2$p.val,8)
    testres3 = cor.test(seu$DC3, seu[[paste0(names(azizi_signatures)[j],"1")]][[1]],method="spearman")
    test_pval3 = round(testres3$p.val,8)
    testres4 = cor.test(seu$DC4, seu[[paste0(names(azizi_signatures)[j],"1")]][[1]],method="spearman")
    test_pval4 = round(testres4$p.val,8)
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
    #maintitle_third = maintitle
    #maintitle_third = paste0(maintitle,"\nSpearman p-val with DC1: ",test_pval1)
    #maintitle_third = paste0(maintitle_third,"\nSpearman p-val with DC2: ",test_pval2)
    #maintitle_third = paste0(maintitle_third,"\nSpearman p-val with DC3: ",test_pval3)

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
    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(length(unique(seu$clonality_group))))
    #maintitle = 'clonality_group'
    #maintitle = paste0(prefix_arr[i],"_",maintitle)
    maintitle = "Clonal expansion"
    color_factor = as.factor(es@phenoData@data$clonality_group)
    print(plot.DiffusionMap(dm, c(1,2,3), col = color_factor, pch = 20, main = maintitle, pal = clonality_palette[levels(color_factor)]))

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

    # par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    # palette(colDC)
    # maintitle = 'treatment_group'
    # maintitle = paste0(prefix_arr[i],"_",maintitle)
    # print(plot.DiffusionMap(dm, c(1,2,3), col = as.factor(es@phenoData@data$treatment_group), pch = 20, main = maintitle))

    # par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    # palette(viridis(length(unique(seu$celltype_bped_main))))
    # maintitle = 'celltype_bped_main'
    # maintitle = paste0(prefix_arr[i],"_",maintitle)
    # print(plot.DiffusionMap(dm, c(1,2,3), col = as.factor(es@phenoData@data$celltype_bped_main), pch = 20, main = maintitle))

    # if (length(unique(seu$original_seurat_clusters))>1)
    # {
    #   par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    #   palette(viridis(length(unique(seu$original_seurat_clusters))))
    #   maintitle = 'original_seurat_cluster'
    #   maintitle = paste0(prefix_arr[i],"_",maintitle)
    #   print(plot.DiffusionMap(dm, c(1,2,3), col = as.factor(es@phenoData@data$original_seurat_clusters), pch = 20, main = maintitle))
    # }
  }
  else
  {
    rownum = rownum+1
    colnum = 1
    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(colDC)
    maintitle = 'treatment_group'
    maintitle = paste0(prefix_arr[i],"_",maintitle)
    print(plot.DiffusionMap(dm, c(1,2), col = as.factor(es@phenoData@data$treatment_group), pch = 20, main = maintitle), vp=vplayout(rownum, colnum))
    figurecount = figurecount+1

    colnum = 2
    par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
    palette(viridis(length(unique(seu$celltype_bped_main))))
    maintitle = 'celltype_bped_main'
    maintitle = paste0(prefix_arr[i],"_",maintitle)
    print(plot.DiffusionMap(dm, c(1,2), col = as.factor(es@phenoData@data$celltype_bped_main), pch = 20, main = maintitle), vp=vplayout(rownum, colnum))
    figurecount = figurecount+1

    if (length(unique(seu$original_seurat_clusters))>1)
    {
      colnum = 3
      par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
      palette(viridis(length(unique(seu$original_seurat_clusters))))
      maintitle = 'original_seurat_cluster'
      maintitle = paste0(prefix_arr[i],"_",maintitle)
      print(plot.DiffusionMap(dm, c(1,2), col = as.factor(es@phenoData@data$original_seurat_clusters), pch = 20, main = maintitle), vp=vplayout(rownum, colnum))
      figurecount = figurecount+1
    }
  }
  dev.off()
}
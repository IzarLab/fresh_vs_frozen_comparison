library(rlist)
library(ggplot2)
library(Seurat)
library(infercnv)

foldersList = c("",
  "",
  "",
  "")
integrated_name_arr = c("reannotate_ribas_melanoma_merged_tcells","reannotate_uveal_melanoma_tcells_reintegrated_with_1_subclustered_dim_num_25_then_15","fresh_vs_frozen_tcells_reannotate_BI5","fresh_vs_frozen_tcells_reannotate_cpoi-uvealprimarydata")
for (i in 3:4) {
  seu = readRDS(paste0("/data/",integrated_name_arr[i],".rds"))
  DefaultAssay(seu) = "RNA"

  if (i==1)
  {
    seu = subset(seu, seurat_clusters %in% c(0,1,2,3,4,6,7,9,10,11,12,14))
    seu$orig.ident[seu$orig.ident=="ribas_310_on"] = "ribas1_on_tcr_S36_L004"
    seu$orig.ident[seu$orig.ident=="ribas_310_on_later"] = "ribas_310_on_later_previd_3_TCR"
    seu$orig.ident[seu$orig.ident=="ribas_310_pre"] = "ribas1_pre_tcr_S35_L004"
  }
  else if (i==2)
  {
    seu = subset(seu, seurat_clusters %in% c("0","1_CD8","4","6","8"))
  }
  else if (i==3)
  {
    seu = subset(seu, orig.ident %in% c("CD45pos","5snseq"))
    seu$orig.ident[seu$orig.ident=="CD45pos"] = "TCRBI5_S1_L001"
    seu$orig.ident[seu$orig.ident=="5snseq"] = "bi005-skcm-5snseq-TCR"
  }
  else
  {
    seu$orig.ident[seu$orig.ident=="SCRNA-5P-NA-E12"] = "UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F2"
    seu$orig.ident[seu$orig.ident=="SCRNA-5P-NA-F1"] = "UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F3"
    seu$orig.ident[seu$orig.ident=="SNRNA-5P-WI-F12"] = "UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-TCR-F9"
  }

  if (i==1 || i==2)
  {
    dim_num = 15

    seu <- ScaleData(object = seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    seu <- RunPCA(object = seu)
    seu <- FindNeighbors(seu, dims = 1:dim_num)
    seu <- FindClusters(seu)
    seu <- RunUMAP(object = seu, dims = 1:dim_num)
  }

  seu$barebarcodes = unlist(lapply(strsplit(colnames(seu),"_"), function(x) {x[1]}))

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

  if (i==2 || i==3)
  {
    # if (i==2)
    # {
    #   seu_reannotate = readRDS("/data/fresh_vs_frozen_tcells_reannotate_BI5.rds")
    # }
    # if (i==3)
    # {
    #   seu_reannotate = readRDS("/data/fresh_vs_frozen_tcells_reannotate_cpoi-uvealprimarydata.rds")
    # }
    
    # seu_reannotate$barebarcodes = unlist(lapply(strsplit(colnames(seu_reannotate),"_"), function(x) {x[1]}))
    # seu_reannotate$manual_annotation_label_tcell
  }

  if (i==1 || i==2)
  {
    seu$clonality_group = "Unmatched with TCR sequencing"
  }
  if (i==3 || i==4)
  {
    seu$clonality_group = "Unmatched with TCR sequencing"
    seu$clonality_group[seu$clonality==0 & seu$manual_annotation_label_tcell=="CD4+ T-cells"] = "CD4+ T-cells unmatched with TCR sequencing"
    seu$clonality_group[seu$clonality==0 & seu$manual_annotation_label_tcell=="CD8+ T-cells"] = "CD8+ T-cells unmatched with TCR sequencing"
  }
  seu$clonality_group[seu$clonality==1] = "Unexpanded clones"
  seu$clonality_group[seu$clonality==2] = "Expanded clones with clonality 2"
  seu$clonality_group[seu$clonality>2 & seu$clonality<=5] = "Expanded clones with clonality > 2 and <= 5"
  seu$clonality_group[seu$clonality>5 & seu$clonality<=20] = "Expanded clones with clonality > 5 and <= 20"
  seu$clonality_group[seu$clonality>20] = "Expanded clones with clonality > 20"

  unique_clonality_groups = unique(seu$clonality_group)
  highlight_list = list()
  for (agroup in unique_clonality_groups)
  {
    highlight_list = list.append(highlight_list, colnames(seu)[seu$clonality_group==agroup])
  }
  names(highlight_list) = unique_clonality_groups

  if (i==1 || i==2)
  {
    clonality_palette = c("gray","blue",colorRampPalette(c("yellow","red"))(4))
    point_scale = c(0.2,1,1,1,1,1)
    names(clonality_palette) = c("Unmatched with TCR sequencing","Unexpanded clones","Expanded clones with clonality 2","Expanded clones with clonality > 2 and <= 5","Expanded clones with clonality > 5 and <= 20","Expanded clones with clonality > 20")
    names(point_scale) = names(clonality_palette)
  }
  else if (i==3 || i==4)
  {
    clonality_palette = c("gray","black","purple","blue",colorRampPalette(c("yellow","red"))(4))
    point_scale = c(0.5,0.5,0.5,1,1,1,1,1)
    names(clonality_palette) = c("Unmatched with TCR sequencing","CD4+ T-cells unmatched with TCR sequencing","CD8+ T-cells unmatched with TCR sequencing","Unexpanded clones","Expanded clones with clonality 2","Expanded clones with clonality > 2 and <= 5","Expanded clones with clonality > 5 and <= 20","Expanded clones with clonality > 20")
    names(point_scale) = names(clonality_palette)
  }

  pdf(paste0(integrated_name_arr[i],"_clonality_umap.pdf"),height=7,width=12)
  umap_df = data.frame(UMAP_1=seu@reductions$umap[[,1]], UMAP_2=seu@reductions$umap[[,2]], clonality_group=seu$clonality_group)
  theme_set(theme_bw())
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "orig.ident", repel = T, label.size = 3, shuffle = T) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)) + ggtitle("orig.ident"))
  print(ggplot(umap_df) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=clonality_group, size=clonality_group)) + scale_color_manual(values=clonality_palette, breaks=names(clonality_palette)) + scale_size_manual(values=point_scale, breaks=names(point_scale)) + xlab("UMAP_1") + ylab("UMAP_2") + guides(color = guide_legend(title="Clonality"), size = "none"))

  if (integrated_name_arr[i]=="reannotate_ribas_melanoma_merged_tcells")
  {
    print(FeaturePlot(seu, features = c("TOP2A","MKI67","TOX","TCF7"), min.cutoff = "1", max.cutoff = "4"))
  }
  else
  {
    print(FeaturePlot(seu, features = c("TOP2A","MKI67","TOX","TCF7"), min.cutoff = "1", max.cutoff = "3"))
  }
  dev.off()
}
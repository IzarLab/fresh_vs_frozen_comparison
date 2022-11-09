library(rlist)
library(ggplot2)
library(Seurat)
library(infercnv)
library(stringr)

### title: Create umaps of T cells with TCR clonality plotted in different colors according to frequency, for cutaneous melanoma, uveal melanoma primary, uveal melanoma metastasis, and sequential cutaneous melanoma samples under pembrolizumab therapy
### author: Yiping Wang date: 11/08/2022

#load list of rds objects corresponding to each dataset
foldersList = c("",
  "",
  "",
  "")
integrated_name_arr = c("reannotate_ribas_melanoma_merged_tcells","reannotate_uveal_melanoma_tcells_reintegrated_with_1_subclustered_dim_num_25_then_15","fresh_vs_frozen_tcells_reannotate_BI5","fresh_vs_frozen_tcells_reannotate_cpoi-uvealprimarydata")
dataset_names = c("ribas","UMEL","BI5","UM")
display_individual_clonalities = FALSE
for (i in 3:4) {
  seu = readRDS(paste0("/data/",integrated_name_arr[i],".rds"))
  DefaultAssay(seu) = "RNA"

  #extract specific seurat clusters or cell types corresponding to t-cells in each dataset
  if (dataset_names[i]=="ribas")
  {
    seu = subset(seu, seurat_clusters %in% c(0,1,2,3,4,6,7,9,10,11,12,14))
  }
  else if (dataset_names[i]=="UMEL")
  {
    seu = subset(seu, seurat_clusters %in% c("0","1_CD8","4","6","8"))
  }
  else if (dataset_names[i]=="BI5")
  {
    seu = subset(seu, orig.ident %in% c("CD45pos","5snseq"))
  }
  else if (dataset_names[i]=="UM")
  {
  }

  #if specific seurat clusters were selected, recluster and rerun umap
  if (dataset_names[i]=="ribas" || dataset_names[i]=="UMEL")
  {
    dim_num = 15

    seu <- ScaleData(object = seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    seu <- RunPCA(object = seu)
    seu <- FindNeighbors(seu, dims = 1:dim_num)
    seu <- FindClusters(seu)
    seu <- RunUMAP(object = seu, dims = 1:dim_num)
  }

  #assign TCR clonality groups to rds object using function
  source("fresh_vs_frozen_comparison/TCR_analyses/assign_TCR_clonality.R")
  seu = assign_TCR_clonality(seu, dataset_names[i])

  #for uveal melanoma metastasis data only, write out csv file containing clonality and clonality_group information
  if (dataset_names[i]=="UM")
  {
    writetable = data.frame(barcode = seu$barebarcodes, orig.ident = seu$orig.ident, clonality = seu$clonality, clonality_group = seu$clonality_group)
    write.table(writetable, "uveal_melanoma_tcell_clonality.csv", sep=",", quote=F, col.names=T, row.names=F)
  }

  #if display_individual_clonalities is TRUE, create list of cell barcodes corresponding to each clonality value
  #label list using human-friendly clonality labels
  if (display_individual_clonalities)
  {
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
  }
  else
  {
    #otherwise, create list of cell barcodes corresponding to each clonality_group value
    unique_clonality_groups = unique(seu$clonality_group)
    highlight_list = list()
    for (agroup in unique_clonality_groups)
    {
      highlight_list = list.append(highlight_list, colnames(seu)[seu$clonality_group==agroup])
    }
    names(highlight_list) = unique_clonality_groups
  }

  #use gray for cell barcodes not matched with TCR sequencing, yellow-red spectrum for expanded clonotypes, and either blue for all unexpanded clones, or additionally black and purple for CD4 and CD8 unexpanded clones
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

  #print umaps of tcell by fresh/frozen status, original sample identity, and clonality group, as well as expression of cell cycle and exhaustion markers
  pdf(paste0(integrated_name_arr[i],"_clonality_umap.pdf"),height=7,width=12)
  umap_df = data.frame(UMAP_1=seu@reductions$umap[[,1]], UMAP_2=seu@reductions$umap[[,2]], clonality_group=seu$clonality_group)
  theme_set(theme_bw())
  print(DimPlot(seu, reduction = "umap", label = T, group.by = "fresh_frozen", repel = T, label.size = 3, shuffle = T) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)) + ggtitle("fresh_frozen"))
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
  #nonsense = nonsense+1
}
library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(grid)

# colBP <- c('#A80D11', '#008DB8')
# colSCSN <- c('#E1AC24', '#288F56')
colDC <- c('#DE8C00', '#F564E3', '#7CAE00', '#00B4F0', '#00C08B')

# system("aws s3 cp s3://uveal-melanoma/Seurat/integrated/um_all_integrated.rds /data/um_all_integrated.rds")

#load in signatures of t-cell exhaustion
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

diff_sigs = read.table("/mnt/vdb/home/ubuntu2/uveal-melanoma-code/azizi_differentiation_state_markers.txt",sep="\t",header=T,quote=NULL)
for (j in 1:length(names(diff_sigs)))
{
  #determine genes that are negative markers in signatures, take negative of gene expression in data for these markers
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

  #calculate signature strengths for each signature, reverse negative gene expression
  seu = AddModuleScore(seu, features = list(na.omit(sig_genenames)), name = names(diff_sigs)[j], assay = "RNA", search = T)

  for (k in 1:length(rev_genenames))
  {
    seu@assays$RNA@counts[rev_genenames[k],] = -seu@assays$RNA@counts[rev_genenames[k],]
    seu@assays$RNA@data[rev_genenames[k],] = -seu@assays$RNA@data[rev_genenames[k],]
    #seu@assays$RNA@scale.data[rev_genenames[k],] = -seu@assays$RNA@scale.data[rev_genenames[k],]
  }

  azizi_signatures[[names(diff_sigs[j])]] = c(sig_genenames,rep("",dim(azizi_signatures)[1]-length(sig_genenames)))
}

seu = readRDS("/data/reannotate_uveal_melanoma_tcells_reintegrated_with_1_subclustered_dim_num_25_then_15.rds")
for (j in 1:length(names(azizi_signatures)))
{
  seu = AddModuleScore(seu, features = list(na.omit(azizi_signatures[[names(azizi_signatures)[j]]])), name = names(azizi_signatures)[j], assay = "RNA", search = T)
}

#print vln plots of signature strength
pdf("uveal_melanoma_tcells_reintegrated_vln_plots.pdf",width=12,height=12)
print(VlnPlot(seu, features = paste0(names(azizi_signatures),"1"), group.by = 'seurat_clusters_with_1_subclustered', pt.size = 0, assay = 'RNA', stack = T, flip = T) + NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)))
dev.off()

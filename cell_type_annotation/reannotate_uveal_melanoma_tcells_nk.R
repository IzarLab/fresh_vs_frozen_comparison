library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(grid)
library(stringr)

seu = readRDS(paste0("/data/reannotate_uveal_melanoma_tcells_reintegrated_with_1_subclustered_dim_num_25_then_15.rds"))

seu_nk = subset(seu, (seurat_clusters %in% c("5","9")))
seu_nk <- ScaleData(object = seu_nk)
seu_nk <- RunPCA(object = seu_nk)
seu_nk <- FindNeighbors(seu_nk, dims = 1:5)
seu_nk <- FindClusters(seu_nk, resolution = 0.1)
seu_nk <- RunUMAP(object = seu_nk, dims = 1:5)

pdf("reannotate_uveal_melanoma_tcells_nk.pdf")
print(DimPlot(seu_nk, reduction = "umap", label = T, group.by = "seurat_clusters", repel = T, label.size = 3) + ggtitle(paste0('nk seurat_clusters')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))
print(DimPlot(seu_nk, reduction = "umap", label = T, group.by = "orig.ident", repel = T, label.size = 3) + ggtitle(paste0('nk orig.ident')) + guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) + theme(legend.text=element_text(size=10)))

DefaultAssay(seu_nk) = "RNA"
canonical_markers = c("NCAM1","KLRF1","FCGR3A","B3GAT1")
for (marker in canonical_markers)
{
  print(FeaturePlot(seu_nk, features = c(marker), min.cutoff = "q9", max.cutoff = "q90") + ggtitle(paste0(marker)))
}
DefaultAssay(seu_nk) = "integrated"

allmarkers = FindAllMarkers(seu_nk, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
print(DoHeatmap(seu_nk, features = top10$gene) + NoLegend() + ggtitle(paste0("DEG Heatmap")))
write.table(allmarkers,paste0("reannotate_uveal_melanoma_tcells_nk_markers.csv"),sep=",",quote=F,row.names=F,col.names=T)
dev.off()

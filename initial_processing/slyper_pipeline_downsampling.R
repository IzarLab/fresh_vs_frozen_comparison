library(Seurat)
library(scuttle)
library(rlist)

foldersList = c("s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline","s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
integrated_name_arr = c("BI5","NR1")

for (i in 1:length(foldersList)) {
  system(paste0("aws s3 cp ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
  seu = readRDS(paste0("/data/",integrated_name_arr[i],"_integrated.rds"))

  unique_idents = unique(seu$orig.ident)
  counts_arr = list()
  for (an_ident in unique_idents)
  {
    subset_obj = subset(seu, orig.ident==an_ident)
    counts_arr = list.append(counts_arr, subset_obj@assays$RNA@counts)
  }
  downsampled_counts_arr = downsampleBatches(counts_arr)
  names(downsampled_counts_arr) = unique_idents
  object.list = c()
  for (an_ident in unique_idents)
  {
    aseu = CreateSeuratObject(downsampled_counts_arr[[an_ident]], project = "a", min.cells = 1, min.features = 1)
    print(an_ident)
    print(median(aseu$nFeature_RNA))
    aseu[["percent.mt"]] = PercentageFeatureSet(aseu, pattern = "^MT-")
    print(median(aseu$percent.mt))
    print(dim(aseu@assays$RNA@counts))
    object.list = c(object.list, aseu)
  }

  # find anchors
  anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:20)

  # integrate data sets
  seu <- IntegrateData(anchorset = anchors, dims = 1:20)

  # normal workflow
  seu <- ScaleData(object = seu)
  seu <- RunPCA(object = seu)
  seu <- FindNeighbors(seu, dims = 1:15)
  seu <- FindClusters(seu)
  seu <- RunUMAP(object = seu, dims = 1:20)

  DefaultAssay(seu) <- "RNA"

  outputname = integrated_name_arr[i]

  saveRDS(seu,paste0("/data/",outputname,"_downsampled_integrated.rds"))
  system(paste0("aws s3 cp /data/",outputname,"_downsampled_integrated.rds ",foldersList[i],"/Seurat/integrated/",outputname,"_downsampled_integrated.rds"))
  system(paste0("rm /data/",outputname,"_downsampled_integrated.rds"))
}

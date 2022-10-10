library(DropletUtils)
library(cellranger)
library(rlist)
library(Seurat)

patsList = list(c("CD45negGEXBI5_S1_L001","CD45posGEXBI5_S1_L001"),
  c("bi005-skcm-5snseq","bi005-skcm","skcm-bi005-5pv2-snseq"),
  c("UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-E12","UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-GEX-F1","UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-GEX-F12"),
  c("NSCL_NR001_SCRNA_5P_NA_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_3P_WI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_NI_BRAIN_GEX","NSCL_NR001_SNSEQ_5P_WI_BRAIN_GEX","NSCL-NR001-5pv2-snseq"),
  c("BI5CST","BI5TST","NR1CST","NR1TST"))

foldersList = c("s3://fresh-vs-frozen-comparison-ohio/BI5/scrna-seq",
  "s3://fresh-vs-frozen-comparison-ohio/BI5/snrna-seq",
  "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata",
  "s3://fresh-vs-frozen-comparison-ohio/nsclc",
  "s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")

origreadcountslist = list(c(64488,46254),c(33369,55664,14196),c(26345,45174,33034),c(68782,45384,28213,48042,47685,14365),c(284834,120278,265922,157799))
downsamplereadcountslist = list(c(64488,23127),c(10000,11132,14196),c(26345,9034,9910),c(6878,45384,28213,9608,9536,8517),c(10000,10000,10000,10000))

do_manual_downsampling = FALSE

saturation_arr = c()
downsample_saturation_arr = c()
median_counts_arr = c()
median_genes_arr = c()
filtered_median_counts_arr = c()
filtered_median_genes_arr = c()
double_filtered_median_counts_arr = c()
double_filtered_median_genes_arr = c()
double_filtered_genes_dist_arr = list()
double_filtered_seurat_arr = list()
processed_pats_arr = c()

for (z in 1:length(foldersList)) {
  pats = patsList[[z]]
  origreadcounts = origreadcountslist[[z]]
  downsamplereadcounts = downsamplereadcountslist[[z]]

  #should be same as version in fresh_vs_frozen_comparison
  Rcpp::sourceCpp("~/R/x86_64-pc-linux-gnu-library/4.0/DropletUtils/src/downsample_run_from_DropletUtils.cpp")
  for (z1 in 1:length(pats)) {
    system(paste0("aws s3 cp ",foldersList[z],"/cellranger/",pats[z1],"/outs/molecule_info.h5 /data/molecule_info.h5"))
    mol.info = read10xMolInfo("/data/molecule_info.h5")

    if (do_manual_downsampling)
    {
      #downsample.count.matrix <- downsampleReads("/data/molecule_info.h5",prop=0.5)
      saturation = 1-(length(mol.info$data$reads))/(sum(mol.info$data$reads))
      #downsample_arr = sample(x=seq(1,length(mol.info$data$reads),1), size=sum(mol.info$data$reads)/2, prob=mol.info$data$reads, replace=T)
      downsample_arr = sample(x=seq(1,length(mol.info$data$reads),1), size=downsamplereadcounts[z1]/origreadcounts[z1]*sum(mol.info$data$reads), prob=mol.info$data$reads, replace=T)
      downsample_arr_table = table(downsample_arr)
      downsample_saturation = 1-(length(downsample_arr_table))/(sum(downsample_arr_table))
      saturation_arr = c(saturation_arr, saturation)
      downsample_saturation_arr = c(downsample_saturation_arr, downsample_saturation)
      processed_pats_arr = c(processed_pats_arr, pats[z1])
    }
    else
    {
      testmat = downsampleReads("/data/molecule_info.h5",downsamplereadcounts[z1]/origreadcounts[z1])

      saturation = 1-(length(mol.info$data$reads))/(sum(mol.info$data$reads))
      new.read.counts = downsample_run(mol.info$data$reads, downsamplereadcounts[z1]/origreadcounts[z1])
      new.read.counts = new.read.counts[new.read.counts!=0]
      downsample_saturation = 1-(length(new.read.counts))/(sum(new.read.counts))

      saturation_arr = c(saturation_arr, saturation)
      downsample_saturation_arr = c(downsample_saturation_arr, downsample_saturation)
      processed_pats_arr = c(processed_pats_arr, pats[z1])
      
      out <- emptyDropsCellRanger(testmat)
      is.cell <- out$FDR <= 0.01
      testmat_filtered <- testmat[,which(is.cell),drop=FALSE]
      testmat_double_filtered = testmat_filtered[,colSums(testmat_filtered!=0)>=300]

      median_counts_arr = c(median_counts_arr, median(colSums(testmat)))
      median_genes_arr = c(median_genes_arr, median(colSums(testmat!=0)))
      filtered_median_counts_arr = c(filtered_median_counts_arr, median(colSums(testmat_filtered)))
      filtered_median_genes_arr = c(filtered_median_genes_arr, median(colSums(testmat_filtered!=0)))
      double_filtered_median_counts_arr = c(double_filtered_median_counts_arr, median(colSums(testmat_double_filtered)))
      double_filtered_median_genes_arr = c(double_filtered_median_genes_arr, median(colSums(testmat_double_filtered!=0)))
      double_filtered_genes_dist_arr = list.append(double_filtered_genes_dist_arr, colSums(testmat_double_filtered!=0))
      double_filtered_seurat_arr = list.append(double_filtered_seurat_arr, CreateSeuratObject(testmat_filtered, project = "sampling_to_uniform", assay = "RNA", min.cells = 0, min.features = 0, names.field = 1, names.delim = "_", meta.data = NULL))

      write10xCounts(paste0("/data/",pats[z1],"_downsampled_raw_feature_bc_matrix.h5"), testmat, gene.id=rownames(testmat), gene.symbol=rownames(testmat), barcodes=colnames(testmat), type="HDF5")
      system("rm /data/molecule_info.h5")
      system(paste0("aws s3 cp /data/",pats[z1],"_downsampled_raw_feature_bc_matrix.h5 ",foldersList[z],"/cellranger_downsampled/",pats[z1],"/",pats[z1],"_downsampled_raw_feature_bc_matrix.h5"))
      system(paste0("rm /data/",pats[z1],"_downsampled_raw_feature_bc_matrix.h5"))
    }
  }
}

names(saturation_arr) = processed_pats_arr
names(downsample_saturation_arr) = processed_pats_arr
names(median_counts_arr) = processed_pats_arr
names(median_genes_arr) = processed_pats_arr
names(filtered_median_counts_arr) = processed_pats_arr
names(filtered_median_genes_arr) = processed_pats_arr
names(double_filtered_median_counts_arr) = processed_pats_arr
names(double_filtered_median_genes_arr) = processed_pats_arr
names(double_filtered_genes_dist_arr) = processed_pats_arr
names(double_filtered_seurat_arr) = processed_pats_arr
for (i in 1:length(double_filtered_seurat_arr))
{
  double_filtered_seurat_arr[[i]]$orig.ident = processed_pats_arr[i]
}
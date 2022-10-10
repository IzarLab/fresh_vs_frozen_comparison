library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(ggpubr)
library(ggrastr)
library(Seurat)
library(infercnv)
library(ggpubr)
library(rlist)

#download object containing all ribas data from s3, subset to data for 310 sample
system("aws s3 cp s3://melanoma-ribas/ribas1/Seurat/integrated/ribas_integrated_titrate_thresh_integrated.rds /data/ribas_integrated_titrate_thresh_integrated.rds")
integrated_rds = readRDS("/data/ribas_integrated_titrate_thresh_integrated.rds")
system("rm /data/ribas_integrated_titrate_thresh_integrated.rds")
DefaultAssay(integrated_rds) = "RNA"
integrated_rds$placeholder = FALSE
integrated_rds$placeholder[grep("310",integrated_rds$orig.ident)] = TRUE
integrated_rds = subset(integrated_rds, placeholder)

#read in up and downregulated genes in immune resistance signature
#initialize arrays for average expression of up and downregulated genes in each cell, and total number of up and downregulated genes found in each cell
infercnv_pats = unique(integrated_rds$orig.ident)
up_genes = read.table("fresh_vs_frozen_comparison/jerbyarnon_TableS4A_immune_resistance_up_genes.txt",sep="\t",quote=NULL,header=T)$gene
integrated_rds$up_genes_average_infercnv = 0
integrated_rds$num_up_genes_found = 0
down_genes = read.table("fresh_vs_frozen_comparison/jerbyarnon_TableS4A_immune_resistance_down_genes.txt",sep="\t",quote=NULL,header=T)$gene
integrated_rds$down_genes_average_infercnv = 0
integrated_rds$num_down_genes_found = 0
fisher_pval_list = list()

for (infercnv_pat in infercnv_pats)
{
  #download infercnv data for each sample in ribas data, add infercnv copy number to rds object
  integrated_rds$barebarcodes = unlist(lapply(strsplit(colnames(integrated_rds),"_"), function(x) x[1]))
  integrated_rds$barebarcodes[integrated_rds$orig.ident!=infercnv_pat] = ""
  system(paste0("aws s3 cp s3://melanoma-ribas/ribas1/inferCNV/without_plots/inferCNV_cellbender_subcluster_",infercnv_pat,"_titrate_thresh/run.final.infercnv_obj /data/",infercnv_pat,"_infercnv.rds"))
  infercnv_obj = readRDS(paste0("/data/",infercnv_pat,"_infercnv.rds"))
  infercnv_obj@expr.data = infercnv_obj@expr.data[,infercnv_obj@observation_grouped_cell_indices$Melanocytes]

  #for each up and downregulated gene in immune resistance signature, check if gene is found in rds object
  #if so, increment number of up and downregulated genes found in each cell, and add expression of that gene in each cell
  for (up_gene in up_genes)
  {
    if (sum(rownames(infercnv_obj@expr.data)==up_gene)!=0 && sum(rownames(integrated_rds)==up_gene)!=0)
    {
      match_idxs = match(colnames(infercnv_obj@expr.data),integrated_rds$barebarcodes)
      integrated_rds$up_genes_average_infercnv[match_idxs] = integrated_rds$up_genes_average_infercnv[match_idxs] + infercnv_obj@expr.data[up_gene,]
      integrated_rds$num_up_genes_found[match_idxs] = integrated_rds$num_up_genes_found[match_idxs] + 1
      print(up_gene)
    }
  }

  for (down_gene in down_genes)
  {
    if (sum(rownames(infercnv_obj@expr.data)==down_gene)!=0 && sum(rownames(integrated_rds)==down_gene)!=0)
    {
      match_idxs = match(colnames(infercnv_obj@expr.data),integrated_rds$barebarcodes)
      integrated_rds$down_genes_average_infercnv[match_idxs] = integrated_rds$down_genes_average_infercnv[match_idxs] + infercnv_obj@expr.data[down_gene,]
      integrated_rds$num_down_genes_found[match_idxs] = integrated_rds$num_down_genes_found[match_idxs] + 1
      print(down_gene)
    }
  }

  #determine genes which have recurrent amplifications or deletions (greater than 75th percentile of copy number, or less than 25th percentile, in at least 10% of cells)
  fisher_pval_arr = c()
  linear_arr = infercnv_obj@expr.data[1:(length(rownames(infercnv_obj@expr.data))*length(colnames(infercnv_obj@expr.data)))]
  over_thresh = quantile(linear_arr,0.25)+0.01
  under_thresh = quantile(linear_arr,0.75)-0.01
  is_recurrent_over = rep(FALSE,length(rownames(infercnv_obj@expr.data)))
  is_recurrent_under = rep(FALSE,length(rownames(infercnv_obj@expr.data)))
  for (z in 1:length(rownames(infercnv_obj@expr.data)))
  {
    agene = rownames(infercnv_obj@expr.data)[z]
    num_over = sum(infercnv_obj@expr.data[agene,]>over_thresh)
    num_under = sum(infercnv_obj@expr.data[agene,]<under_thresh)
    if (num_over>length(colnames(infercnv_obj@expr.data))/10)
    {
      is_recurrent_over[z] = TRUE
    }
    if (num_under>length(colnames(infercnv_obj@expr.data))/10)
    {
      is_recurrent_under[z] = TRUE
    }
  }
  run_arr_breakpoints_list = list()
  for (abarcode in colnames(infercnv_obj@expr.data))
  {
    #determine genes which are part of contiguous runs of at least 10 recurrently amplified or deleted genes
    rle_res = rle(infercnv_obj@expr.data[,abarcode]>over_thresh & is_recurrent_over)
    is_in_run_arr_over = c()
    for (z in 1:length(rle_res$lengths))
    {
      if (rle_res$lengths[z]>10)
      {
        is_in_run_arr_over = c(is_in_run_arr_over, rep(rle_res$values[z],rle_res$lengths[z]))
      }
      else
      {
        is_in_run_arr_over = c(is_in_run_arr_over, rep(FALSE,rle_res$lengths[z]))
      }
    }
    names(is_in_run_arr_over) = rownames(infercnv_obj@expr.data)
    rle_res = rle(infercnv_obj@expr.data[,abarcode]<under_thresh & is_recurrent_under)
    is_in_run_arr_under = c()
    for (z in 1:length(rle_res$lengths))
    {
      if (rle_res$lengths[z]>10)
      {
        is_in_run_arr_under = c(is_in_run_arr_under, rep(rle_res$values[z],rle_res$lengths[z]))
      }
      else
      {
        is_in_run_arr_under = c(is_in_run_arr_under, rep(FALSE,rle_res$lengths[z]))
      }
    }
    names(is_in_run_arr_under) = rownames(infercnv_obj@expr.data)

    #store locations where runs of recurrently amplified or deleted genes begin and end run_arr_breakpoints variable
    run_arr_breakpoints = c()
    run_arr_over_rle = rle(is_in_run_arr_over)
    runninglength = 0
    for (z in 1:length(run_arr_over_rle$lengths))
    {
      runninglength = runninglength + run_arr_over_rle$lengths[z]
      run_arr_breakpoints = c(run_arr_breakpoints, runninglength)
    }
    run_arr_under_rle = rle(is_in_run_arr_under)
    runninglength = 0
    for (z in 1:length(run_arr_under_rle$lengths))
    {
      runninglength = runninglength + run_arr_under_rle$lengths[z]
      run_arr_breakpoints = c(run_arr_breakpoints, runninglength)
    }
    run_arr_breakpoints_list = list.append(run_arr_breakpoints_list, sort(unique(run_arr_breakpoints)))

    #determine numbers of up and downregulated immune resistance genes that are also part of recurrently amplified or deleted runs
    match_up_genes = intersect(rownames(infercnv_obj@expr.data), up_genes)
    match_up_genes_cnv = infercnv_obj@expr.data[match_up_genes,abarcode]
    match_up_genes_cnv_over_thresh = sum(is_in_run_arr_over[match_up_genes])#sum(match_up_genes_cnv>1.1)
    match_up_genes_cnv_under_thresh = sum(is_in_run_arr_under[match_up_genes])#sum(match_up_genes_cnv<0.9)
    
    match_down_genes = intersect(rownames(infercnv_obj@expr.data), down_genes)
    match_down_genes_cnv = infercnv_obj@expr.data[match_down_genes,abarcode]
    match_down_genes_cnv_over_thresh = sum(is_in_run_arr_over[match_down_genes])#sum(match_down_genes_cnv>1.1)
    match_down_genes_cnv_under_thresh = sum(is_in_run_arr_under[match_down_genes])#sum(match_down_genes_cnv<0.9)

    #perform fisher's exact test among up or downregulated genes located in amplified or deleted regions
    test_mat = matrix(0,2,2)
    rownames(test_mat) = c("Amplified Region","Deleted Region")
    colnames(test_mat) = c("Up Gene","Down Gene")
    test_mat[1,1] = match_up_genes_cnv_over_thresh
    test_mat[1,2] = match_down_genes_cnv_over_thresh
    test_mat[2,1] = match_up_genes_cnv_under_thresh
    test_mat[2,2] = match_down_genes_cnv_under_thresh
    fisherres = fisher.test(test_mat)
    fisher_pval_arr = c(fisher_pval_arr, fisherres$p.value)
  }
  fisher_pval_list[[infercnv_pat]] = fisher_pval_arr

  #begin heatmap code
  if (infercnv_pat=="ribas_310_on_GEX_5pv2_S27_L004")
  {
    shortname = "on"
  }
  if (infercnv_pat=="ribas_310_on_later_previd_3_GEX")
  {
    shortname = "on_later"
  }
  if (infercnv_pat=="ribas_310_pre_GEX_5pv2_S26_L004")
  {
    shortname = "pre"
  }

  #reformat copy number matrix from infercnv object, using melt function
  expr_data = infercnv_obj@expr.data
  rownames(expr_data) = (1:length(rownames(expr_data)))
  colnames(expr_data) = (1:length(colnames(expr_data)))
  df_exprdata = melt(infercnv_obj@expr.data)
  names(df_exprdata) = c("gene_num","cell_num","infercnv")
  df_exprdata$gene_num = as.integer(df_exprdata$gene_num)
  df_exprdata$cell_num = as.integer(df_exprdata$cell_num)

  #clip copy number values to 1.15 and .93, imitating procedure in infercnv figure generation
  df_exprdata$infercnv[df_exprdata$infercnv>1.15] = 1.15
  df_exprdata$infercnv[df_exprdata$infercnv<.93] = .93

  #calculate mindist as the smallest difference between the median of the copy number matrix and any non-median copy number value, divided by 10
  absdist = abs(df_exprdata$infercnv-median(df_exprdata$infercnv))
  mindist = min(absdist[absdist!=0])/10
  #df_exprdata$infercnv[absdist==0] = median(df_exprdata$infercnv)
  #df_exprdata$infercnv[absdist<=.01] = median(df_exprdata$infercnv)

  #plot heat map of copy number matrix
  hm <- ggplot(data = df_exprdata, aes(x = factor(gene_num), y = cell_num, fill = infercnv)) + rasterise(geom_tile()) + 
      #scale_fill_distiller(name = "InferCNV infercnv", palette = "RdBu", direction = 1, na.value = "transparent") +
      scale_fill_gradientn(breaks = c(min(df_exprdata$infercnv), median(df_exprdata$infercnv)-.02, median(df_exprdata$infercnv), median(df_exprdata$infercnv)+.02, max(df_exprdata$infercnv)), colors = c("blue","white","white","white","red")) + 
      scale_x_discrete(breaks = unique(df_exprdata$gene_num), labels = unique(df_exprdata$gene_num)) + theme_gray() +
      theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_text(size = 30), legend.key.size = unit(1,"cm"), legend.text = element_text(size = 20), axis.text = element_text(size = 20), axis.title = element_text(size = 30), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.key.width= unit(8, 'cm')) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

  tmp <- ggplot_gtable(ggplot_build(hm))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]

  hm.clean <- hm# +theme(axis.title.y = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(), axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position="none")

  #add vertical lines to heatmap for chromosome lengths
  chrlengths = rle(as.vector(infercnv_obj@gene_order$chr))$lengths
  runninglength = 0
  #nonsense = nonsense+1
  for (z in 1:length(chrlengths))
  {
    runninglength = runninglength + chrlengths[z]
    hm.clean <- hm.clean + geom_vline(xintercept = runninglength)
    if (z<=2)
    {
      #hm.clean <- hm.clean + geom_text(aes(label=paste0("chr",chrlengths[z]), y=-600))
    }
  }

  tempdf2 = data.frame(location = integer(), cellnum = integer())
  for (z in 1:length(run_arr_breakpoints_list))
  {
    tempdf = data.frame(location = run_arr_breakpoints_list[[z]], cellnum = z)
    tempdf2 = rbind(tempdf2, tempdf)
  }
  # don't create boxes for run arrs
  # hm.clean <- hm.clean + geom_point(data = tempdf2, mapping = aes(x=location, y=cellnum), size = 0.5, inherit.aes = FALSE)

  #create dataframe with immune resistance genes, and corresponding locations in copy number matrix
  up_genes = read.table("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_comparison/jerbyarnon_TableS4A_immune_resistance_up_genes.txt",sep="\t",quote=NULL,header=T)$gene
  matchidxs = match(up_genes,rownames(infercnv_obj@gene_order))
  df_icr = data.frame(gene_name = up_genes, gene_loc = matchidxs, type = "icr_up_genes")
  down_genes = read.table("/mnt/vdb/home/ubuntu2/fresh_vs_frozen_comparison/jerbyarnon_TableS4A_immune_resistance_down_genes.txt",sep="\t",quote=NULL,header=T)$gene
  matchidxs = match(down_genes,rownames(infercnv_obj@gene_order))
  df_icr2 = data.frame(gene_name = down_genes, gene_loc = matchidxs, type = "icr_down_genes")
  df_icr = rbind(df_icr, df_icr2)
  df_icr = df_icr[!is.na(df_icr$gene_loc),]
  df_icr = df_icr[order(df_icr$gene_loc),]

  theme_set(theme_bw())
  #plot histogram of up and down regulated genes in immune resistance signature
  dens_figure <- ggplot(df_icr, aes(x=gene_loc, fill=type)) + scale_fill_manual(breaks=c("icr_down_genes","icr_up_genes"), values=c("blue","red")) + geom_density(alpha=0.5) + scale_x_continuous(expand = c(0, 0)) + theme(legend.position = "none", plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + labs(x = "ICR Gene Density")
  hist_figure <- ggplot(df_icr, aes(x=gene_loc, fill=type)) + scale_fill_manual(breaks=c("icr_down_genes","icr_up_genes"), values=c("blue","red")) + geom_histogram(alpha=0.2, position = "identity") + scale_x_continuous(expand = c(0, 0)) + theme(legend.position = "none", plot.margin = margin(t = 0, r = 0, b = 0, l = 0), axis.text = element_text(size = 20), axis.title = element_text(size = 30)) + labs(x = "ICR Gene Density")

  #add vertical bars for chromosome lengths
  runninglength = 0
  for (z in 1:length(chrlengths))
  {
    runninglength = runninglength + chrlengths[z]
    dens_figure <- dens_figure + geom_vline(xintercept = runninglength)
    hist_figure <- hist_figure + geom_vline(xintercept = runninglength)
  }

  pdf(paste0("r310_infercnv_and_resistance_sig_heatmap_",shortname,".pdf"),width=20,height=20)
  #print(grid.arrange(hist_figure, dens_figure, hm.clean, nrow = 3, ncol = 1, widths = c(30), heights = c(4, 4, 15)))
  print(grid.arrange(hist_figure, hm.clean, nrow = 2, ncol = 1, widths = c(30), heights = c(5, 15)))
  dev.off()
}

fisher_pval_df = data.frame(sample = character(), fisher_pval = double())
for (infercnv_pat in names(fisher_pval_list))
{
  if (infercnv_pat=="ribas_310_on_GEX_5pv2_S27_L004")
  {
    shortname = "on"
  }
  if (infercnv_pat=="ribas_310_on_later_previd_3_GEX")
  {
    shortname = "on_later"
  }
  if (infercnv_pat=="ribas_310_pre_GEX_5pv2_S26_L004")
  {
    shortname = "pre"
  }
  for (z in 1:length(fisher_pval_list[[infercnv_pat]]))
  {
    tempdf = data.frame(sample = shortname, fisher_pval = fisher_pval_list[[infercnv_pat]][z])
    fisher_pval_df = rbind(fisher_pval_df, tempdf)
  }
}

fisher_pval_df$sample[fisher_pval_df$sample=="pre"] = "Pre-treatment"
fisher_pval_df$sample[fisher_pval_df$sample=="on"] = "On-treatment"
fisher_pval_df$sample[fisher_pval_df$sample=="on_later"] = "On-treatment (later)"

fisher_pval_df$sample = factor(fisher_pval_df$sample, levels=c("Pre-treatment","On-treatment","On-treatment (later)"))

#ggboxplot(fisher_pval_df, x="sample", y="fisher_pval", color="sample") + theme(axis.text.x = element_text(angle = 0))
ggplot(fisher_pval_df, aes(x=sample, y=fisher_pval, fill=sample))+ geom_violin() + theme(axis.text.x = element_text(angle = 0), legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text = element_text(size=20), axis.title = element_text(size=20))# + scale_x_discrete(breaks = c("pre","on","on_later"))
ggsave("r310_infercnv_and_resistance_sig_fisher_vlnplot.pdf",width=12,height=10)

for (infercnv_pat in infercnv_pats)
{
  system(paste0("aws s3 cp s3://melanoma-ribas/ribas1/inferCNV/without_plots/inferCNV_cellbender_subcluster_",infercnv_pat,"_titrate_thresh/run.final.infercnv_obj /data/",infercnv_pat,"_infercnv.rds"))
  infercnv_obj = readRDS(paste0("/data/",infercnv_pat,"_infercnv.rds"))
  infercnv_obj@expr.data = infercnv_obj@expr.data[,infercnv_obj@observation_grouped_cell_indices$Melanocytes]

}
nonsense = nonsense+1

integrated_rds$up_genes_average_infercnv = integrated_rds$up_genes_average_infercnv/integrated_rds$num_up_genes_found
integrated_rds$down_genes_average_infercnv = integrated_rds$down_genes_average_infercnv/integrated_rds$num_down_genes_found

integrated_rds = subset(integrated_rds, celltype_bped_main=="Melanocytes")
# object.list = c(subset(integrated_rds, orig.ident==infercnv_pats[1]), subset(integrated_rds, orig.ident==infercnv_pats[2]), subset(integrated_rds, orig.ident==infercnv_pats[3]))
# anchors = FindIntegrationAnchors(object.list = object.list, dims = 1:20)
# integrated_rds = IntegrateData(anchorset = anchors, dims = 1:20)
# integrated_rds = ScaleData(object = integrated_rds)
# integrated_rds = RunPCA(object = integrated_rds)
# integrated_rds = FindNeighbors(integrated_rds, dims = 1:15)
# integrated_rds = FindClusters(integrated_rds)
# integrated_rds = RunUMAP(object = integrated_rds, dims = 1:20)

integrated_rds = AddModuleScore(integrated_rds, features = list(na.omit(up_genes)), name = "immune_resistance_up_genes", assay = "RNA", search = T)
integrated_rds = AddModuleScore(integrated_rds, features = list(na.omit(down_genes)), name = "immune_resistance_down_genes", assay = "RNA", search = T)
# pdf("r310_infercnv_and_resistance_sig.pdf")
# print(FeaturePlot(integrated_rds, features=c("up_genes_average_infercnv"), min.cutoff = "q0", max.cutoff = "q100"))
# print(FeaturePlot(integrated_rds, features=c("immune_resistance_up_genes1"), min.cutoff = "q0", max.cutoff = "q100"))
# print(FeaturePlot(integrated_rds, features=c("down_genes_average_infercnv"), min.cutoff = "q0", max.cutoff = "q100"))
# print(FeaturePlot(integrated_rds, features=c("immune_resistance_down_genes1"), min.cutoff = "q0", max.cutoff = "q100"))

print(cor.test(integrated_rds$immune_resistance_up_genes1, integrated_rds$up_genes_average_infercnv, method="spearman"))
print(cor.test(integrated_rds$immune_resistance_down_genes1, integrated_rds$down_genes_average_infercnv, method="spearman"))
dev.off()
library(Seurat)
library(ggplot2)
library(stringr)
library(vegan)
library(viridis)
library(seriation)
library(RColorBrewer)

### title: Create barplots of cell type frequency, immune and non-immune simpson diversity, and malignant cell fractions
### author: Yiping Wang date: 11/08/2022

#load in folder names for each dataset
useRibas = FALSE
if (useRibas)
{
  foldersList = c("")
  integrated_name_arr = c("ribas_310")
  integrated_name_arr_underscore = c("ribas_310")
  outputname = "ribas_310"
} else {
  foldersList = c("",
    "s3://fresh-vs-frozen-comparison-ohio/cpoi-uvealprimarydata",
    "s3://fresh-vs-frozen-comparison-ohio/nsclc",
    "s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline",
    "s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
  integrated_name_arr = c("BI5","cpoi-uvealprimarydata","nsclc","BI5","NR1")
  integrated_name_arr_underscore = c("BI5","cpoi_uvealprimarydata","nsclc","BI5","NR1")
  outputname = "fresh_vs_frozen"
}

immune_simpson_list = c()
nonimmune_simpson_list = c()

immune_celltypes = c("B cells/Plasma cells","Myeloid","T cells/NK cells","NK cells","Dendritic cells")

bardf = data.frame(sample=character(), celltype=character(), frequency=integer())

for (i in 1:length(foldersList))
{
  #load in rds files for each dataset
  using_slyper = FALSE
  if (useRibas)
  {
    integrated_rds = readRDS("/data/fresh_vs_frozen_all_reannotate_ribas_310.rds")
  }
  else
  {
    if (i==1)
    {
      integrated_rds = readRDS("/data/fresh_vs_frozen_all_reannotate_BI5.rds")
    }
    else
    {
      if (foldersList[i]=="s3://fresh-vs-frozen-comparison-ohio/slyper_pipeline")
      {
        using_slyper = TRUE
      }
      system(paste0("aws s3 cp ",foldersList[i],"/Seurat/integrated/",integrated_name_arr[i],"_integrated.rds /data/",integrated_name_arr[i],"_integrated.rds"))
      integrated_rds = readRDS(paste0("/data/",integrated_name_arr[i],"_integrated.rds"))
      system(paste0("rm /data/",integrated_name_arr[i],"_integrated.rds"))
    }
  }

  if (using_slyper)
  {
    integrated_rds = subset(integrated_rds, orig.ident!="3snseq")
    integrated_rds = subset(integrated_rds, orig.ident!="SNSEQ_3P_NI")
    integrated_rds = subset(integrated_rds, orig.ident!="SNSEQ_3P_WI")
  }

  #calculate Simpson diversity indeces for immune and non-immune cell types
  integrated_rds$orig.ident = paste0(integrated_name_arr_underscore[i],"_",integrated_rds$orig.ident)
  unique_idents = unique(integrated_rds$orig.ident)
  for (i1 in 1:length(unique_idents))
  {
    diversity_table = table(integrated_rds$manual_annotation_label[integrated_rds$orig.ident==unique_idents[i1]])
    immune_simpson_list = c(immune_simpson_list,diversity(diversity_table[names(diversity_table) %in% immune_celltypes],index="simpson"))
    nonimmune_simpson_list = c(nonimmune_simpson_list,diversity(diversity_table[!(names(diversity_table) %in% immune_celltypes)],index="simpson"))
  }
  names(immune_simpson_list)[(length(immune_simpson_list)-length(unique_idents)+1):length(immune_simpson_list)] = unique_idents
  names(nonimmune_simpson_list)[(length(nonimmune_simpson_list)-length(unique_idents)+1):length(nonimmune_simpson_list)] = unique_idents

  #store frequency of each cell type for each sample in a dataframe
  unique_idents = unique(integrated_rds$orig.ident)
  unique_celltypes = unique(integrated_rds$manual_annotation_label)
  for (i1 in 1:length(unique_idents))
  {
    for (j1 in 1:length(unique_celltypes))
    {
      tempdf = data.frame(sample=unique_idents[i1], celltype=unique_celltypes[j1], frequency=sum(integrated_rds$orig.ident==unique_idents[i1] & integrated_rds$manual_annotation_label==unique_celltypes[j1]))
      bardf = rbind(bardf, tempdf)
    }
  }
}

bardf$celltype[bardf$celltype=="Melanocytes"] = "Melanoma cells"
bardf$celltype[bardf$celltype=="Lung cells"] = "Lung cancer cells"

if (useRibas)
{
  #generate random color palette to color each cell type in bar plot
  pdf(paste0("fresh_vs_frozen_cell_comp/",outputname,"_cell_comp.pdf"),height=7,width=5)
  coul = brewer.pal(9,"Set1")
  scatter_colors = colorRampPalette(coul)(length(unique(bardf$celltype)))
  scatter_colors = scatter_colors[c(seq(1,length(scatter_colors),3),seq(2,length(scatter_colors),3),seq(3,length(scatter_colors),3))]
  names(scatter_colors) = unique(bardf$celltype)
  scatter_colors = scatter_colors[order(names(scatter_colors))]

  #make new sample labels to plot on barplot
  unique_samples = sort(unique(bardf$sample))
  unique_samples = unique_samples[c(3,1,2)]
  label_samples = unique_samples
  label_samples[grep("on_GEX",label_samples)] = "on"
  label_samples[grep("pre_GEX",label_samples)] = "pre"
  label_samples[grep("on_later",label_samples)] = "on_later"

  #plot barplot of cell type frequency
  print(ggplot(bardf, aes(fill=celltype, y=frequency, x=sample, color="black")) + scale_fill_manual(breaks = names(scatter_colors), values = scatter_colors) + scale_colour_manual(breaks = c("black"), values = c("black")) + scale_x_discrete(limits = unique_samples, labels=label_samples) + geom_bar(position="fill", stat="identity") + guides(color = "none") + xlab("sample") + theme(axis.text.x = element_text(angle = 90, hjust = 1)))

  dev.off()

  #calculate frequencies of tumor cells only versus other cell types, plot barplot
  malig_df = data.frame(sample=character(), malignant=character(), frequency=integer())
  for (unique_sample in unique_samples)
  {
    tempdf = data.frame(sample=unique_sample, malignant="Yes", frequency=sum(bardf$frequency[bardf$sample==unique_sample & (bardf$celltype %in% c("Melanoma cells","Lung cancer cells"))]))
    malig_df = rbind(malig_df, tempdf)
    tempdf = data.frame(sample=unique_sample, malignant="No", frequency=sum(bardf$frequency[bardf$sample==unique_sample & !(bardf$celltype %in% c("Melanoma cells","Lung cancer cells"))]))
    malig_df = rbind(malig_df, tempdf)
  }
  pdf(paste0("fresh_vs_frozen_cell_comp/",outputname,"_malignant.pdf"),height=7,width=5)
  print(ggplot(malig_df, aes(fill=malignant, y=frequency, x=sample, color="black")) + scale_fill_manual(breaks = c("Yes", "No"), values = c("red","blue")) + scale_colour_manual(breaks = c("black"), values = c("black")) + scale_x_discrete(limits = unique_samples, labels=label_samples) + geom_bar(position="fill", stat="identity") + guides(color = "none") + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  dev.off()

  #print barplot of simpson diversity for immune cells
  simpsondf = data.frame(simpson_diversity=immune_simpson_list, sample=names(immune_simpson_list))
  simpsondf = simpsondf[order(simpsondf$sample),]
  pdf(paste0("fresh_vs_frozen_cell_comp/",outputname,"_immune_simpson_diversity.pdf"),height=7,width=3)
  theme_set(theme_bw())
  print(ggplot(simpsondf, aes(y=simpson_diversity, x=sample)) + geom_bar(stat="identity") + scale_x_discrete(limits = unique_samples, labels=label_samples) + ylab("Immune_Simpson Diversity") + xlab("sample") + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  dev.off()

  #print barplot of simpson immunity for non-immune cells
  simpsondf = data.frame(simpson_diversity=nonimmune_simpson_list, sample=names(nonimmune_simpson_list))
  simpsondf = simpsondf[order(simpsondf$sample),]
  pdf(paste0("fresh_vs_frozen_cell_comp/",outputname,"_nonimmune_simpson_diversity.pdf"),height=7,width=3)
  theme_set(theme_bw())
  print(ggplot(simpsondf, aes(y=simpson_diversity, x=sample)) + geom_bar(stat="identity") + scale_x_discrete(limits = unique_samples, labels=label_samples) + ylab("Nonimmune_Simpson Diversity") + ylim(c(0,1)) + xlab("sample") + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  dev.off()
} else {
  pdf(paste0("fresh_vs_frozen_cell_comp/",outputname,"_cell_comp.pdf"),height=7,width=12)

  #rename sample IDs to more human-readable format
  rename_arr_list = list(bardf$sample, names(immune_simpson_list), names(nonimmune_simpson_list))
  for (i in 1:length(rename_arr_list))
  {
    rename_arr = rename_arr_list[[i]]
    rename_arr[rename_arr=="BI5_CD45neg"] = "Mel_sc_5_CD45-"
    rename_arr[rename_arr=="BI5_CD45pos"] = "Mel_sc_5_CD45+"
    rename_arr[rename_arr=="BI5_3snseq"] = "Mel_sn_3"
    rename_arr[rename_arr=="BI5_5pv2-snseq"] = "Mel_sn_5v2"
    rename_arr[rename_arr=="BI5_5snseq"] = "Mel_sn_5"
    rename_arr[rename_arr=="cpoi_uvealprimarydata_SCRNA-5P-NA-E12"] = "UM_sc_5"
    rename_arr[rename_arr=="cpoi_uvealprimarydata_SCRNA-5P-NA-F1"] = "UM_sc_5_CD45+"
    rename_arr[rename_arr=="cpoi_uvealprimarydata_SNRNA-5P-WI-F12"] = "UM_sn_5_inhib"
    rename_arr[rename_arr=="nsclc_SCRNA_5P_NA"] = "NSCLC_sc_5"
    rename_arr[rename_arr=="nsclc_5pv2-snseq"] = "NSCLC_sn_5v2"
    rename_arr[rename_arr=="nsclc_SNSEQ_3P_NI"] = "NSCLC_sn_3"
    rename_arr[rename_arr=="nsclc_SNSEQ_3P_WI"] = "NSCLC_sn_3_inhib"
    rename_arr[rename_arr=="nsclc_SNSEQ_5P_NI"] = "NSCLC_sn_5"
    rename_arr[rename_arr=="nsclc_SNSEQ_5P_WI"] = "NSCLC_sn_5_inhib"
    if (i==1)
    {
      bardf$sample = rename_arr
    }
    else if (i==2)
    {
      names(immune_simpson_list) = rename_arr
    }
    else if (i==3)
    {
      names(nonimmune_simpson_list) = rename_arr
    }
  }

  #add dummy entries to barplot dataframe, in order to create spaces between datasets on final barplot
  bardf2 = bardf
  bardf_temp = data.frame(sample="Bz",celltype="Mitochondrial",frequency=0)
  bardf2 = rbind(bardf2, bardf_temp)
  bardf_temp = data.frame(sample="Mz",celltype="Mitochondrial",frequency=0)
  bardf2 = rbind(bardf2, bardf_temp)
  bardf_temp = data.frame(sample="NSCLCz",celltype="Mitochondrial",frequency=0)
  bardf2 = rbind(bardf2, bardf_temp)
  bardf_temp = data.frame(sample="Uz",celltype="Mitochondrial",frequency=0)
  bardf2 = rbind(bardf2, bardf_temp)
  bardf2 = bardf2[order(bardf2$sample),]

  #reorder sample names, so that labels are in logical order in final plot
  unique_samples = sort(unique(bardf2$sample))

  #unique_samples = unique_samples[c(1,2,3)]
  #unique_samples = unique_samples[c(14,15,16,17,18,19,20,4,5,6,7,8,9,21,22,23,24,1,2,3,10,11)]
  unique_samples = unique_samples[c(12,17,13,14,15,16,18,4,5,6,7,8,9,19,20,21,22,1,2,3,10,11)]

  #set array of colors for fresh/frozen status in bar plot
  x_colors = rep("red",length(unique_samples))
  x_colors[unique_samples=="Mel_sc_5_CD45-"] = "blue"
  x_colors[unique_samples=="Mel_sc_5_CD45+"] = "blue"
  x_colors[unique_samples=="UM_sc_5"] = "blue"
  x_colors[unique_samples=="UM_sc_5_CD45+"] = "blue"
  x_colors[unique_samples=="NSCLC_sc_5"] = "blue"
  names(x_colors) = unique_samples

  #modify sample labels on barplot to make them look nicer
  label_samples = unique_samples
  label_samples = str_replace(label_samples,"nsclc_","")
  label_samples = str_replace(label_samples,"BI5_","")
  label_samples = str_replace(label_samples,"cpoi_uvealprimarydata_","")
  label_samples[label_samples=="Bz"] = ""
  label_samples[label_samples=="Mz"] = ""
  label_samples[label_samples=="NSCLCz"] = ""
  label_samples[label_samples=="Uz"] = ""

  #create random array of colors for each cell type
  coul = brewer.pal(9,"Set1")
  scatter_colors = colorRampPalette(coul)(length(unique(bardf$celltype)))
  scatter_colors = scatter_colors[c(seq(1,length(scatter_colors),3),seq(2,length(scatter_colors),3),seq(3,length(scatter_colors),3))]
  names(scatter_colors) = unique(bardf$celltype)
  scatter_colors = scatter_colors[order(names(scatter_colors))]

  #create one large string that will label each dataset on lower part of x-axis label
  #use spaces to position dataset labels beneath barplot for each dataset
  xlab_string = paste0(strrep(" ",30),"NSCLC",strrep(" ",27),strrep(" ",30),"Mel",strrep(" ",27),strrep(" ",16),"UM",strrep(" ",13),strrep(" ",10),"BI5 slyper",strrep(" ",10),strrep(" ",10),"NR1 slyper",strrep(" ",10))

  #plot barplot of cell type frequency
  print(ggplot(bardf2, aes(fill=celltype, y=frequency, x=sample, color="black")) + scale_fill_manual(breaks = names(scatter_colors), values = scatter_colors) + scale_colour_manual(breaks = c("black"), values = c("black")) + scale_x_discrete(limits = unique_samples, labels=label_samples) + geom_bar(position="fill", stat="identity") + guides(color = "none") + xlab(xlab_string) + theme(axis.text.x = element_text(angle = 90, hjust = 1, color = x_colors), axis.title.x = element_text(angle = 0, hjust = 0)))
  dev.off()

  #calculate frequencies of tumor cells only versus other cell types, plot barplot
  malig_df = data.frame(sample=character(), malignant=character(), frequency=integer())
  for (unique_sample in unique_samples)
  {
    tempdf = data.frame(sample=unique_sample, malignant="Yes", frequency=sum(bardf2$frequency[bardf2$sample==unique_sample & (bardf2$celltype %in% c("Melanoma cells","Lung cancer cells"))]))
    malig_df = rbind(malig_df, tempdf)
    tempdf = data.frame(sample=unique_sample, malignant="No", frequency=sum(bardf2$frequency[bardf2$sample==unique_sample & !(bardf2$celltype %in% c("Melanoma cells","Lung cancer cells"))]))
    malig_df = rbind(malig_df, tempdf)
  }
  pdf(paste0("fresh_vs_frozen_cell_comp/",outputname,"_malignant.pdf"),height=7,width=12)
  print(ggplot(malig_df, aes(fill=malignant, y=frequency, x=sample, color="black")) + scale_fill_manual(breaks = c("Yes", "No"), values = c("red","blue")) + scale_colour_manual(breaks = c("black"), values = c("black")) + scale_x_discrete(limits = unique_samples, labels=label_samples) + geom_bar(position="fill", stat="identity") + guides(color = "none") + xlab(xlab_string) + theme(axis.text.x = element_text(angle = 90, hjust = 1, color = x_colors), axis.title.x = element_text(angle = 0, hjust = 0)))
  dev.off()

  xlab_string_simpson = paste0(strrep(" ",35),"NSCLC",strrep(" ",35),strrep(" ",27),"Mel",strrep(" ",27),strrep(" ",16),"UM",strrep(" ",16))

  #print barplot of simpson immunity for immune cells
  simpsondf = data.frame(simpson_diversity=immune_simpson_list, sample=names(immune_simpson_list))
  simpsondf_temp = data.frame(sample="Mz",simpson_diversity=0)
  simpsondf = rbind(simpsondf, simpsondf_temp)
  simpsondf_temp = data.frame(sample="Nz",simpson_diversity=0)
  simpsondf = rbind(simpsondf, simpsondf_temp)
  simpsondf = simpsondf[order(simpsondf$sample),]
  pdf(paste0("fresh_vs_frozen_cell_comp/",outputname,"_immune_simpson_diversity.pdf"),height=7,width=12)
  theme_set(theme_bw())
  print(ggplot(simpsondf, aes(y=simpson_diversity, x=sample)) + geom_bar(stat="identity") + scale_x_discrete(limits = unique_samples, labels=label_samples) + ylab("Immune_Simpson Diversity") + xlab(xlab_string_simpson) + theme(axis.text.x = element_text(angle = 90, hjust = 1, color = x_colors), axis.title.x = element_text(angle = 0, hjust = 0)))
  dev.off()

  #print barplot of simpson immunity for non-immune cells
  simpsondf = data.frame(simpson_diversity=nonimmune_simpson_list, sample=names(nonimmune_simpson_list))
  simpsondf_temp = data.frame(sample="Mz",simpson_diversity=0)
  simpsondf = rbind(simpsondf, simpsondf_temp)
  simpsondf_temp = data.frame(sample="Nz",simpson_diversity=0)
  simpsondf = rbind(simpsondf, simpsondf_temp)
  simpsondf = simpsondf[order(simpsondf$sample),]
  pdf(paste0("fresh_vs_frozen_cell_comp/",outputname,"_nonimmune_simpson_diversity.pdf"),height=7,width=12)
  theme_set(theme_bw())
  print(ggplot(simpsondf, aes(y=simpson_diversity, x=sample)) + geom_bar(stat="identity") + scale_x_discrete(limits = unique_samples, labels=label_samples) + ylab("Nonimmune_Simpson Diversity") + ylim(c(0,1)) + xlab(xlab_string_simpson) + theme(axis.text.x = element_text(angle = 90, hjust = 1, color = x_colors), axis.title.x = element_text(angle = 0, hjust = 0)))
  dev.off()

}
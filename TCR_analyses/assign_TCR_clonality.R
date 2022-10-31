assign_TCR_clonality <- function(seu, dataset_name) {
  #rename orig.ident to more human-readable names
  seu$orig.ident[seu$orig.ident=="ribas_310_on"] = "ribas1_on_tcr_S36_L004"
  seu$orig.ident[seu$orig.ident=="ribas_310_on_later"] = "ribas_310_on_later_previd_3_TCR"
  seu$orig.ident[seu$orig.ident=="ribas_310_pre"] = "ribas1_pre_tcr_S35_L004"
  seu$orig.ident[seu$orig.ident=="CD45pos"] = "TCRBI5_S1_L001"
  seu$orig.ident[seu$orig.ident=="5snseq"] = "bi005-skcm-5snseq-TCR"
  seu$orig.ident[seu$orig.ident=="SCRNA-5P-NA-E12"] = "UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F2"
  seu$orig.ident[seu$orig.ident=="SCRNA-5P-NA-F1"] = "UMEL-CUUM1-SCRNA-5P-NA-PRIMARY-TCR-F3"
  seu$orig.ident[seu$orig.ident=="SNRNA-5P-WI-F12"] = "UMEL-CUUM1-SNRNA-5P-WI-PRIMARY-TCR-F9"

  seu$barebarcodes = unlist(lapply(strsplit(colnames(seu),"_"), function(x) {x[1]}))

  #read in TR_frequency.csv files for each sample, which contain barcodes associated with each clonotype
  #add in clonotype lengths into clonality field
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

  #annotate barcodes that are unmatched with TCR sequencing, and therefore not found in TR_frequency.csv files
  #further annotate as CD4 or 8 T-cells for cutaneous melanoma and uveal melanoma primary datasets
  #annotate size range of clonotypes in clonality_group field
  if (dataset_name=="ribas" || dataset_name=="UMEL")
  {
    seu$clonality_group = "Unmatched with TCR sequencing"
  }
  if (dataset_name=="BI5" || dataset_name=="UM")
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

  return(seu)
}
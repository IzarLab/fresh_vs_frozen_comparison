### title: Print out csv files of count matrices for sequential cutaneous melanoma samples
### author: Yiping Wang date: 11/08/2022

system("aws s3 cp s3://melanoma-ribas/ribas1/Seurat/integrated/ribas_integrated_titrate_thresh_integrated.rds /data/ribas_integrated_titrate_thresh_integrated.rds")
integrated_rds = readRDS("/data/ribas_integrated_titrate_thresh_integrated.rds")
system("rm /data/ribas_integrated_titrate_thresh_integrated.rds")

r310_pre = subset(integrated_rds, orig.ident=="ribas_310_pre_GEX_5pv2_S26_L004")
r310_pre = subset(r310_pre, celltype_bped_main=="Melanocytes")
r310_on = subset(integrated_rds, orig.ident=="ribas_310_on_GEX_5pv2_S27_L004")
r310_on = subset(r310_on, celltype_bped_main=="Melanocytes")
r310_on_later = subset(integrated_rds, orig.ident=="ribas_310_on_later_previd_3_GEX")
r310_on_later = subset(r310_on_later, celltype_bped_main=="Melanocytes")

write.table(as.matrix(r310_pre@assays$RNA@counts), file=paste0("/data/r310_pre_counts.csv"), row.names=T, col.names=T, quote=FALSE, sep=",")
write.table(as.matrix(r310_on@assays$RNA@counts), file=paste0("/data/r310_on_counts.csv"), row.names=T, col.names=T, quote=FALSE, sep=",")
write.table(as.matrix(r310_on_later@assays$RNA@counts), file=paste0("/data/r310_on_later_counts.csv"), row.names=T, col.names=T, quote=FALSE, sep=",")

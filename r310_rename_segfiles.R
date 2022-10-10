wes_table = read.table("r310_pre_rg.sorted.bam_CNVs_manual.seg",quote=NULL,header=T,sep="\t")
wes_table$sample = "Pre-treatment WES"
write.table(wes_table,"Pre-treatment WES.seg",sep="\t",col.names=T,row.names=F,quote=FALSE)
#system("cp r310_pre_rg.sorted.bam_CNVs_manual.seg Pre-treatment\\\ WES.seg")

wgs_table = read.table("Ribas-1Pre_S102_L001_tumor.seg",quote=NULL,header=T,sep="\t")
wgs_table$sample = "Pre-treatment US-WGS"
write.table(wgs_table,"Pre-treatment US-WGS.seg",sep="\t",col.names=T,row.names=F,quote=FALSE)
#system("cp Ribas-1Pre_S102_L001_tumor.seg Pre-treatment\\\ US-WGS.seg")
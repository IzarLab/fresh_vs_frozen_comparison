gen_pos_gene_name_strip = read.table("/data/gen_pos_gene_name_strip.txt",header=F,sep="\t",quote=NULL)
gen_pos_gene_name_strip$V4 = "+"
gen_pos_gene_name_strip$V5 = "-"
write.table(gen_pos_gene_name_strip,"/data/gen_pos_gene_name_strip_add_cols.txt",col.names=F,row.names=F,sep="\t",quote=FALSE)

KAPA = read.table("/data/KAPA_HyperExome_hg19_capture_targets.bed",header=F,sep="\t",quote=NULL)
KAPA$V4 = "+"
KAPA$V5 = "-"
write.table(KAPA,"/data/KAPA_HyperExome_hg19_capture_targets_strip_add_cols.txt",col.names=F,row.names=F,sep="\t",quote=FALSE)

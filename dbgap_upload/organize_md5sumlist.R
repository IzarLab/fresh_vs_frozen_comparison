library(stringr)

md5sumlist = read.table("md5sumlist_edited.txt",sep=" ",header=F,quote=NULL)

seenfilenames = unlist(lapply(str_split(md5sumlist$V2,"/"), function(x) {x[length(x)]}))
md5sums = md5sumlist$V1

source("rename_md5sumlist_function.R")

filenamedf = data.frame(filename = seenfilenames, md5sum = md5sums)
filenamedf = filenamedf[grep("(_R1_|_R2_|HM)",filenamedf$filename),]

filenamedf$shortname = ""
for (i in 1:length(filenamedf$filename))
{
  afilename = filenamedf$filename[i]
  #if (sum(grep("fastq",afilename))!=0 && (sum(grep("_R1_",afilename))!=0 || sum(grep("_R2_",afilename))!=0))
  #{
    ashortname = rename_md5sumlist_function(afilename)
    if (ashortname=="")
    {
      print(afilename)
    }
    filenamedf$shortname[filenamedf$filename==afilename] = ashortname
  #}
}

filenamedf2 = data.frame(common_name = character(), R1_file = character(), R1_md5sum = character(), R2_file = character(), R2_md5sum = character())
prevfilename = ""
prevmd5sum = ""
prevshortname = ""
for (i in 1:length(filenamedf$filename)) {
  afilename = filenamedf$filename[i]
  anmd5sum = filenamedf$md5sum[i]
  ashortname = filenamedf$shortname[i]
  if (prevfilename=="")
  {
    prevfilename = afilename
    prevmd5sum = anmd5sum
    prevshortname = ashortname
  }
  else
  {
    if (str_replace(prevfilename,"_R1_","_R2_")==afilename || (prevfilename=="HM-baseline_TGACCA_R1.fastq.gz" && afilename=="HM-baseline_TGACCA_R2.fastq.gz"))
    {
      tempdf = data.frame(common_name = ashortname, R1_file = prevfilename, R1_md5sum = prevmd5sum, R2_file = afilename, R2_md5sum = anmd5sum)
      filenamedf2 = rbind(filenamedf2, tempdf)
      prevfilename = ""
      prevmd5sum = ""
      prevshortname = ""
    }
    else
    {
      tempdf = data.frame(common_name = prevshortname, R1_file = prevfilename, R1_md5sum = prevmd5sum, R2_file = "", R2_md5sum = "")
      filenamedf2 = rbind(filenamedf2, tempdf)
      prevfilename = afilename
      prevmd5sum = anmd5sum
      prevshortname = ashortname
    }
  }
}

columndf = read.table("dbGaP_sample_ID_column.txt",sep=",",header=F,quote=NULL)

printdf = filenamedf2[match(columndf$V1, filenamedf2$common_name),]

write.table(printdf, "organize_md5sumlist.csv", sep=",", row.names=F, col.names=T, quote=FALSE)
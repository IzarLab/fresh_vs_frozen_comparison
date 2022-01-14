library(rgl)
clonality_table = read.table("ribas_melanoma_cd8_clonality_over_time.csv",sep=",",header=T,quote=NULL)
clonality_names = names(clonality_table)
for (i in 4:length(clonality_names))
{
  clonality_table_temp = data.frame(DC1=clonality_table$DC1, DC2=clonality_table$DC2, DC3=clonality_table$DC3, clonality_over_time=clonality_table[[clonality_names[i]]])
  clonality_table_temp$color = "gray"
  clonality_table_temp$color[clonality_table_temp$clonality_over_time=="ribas_310_pre"] = "blue"
  clonality_table_temp$color[clonality_table_temp$clonality_over_time=="ribas_310_on"] = "green"
  clonality_table_temp$color[clonality_table_temp$clonality_over_time=="ribas_310_on_later"] = "red"
  clonality_table_temp$size = 0.00025
  clonality_table_temp$size[clonality_table_temp$clonality_over_time!=""] = 0.0015
  #clonality_table_temp = clonality_table_temp[1:1000,]
  rgl.viewpoint(theta=160, phi=15, fov=0)
  plot3d(x=clonality_table_temp$DC1, y=clonality_table_temp$DC3, z=clonality_table_temp$DC2, xlab="DC1", ylab="DC3", zlab="DC2", type="s", radius=clonality_table_temp$size, col=clonality_table_temp$color, lit = FALSE)
  legend3d("top", legend = clonality_names[i])
  #rgl.postscript("test3.pdf", fmt="pdf")
  rgl.snapshot(paste0("rglscript_output/",clonality_names[i],".png"), fmt="png")
}

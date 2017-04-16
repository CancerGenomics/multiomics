#OJO QUE QUEDARON LOS PATH CLAVADOS DE MI MAQUINA.
meth.platform<-read.table("D:\\desarrollo\\workspaces\\R\\multiomics\\examples\\methylation_X_mrnas\\cg-gene.txt", header = TRUE)
mrna<-read.table("D:\\desarrollo\\workspaces\\R\\multiomics\\examples\\methylation_X_mrnas\\mrnas.csv", header = TRUE)
meth<-read.table("D:\\desarrollo\\workspaces\\R\\multiomics\\examples\\methylation_X_mrnas\\meth.csv", header = TRUE)
methXMrnas(mrna, meth, meth.platform, output.path="D:\\desarrollo\\workspaces\\R\\multiomics\\examples\\methylation_X_mrnas\\", 
           output.file.name="methylationXMrna.csv",
           r.minimium=0.4, 
           pearsons.method = "pearson", 
           inc.progress = F)


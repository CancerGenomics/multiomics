# TODO: Add comment
# 
# Author: Matias
###############################################################################


#paramDir<- "D:\\desarrollo\\bioclassifier\\bioclassifier_R" # the location of unchanging files such as the function library and main program

#source("D:/desarrollo/workspaces/R/multiomics/src/3-do/Private/bioclassifier/subtypePrediction_distributed.R", echo=FALSE, encoding="Cp1252")
#source("D:/desarrollo/workspaces/R/multiomics/src/3-do/API/API_Pam_50.R", echo=FALSE, encoding="Cp1252")
#source("D:/desarrollo/workspaces/R/multiomics/src/3-do/Private/bioclassifier/subtypePrediction_functions.R", echo=FALSE, encoding="Cp1252")

input.folder<- "D:\\datosGenomica\\multiomics\\pam50"  
input.file.name <- "PAM50fromNormals&DCIS_normalizedBysample&genes.txt"
prefix.name.for.output.files<-"PAM50_median_centered.txt"
config.dir<-"D:\\desarrollo\\workspaces\\R\\multiomics\\resources\\bioclassifier"
compare.with.pam.50 (config.dir, input.folder, input.file.name, prefix.name.for.output.files)

#source("D:/desarrollo/workspaces/R/multiomics/src/3-do/Private/bioclassifier/subtypePrediction_distributed.R", echo=FALSE, encoding="Cp1252")
#source("D:/desarrollo/workspaces/R/multiomics/src/3-do/API/API_Pam_50.R", echo=FALSE, encoding="Cp1252")
#source("D:/desarrollo/workspaces/R/multiomics/src/3-do/Private/bioclassifier/subtypePrediction_functions.R", echo=FALSE, encoding="Cp1252")

input.folder<- "D:\\matias\\academia\\investigacion\\medicina personalizada\\1-analisis realizados\\2015-08-17---ratas de aldaz\\input"  
input.file.name <- "MPAstudy_pam50genes.csv"
prefix.name.for.output.files<-"MPAstudy_pam50genes_multiomics.txt"
config.dir<-"D:\\desarrollo\\workspaces\\R\\multiomics\\resources\\bioclassifier"
compare.with.pam.50 (config.dir, input.folder, input.file.name, prefix.name.for.output.files)

# TODO: Add comment
# 
# Author: Matias
###############################################################################


install.packages("XML")
install.packages("RCurl")
multimirFile <- "https://gitlab.com/cancergenomics/multiomics/raw/master/src/1-installation/dependencies/multiMiR_1.0.1.tar.gz"    
# for local instalnation, require cloning the git repo 
# multimirFile <- paste(getwd(),"/src/1-installation/dependencies/multiMiR_1.0.1.tar.gz", sep="")
install.packages(multimirFile, repos=NULL, type="source")
install.packages("GGally")
install.packages("Hmisc")
install.packages("ggplot2")

install.packages("shiny")
install.packages("DT")
install.packages("shinyBS")
install.packages("shinyjs")
install.packages("rclipboard")
install.packages("clipr")

install.packages("heatmap.plus")
source("http://bioconductor.org/biocLite.R")
biocLite("ctc")
biocLite("impute")
biocLite("survcomp")
biocLite("org.Hs.eg.db")
biocLite("devtools")
# if you're behind proxy see: http://stackoverflow.com/questions/16740256/error-installing-packages-from-github
biocLite("mtmorgan/xenar")
install.packages("reshape2")
source("http://bioconductor.org/biocLite.R") 
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
install.packages("WGCNA")

#MDV: 19/5/2018
install.packages("Matrix")

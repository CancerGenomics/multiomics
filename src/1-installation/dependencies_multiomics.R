# TODO: Add comment
# 
# Author: Matias
###############################################################################


install.packages("XML")
install.packages("RCurl")
install.packages("/home/hernan/Escritorio/multiMiR_1.0.1.tar.gz", repos=NULL, type="source")

install.packages("heatmap.plus")
source("http://bioconductor.org/biocLite.R")
biocLite("ctc")
biocLite("impute")
biocLite("survcomp")



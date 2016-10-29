# TODO: Add comment
# 
# Author: Matias
###############################################################################


install.packages("XML")
install.packages("RCurl")
multimirFile <- paste(getwd(),"/src/1-installation/dependencies/multiMiR_1.0.1.tar.gz", sep="")
install.packages(multimirFile, repos=NULL, type="source")
install.packages("GGally")
install.packages("Hmisc")

install.packages("heatmap.plus")
source("http://bioconductor.org/biocLite.R")
biocLite("ctc")
biocLite("impute")
biocLite("survcomp")
biocLite("org.Hs.eg.db")
biocLite("devtools")
biocLite("mtmorgan/xenar")


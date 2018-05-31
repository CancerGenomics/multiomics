#   cnv.path: it is the path of a file with the following format
#      -Row 1: It has the sample labels
#	   -Column1: It has the gene symbol (for example: A1BG, A2M)
#	   -The cells has got the copy number variations
#
#The row.names is null because the file can have more than once the same gene.
readCNVFile <- function(cnv.path, ncol.for.expression.id=1) {

  print("Reading cnv file...")
  cnv <- na.omit(read.table(cnv.path, header=TRUE,fill=TRUE, row.names = NULL, sep="\t"))
  cnv <-SortMatrixByColumnName(cnv, 1)
  
  return (cnv)
}


#   expression.file: it is the path of a file with the following format
#      -Row 1: It has the sample labels
#	   -Column1: It has the gene symbol (for example: A1BG, A2M)
#	   -The cells has got the expression level of each gene for each sample
readMrnaExpressionFile <- function(expression.file, ncol.for.expression.id=1) {
  
  print("Reading the mrna file...")
  expression <- na.omit(read.table(expression.file, header=TRUE,fill=TRUE, row.names = NULL, sep = "\t"))
  unique.names<-make.unique(as.character(expression[,1]),  sep = "_")
  rownames(expression)<-unique.names
  #expression<-expression[,2:ncol(expression)]
  print("Sorting the mrna data...")
  expression <-SortMatrixByColumnName(expression, 1)
  #rownames(expression) <- expression[,1]
  return (expression)
}

keepSameColumns <- function(matrix1, matrix2) {
  cols_to_keep <- intersect(colnames(matrix1)[-1],colnames(matrix2)[-1])
  matrix1 <- matrix1[,c(colnames(matrix1)[1],cols_to_keep), drop=FALSE]
  matrix2 <- matrix2[,c(colnames(matrix2)[1],cols_to_keep), drop=FALSE]
  return (list(matrix1,matrix2))
  
}



#Ejemplo de cohort: TCGA Uterine Carcinosarcoma (UCS)
#Ejemplo de dataset: TCGA/TCGA.UCS.sampleMap/RPPA
#Prefijo del xenahub: https://tcga.xenahubs.net/download/
getUrlFromTCGAXenaHub <- function(dataset){
  prefix<-"https://tcga.xenahubs.net/download/"
  url<-paste0(prefix, dataset, ".gz")
}


generateAllURLSFromXenaHub <- function(){
  cohorts.names<-cohorts(XenaHub(hosts = "https://tcga.xenahubs.net"))
  for (i in 1:length(cohorts.names)) {
    datasets<-datasets(XenaHub(cohorts=cohorts.names[i]))
    url.parts<-unlist(strsplit(datasets[1], '/'))
    url.sufix<-url.parts[2]
    if (length(url.parts)>=3){
      for (x in 3:length(url.parts)) {
        url.sufix<-paste(url.sufix, url.parts[x],sep="/")
      }
    }
    print(paste("https://tcga.xenahubs.net/download/",  url.sufix, sep=""))
  }
}


generateAllURLSFromXenaHub2 <- function(){
	cohorts.names<-cohorts(XenaHub(hosts = "https://tcga.xenahubs.net"))
	for (i in 1:length(cohorts.names)) {
		ds<-datasets(XenaHub(cohorts=cohorts.names[i]))
		for (j in 1:length(ds)) {
		  url <- getUrlFromTCGAXenaHub(ds[j])
		  print(url)
		}
	}
}


readMethylationFile <- function(meth.path, ncol.for.expression.id=1) {
  
  print("Reading methylation file...")
  meth <- na.omit(read.table(meth.path, header=TRUE,fill=TRUE, row.names = NULL, sep="\t"))
  meth <-SortMatrixByColumnName(meth, 1)
  return (meth)
}
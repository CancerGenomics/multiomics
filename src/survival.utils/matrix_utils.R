# TODO: Add comment
# 
# Author: Matias
###############################################################################

#It sorts the received matrix by columnName
#   x: The matrix
#   colsToExclude: you can exclude from the sorting, some columns. For example, if the first column is the row name
SortMatrixByColumnName <- function(x, colsToExclude=0){
	return (x[,union(c(1), order(colnames(x)[2:length(colnames(x))])+1)])	
}


#It removes rows containing at least on NA
#   dataframe: The matrix
removeRowsWithNA <- function(dataframe){
	for(i in 1:nrow(dataframe)){
		filter[i]<-!any(is.na(dataframe[i,]))
	}
	return (dataframe[filter,])
}

writeInFile <- function(matrix, output.file.path){
	write.table(matrix, output.file.path, sep="\t",row.names=FALSE)	
	
}

orderMatrixRows <- function(mat, decr = F, cols = NULL){

	return (mat[order(as.numeric(mat[,cols]),decreasing = FALSE, na.last = TRUE), ])

}

convertVectorToMatrix <- function(suspectedVector) {
  if(is.vector(suspectedVector)){
    suspectedVector  <- t(as.matrix(suspectedVector))
  }  
  return (suspectedVector)
}

############################################################################### 
# Author: Matias
###############################################################################

# Chequea si la matrix de expresi?n y la matriz de mirnas, tienen los mismos sample names y en el mismo orden.
checkSamplesFormIRNArnaCorrelation <- function(expressionMatrix, mirnaMatrix, ncolForExpressionId){
	firstIndexForExpression<-ncolForExpressionId+1
	samplesOfDataMatrix1 <- gsub("\\s","",toupper(colnames(expressionMatrix)))[firstIndexForExpression:ncol(expressionMatrix)]
	samplesOfDataMatrix2<- gsub("\\s","",toupper(colnames(mirnaMatrix)))[2:ncol(mirnaMatrix)]
	if (!(all(samplesOfDataMatrix1==samplesOfDataMatrix2))) stop("different samples")
}




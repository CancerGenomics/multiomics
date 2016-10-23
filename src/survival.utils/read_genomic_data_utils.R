#   cnv.path: it is the path of a file with the following format
#      -Row 1: It has the sample labels
#	   -Column1: It has the gene symbol (for example: A1BG, A2M)
#	   -The cells has got the copy number variations
#
#The row.names is null because the file can have more than once the same gene.
readCNVFile <- function(cnv.path, ncol.for.expression.id=1) {
  print("Reading the mrna file...")
  cnv <- na.omit(read.table(cnv.path, header=TRUE,fill=TRUE, row.names = NULL))
  print("Sorting the mrna data...")
  cnv <-SortMatrixByColumnName(cnv, 1)
  return (cnv)
}


#   expression.file: it is the path of a file with the following format
#      -Row 1: It has the sample labels
#	   -Column1: It has the gene symbol (for example: A1BG, A2M)
#	   -The cells has got the expression level of each gene for each sample
readMrnaExpressionFile <- function(expression.file, ncol.for.expression.id=1) {
  print("Reading the mrna file...")
  expression <- na.omit(read.table(expression.file, header=TRUE,fill=TRUE))
  print("Sorting the mrna data...")
  expression <-SortMatrixByColumnName(expression, 1)
  rownames(expression) <- expression[,1]
  return (expression)
}
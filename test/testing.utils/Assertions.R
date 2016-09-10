# TODO: Add comment
# 
# Author: Matias
###############################################################################

areEquals <-function(number1, number2){
	return (number1==number2)
}

assertEquals <- function(result, expected, testName){
	if (areEquals(result,expected)) print(paste(testName, "OK")) else {print(paste(testName, " failed. ", "Expected value: ", expected, ". Result: ", result))}
	
}


assertAllEquals <- function(result, expected, testName){
	if (all(result==expected)) print(paste(testName, "OK")) else {print(paste(testName, " failed. ", "Expected value: ", expected, ". Result: ", result))}
	
}


debug <-function(what.to.print, path){
	write.table(what.to.print, path,  row.names=FALSE, sep="\t")
}
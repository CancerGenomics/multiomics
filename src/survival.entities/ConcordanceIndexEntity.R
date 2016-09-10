##############################################################################
# It Represents the result of concordance index 
# Function to compute the concordance index for a risk prediction, i.e. the probability that, for a pair of randomly chosen comparable samples, the sample with the higher risk prediction will experience an event before the other sample or belongs to a higher binary class. 
# Author: Matias
###############################################################################

setClass("ConcordanceIndex",
		representation(concordance.index="numeric"))

getConcordanceIndex <-function(time.vector, event.vector, groups.vector){
	library("survcomp")
	cluster.index.result<-concordance.index(groups.vector, time.vector, event.vector, method="noether")
	result<-new("ConcordanceIndex", concordance.index=cluster.index.result$c.index)
	return (result)
}

concordnaceIndextoString <- function(concordance.index.object){
	elements<-c(paste("Concordance index: ", concordance.index.object@concordance.index))
	return (elements)
	
}

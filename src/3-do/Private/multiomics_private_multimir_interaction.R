# TODO: Add comment
# 
# Author: Matias
###############################################################################


getPredictedFromMulimir <- function(mirna){

	multimir = get.multimir(org="hsa", mirna=mirna, table="all",summary=TRUE, predicted.cutoff.type="p", predicted.cutoff=30)
	
	multimirPredicted<-multimir$predicted
	varsToKeep <- names(multimirPredicted) %in% c("database", "mature_mirna_id", "target_symbol", "score")
	multimirPredictedFiltered <- multimirPredicted[varsToKeep]
	
	multimirValidated<-multimir$validated
	varsToKeep <- names(multimirValidated) %in% c("database", "mature_mirna_id", "target_symbol", "pubmed_id")
	multimirValidatedFiltered <- multimirValidated[varsToKeep]
	
	if (is.null(multimirValidatedFiltered)) multimirValidatedFiltered<-data.frame(database= character(0), target_symbol = character(0), pubmed_id=character(0))
	
	
	if (is.null(multimirPredictedFiltered)) multimirPredictedFiltered<-data.frame(database= character(0), target_symbol = character(0), score=numeric(0))
	
	multimir<-merge(multimirPredictedFiltered, multimirValidatedFiltered, all.x=TRUE, all.y=TRUE)     
	
	return (multimir)
}


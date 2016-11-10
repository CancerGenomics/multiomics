# TODO: Add comment
# 
# Author: Matias
###############################################################################


getPredictedFromMulimir <- function(mirna){

  print("get.multimir")
  mirnaDBs <- "all"
  #mirnaDBs <- "diana_microt"
	multimir = get.multimir(org="hsa", mirna=mirna, table=mirnaDBs,summary=TRUE, 
	                        predicted.cutoff.type="p", predicted.cutoff=30)
	
	print("multimir predicted")
	multimirPredicted<-multimir$predicted
	varsToKeep <- names(multimirPredicted) %in% c("database", "mature_mirna_id", "target_symbol", "score")
	multimirPredictedFiltered <- multimirPredicted[varsToKeep]
	
	print("multimir validated")
	multimirValidated<-multimir$validated
	varsToKeep <- names(multimirValidated) %in% c("database", "mature_mirna_id", "target_symbol", "pubmed_id")
	multimirValidatedFiltered <- multimirValidated[varsToKeep]
	
	if (is.null(multimirValidatedFiltered)) {
	  print("no null validated")
	  multimirValidatedFiltered<-data.frame(database= character(0), 
	                                        target_symbol = character(0), 
	                                        pubmed_id=character(0))
	}
	
	
	if (is.null(multimirPredictedFiltered)) {
	  print("no null predicted")
	  multimirPredictedFiltered<-data.frame(database= character(0), 
	                                        target_symbol = character(0), 
	                                        score=numeric(0))
	}
	
	print("merging")
	multimir<-merge(multimirPredictedFiltered, multimirValidatedFiltered, all.x=TRUE, all.y=TRUE)     
	print("end merge")
	
	return (multimir)
}


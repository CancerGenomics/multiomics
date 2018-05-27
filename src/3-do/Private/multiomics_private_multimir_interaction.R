# TODO: Add comment
# 
# Author: Matias
###############################################################################


getPredictedFromMulimir <- function(mirna){

  print("get.multimir")
  mirnaDBs <- "all"
  #mirnaDBs <- "validated"
  #mirnaDBs <- "diana_microt"
	multimir = get_multimir(org="hsa", mirna=mirna, table=mirnaDBs,summary=FALSE, 
	                        predicted.cutoff.type="p", predicted.cutoff=30)

		#write.table(multimir,'c:\\temp\\all' )
  if (length(multimir)==0) return (NULL)

	if (mirnaDBs=="all") {
	  multimir<-multimir@data
	} else {
	  print("NO ES ALL")
	  #multimirPredicted<-get(mirnaDBs, multimir)
	  multimir<-get(mirnaDBs, multimir)
	}

	
	
	#varsToKeep <- names(multimirPredicted) %in% c("database", "mature_mirna_id", "target_symbol", "score")
	#	varsToKeep <- c("database", "mature_mirna_id", "target_symbol", "score")
	#multimirPredictedFiltered <- multimirPredicted[,varsToKeep]
	#	multimir<-get(mirnaDBs, multimir)
	
	
	#print("multimir validated")
	#multimirValidated<-multimir$validated
	#varsToKeep <- names(multimirValidated) %in% c("database", "mature_mirna_id", "target_symbol", "pubmed_id")
	#multimirValidatedFiltered <- multimirValidated[varsToKeep]
	
	#	if (is.null(multimirValidatedFiltered)) {
	#	  print("no null validated")
	#	  multimirValidatedFiltered<-data.frame(database= character(0), target_symbol = character(0), pubmed_id=character(0))
	#	}
	
	
	#	if (is.null(multimirPredictedFiltered)) {
	#	  print("no null predicted")
	#	  multimirPredictedFiltered<-data.frame(database= character(0), target_symbol = character(0),score=numeric(0))
	#	}
	
	#print("merging")
	#multimir<-merge(multimirPredictedFiltered, multimirValidatedFiltered, all.x=TRUE, all.y=TRUE)     
	#print("end merge")
	
	
	return (multimir)
}


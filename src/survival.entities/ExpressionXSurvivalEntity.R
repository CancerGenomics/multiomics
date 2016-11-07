# TODO: Add comment
# 
# Author: Matias
###############################################################################

setClass("ExpressionXSurvival",
		representation(expression="data.frame", time="numeric", event="numeric", sample.ids="character", nsamples="numeric", ngenes="numeric"))

#It reads expressionXSurvival from a file with the following format
#      -Row 1: Sample labels
#	   -Row 2: Surival time (it could be for example relapse free survival time or overall suvival)
#	   -Row 3: Survival event. 0 means the event did not occur. 1 the event occurred.
#	   -From Row 4: expression data 
#	   -Column 1: Gene symbol (ejemplo: A1BG, A2M) except in row 1, 2 and 3 that are the labels for the rows descripted below.
readExpressionXSurvivalFromTypicalFile <- function(expression.file.path){

#ptm <- proc.time()
	the.original.data.frame<-read.table(expression.file.path, header=FALSE)

#print(paste("reading table: ", (proc.time() - ptm)["elapsed"]))	

#ptm <- proc.time()
	the.original.data.frame<-the.original.data.frame[,colSums(is.na(the.original.data.frame)) == 0]
#print(paste("Eliminando los nas: ", (proc.time() - ptm)["elapsed"]))
	
#ptm <- proc.time()
	event.vector<-the.original.data.frame[3,2:ncol(the.original.data.frame)]
	event.vector<-as.numeric(as.character(unlist(event.vector[1,])))
#print(paste("Creando el vector de eventos: ", (proc.time() - ptm)["elapsed"]))

#ptm <- proc.time()
	time.vector<-the.original.data.frame[2,2:ncol(the.original.data.frame)]
	time.vector<-as.numeric(as.character(unlist(time.vector[1,])))		
#print(paste("Creando el vector de time: ", (proc.time() - ptm)["elapsed"]))
	
#ptm <- proc.time()
	sample.ids.vector<-the.original.data.frame[1,2:ncol(the.original.data.frame)]
	sample.ids.vector<-as.character(unlist(sample.ids.vector[1,]))
#print(paste("Creando el vector de samplesIds: ", (proc.time() - ptm)["elapsed"]))	

	the.ngenes<-nrow(the.original.data.frame)-3
	the.nsamples<-ncol(the.original.data.frame)-1

#ptm <- proc.time()	
	expression.data.frame <-the.original.data.frame[4:nrow(the.original.data.frame),2:ncol(the.original.data.frame)]
	row.names(expression.data.frame)<-as.character(unlist(the.original.data.frame[4:nrow(the.original.data.frame),1]))
	colnames(expression.data.frame)<-as.character(unlist(the.original.data.frame[1, 2:ncol(the.original.data.frame)]))
#print(paste("Creando el vector de expresion: ", (proc.time() - ptm)["elapsed"]))
	
	rm(the.original.data.frame)

#ptm <- proc.time()
	expressionXSurvival<-new("ExpressionXSurvival", time=time.vector, event=event.vector, expression=expression.data.frame, sample.ids=sample.ids.vector, nsamples=the.nsamples, ngenes=the.ngenes) 
#print(paste("Creando el objeto expressionXSurvival: ", (proc.time() - ptm)["elapsed"]))
	
	return (expressionXSurvival)
}


expressionVectorByPosition <- function (expression.x.survival.object, gene.position){
	#expression<-expression.x.survival@expression[i,]
	#expression<-gsub(",", ".", as.character(unlist(expression[1,])))
	#expression<-as.numeric(expression)
	expression.of.partner.gene<-expression.x.survival.object@expression[gene.position,]
	expression.of.partner.gene<-as.numeric(as.character((t(expression.of.partner.gene))[,1]))
	
}


expressionVectorByGeneName <- function (expression.x.survival.object, gene.name){
	#expression<-expression.x.survival@expression[gene.name,]
	#expression<-gsub(",", ".", as.character(unlist(expression[1,])))
	#expression<-as.numeric(expression)
	
	expression.of.partner.gene<-expression.x.survival.object@expression[gene.name,]
	expression.of.partner.gene<-as.numeric(as.character((t(expression.of.partner.gene))[,1]))
	
}

geneNameInPosition <- function(expression.x.survival, i){
	return (row.names(expression.x.survival@expression)[i])
}

getNumericMatrixWithTheExpressionOfPairOfGenesFaster<-function (expression.x.survival.object, expression.of.candidate.gene, gene.partner.name.position){
	
#ptm<-proc.time()
	expression.of.partner.gene<-expression.x.survival.object@expression[gene.partner.name.position,]
	expression.of.partner.gene<-as.numeric(as.character((t(expression.of.partner.gene))[,1]))
#print(paste("Conversion: ", (proc.time() - ptm)["elapsed"]))


	
	
	##########funca pero lento
	#expression.of.partner.gene.1 <- expression.of.partner.gene[1,]
	#expression.of.partner.gene.1.unlisted<-unlist(expression.of.partner.gene.1)
	#expression.of.partner.gene.1.unlisted.as.character<-as.character(expression.of.partner.gene.1.unlisted)
	#expression.of.partner.gene<-gsub(",", ".", expression.of.partner.gene.1.unlisted.as.character)
	#expression.of.partner.gene<-as.numeric(expression.of.partner.gene)
	
	
	

	mat <- matrix(, 2, ncol = expression.x.survival.object@nsamples)
	
	mat[1,] <- expression.of.candidate.gene

	mat[2,] <- expression.of.partner.gene

	return (mat)
}



getNumericMatrixWithTheExpressionOfPairOfGenes<-function (expression.x.survival.object, gene.candidate.name, gene.partner.name){
	
	
	expression.of.candidate.gene<-expression.x.survival.object@expression[gene.candidate.name,]
	expression.of.candidate.gene<-gsub(",", ".", as.character(unlist(expression.of.candidate.gene[1,])))
	expression.of.candidate.gene<-as.numeric(expression.of.candidate.gene)
	
	expression.of.partner.gene<-expression.x.survival.object@expression[gene.partner.name,]
	expression.of.partner.gene<-gsub(",", ".", as.character(unlist(expression.of.partner.gene[1,])))
	expression.of.partner.gene<-as.numeric(expression.of.partner.gene)

	mat <- matrix(, 2, ncol = expression.x.survival.object@nsamples)
	mat[1,] <- expression.of.candidate.gene
	mat[2,] <- expression.of.partner.gene
	
	return (mat)
}



filter.genes <- function(expression.x.survival.object, genes.to.keep){
	row.names(expression.x.survival.object@expression)
	expression.x.survival.object@expression<-expression.x.survival.object@expression[genes.to.keep,]
	expression.x.survival.object@ngenes<-nrow(expression.x.survival.object@expression)
	return (expression.x.survival.object)
	
}


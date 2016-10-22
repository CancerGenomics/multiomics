###############################################################################
# Author: Matias
###############################################################################


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

# mirna.file: it is the path of a file with the following format
#  -Row 1: It has the sample labels
#  -Column 1: Mature mirna ID (for example: hsa-miR-22-3p). Note that R is un upper case. 
#             This R in uppercase is the way to identify that it is a mature-mirna and not a pre-mirna. 
#             The pre-mirna for this mature-mirna is hsa-mir-22. 
#             Multimir doesnt work with premirna because the databases it queries, doesnt work with pre-mirna. 
#             So, if you pass the pre-mirna (for example hsa-miR-22) multimir does not return anything. 
#             So, it is importante that this file contains mature mirna ids. 
#             With the accession it also works (MIMAT0004982). 
#             These IDs could be find in http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0005761 
#  -Cells has the expression levels of each mirna for each sample. 
readMirnaExpressionFile <- function(mirna.file, ncol.for.expression.id=1) {
  print("Reading the mirna file...")
  mirna <- na.omit(read.table(mirna.file, header=TRUE, fill=TRUE, sep="\t"))  
  print("Sorting the mirna data...")
  mirna <-SortMatrixByColumnName(mirna, 1)
  return (mirna)
}  



#Input
#  expression: matrix obtained with readMirnaExpressionFile function
#  mirna: matrix obtained with readMrnaExpressionFile function
#  output.path: The folder for output file. 
#  ncol.for.expression.id: It defines hoy many ID columns has got the file. 
#                          First one will be take into account to identify the gene or mirna; 
#                          others will be discarded. The default value is 1.
#  r.minimium: The minimium value for pearson coefficient for considering the 
#              correlation between a gene and a mirna.
#  pearsons.method: Pearson's product moment correlation coefficient.
#                   Possible values c("pearson", "kendall", "spearman"). Default "pearsons"
#  inc.progress: when running from shiny, increments progress bar during calculation.
#
#  Output
#	It generates a file with 4 columns: Gen - Mirna - Pearson correlation - p-value 
#	Pearson correlation indicates if there is or not a correlation between mrna 
#   expression and mirna expression taking into account all samples. 
#	If there is a correlation (both expression has similar values in all samples 
#   (positive correlation) or so different values in all samples (negative values)). 
#	But how we define if there are enough similar or enough different to say 
#   there are or there are not (positive or negative) correlation? Using statistics. 
#   In particular the Pearson method that works exactly for this case: compare 
#   two vector of values of the same size to check if there is correlation with 
#   statistical significance.	
#	It will be considered that both vectors correlates if (r>0.7 with p-value<0.05).
#   In this case it would be indicate that this mirna participates in the 
#   regulation of this gene.
#	In general, the correlation will be negative because mirnas inhibes the 
#   arnM translation (yes translation).
#   
CalculateCorrelationsMirnaMrna <- function(expression, mirna, output.path="~/", 
                                           output.file.name="inputStep2-matureMirnaXmrna.csv",
                                           r.minimium=0.7, 
										   pearsons.method = "pearson", 
                                           inc.progress = F){
	ptm <- proc.time()
	total.rows=nrow(expression)*nrow(mirna)
	print(paste("Running pipeline with", r.minimium, 
				"threshold and pearson's method:", pearsons.method, sep=" "))
	
	# The result matix is created
	res <- matrix(nrow=total.rows,ncol=4)
	colnames(res)<-(c("Gene symbol","Marure mirna id", "Mirna-Mrna correlation", "p_value of correlation"))

	actual<-1
	actual.n.correlated<-1
	print("Start process!")
	for (i in 1:nrow(mirna)) {
		actual.mirna<-mirna[i,1]
		mirna.para.ese.gen<-mirna[i,2:ncol(mirna)]
		for (j in 1:nrow(expression)) {
		  actual<-actual+1
			actual.gen<-expression[j,1]
			expression.para.ese.gen<-expression[j,2:ncol(expression)]
			if ((actual)%%500==0)print(paste("analised ", actual, " from ", total.rows))
			if ((actual)%%1000==0) {
			  elapsedTime <- (proc.time() - ptm)[3]
			  print(paste(
			    "elapsed time: (seconds)", format2Print(elapsedTime), 
			    " - (minutes)", format2Print(elapsedTime/60), 
			    " - (hours)", format2Print(elapsedTime/60/60)
			  ))
			  remainingTime <- ((total.rows*elapsedTime)/actual) - elapsedTime
			  print(paste("estimated remaining time (seconds)", format2Print(remainingTime),
			              " - (minutes)", format2Print(remainingTime/60), 
			              " - (hours)", format2Print(remainingTime/60/60)
			  ))
			}
			resultado.pearson<-cor.test(as.numeric(expression.para.ese.gen),
					                    as.numeric(mirna.para.ese.gen), 
										method = pearsons.method)
			if (!is.na(abs(resultado.pearson$estimate))) {
				if (abs(resultado.pearson$estimate) > r.minimium) {
				  newValue<-c(as.character(actual.gen), as.character(actual.mirna), 
							  resultado.pearson$estimate, resultado.pearson$p.value)
				  res[actual.n.correlated,1:4] <- newValue
				  actual.n.correlated<-actual.n.correlated+1
				}
			}
		}
		
		if(inc.progress) {
			incProgress(1/nrow(mirna));
		}	
	}
	# deleting useless and unused rows
	res <- res[c(1:actual.n.correlated-1),c(1:4)]

	#if (!(folder.exists(output.path))) {dir.create(output.path)}
	file.path<-paste(output.path, output.file.name, sep="")
	write.table(res, file.path, sep="\t", row.names=FALSE, 
			    col.names=TRUE, quote=FALSE)
	print(proc.time() - ptm)
	
	return (res)
}


#Input
#	genes.x.mirnas.path: es el path de un file que contiene los pares mrna-mirna con alta correlaci?n (en general deber?a ser negativa porque el mirna y el mrna suelen correlacionar negativamente. Cuando hay mirna, se suele inhibir el mrna). El file tiene las siguientes caracter?sticas
#		-Tiene 4 columnas: Gen_symbol	mature_mirna_id	Mirna_Mrna_Correlation	p_value_Of_Mirna_Mrna_Correlation.
#		-Los ?ltimos dos campos son la correlaci?n de expresi?n entre el gene y el mirna
#		-Ejmplo: ABCA2  hsa-miR-1292-5p	-0.701082209957814	2.72363288834111e-06
#	predicted.cut.off=30: A correlation registered in genes.x.mirnas.path will be considered in the final result if it is in the best "predcited.cut.off" % of genes that correlates for a mirna. 
#
#
#Output
# 	It returns a table with the 30% (configurable using predicted.cut.off parameter) mirna X mrna correlations which has got better predicted and/or validated score in mrna databases. This step also adds predicition and validation information for those correlations. The other 70% will be discarded.
#	Each row in the result represents a predicted or validated correlated in the database. So the same pair Gen_symbol - mature_mirna_id can have more than one row for registering the predicction/validation information in a particular database. Obviously if it is in the final result this correlation in this database is in the best of the correlations for this mirna.   
#		Gen_symbol: gene sybol identifying ARNm
#	    mature_mirna_id
#		Mirna_Mrna_Correlation: Pearson value for this correlation in a particular database
#		p_value_Of_Mirna_Mrna_Correlation: p value
#		mirna_database: mirna database. Possible values are: diana_microt, elmmo,  microcosm,  miranda,  mirdb,  pictar,  pita,  targetscan ,  mirecords,  mirtarbase,  tarbase, mir2disease,  pharmaco_mir,  phenomir
#		database_predicted_score: it is an socre for the prediction but take into account that it is not standarized. So, each database uses differente ranges.
#		validation_pubmed_id: it will be completed just if this correlation was validated (that is that was experimentally validated and was described in a paper).
keepBestGeneXMirnaAccordingCorrelationAndAddMirnaDbInfo <- function(genes.x.mirnas.path, output.path, output.file="inputStep3-mirnaXmRNAWithPredictedAndValidated.csv", predicted.cut.off=30){
	genes.x.mirnas<-read.table(genes.x.mirnas.path, header=TRUE)
	
	mirnas<-genes.x.mirnas[,2]
	mirnas<-unique(as.character(mirnas))
	
	genes<-genes.x.mirnas[,1]
	genes<-unique(as.character(genes))
	
	result <- data.frame(Gene_Symbol=character(0), mature_mirna_id=character(0), Mirna_Mrna_Correlation=numeric(0), p_value_Of_Mirna_Mrna_Correlation=numeric(0),Database=character(0),Database_Predicted_Score=numeric(0),pubMedID=character(0))  
	for (i in 1:length(mirnas)) {
		
		print(paste("mirna",i, "/", length(mirnas), ": ", mirnas[i]), sep="")
		multimir<- getPredictedFromMulimir(mirnas[i])
	
		if (!is.null(multimir) && nrow(multimir)>0)
		{resultTemp<-merge(genes.x.mirnas, multimir,  by.x=c(colnames(genes.x.mirnas)[1],colnames(genes.x.mirnas)[2]),by.y=c("target_symbol", "mature_mirna_id"))
			result<-rbind(result, resultTemp)
		}
	}
	colnames(result)<-c("Gen_symbol","mature_mirna_id","Mirna_Mrna_Correlation","p_value_Of_Mirna_Mrna_Correlation","mirna_database","database_predicted_score", "validation_pubmed_id")
	csvOutputFile<-paste(output.path, output.file, sep="")
	write.table(result, csvOutputFile, sep="\t",row.names=FALSE)
}






#Input
#	genes.x.mirnas.path: It is the file path es el path with mrna-mirna correlations enriched with predicted and validated data get from mirna databases. The same mirna-mrna could be more than once if it was predicted by multiple databases and all prediction were in the best 30% correlation score  
#	The file must have the following fields
#			Gen_Symbol: Gene Symbol, for example SAMD4A
#			mirnaId: It must be the mature mirna id (hsa-miR-1292-5p) or the accession (MIMAT0005943)
#			Mirna_Mrna_Correlation: Statistic of the correlation between this mrna and mirna according to the input files (it would be Pearson)
#			p_value_Of_Mirna_Mrna_Correlation: Pearson p_value
#			mirna_database
#			database_predicted_score: Take into account that this prediction score is nos standarized. Each database has got different scoring.
#			validation_pubmed_id: If this mirna x mrna was validated this field have got the pubmed id of the paper describing the value . 
#
#
#	output.path: folder in which the result will be stored.
#
#Output
#   It generates a file with the same columns but collapsing all rows referencing the same pair (mirna, mrna) in one row. The databases will be just one column with all databases separated by comma. 
#	Some considerations
#		In the databases column will be all the databases that was in the input file, separated by comma.
#		Database_predicted_score and validation_pubmed_id are also colapsed and they keep the same order than databases. So it is easy to match the score with the corresponding database and the validation_pubmed_id with the corresponding database.
#		The value "NO" in the field validation_pubmed_id means that this database hasnt got validated the pair (mirna, mrna)
#		Sometimes for the same pair (minra, mrna) the score is different for the same database. But this is the way, multimir returns the results. The difference is in the ensemblId (Check)
#
#
#Example
#   Input
#		MIR1292	hsa-miR-1292-5p	MIMAT0005943	ESRRG	-0.747206480568877	2,52E+07	NO_PUB	diana_microt
#		MIR1292	hsa-miR-1292-5p	MIMAT0005943	ESRRG	-0.747206480568877	2,52E+07	NO_PUB	pita
#   will be colapsed into 
#		MIR1292	hsa-miR-1292-5p	MIMAT0005943	ESRRG	-0.747206480568877	2,52E+07	NO_PUB	diana_microt, pita
#
ColapseMirnaXMrna <- function(mirna.X.mRNA.with.predicted.and.validated.path, output.path){
	mirnaXmrna<-read.table(mirna.X.mRNA.with.predicted.and.validated.path, header=T)

	#It keeps just the ones correlating negatively
	mirnaXmrnaNegative=mirnaXmrna[mirnaXmrna$Mirna_Mrna_Correlation<0,]

	#It converts NA in the validation_pubmed_id into NO
	mirnaXmrnaNegative$validation_pubmed_id <- as.character(mirnaXmrnaNegative$validation_pubmed_id)
	mirnaXmrnaNegative$mirna_database <- as.character(mirnaXmrnaNegative$mirna_database)
	mirnaXmrnaNegative$validation_pubmed_id[is.na(mirnaXmrnaNegative$validation_pubmed_id)] <- "NO"
	
	#It collapses all the rows with the same mirna and mrna.
	result <- aggregate(cbind(mirna_database,database_predicted_score, validation_pubmed_id)~Gen_symbol+mature_mirna_id+Mirna_Mrna_Correlation+p_value_Of_Mirna_Mrna_Correlation,paste,collapse=",",data=mirnaXmrnaNegative, na.action=na.pass)
	    
	
	#It ordes by mrna and then by mirna
	result  <- result[order(result$Gen_symbol, result$mature_mirna_id),] 

	#Write output
	csvOutputFile<-paste(output.path, "\\pipelineOutput-mirnaXmRNAWithPredictedAndValidatedColapsed.csv", sep="")
	write.table(result, csvOutputFile, sep="\t",row.names=FALSE,quote=FALSE)
}

# Formats a number with 2 decimal places
format2Print <- function(number){
  return (format(round(number, 2), nsmall = 2))
}
# TODO: Add comment
# 
# Author: Matias
###############################################################################




#It gets as input the mrna file (in tcga format), the clincal data file (in tcga format, that is clinical data as columns and samples as rows), the survival column name and the event column name and it generates a file combining expression data and clincal data.
#Some considerations
#	The samples that has not got survival and/or event will be discarded; it will not be in the final result.
#	The difference with the function joinExpressionAndAllClinicalDataFromTCGA is that this function extracts just the survival and event column. The joinExpressionAndAllClinicalDataFromTCGA function writes all the attributes.
#
#
#Input
#   expression.file.path: It is the file path of the expression data file in TCGA format. It must respect the following format.
#      -Row 1: Sample labels
#	   -From Row 2: Expression data
#	   -Column 1: Gene symbol (ejemplo: A1BG, A2M) except in row 1, 2 and 3 that are the labels for the rows descripted below.
#	   -Example: 
#				sample	TCGA-A8-A092-01	TCGA-A7-
#				ARHGEF10L	8.8784	9.9980	
#				HIF3A	0.3722	6.1594	2.4062	
#				RNF17	0.0000	4.1813	0.0000	
#		
#   clinical.file.path: It is the file in TCGA format. It must respect the following format.
#      -Row 1: Clinical data column name
#	   -Column 1: Sample labels (it should be the same than the row 1 of expression.file
#	   -Column 1: The cells will have the values of each sample in each clinical data.
#	   -Example: 
#				sampleID	AJCC_Stage_nature2012	OS_Time_nature2012	OS_event_nature2012
#				TCGA-A1-A0SB-01	Stage I	259	0
#				TCGA-A1-A0SD-01	Stage IIA	437	0
#				TCGA-A1-A0SE-01	Stage I	1320	0
#				TCGA-A1-A0SF-01	Stage IIA	1463	0
#				TCGA-A1-A0SG-01	Stage IIB	433	0
#
#   clinical.survival.column.name: The name of the column that have the survival data.
#
#	clinical.event.column.name: The name of the column that have the event data.
#
#	output.path: It is the folder in which the outfile will be stored (ExpressionWithSurvival.csv)
#
#
#Output
#   It generates a file with the following format
#      -Row 1: Sample labels
#	   -Row 2: Surival time (it could be for example relapse free survival time or overall suvival)
#	   -Row 3: Survival event. 0 means the event did not occur. 1 the event occurred.
#	   -From Row 4: expression data 
#	   -Column 1: Gene symbol (ejemplo: A1BG, A2M) except in row 1, 2 and 3 that are the labels for the rows descripted below.
#	All the samples which has NA for the columns event or survival will be eliminated and they will not be i
#   
generateExpressionAndSurvivalDataFromTCGA <- function(expression.file.path, clinical.file.path, clinical.survival.column.name, clinical.event.column.name, output.path){
  
  #It reads expression file. When reading the "-" are changed by "." that is why it is changed again to "-".
  expression <- (read.table(expression.file.path, sep="\t", header=TRUE, na.strings=c("", "NA"),stringsAsFactors=FALSE))
  colnames(expression)<-gsub("\\.", "-", colnames(expression))
  rownames(expression) <- expression[,1]
  
  #It reads the clincal data file, it transposse because samples are in rows in this file, 
  clinical <- (read.table(clinical.file.path, sep="\t", header=TRUE, na.strings=c("", "NA")))
  rownames(clinical) <- clinical[,1]
  timeAndEvent<-clinical[,c(colnames(clinical)[1],clinical.survival.column.name, clinical.event.column.name)]
  timeAndEvent<-t(timeAndEvent)
  timeAndEvent<-cbind(c("sampleID", clinical.survival.column.name,clinical.event.column.name),timeAndEvent )
  colnames(timeAndEvent) <- timeAndEvent[1,]
  
  #It gets the samples in common between both files
  samplesInCommon <- intersect(colnames(expression), colnames(timeAndEvent))
  
  #It filters the expressionMatrix for getting just the samples that are in both in expression file and in clincal file
  expressionThatHasCorrespondingClinicalDataInTheOtherFile<-expression[,c(colnames(expression)[1], samplesInCommon)]
  expressionThatHasCorrespondingClinicalDataInTheOtherFile<-SortMatrixByColumnName(expressionThatHasCorrespondingClinicalDataInTheOtherFile, 1)
  
  #It filters the clincalDataMatrix for getting just the samples that are in both in expression file and in clincal file
  clinicalThatHasCorrespondingexpressionDataInTheOtherFile<-timeAndEvent[,c(colnames(timeAndEvent)[1], samplesInCommon)]
  clinicalThatHasCorrespondingexpressionDataInTheOtherFile <-SortMatrixByColumnName(clinicalThatHasCorrespondingexpressionDataInTheOtherFile, 1)
  
  #It gets the sample names that will be in the final result.
  sampleNames<-data.frame(clinicalThatHasCorrespondingexpressionDataInTheOtherFile[1,], stringsAsFactors = FALSE)
  sampleNames<-t(sampleNames)
  sampleNames<-clinicalThatHasCorrespondingexpressionDataInTheOtherFile[1,]
  
  #It gets the survival data that will be in the final result.
  survivalData<-data.frame(clinicalThatHasCorrespondingexpressionDataInTheOtherFile[2:3,],stringsAsFactors = FALSE)
  
  #It gets the expression data that will be in the final result.
  expressionData<-data.frame(expressionThatHasCorrespondingClinicalDataInTheOtherFile[1:nrow(expressionThatHasCorrespondingClinicalDataInTheOtherFile),], stringsAsFactors = FALSE)
  
  #It sets the same colnames for allow the rbind
  colnames(expressionData)<-colnames(clinicalThatHasCorrespondingexpressionDataInTheOtherFile)
  colnames(survivalData)<-colnames(clinicalThatHasCorrespondingexpressionDataInTheOtherFile)
  
  #It generates the resutl
  result<-rbind(sampleNames,survivalData,expressionData)
  result<-result[,colSums(is.na(result)) == 0]
  
  #It writes the result
  write.table(result, output.path, col.names=FALSE, quote=FALSE, row.names=FALSE, sep="\t")
}







#It gets as input the mrna file (in tcga format), the clincal data file (in tcga format, that is clinical data as columns and samples as rows), expression data and all clincal data.
#Some considerations
#	The samples that are not in both files will be discarded 
#	
#
#Input
#   expression.file.path: It is the file path of the expression data file in TCGA format. It must respect the following format.
#      -Row 1: Sample labels
#	   -From Row 2: Expression data
#	   -Column 1: Gene symbol (ejemplo: A1BG, A2M) except in row 1, 2 and 3 that are the labels for the rows descripted below.
#	   -Example: 
#				sample	TCGA-A8-A092-01	TCGA-A7-
#				ARHGEF10L	8.8784	9.9980	
#				HIF3A	0.3722	6.1594	2.4062	
#				RNF17	0.0000	4.1813	0.0000	
#		
#   clinical.file.path: It is the file in TCGA format. It must respect the following format.
#      -Row 1: Clinical data column name
#	   -Column 1: Sample labels (it should be the same than the row 1 of expression.file
#	   -Column 1: The cells will have the values of each sample in each clinical data.
#	   -Example: 
#				sampleID	AJCC_Stage_nature2012	OS_Time_nature2012	OS_event_nature2012
#				TCGA-A1-A0SB-01	Stage I	259	0
#				TCGA-A1-A0SD-01	Stage IIA	437	0
#				TCGA-A1-A0SE-01	Stage I	1320	0
#				TCGA-A1-A0SF-01	Stage IIA	1463	0
#				TCGA-A1-A0SG-01	Stage IIB	433	0
#
#	output.path: It is the folder in which the outfile will be stored (ExpressionWithSurvival.csv)
#
#
#Output
#   It generates a file with the following format
#      -Row 1: Sample labels
#	   -Row2 - Row n: clinical data (assuming there are n-1 clincal attributes) 
#	   -From Row n+1: expression data 
#	   -Column 1: Gene symbol (ejemplo: A1BG, A2M) except in row 1, 2 and 3 that are the labels for the rows descripted below.
#   
joinExpressionAndAllClinicalDataFromTCGA <- function(expression.file.path, clinical.file.path, output.path){
  
  #It reads expression file. When reading the "-" are changed by "." that is why it is changed again to "-".
  expression <- (read.table(expression.file.path, sep="\t", header=TRUE, na.strings=c("", "NA"),stringsAsFactors=FALSE))
  colnames(expression)<-gsub("\\.", "-", colnames(expression))
  rownames(expression) <- expression[,1]
  
  #It reads the clincal data file, it transposse because samples are in rows in this file, 
  clinical <- (read.table(clinical.file.path, sep="\t", header=TRUE, na.strings=c("", "NA")))
  clinical<-t(clinical)
  clinical<-cbind(rownames(clinical),clinical)
  colnames(clinical)<-gsub("\\.", "-", clinical[1,])
  
  #It gets the samples in common between both files
  samplesInCommon <- intersect(colnames(expression), colnames(clinical))
  
  #It filters the expressionMatrix for getting just the samples that are in both in expression file and in clincal file
  expressionThatHasCorrespondingClinicalDataInTheOtherFile<-expression[,c(colnames(expression)[1], samplesInCommon)]
  expressionThatHasCorrespondingClinicalDataInTheOtherFile<-SortMatrixByColumnName(expressionThatHasCorrespondingClinicalDataInTheOtherFile, 1)
  
  #It filters the clincalDataMatrix for getting just the samples that are in both in expression file and in clincal file
  clinicalThatHasCorrespondingexpressionDataInTheOtherFile<-clinical[,c(colnames(clinical)[1], samplesInCommon)]
  clinicalThatHasCorrespondingexpressionDataInTheOtherFile <-SortMatrixByColumnName(clinicalThatHasCorrespondingexpressionDataInTheOtherFile, 1)
  
  #It gets the sample names that will be in the final result.
  sampleNames<-data.frame(clinicalThatHasCorrespondingexpressionDataInTheOtherFile[1,], stringsAsFactors = FALSE)
  sampleNames<-t(sampleNames)
  sampleNames<-clinicalThatHasCorrespondingexpressionDataInTheOtherFile[1,]
  
  #It gets the survival data that will be in the final result.
  #survivalData<-data.frame(clinicalThatHasCorrespondingexpressionDataInTheOtherFile[2:3,],stringsAsFactors = FALSE)
  survivalData<-data.frame(clinicalThatHasCorrespondingexpressionDataInTheOtherFile,stringsAsFactors = FALSE)
  
  #It gets the expression data that will be in the final result.
  expressionData<-data.frame(expressionThatHasCorrespondingClinicalDataInTheOtherFile[1:nrow(expressionThatHasCorrespondingClinicalDataInTheOtherFile),], stringsAsFactors = FALSE)
  
  #It sets the same colnames for allow the rbind
  colnames(expressionData)<-colnames(clinicalThatHasCorrespondingexpressionDataInTheOtherFile)
  colnames(survivalData)<-colnames(clinicalThatHasCorrespondingexpressionDataInTheOtherFile)
  
  #It generates the resutl
  result<-rbind(survivalData,expressionData)
  #result<-result[,colSums(is.na(result)) == 0]
  
  #It writes the result
  write.table(result, output.path, col.names=FALSE, quote=FALSE, row.names=FALSE, sep="\t")
}
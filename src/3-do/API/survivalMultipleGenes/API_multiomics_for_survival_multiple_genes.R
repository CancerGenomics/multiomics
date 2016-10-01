###############################################################################
# Author: Matias
###############################################################################

#Source
#source("D:/desarrollo/workspaces/R/multiomics/src/survival.entities/CoxphEntity.R", echo=FALSE, encoding="Cp1252")


######################################API#############################################
######################################API#############################################
######################################API#############################################
#Some considerations on the results
#	-Correlation between expression and survival data
#		-coxph and wald: It measures the correlations between the expressionData and the survivalData without creating gorups.
#
#	-Correlation between assigned group to sample and survival data (the expression data is not used by this tests, but in general you will use multiomics.G1PairUpG2PairDown function that assign groups depending on the expression)
#		-survDiff G-Rho p-value: It is the p-value of a test that measures the difference between the kaplanMeier curves. So if the p-value is lower than 0.05 there is saome statistics significance on this test (i.e for this dataset the kaplan-Meier curves has got some differences, so in this dataset this pair of genes divides the samples of this datasets into groups that correlates with survival data)
#		-Concordance.index: It is another index to measure the correlation between assigned groups and survival data .  







#First, it will find the best (number.of.genes.to.keep.during.training=10) genePartners which maximize correlation between the group (based on gene expression) and the clinical data in the expression.with.survival.file.path.training dataset . Todays, the index to measure the correlation is survDiff-g-rho-p-value.
#Then, it will generate 
#	1-A PDF (output.file.path with the suffix _TRAINING) with the kaplanMeier+logRankTest+coxph+wald.test in training test, for the best found in previous step. 
#	2-Another PDF (output.file.path with the suffix _TESTING) with the kaplanMeier+logRankTest+coxph+wald.test in testing test, for the best found in previous step.
#
#How to interpret the result: Look at the "output.file.path with the suffix _TESTING" file and look if there are some still goodIndexes. If so, it is probably a good marker because it correlates in both datasets (training and testing).  
#
#Input
#   expression.with.survival.file.path.training: It is the file path of the expression and clinical data file for finding the best 10. It must respect the following format.
#      -Row 1: Sample labels
#	   -Row 2: Surival time (it could be for example relapse free survival time or overall suvival)
#	   -Row 3: Survival event. 0 means the event did not occur. 1 the event occurred.
#	   -From Row 4: expression data 
#	   -Column 1: Gene symbol (ejemplo: A1BG, A2M) except in row 1, 2 and 3 that are the labels for the rows descripted below.
#
#   expression.with.survival.file.path.testing: It is the file path of the expression and clinical data file for testing the best 10. It must respect the following format.
#      -Row 1: Sample labels
#	   -Row 2: Surival time (it could be for example relapse free survival time or overall suvival)
#	   -Row 3: Survival event. 0 means the event did not occur. 1 the event occurred.
#	   -From Row 4: expression data 
#	   -Column 1: Gene symbol (ejemplo: A1BG, A2M) except in row 1, 2 and 3 that are the labels for the rows descripted below.
#
#	number.of.clusters: Number of clusters for measuring correlation
#
#	output.file.path: The file name of the otuput (without .pdf) 
#	
#	x.lab: the label for the kaplanMeier plots
#
#	y.lab: the label for the kaplanMeier plots
#
# 	gene.name: Gene name of the candidate. This function will evaluate all the pairs (gene.name, X) where X is each gene in the input file.
#
#	grouping.FUN: The function for assigning cluster id. See expression_grouping_functions.R
#
#	print.surv.diff: TRUE if you want to print survDiff in the KaplanMeier PDF
#
#	print.coxph: TRUE if you want to print coxph and wald.test in the KaplanMeier PDF
#
#	partner.gene.names.to.evaluate: You can pass a vector with the gene names to evaluate if you dont want to evaluate each gene in the input file.
#	
#	minimium.number.of.samples.in.a.group: If the grouping function creates a group with less than "minimium.number.of.samples.in.a.group", the pair will be discarded.   
#
#	number.of.genes.to.keep.during.training: How many bests will keep the first step on training.
#
#Output
#	1-A PDF (output.file.path with the suffix _TRAINING) with the kaplanMeier+logRankTest+coxph+wald.test in training test, for the best found in previous step. 
#	2-Another PDF (output.file.path with the suffix _TESTING) with the kaplanMeier+logRankTest+coxph+wald.test in testing test, for the best found in previous step.
#
bestKaplanMeierForEachGenePartner <- function(expression.with.survival.file.path.training, expression.with.survival.file.path.testing, number.of.clusters, output.file.path, x.lab="Time", y.lab="survival", gene.name, grouping.FUN, print.surv.diff=TRUE, print.concordance.index=TRUE, print.coxph=TRUE, partner.gene.names.to.evaluate=NULL, minimium.number.of.samples.in.a.group=10, number.of.genes.to.keep.during.training=10){
	betters.genes<-betters.partners(expression.with.survival.file.path.training, number.of.clusters, gene.name, grouping.FUN, partner.gene.names.to.evaluate=partner.gene.names.to.evaluate, minimium.number.of.samples.in.a.group=10, number.of.genes.to.keep.during.training)
	output.for.training=paste(output.file.path, "_TRAINING")
	KaplanMeierForEachGenePartner(expression.with.survival.file.path.training, number.of.clusters, output.for.training, x.lab="Time", y.lab="survival", maximum.p.value.accepted=1, gene.name, grouping.FUN, print.surv.diff=TRUE, print.concordance.index=TRUE, print.coxph=TRUE, partner.gene.names.to.evaluate=betters.genes, minimium.number.of.samples.in.a.group=10)
	output.for.testing=paste(output.file.path, "_TESTING")
	KaplanMeierForEachGenePartner(expression.with.survival.file.path.testing, number.of.clusters, output.for.testing, x.lab="Time", y.lab="survival", maximum.p.value.accepted=1, gene.name, grouping.FUN, print.surv.diff=TRUE, print.concordance.index=TRUE, print.coxph=TRUE, partner.gene.names.to.evaluate=betters.genes, minimium.number.of.samples.in.a.group=10)
}




#It will generate a PDF (rolled in 50 pages) with the  kaplanMeier+logRankTest+coxph+wald.test for each gene pair with survDifpValue<maximum.p.value.accepted. The gene pairs are formed of (gene.name+geneX) where geneX is each gene in the inputFile (if partner.gene.names.to.evaluate is set, then all genes that are not in this list will not be evaluated). 
#This method should be used for evaluating the kaplanMeier+logRankTest+coxph+wald.test for a set of genes previously picked up. 
#IT SHOULD NOT BE USED FOR PICKING UP GENES BECAUSE IT DOESNT DO TESTING. IF YOU WANT TO PICK UP USE bestKaplanMeierForEachGenePartner. That is why the parameter partner.gene.names.to.evaluate should be always be used, if the input file has got more than the genes that you need to evaluate.  
#
#Input
#   expression.with.survival.file.path: It is the file path of the expression and clinical data file for finding the best 10. It must respect the following format.
#      -Row 1: Sample labels
#	   -Row 2: Surival time (it could be for example relapse free survival time or overall suvival)
#	   -Row 3: Survival event. 0 means the event did not occur. 1 the event occurred.
#	   -From Row 4: expression data 
#	   -Column 1: Gene symbol (ejemplo: A1BG, A2M) except in row 1, 2 and 3 that are the labels for the rows descripted below.
#
#	number.of.clusters: Number of clusters for measuring correlation
#
#	output.file.path: The file name of the otuput (without .pdf) 
#	
#	x.lab: the label for the kaplanMeier plots
#
#	y.lab: the label for the kaplanMeier plots
#
#	maximum.p.value.accepted: If a survDif-pvalue is higher than this value, then the pair will be discardes (i.e. it will not appear in the PDF)
#
# 	gene.name: Gene name of the candidate. This function will evaluate all the pairs (gene.name, X) where X is each gene in the input file.
#
#	grouping.FUN: The function for assigning cluster id. See expression_grouping_functions.R
#
#	print.surv.diff: TRUE if you want to print survDiff in the KaplanMeier PDF
#
#	print.coxph: TRUE if you want to print coxph and wald.test in the KaplanMeier PDF
#
#	partner.gene.names.to.evaluate: You can pass a vector with the gene names to evaluate if you dont want to evaluate each gene in the input file.
#	
#	minimium.number.of.samples.in.a.group: If the grouping function creates a group with less than "minimium.number.of.samples.in.a.group", the pair will be discarded.   
#
#
#Output
#	PDF (rolled in 50 pages) with the  kaplanMeier+logRankTest+coxph+wald.test for each gene pair with survDifpValue<maximum.p.value.accepted. The gene pairs are formed of (gene.name, geneX) where geneX is each gene in the inputFile (if partner.gene.names.to.evaluate is set, then all genes that are not in this list will not be evaluated). 
KaplanMeierForEachGenePartner <- function(expression.with.survival.file.path, number.of.clusters, output.file.path, x.lab="Time", y.lab="survival", maximum.p.value.accepted=0.05, gene.name, grouping.FUN, print.surv.diff=TRUE, print.concordance.index=TRUE, print.coxph=TRUE, partner.gene.names.to.evaluate=NULL, minimium.number.of.samples.in.a.group=10){
	doKaplanMeierForEachGenePartner(expression.with.survival.file.path, number.of.clusters, output.file.path, x.lab, y.lab, maximum.p.value.accepted, gene.name, grouping.FUN, print.surv.diff, print.concordance.index, print.coxph, partner.gene.names.to.evaluate, minimium.number.of.samples.in.a.group)
}


#It will generate an spreadsheet with the  logRankTest+coxph+wald.test for each gene pair with survDifpValue<maximum.p.value.accepted. The gene pairs are formed of (gene.name+geneX) where geneX is each gene in the inputFile (if partner.gene.names.to.evaluate is set, then all genes that are not in this list will not be evaluated). 
##This method should be used for evaluating the kaplanMeier+logRankTest+coxph+wald.test for a set of genes previously picked up.
#THIS METHOD SHOULD NOT BE USED FOR PICKING UP GENES BECAUSE IT DOESNT DO TESTING. IF YOU WANT TO PICK UP USE bestKaplanMeierForEachGenePartner. It would be used for generate    
#
#Input
#   expression.with.survival.file.path: It is the file path of the expression and clinical data file for finding the best 10. It must respect the following format.
#      -Row 1: Sample labels
#	   -Row 2: Surival time (it could be for example relapse free survival time or overall suvival)
#	   -Row 3: Survival event. 0 means the event did not occur. 1 the event occurred.
#	   -From Row 4: expression data 
#	   -Column 1: Gene symbol (ejemplo: A1BG, A2M) except in row 1, 2 and 3 that are the labels for the rows descripted below.
#
#	number.of.clusters: Number of clusters for measuring correlation
#
#	output.file.path: The file name of the otuput 
#	
#	.xlab: the label for the kaplanMeier plots
#
#	y.lab: the label for the kaplanMeier plots
#
#	maximum.p.value.accepted: If a survDif-pvalue is higher than this value, then the pair will be discardes (i.e. it will not appear in the PDF)
#
# 	gene.name: Gene name of the candidate. This function will evaluate all the pairs (gene.name, X) where X is each gene in the input file.
#
#	grouping.FUN: The function for assigning cluster id. See expression_grouping_functions.R
#
#	print.surv.diff: TRUE if you want to print survDiff in the KaplanMeier PDF
#
#	print.coxph: TRUE if you want to print coxph and wald.test in the KaplanMeier PDF
#
#	partner.gene.names.to.evaluate: You can pass a vector with the gene names to evaluate if you dont want to evaluate each gene in the input file.
#	
#	minimium.number.of.samples.in.a.group: If the grouping function creates a group with less than "minimium.number.of.samples.in.a.group", the pair will be discarded.   
#
#
#Output
#	spreadsheet with the  logRankTest+coxph+wald.test for each gene pair with survDifpValue<maximum.p.value.accepted. The gene pairs are formed of (gene.name+geneX) where geneX is each gene in the inputFile (if partner.gene.names.to.evaluate is set, then all genes that are not in this list will not be evaluated). 
multiomics.to.spreadsheet <- function(expression.with.survival.file.path, number.of.clusters, output.file.path, maximum.p.value.accepted=0.05, gene.name, grouping.FUN, print.surv.diff=TRUE, print.concordance.index=TRUE, print.coxph=TRUE, partner.gene.names.to.evaluate=NULL, minimium.number.of.samples.in.a.group=10){
	do.multiomics.to.spreadsheet(expression.with.survival.file.path, number.of.clusters, output.file.path, maximum.p.value.accepted, gene.name, grouping.FUN, print.surv.diff, print.concordance.index, print.coxph, partner.gene.names.to.evaluate, minimium.number.of.samples.in.a.group)
}
	


######################################END API#############################################
######################################END API#############################################
######################################END API#############################################
















######################################PRIVATE#############################################
######################################PRIVATE#############################################
######################################PRIVATE#############################################


doKaplanMeierForEachGenePartner <- function(expression.with.survival.file.path, number.of.clusters, output.file.path, x.lab="Time", y.lab="survival", maximum.p.value.accepted=0.05, gene.name, grouping.FUN, print.surv.diff=TRUE, print.concordance.index=TRUE, print.coxph=TRUE, partner.gene.names.to.evaluate=NULL, minimium.number.of.samples.in.a.group=10){
	library(survival)
	print("doKaplanMeierForEachGenePartner")
	
	#Constants
	col.position.gene.pair=1
	col.name.observations="observations"
	col.position.observations=2
	col.name.gene.pair="gene.pair"
	ERROR.GROUPING="ERROR.GROUPING"
	ERROR.SURVDIFF="ERROR.EXECUTING.SURVDIFF"
	ERROR.EXECUTING.SURV.FIT.FOR.PLOTTING="ERROR.EXECUTING.SURVFIT.FOR.PLOTTING"
	ERROR.CONCORDANCE.INDEX="ERROR.EXECUTING.CONCORDANCE.INDEX"
	ERROR.COXPH="ERROR.EXECUTING.COXPH"
	
	#Create the output folder for the csv file, if the folder doesnt exists	
	pdf.file.path<-paste(output.file.path, "_1.pdf", sep="")
	createAndOpenPDF(pdf.file.path)
	
	#It reads the file and creates an expressionXSruvivalObject containing the expression and the clinical data
	expression.x.survival.object<-readExpressionXSurvivalFromTypicalFile(expression.with.survival.file.path)
	print(paste("The Input file: ", expression.with.survival.file.path, "was read", sep=""))
	
	#It gets the expression vector for the candidate gene.	
	expressionOfCandidateGene<-expressionVectorByGeneName(expression.x.survival.object, gene.name)
	
	#Filter genes keeping just the ones the user pass as parameters
	if (!is.null(partner.gene.names.to.evaluate)) {
		genes.to.keep<-partner.gene.names.to.evaluate
		if (!is.element(gene.name, partner.gene.names.to.evaluate)) genes.to.keep<-append (genes.to.keep, gene.name)
		expression.x.survival.object<-filter.genes(expression.x.survival.object, genes.to.keep)
	}
	
	cant<-0
	npdf<-2
	
	for (i in 1:expression.x.survival.object@ngenes) {
		
		if (i%%10 == 0)print(paste("Fila del archivo original: ", i))
		
		expression.matrix.for.the.gene.pair<-getNumericMatrixWithTheExpressionOfPairOfGenesFaster(expression.x.survival.object, expressionOfCandidateGene, i)
		
		gene.partner.name<-geneNameInPosition(expression.x.survival.object, i)
		
		gene.pair.name<-paste("(", gene.name, ",", gene.partner.name, ")", sep="")
		
		theLegend<-c("")
		
		tryCatch({
					
					#Grouping
					tryCatch(
							{	groups<-grouping.FUN(expression.matrix.for.the.gene.pair, number.of.clusters)
								#Validation for checking if all groups have got well formed groups
								result.check.groups.are.well.formed<-checkGroupsAreWellFormed(groups, gene.pair.name, minimium.number.of.samples.in.a.group)	
								if (!result.check.groups.are.well.formed@OK) stop(result.check.groups.are.well.formed@message)
								
							},error=function(e){stop(formatErrorMessage(error.type=ERROR.GROUPING, error.detail=e$message))})
					
					n.groups <- length(unique(groups))
					
					#Survdiff
					tryCatch({	
								surv.diff.object<-getSurvDiff(expression.x.survival.object@time, expression.x.survival.object@event, groups)
							},error=function(e){stop(formatErrorMessage(error.type=ERROR.SURVDIFF, error.detail=e$message))})
					
					
					
					
					if ((gene.name!=gene.partner.name) && (surv.diff.object@surv.diff.p.value<maximum.p.value.accepted)){
						
						#coxph
						tryCatch({	
									if (print.coxph){ 
										multiomics.coxph.object<-multiomics.coxph.for.two.genes(gen=gene.pair.name, time.vector=expression.x.survival.object@time, event.vector=expression.x.survival.object@event, expression.of.candidate.gene=expressionOfCandidateGene, expression.of.potential.partner.gene=expression.matrix.for.the.gene.pair[2,])
										theLegend<-append(theLegend, coxphToString(multiomics.coxph.object))
									}
								},error=function(e){stop(formatErrorMessage(error.type=ERROR.COXPH, error.detail=e$message))})
						
						
						#Survdiff
						tryCatch({
									if (print.surv.diff){
										theLegend<-append(theLegend, survDifftoString(surv.diff.object))
									}
								},error=function(e){stop(formatErrorMessage(error.type=ERROR.SURVFIT, error.detail=e$message))})						
						
						
						
						#concordance index
						tryCatch({
									if (print.concordance.index){
										numeric.groups=replace(replace(groups, groups=="up", 1),groups=="down", 2) 
										concordance.index.object<-getConcordanceIndex(expression.x.survival.object@time, expression.x.survival.object@event, numeric.groups)
										theLegend<-append(theLegend, concordnaceIndextoString(concordance.index.object))
									}
								},error=function(e){stop(formatErrorMessage(error.type=ERROR.CONCORDANCE.INDEX, error.detail=e$message))})
						
						
						#SurvFit for plotting
						tryCatch({
									surv.fit<-survfit(formula = Surv(expression.x.survival.object@time, expression.x.survival.object@event) ~ groups)
									plot(surv.fit,col=c("blue", "red", "black"), xlab=x.lab, ylab=y.lab)
									title<-gene.pair.name
									drawLegends(groups, theLegend, title)
								},error=function(e){stop(formatErrorMessage(error.type=ERROR.EXECUTING.SURV.FIT.FOR.PLOTTING, error.detail=e$message))})
						
						
						cant=cant+1
						if (cant%%50==0) {
							dev.off()
							createAndOpenPDF(paste(output.file.path, "_", npdf ,".pdf"))
							npdf<-npdf+1
						}
						
						
					}
				}, error=function(e){
					if (grepl(ERROR.CONCORDANCE.INDEX, e$message)) e$message<-paste(ERROR.CONCORDANCE.INDEX, ". Possible causes: groups are not enough balanced", sep="")
					plot(0,type='n',axes=FALSE,ann=FALSE)
					legend("top", legend=gene.pair.name, cex = 2, bty='n')
					legend("center", e$message)
					print(paste("Error on pair ", gene.pair.name, ". Message: ", e$message, sep=""))
				})	
	}  
	dev.off()
	
}

do.multiomics.to.spreadsheet <- function(expression.with.survival.file.path, number.of.clusters, output.file.path, maximum.p.value.accepted=0.05, gene.name, grouping.FUN, print.surv.diff=TRUE, print.concordance.index=TRUE, print.coxph=TRUE, partner.gene.names.to.evaluate=NULL, minimium.number.of.samples.in.a.group=10){
	library(survival)
	print("multiomics.to.spreadsheet")
	
	#Constants
	col.position.gene.pair=1
	col.name.observations="observations"
	col.position.observations=2
	col.name.gene.pair="gene.pair"
	ERROR.GROUPING="ERROR.GROUPING"
	ERROR.SURVDIFF="ERROR.EXECUTING.SURVDIFF"
	ERROR.SURVDIFF="ERROR.EXECUTING.SURVFIT"
	ERROR.CONCORDANCE.INDEX="ERROR.EXECUTING.CONCORDANCE.INDEX"
	ERROR.COXPH="ERROR.EXECUTING.COXPH"
	
	COL_SURVDIFF_P_VALUE = "survDiff G-RHO p-value"
	COL_SURVDIFF = "survDiff G-RHO score"
	
	
	
	#Shared Environment with error functions
	res.env<-new.env()
	
	#Create the output folder for the csv file, if the folder doesnt exists	
	createFolderIfDoesntExist(output.file.path)
	
	#It reads the file and creates an expressionXSruvivalObject containing the expression and the clinical data
	expression.x.survival.object<-readExpressionXSurvivalFromTypicalFile(expression.with.survival.file.path)
	print(paste("The Input file: ", expression.with.survival.file.path, "was read", sep=""))
	
	#It gets the expression vector for the candidate gene.	
	expressionOfCandidateGene<-expressionVectorByGeneName(expression.x.survival.object, gene.name)
	
	#Filter genes keeping just the ones the user pass as parameters
	if (!is.null(partner.gene.names.to.evaluate)) {
		genes.to.keep<-partner.gene.names.to.evaluate
		if (!is.element(gene.name, partner.gene.names.to.evaluate)) genes.to.keep<-append (genes.to.keep, gene.name)
		expression.x.survival.object<-filter.genes(expression.x.survival.object, genes.to.keep)
	}
	
	#It defines the columns for the output file, considering the user parameters
	colnames<-c("gene.pair")
	if (print.coxph)colnames<-append(colnames, c("coxph.coef", "coxph.exp.coef", "coxph.Rsquare", "coxph.concordance", "coxph.log.rank.p.value","wald.test.p.value")) 
	if (print.surv.diff) colnames<-append(colnames, c(COL_SURVDIFF, COL_SURVDIFF_P_VALUE))
	if (print.concordance.index) colnames<-append(colnames, c("Concordance index"))
	colnames<-append(colnames, c("result"))
	
	
	
	
	#It creates the matrix
	res.env$mat <- matrix(, nrow = expression.x.survival.object@ngenes, ncol = length(colnames))
	colnames(res.env$mat)<-colnames
	
	#It calculates the number of fixed and dynamic columns
	a.number.of.fixed.columns=2
	a.number.of.metrics=length(colnames)-a.number.of.fixed.columns
	
	
	for (i in 1:expression.x.survival.object@ngenes) {
		
		if (i%%10 == 0)print(paste("Fila del archivo original: ", i))
		
		expression.vector.of.potential.partner <- expressionVectorByPosition(expression.x.survival.object, i)
		
		expression.matrix.for.the.gene.pair<-getNumericMatrixWithTheExpressionOfPairOfGenesFaster(expression.x.survival.object, expressionOfCandidateGene, i)
		
		gene.partner.name<-geneNameInPosition(expression.x.survival.object, i)
		
		gene.pair.name<-paste("(", gene.name, ",", gene.partner.name, ")", sep="")
		
		tryCatch({
					
					newrow<-c(gene.pair.name)
					
					#if ((gene.name!=gene.partner.name) && (surv.diff.object@surv.diff.p.value<maximum.p.value.accepted)){
					if (gene.name!=gene.partner.name) {
					
					
						#coxph
						tryCatch({	
									if (print.coxph){ 
										multiomics.coxph.object<-multiomics.coxph.for.two.genes(gen=gene.pair.name, time.vector=expression.x.survival.object@time, event.vector=expression.x.survival.object@event, expression.of.candidate.gene=expressionOfCandidateGene, expression.of.potential.partner.gene=expression.vector.of.potential.partner)
										newrow<-append(newrow, c(multiomics.coxph.object@coxph.coef, multiomics.coxph.object@coxph.exp.coef, multiomics.coxph.object@coxph.Rsquare, multiomics.coxph.object@coxph.concordance, multiomics.coxph.object@coxph.log.rank.p.value,multiomics.coxph.object@wald.test.p.value))
									}
								},error=function(e){stop(formatErrorMessage(error.type=ERROR.COXPH, error.detail=e$message))})
						
						
					#Grouping
					tryCatch(
							{	groups<-grouping.FUN(expression.matrix.for.the.gene.pair, number.of.clusters)
								#Validation for checking if all groups have got well formed groups
								result.check.groups.are.well.formed<-checkGroupsAreWellFormed(groups, gene.pair.name, minimium.number.of.samples.in.a.group)	
								if (!result.check.groups.are.well.formed@OK) stop(result.check.groups.are.well.formed@message)
							},error=function(e){stop(formatErrorMessage(error.type=ERROR.GROUPING, error.detail=e$message))})
					
					
					
					
					#Survdiff
					tryCatch({	
								surv.diff.object<-getSurvDiff(expression.x.survival.object@time, expression.x.survival.object@event, groups)
							},error=function(e){stop(formatErrorMessage(error.type=ERROR.SURVDIFF, error.detail=e$message))})
					
					
						#Survdiff
						tryCatch({
									if (print.surv.diff){
										#newrow<-append(newrow, surv.diff.object@surv.diff.p.value)
										newrow<-append(newrow, c(surv.diff.object@surv.diff.chi.squared, surv.diff.object@surv.diff.p.value))
									}
								},error=function(e){stop(formatErrorMessage(error.type=ERROR.SURVFIT, error.detail=e$message))})						
						
						
						
						#concordance index
						tryCatch({
									if (print.concordance.index){
										numeric.groups=replace(replace(groups, groups=="up", 1),groups=="down", 2) 
										concordance.index.object<-getConcordanceIndex(expression.x.survival.object@time, expression.x.survival.object@event, numeric.groups)
										newrow<-append(newrow, concordance.index.object@concordance.index)
									}
								},error=function(e){stop(formatErrorMessage(error.type=ERROR.CONCORDANCE.INDEX, error.detail=e$message))})
						
						
						res.env$mat[i,] <- append(newrow, "OK")								
						
					}
#				}, error=function(e){
#					if (grepl(ERROR.CONCORDANCE.INDEX, e$message)) e$message<-paste(ERROR.CONCORDANCE.INDEX, ". Possible causes: groups are not enough balanced", sep="")
#					error.row<-c(gene.pair.name, e$message)#gene.pair.name, paste("Error - Just ", length(which(groups==group.names[k])), " sample(s) in group ", group.names[k]))
#					for(j in a.number.of.fixed.columns+1:a.number.of.metrics+a.number.of.fixed.columns){
#						error.row<-append(error.row, 999)
#					}
#					res.env$mat[i,] <- error.row
#					print(paste("Error on pair ", gene.pair.name, ". Message: ", e$message, sep=""))
#				})	
				
				
				}, error=function(e){
					if (grepl(ERROR.CONCORDANCE.INDEX, e$message)) e$message<-paste(ERROR.CONCORDANCE.INDEX, ". Possible causes: groups are not enough balanced", sep="")
					error.row<-newrow
					if (length(error.row)+1<a.number.of.metrics+a.number.of.fixed.columns){
						for(j in (length(newrow)+1):(a.number.of.metrics+a.number.of.fixed.columns-1)){
							error.row<-append(error.row, NaN)
						}
					}
					error.row<-append(error.row, e$message)
					res.env$mat[i,] <- error.row
					print(paste("Error on gene ", gene.pair.name, ". Message: ", e$message, sep=""))
				})
				
				
	}  
	#remover NAS
	res.env$mat<- orderMatrixRows(res.env$mat, decr = F, cols = COL_SURVDIFF)
	write.table(res.env$mat, output.file.path,  row.names=FALSE, sep="\t", quote=FALSE)
	
}

#It picks up the best genes considering survDiff-pvalue. It returns those better genes as a list.
###revisar luego del refactor que le hice a la busqueda de un solo gen. este no lo modifique. se debe pasar el atributo observatios al final y renombrarlo por result.
betters.partners <- function (expression.with.survival.file.path, number.of.clusters, gene.name, grouping.FUN, partner.gene.names.to.evaluate=NULL, minimium.number.of.samples.in.a.group=10, number.of.genes.to.keep.during.training=10){
	library(survival)
	print("betters.partners")
	
	#Constants
	col.position.gene.pair=1
	col.name.observations="observations"
	col.position.observations=2
	col.name.gene.pair="gene.pair"
	ERROR.GROUPING="ERROR.GROUPING"
	ERROR.SURVDIFF="ERROR.EXECUTING.SURVDIFF"
	ERROR.EXECUTING.SURV.FIT.FOR.PLOTTING="ERROR.EXECUTING.SURVFIT.FOR.PLOTTING"
	ERROR.CONCORDANCE.INDEX="ERROR.EXECUTING.CONCORDANCE.INDEX"
	ERROR.COXPH="ERROR.EXECUTING.COXPH"
	
	better.genes.matrix<-matrix(, nrow = number.of.genes.to.keep.during.training, ncol = 2)
	for (i in 1:number.of.genes.to.keep.during.training){better.genes.matrix[i,]<-c("NA", 9999)}
	
	#It reads the file and creates an expressionXSruvivalObject containing the expression and the clinical data
	expression.x.survival.object<-readExpressionXSurvivalFromTypicalFile(expression.with.survival.file.path)
	print(paste("The Input file: ", expression.with.survival.file.path, "was read", sep=""))
	
	#It gets the expression vector for the candidate gene.	
	expressionOfCandidateGene<-expressionVectorByGeneName(expression.x.survival.object, gene.name)
	
	#Filter genes keeping just the ones the user pass as parameters
	if (!is.null(partner.gene.names.to.evaluate)) {
		genes.to.keep<-partner.gene.names.to.evaluate
		if (!is.element(gene.name, partner.gene.names.to.evaluate)) genes.to.keep<-append (genes.to.keep, gene.name)
		expression.x.survival.object<-filter.genes(expression.x.survival.object, genes.to.keep)
	}
	
	for (i in 1:expression.x.survival.object@ngenes) {
		
		if (i%%10 == 0)print(paste("Fila del archivo original: ", i))
		
		expression.matrix.for.the.gene.pair<-getNumericMatrixWithTheExpressionOfPairOfGenesFaster(expression.x.survival.object, expressionOfCandidateGene, i)
		
		gene.partner.name<-geneNameInPosition(expression.x.survival.object, i)
		gene.pair.name<-paste("(", gene.name, ",", gene.partner.name, ")", sep="")
		
		tryCatch({
					#Grouping
					tryCatch(
							{	groups<-grouping.FUN(expression.matrix.for.the.gene.pair, number.of.clusters)
								#Validation for checking if all groups have got well formed groups
								result.check.groups.are.well.formed<-checkGroupsAreWellFormed(groups, gene.pair.name, minimium.number.of.samples.in.a.group)	
								if (!result.check.groups.are.well.formed@OK) stop(result.check.groups.are.well.formed@message)
							},error=function(e){stop(formatErrorMessage(error.type=ERROR.GROUPING, error.detail=e$message))})
					
					#Survdiff
					tryCatch({	
								surv.diff.object<-getSurvDiff(expression.x.survival.object@time, expression.x.survival.object@event, groups)
							},error=function(e){stop(formatErrorMessage(error.type=ERROR.SURVDIFF, error.detail=e$message))})
					
					if (gene.name!=gene.partner.name){
						
						
						#Si es mas chico que el mas grande hasta el momento, primero se hacce el calculo de concordance y coxph para garantizar que no de error. Si da error se descarta.
						#Para determinar quien es mejor que quien, utilizamos log.rank.test (que es el survdiff con rho=0)
						if (surv.diff.object@surv.diff.p.value<as.double(better.genes.matrix[number.of.genes.to.keep.during.training,2])){
							
							#coxph
							tryCatch({	
										multiomics.coxph.object<-multiomics.coxph.for.two.genes(gen=gene.pair.name, time.vector=expression.x.survival.object@time, event.vector=expression.x.survival.object@event, expression.of.candidate.gene=expressionOfCandidateGene, expression.of.potential.partner.gene=expression.matrix.for.the.gene.pair[2,])
									},error=function(e){stop(formatErrorMessage(error.type=ERROR.COXPH, error.detail=e$message))})
							
							#concordance index
							tryCatch({
										numeric.groups=replace(replace(groups, groups=="up", 1),groups=="down", 2)
										concordance.index.object<-getConcordanceIndex(expression.x.survival.object@time, expression.x.survival.object@event, numeric.groups)
									},error=function(e){stop(formatErrorMessage(error.type=ERROR.CONCORDANCE.INDEX, error.detail=e$message))})
							
							
							
							#Si llega aqui es porque pudo calcular las 3 metricas
							better.genes.matrix<-rbind(better.genes.matrix, c(gene.partner.name, surv.diff.object@surv.diff.p.value))
							better.genes.matrix<-better.genes.matrix[sort.list(better.genes.matrix[,2]), ]
							better.genes.matrix<-better.genes.matrix[1:number.of.genes.to.keep.during.training,]
						}
					}
					
				}, error=function(e){
					print(paste("Error on pair ", gene.pair.name, ". Message: ", e$message, sep=""))
				})	
	}  
	return (better.genes.matrix[,1])
}






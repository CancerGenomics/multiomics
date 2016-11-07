##########################################################################################################################
# #This pipeline calculates, gene by gene in an independente way, the correlation between survival and expression data. 
# Author: Matias
###########################################################################################################################


#This pipeline do 3 steps:
#		Step1: It generates the file (For example 1-TCGAExpressionAndOS.csv) joining the expression and followUpData, keeping just the samples that are in both files and that has the follow up data
#		Step2: It generates the file (For example 2-TCGACorrelationBetweenExpressionAndSurvivalForEachGene) that shows the correlation between each gen expression and its follow up date. The file will have one row per gen showing the coxph.socre and the coxph.p.value 
#		Step3: It generates the pdf (For example 3-TCGAKaplanMeierGeneByGene.pdf) with the kaplanMeier+survdiff p.value (showing if the difference between curves is statistical significance) just with the genes with surv.diff p-value lower than a parameter (now it is set to 0.01) 
#
#If you have the file with the fowllowing format you should avoid first step.
#      -Row 1: Sample labels
#	   -Row 2: Surival time (it could be for example relapse free survival time or overall suvival)
#	   -Row 3: Survival event. 0 means the event did not occur. 1 the event occurred.
#	   -From Row 4: expression data 
#	   -Column 1: Gene symbol (ejemplo: A1BG, A2M) except in row 1, 2 and 3 that are the labels for the rows descripted below.
#
#If you have two different files: one for the expression data and another one with the clinical data, you can use the "generateExpressionAndSurvivalData" function of this package.
#
#This pipeline example uses TCGA data 

###################################################################################################################################
##############################################SOURCE###############################################################################
###################################################################################################################################

#source("D:/desarrollo/workspaces/R/multiomics/3-do/API/API_multiomics_for_survival_gene_by_gene.R", echo=FALSE, encoding="Cp1252")
#source("D:/desarrollo/workspaces/R/BioplatUtils/matrix_utils.R", echo=FALSE, encoding="Cp1252")
# path adonde estan los archivos de entrada
sourceBaseLocation="D:\\desarrollo\\workspaces\\R\\multiomics\\"
source(paste(sourceBaseLocation, "src/3-do/API/survivalGeneByGene/API_multiomics_for_survival_gene_by_gene.R",sep=""), echo=FALSE, encoding="Cp1252")# path for the input files
source(paste(sourceBaseLocation, "src/survival.utils/matrix_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "src/survival.utils/expression_grouping_functions.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "src/survival.utils/survival_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "src/survival.utils/general_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "src/3-do/Private/multiomics_private_tcga.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "src/survival.utils/file_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "src/survival.entities/SurvdiffEntity.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "src/survival.entities/ConcordanceIndexEntity.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "src/survival.entities/CoxphEntity.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "src/survival.entities/ExpressionXSurvivalEntity.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "src/survival.entities/PotentialGene.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "src/survival.entities/ValidationResult.R",sep=""), echo=FALSE, encoding="Cp1252")


##################EXAMPLE OF getPrognosticStatistics################################################################################
sourceBaseLocation="D:\\desarrollo\\workspaces\\R\\multiomics\\"
expression.file.path<-paste(sourceBaseLocation, "examples\\survival_gene_by_gene\\input\\input1_expressionMatrix", sep="")
clinical.file.path<-paste(sourceBaseLocation, "examples\\survival_gene_by_gene\\input\\input2_clinicalData.tsv", sep="")
clinical.survival.column.name="OVERALL_SURVIVAL"
clinical.event.column.name="overall_survival_indicator"
expression <- read.table(expression.file.path, sep="\t", header=TRUE, na.strings=c("", "NA"),stringsAsFactors=FALSE)
result<-getPrognosticStatistic(expression, number.of.clusters, groupin.FUN='multiomics.cut2', clinical.file.path, clinical.survival.column.name, clinical.event.column.name, minimium.number.of.samples.in.a.group=10)
##########################################################################################################################################################




############################################################################################################################################
##########CONFIG: You should change this variable to set your values########################################################################
############################################################################################################################################
#input.path is the folder in which the expression and clinical file are 

#example
input.path="D:\\desarrollo\\workspaces\\R\\multiomics\\examples\\survival_gene_by_gene\\input\\"
input.expression.file.name="input1_expressionMatrix"
input.clinical.file.name="input2_clinicalData.tsv"
input.clinical.survival.column.name="OVERALL_SURVIVAL"
input.clinical.event.column.name="overall_survival_indicator"

input.clinical.file.path=paste(input.path, input.clinical.file.name, sep="")
input.expression.file.path=paste(input.path, input.expression.file.name, sep="")

#output.pat is the folder in which all generated files will be placed.
output.path="D:\\desarrollo\\workspaces\\R\\multiomics\\examples\\survival_gene_by_gene\\output\\"








####################################################################################################################################################################################################################################
############################# Step1 - it Generates the file ExpressionWithOS.csv containing the expression and just the overall suvival data. It takes as input the expression file (in TCGA format) and the clincal data (in TCGA format). 
####################################################################################################################################################################################################################################

#-----------------Output of this step---------------------------------------------------------------------------------
#
#   It generates a file with the following format
#      -Row 1: Sample labels
#	   -Row 2: Surival time (it could be for example relapse free survival time or overall suvival)
#	   -Row 3: Survival event. 0 means the event did not occur. 1 the event occurred.
#	   -From Row 4: expression data 
#	   -Column 1: Gene symbol (ejemplo: A1BG, A2M) except in row 1, 2 and 3 that are the labels for the rows descripted below.
#
#	All the samples which has NA for the columns event or survival will be eliminated and they will not be i
#
#	Example
#		sampleID	TCGA-A1-A0SB-01	TCGA-A1-A0SD-01	TCGA-A1-A0SE-01
#		OS_event_nature2012	259	437	1320
#		OS_Time_nature2012	0	0	0
#		ARHGEF10L	108.075	91.657	97.179
#
#	If you have a file with the described format, you can avoid this step
#-------------------------------------------------------------------------------------------------------------------------------


#---------CONFIG (you dont need to change it)------------------
output.expression.with.survival.file.name="1-TCGAExpressionAndOS.csv"
output.expression.with.survival.path=paste(output.path, output.expression.with.survival.file.name, sep="")

#----------DO IT - it Generates the file ExpressionWithSurvival.csv taking as input the expression file (in TCGA format) and the clincal data (in TCGA format)#########
generateExpressionAndSurvivalDataFromTCGA(input.expression.file.path, input.clinical.file.path, input.clinical.survival.column.name, input.clinical.event.column.name, output.expression.with.survival.path)
##################################################################################################################################################################
############################# END STEP1 ######################################################################################################################### 
##################################################################################################################################################################




####################################################################################################################################################################################################################################
#############################STEP2 (bestKaplanMeierGeneByGene)- It calculates, gene by gene, the correlation between its expression and the follow up data. It uses coxph that is why it is not necessary to do clustering. 
####################################################################################################################################################################################################################################
#---------GENERAl Config-----------------
expression.with.survival.file.path.training="D:\\desarrollo\\workspaces\\R\\multiomics\\examples\\survival_gene_by_gene\\output\\1-TCGAExpressionAndOS.csv"
expression.with.survival.file.path.testing="D:\\desarrollo\\workspaces\\R\\multiomics\\examples\\survival_gene_by_gene\\output\\1-TCGAExpressionAndOS.csv"
#expression.grouping.FUN=multiomics.cut2
expression.grouping.FUN=multiomics.kmeans

#---------SpecficConfig-----------------
a.minimium.number.of.samples.in.a.group=3
number.of.clusters=2
output.file.path=output.path
#---------Do it-----------------
bestKaplanMeierGeneByGene(expression.with.survival.file.path.training, expression.with.survival.file.path.testing, number.of.clusters, output.file.path, x.lab="Time", y.lab="survival", expression.grouping.FUN, print.surv.diff=TRUE, print.concordance.index=TRUE, print.coxph=TRUE, gene.names.to.evaluate=NULL, minimium.number.of.samples.in.a.group=10, number.of.genes.to.keep.during.training=10)

####################################################################################################################################################################################################################################
#############################END STEP2 #############################################################################################################################################################################################  
####################################################################################################################################################################################################################################
  




####################################################################################################################################################################################################################################
#############################STEP3 - It calculates the survival gene by gen 
####################################################################################################################################################################################################################################

#-----------------Output of this step---------------------------------------------------------------------------------
# It will generate a file "TCGA_kaplanMeier_geneByGene.pdf". 
#
#		gene: Gene symbol
input.expression.with.survival.file.name="1-TCGAExpressionAndOS.csv"
input.expression.with.survival.path=paste(output.path, input.expression.with.survival.file.name, sep="")
output.kaplan.meier.file.name="3-TCGAKaplanMeierGeneByGene.pdf"
kaplanMeierForAGeneUsingCut2GeneByGene(input.expression.with.survival.path,  2, output.path, output.kaplanmeier.file.name=output.kaplan.meier.file.name, x.lab="Time", y.lab="survival", maximum.p.value.accepted=0.01)




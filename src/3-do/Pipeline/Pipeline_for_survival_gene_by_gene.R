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


############################################################################################################################################
##########CONFIG: You should change this variable to set your values########################################################################
############################################################################################################################################
#input.path is the folder in which the expression and clinical file are 
input.path="D:\\datosGenomica\\multiomics\\paper\\TCGA\\input\\"
input.expression.file.name="mrnaGenomicMatrix"
input.clinical.file.name="clinical_data"
input.clinical.survival.column.name="OS_Time_nature2012"
input.clinical.event.column.name="OS_event_nature2012"

input.clinical.file.path=paste(input.path, input.clinical.file.name, sep="")
input.expression.file.path=paste(input.path, input.expression.file.name, sep="")

#output.pat is the folder in which all generated files will be placed.
output.path="D:\\datosGenomica\\multiomics\\paper\\TCGA\\output\\"








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

#---------SOURCE-----------------
source("D:/desarrollo/workspaces/R/TCGA/tcga.R", echo=FALSE, encoding="Cp1252")

#---------CONFIG------------------
input.clinical.survival.column.name="OS_Time_nature2012"
input.clinical.event.column.name="OS_event_nature2012"
output.expression.with.survival.file.name="1-TCGAExpressionAndOS.csv"
output.expression.with.survival.path=paste(output.path, output.expression.with.survival.file.name, sep="")

#These variables are used by the step 1 but were set above. You can change a value if you need.
#input.path="D:\\datosGenomica\\multiomics\\paper\\TCGA\\input\\"
#input.expression.file.name="mrnaGenomicMatrix"
#input.clinical.file.name="clinical_data"
#output.pat is the folder in which all generated files will be placed.
#output.path="D:\\datosGenomica\\multiomics\\paper\\TCGA\\output\\"
#expression.file.name is the name of the file containing expression with survival file. The format of this file was desccribed above. 

#----------DO IT - it Generates the file ExpressionWithSurvival.csv taking as input the expression file (in TCGA format) and the clincal data (in TCGA format)#########
generateExpressionAndSurvivalDataFromTCGA(input.expression.file.path, input.clinical.file.path, input.clinical.survival.column.name, input.clinical.event.column.name, output.expression.with.survival.path)
##################################################################################################################################################################
############################# END STEP1 ######################################################################################################################### 
##################################################################################################################################################################









####################################################################################################################################################################################################################################
#############################STEP2 - It calculates, gene by gene, he correlation between its expression and the follow up data. It uses coxph that is why it is not necessary to do clustering. 
####################################################################################################################################################################################################################################

#---------SOURCE-----------------
source("D:/desarrollo/workspaces/R/GenomicAnalysis/ABCC4-PANCREAS/Config.R", echo=FALSE, encoding="Cp1252")

#---------GENERAl Config-----------------
input.expression.with.survival.file.path = output.expression.with.survival.file.path
expression.grouping.FUN=multiomics.G1PairUpG2PairDown #If ABCC4 and the partner are overexpressed, then the sample is in group1; if not is iin group2.

#---------SpecficConfig-----------------
output.multiomics.to.spreadsheet.file.path<-paste(abcc4pancreas.root.path, "ABCC4.Withpartners.OS.G1UG2D0.7.csv", sep="")
a.maximum.p.value.accepted=1
a.candidate.gene.name="ABCC4"
#partner.genes.of.interest=c("BPESC1","MLLT4","OR7G1","PIGS","PRDM10","PRDM10","SEC63","NDUFS1","PDSS2","TSNAX","VPS41","TAS2R38","GPR107","RAI1","LMTK2","KAL1","KIAA0232","STC1","UBXN7","SRCIN1","LASS6","LOC729082","TAF1L","OR2G2","CMTM4","TMEM106B","SHPRH","RGPD8","PPA2","ZNF611","ZNF616","ASH1L","MYEF2","ATF6","SLC38A2","RAB3GAP1","SDCCAG1","MED17","SYN3","ZNF223","PTPDC1","DLG1","PEX19","ZNF28","C20orf117","CCNT1","EARS2","C18orf25","SLC22A25","C11orf42","AKTIP","UHMK1","AHI1","DYNC1H1","GPR160","OR5T2","AGPAT1","MRE11A","PIP4K2B","BCLAF1","KIAA1826","ZNF322A","ZNF609","COX7B2","DDX52","CDC42BPB","OR11L1","TOR1AIP2","NBAS","LACE1","ZFC3H1","TIRAP","CBX5","MKS1","SNX19","SBNO1","EIF5","ASB4","ASB8","PRPF40B","TMEM184C","REST","MAP3K13","CEACAM6","DHX15","IKBKAP","CNNM4","DIP2B","GRHL2","RBM16","KIAA1731","RPL10L","OR5M11","PDE4A","PHF3")
a.minimium.number.of.samples.in.a.group=10

#---------Do it-----------------
multiomics.to.spreadsheet(input.expression.with.survival.file.path, number.of.clusters=2, output.file.path=output.multiomics.to.spreadsheet.file.path, maximum.p.value.accepted=a.maximum.p.value.accepted, gene.name=a.candidate.gene.name, grouping.FUN=expression.grouping.FUN, print.surv.diff=TRUE, print.concordance.index=TRUE, print.coxph=TRUE,partner.gene.names.to.evaluate=NULL,minimium.number.of.samples.in.a.group=a.minimium.number.of.samples.in.a.group)

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




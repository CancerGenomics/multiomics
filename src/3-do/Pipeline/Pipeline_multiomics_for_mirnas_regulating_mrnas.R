###################################################################################
#Author: Matias
#This pipeline take as input mirna expression data and mrna expression data, find 
#all the pairs mirna x mrna which correlates, keep just the best 30% correlations
#according database mirna score, and add information about prediction and publications
#on this correlation. A correlation could suggest that this mirna is involved in the 
#regulation of this gene expression. 
#The final output will have the following columns
#		Gen_symbol	
#		mature_mirna_id	
#		Mirna_Mrna_Correlation	
#		p_value_Of_Mirna_Mrna_Correlation	
#		mirna_database	database_predicted_score	
#		validation_pubmed_id
###################################################################################


###################################################################################################################################
##############################################SOURCE###############################################################################
###################################################################################################################################

# path adonde estan los archivos de entrada, se puede cambiar por un path absoluto
sourceBaseLocation=getwd()

source(paste(sourceBaseLocation, "/src/2-load/load_multiomics.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_private_data_validation.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_private_multimir_interaction.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/mirnasRegulatingMrnas/API_multiomics_for_finding_mirnas_regulating_mrnas.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/survivalGeneByGene/API_multiomics_for_survival_gene_by_gene.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/survivalGeneByGene/API_multiomics_for_survival_gene_by_gene.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/matrix_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/read_genomic_data_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
library("multiMiR")

##########CONFIG: You should change this 3 variables to set your values##################
#Where input files are and where the intermediate and result files will be stored
#working.path="D:\datosGenomica\multiomics\paper\DCIS\5-mirna"  
#working.path="D:/matias/academia/investigacion/medicina personalizada/8-DatosGenomica/2016-09-13---paper multiomics/"

#example
working.path=paste(sourceBaseLocation, "/test/examples/miRnas_regulating_mRnas/",sep="")

# The mrna expression data. It is the path of a file with the following format
#      -Row 1: It has the sample labels
#	   -Column1: It has the gene symbol (for example: A1BG, A2M)
#	   -The cells has got the expression level of each gene for each sample
#mrna.dif.expr.path.file="gastrica-mrna-expression-chico.csv"
#example
mrna.dif.expr.path.file="mrnas-with-extra-samples.csv"

#The mirna expression data. it must have the following format:
#      -Row 1: It has the sample labels
#	   -Column 1: Mature mirna ID (for example: hsa-miR-22-3p). Note that R is un upper case. This R in uppercase is the way to identify that it is a mature-mirna and not a pre-mirna. The pre-mirna for this mature-mirna is hsa-mir-22. Multimir doesnt work with premirna because the databases it queries, doesnt work with pre-mirna. So, if you pass the pre-mirna (for example hsa-miR-22) multimir does not return anything. So, it is important that this file contains mature mirna ids. With the accession it also works (MIMAT0004982). These IDs could be find in http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0005761 
#	   -Cells has the expression levels of each mirna for each sample. 
#mirna.dif.expr.path.file="gastrica-mirnas-expression.csv"
#example
mirna.dif.expr.path.file="mirnas.csv"
##################END CONFIG###########################################################

###################################OPTIONAL CONFIG#########################################################################################################################
#The step2 keeps just the best 30% of mirna X mrna correlations according the coorelation score given by mirna databases. You can change this 30% modifying this parameter
my.predicted.cut.off=10

#Intermediate file names (it is not necessary to change). You can keep this values.
maturemirna.x.mrna.correlation.file.name="inputStep2-matureMirnaXmrna.csv"
just.betters.maturemirna.X.mrna.considering.mirna.databases="inputStep3_justBettersMirnaXmRNA_ConsideringMirnaDatabases.csv"
###################################END OPTIONAL CONFIG#####################################################################################################################



####### Step1 - It looks for mirna x mrna pairs with high correlation ######### 
#For doing it, it evaluates the correlation between each row of mirna file (representing the expression for a particular mature mirna) against each mrna expresion line (representing the expression for a particular gene) using Pearson 
mrna.dif.expr.path<-paste(working.path, mrna.dif.expr.path.file, sep="")
mirna.dif.expr.path<-paste(working.path, mirna.dif.expr.path.file, sep="")
#It reads mrna file and mirna file. It sorts each file by sample name (that is by column)
mrna.dif.expr <- readMrnaExpressionFile(mrna.dif.expr.path)
mirna.dif.expr <- readMirnaExpressionFile(mirna.dif.expr.path)

#Keep columns which are in both databases
intersection<-keepSameColumns(mrna.dif.expr,mirna.dif.expr)
mrna.dif.expr<-(intersection[[1]])
mirna.dif.expr<-(intersection[[2]])


#BIG FILES FOR TESTING
#mrna.dif.expr<-read.table("C:\\desarrollo\\workspaces\\r\\UAI\\multiomics2\\examples\\mRNAvsmiRNA\\mRNA_TCGA_breast.csv", header = T)
#mirna.dif.expr<-read.table("C:\\desarrollo\\workspaces\\r\\UAI\\multiomics2\\examples\\mRNAvsmiRNA\\miRNA_TCGA_breast.csv", header = T)


calculated <- CalculateCorrelationsMirnaMrnaUsingWCGNA(mrna.dif.expr,mirna.dif.expr, working.path, 
                                                       output.file.name=maturemirna.x.mrna.correlation.file.name,
                                                      r.minimium=0.3)



####### STEP2 - For each high correlation mirna X mrna found in step 1, it get just the 30% (configurable) correlations which has got better predicted and/or validated score in mirna databases. This step also adds predicition and validation information for those correlations. The other 70% will be discarded from the analysis .####
#There are many mirna databases for evaluating prediction and validation. A correlation is considered predicted by a mirna database if this correlation mirnaXmrna is registered in the database. A correlation is considered validated if there is a publication showing the correlation experimentally. The multimir packages offers a function for reading all databases in a single call and keeping the ones with highest correlation score.  
#The added predicted information is mirna_database, Mirna_Mrna_Correlation, p_value_Of_Mirna_Mrna_Correlation, database_predicted_score. And for the ones that were also validated it adds the pubmedId.  
mirnaxrna.path<-paste(working.path, maturemirna.x.mrna.correlation.file.name, sep="")
genes.x.mirnas.path <- paste(working.path,maturemirna.x.mrna.correlation.file.name, sep="")
#genes.x.mirnas<-read.table(genes.x.mirnas.path, header=TRUE)
genes.x.mirnas <- calculated
best <- keepBestGeneXMirnaAccordingCorrelationAndAddMirnaDbInfo(genes.x.mirnas, working.path, 
                                                        output.file=just.betters.maturemirna.X.mrna.considering.mirna.databases,
                                                        predicted.cut.off=my.predicted.cut.off)


####### STEP3 - At this point you have multiple lines for the same pair mrna-mirna because of the multiple databases that predict and/or validate this pair. Assuming that both lines were in the best 30% correlation score ####### 
# In this step you colapse the result, for getting just one row for each pair mrna-mirna. The databases will be colapsed in one column separated by commas.
# For example 
#		MIR1292	hsa-miR-1292-5p	MIMAT0005943	ESRRG	-0.747206480568877	2,52E+07	NO_PUB	diana_microt
#		MIR1292	hsa-miR-1292-5p	MIMAT0005943	ESRRG	-0.747206480568877	2,52E+07	NO_PUB	pita
# will be colapsed into 
#		MIR1292	hsa-miR-1292-5p	MIMAT0005943	ESRRG	-0.747206480568877	2,52E+07	NO_PUB	diana_microt, pita
#
just.betters.maturemirna.X.mrna.considering.mirna.databases.path<-paste(working.path, just.betters.maturemirna.X.mrna.considering.mirna.databases, sep="")
mirnaXmrna <- read.table(just.betters.maturemirna.X.mrna.considering.mirna.databases.path, header=T)
mirnaXmrna <- best
collapsed <- ColapseMirnaXMrna(mirnaXmrna, working.path)





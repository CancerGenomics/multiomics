# TODO: Add comment
# 
# Author: Matias
###############################################################################

#input.dir: the location of the data matrix, and where output will be located
#input.file.name: the input data matrix as a tab delimited text file (simil sampleInputFile_200subjects.txt that comes with the package)
#prefix.name.for.output.files<-short name that will be used for output files
#calibration.parameters<-the column of the "mediansPerDataset.txt" file to use for calibration; the column of the "mediansPerDataset.txt" file to use for calibration; NA will force centering within the test set & -1 will not do any adjustment (when adjustment performed by used)
#has.clinical <- may include tumor size as second row, with 'T' as the gene name, and encoded as binary (0 for size <= 2cm or 1 for size > 2cm), set this variable to FALSE if tumor size is not available
#colapse.method <-  can be mean or iqr (probe with max iqr is selected) # typically, mean is preferred for long oligo and # iqr is preferred for short oligo platforms. The default is MEAN
compare.with.pam.50 <- function(config.dir, input.folder, input.file.name, prefix.name.for.output.files, colapse.method="mean", has.clinical=FALSE){
 
	#Loading libraries
	library(ctc)
	library(heatmap.plus)
 
	#Set the variables that bioclassifier uses
	
	paramDir <<- config.dir
	inputDir<<- input.folder
	inputFile<<-input.file.name
	short<<-prefix.name.for.output.files
	calibrationParameters<<- NA
	hasClinical<<-has.clinical

	####
	# run the assignment algorithm
	####
	subtypePrediction_distributed()
	#source(paste(paramDir,"subtypePrediction_functions.R",sep="\\"))
	#source(paste(paramDir,"subtypePrediction_distributed.R",sep="\\"))
	

}

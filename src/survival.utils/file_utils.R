###############################################################################
# Author: Matias
###############################################################################



createFolderIfDoesntExist <- function(file.path){
	dir<-dirname(file.path)
	if (file.exists(dir)==FALSE){
		dir.create(dir, recursive=TRUE)
	}
}


createAndOpenPDF<-function(output.file.path){
	
	createFolderIfDoesntExist(output.file.path)
	pdf(output.file.path,width=13,height=8)
	
}

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

write.to.file <- function(table, output.path, output.file.name) {
  # Write the result to a file
  write.start <- proc.time()
  file.path<-paste(output.path, paste(output.file.name), sep="")
  fwrite(table,file.path,sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  cat("WRITE TIME: : ", (proc.time() - write.start)["elapsed"], "\n")
}
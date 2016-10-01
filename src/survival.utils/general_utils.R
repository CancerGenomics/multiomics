# TODO: Add comment
# 
# Author: Matias
###############################################################################

formatErrorMessage<-function(error.type, error.detail){
  error.message<-paste(error.type, ". ", error.detail, sep="")
  return (error.message)
}




# TODO: Add comment
# 
# Author: Matias
###############################################################################

checkLengthOfVectorsForKaplanMeier <- function(groups, event, time){
	if (length(groups)!=length(event)) stop(paste("The length of the vector groups is ",length(groups)," while the length of the event group is ", length(event),". Both vector must have the same length. ",sep=""))
	if (length(time)!=length(event)) stop(paste("The length of the vector groups is ",length(time)," while the length of the event group is ", length(event),". Both vector must have the same length. ",sep=""))

}


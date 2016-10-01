# TODO: Add comment
# 
# Author: Matias
###############################################################################

checkLengthOfVectorsForKaplanMeier <- function(groups, event, time){
	if (length(groups)!=length(event)) stop(paste("The length of the vector groups is ",length(groups)," while the length of the event group is ", length(event),". Both vector must have the same length. ",sep=""))
	if (length(time)!=length(event)) stop(paste("The length of the vector groups is ",length(time)," while the length of the event group is ", length(event),". Both vector must have the same length. ",sep=""))

}


checkGroupsAreWellFormed <- function(groups, gene.pair.name, minimium.number.of.samples.in.a.group){
  group.names<-unique(groups)
  for(k in 1:length(group.names)){
    if (length(which(groups==group.names[k])) < minimium.number.of.samples.in.a.group){
      a.message<-paste("Just ", length(which(groups==group.names[k])), " sample(s) in group ", group.names[k])
      result<-new("ValidationResult", OK=FALSE, message=a.message)
      return (result)
    }
  }
  return (new("ValidationResult", OK=TRUE, message=""))
}

# TODO: Add comment
# Represents a class containing the result of Testing if there is a difference between two or more survival curves using the G-rho family of tests, or for a single curve against a known alternative. 
# Author: Matias
###############################################################################


setClass("SurvDiff",
		representation(surv.diff.chi.squared="numeric", surv.diff.p.value="numeric"))

getSurvDiff <-function(time.vector, event.vector, groups.vector){
	library("survival")
	#It is used the parameter rho=0 (default value), so this test is the log.rank test
	surv.diff<-survdiff(formula = Surv(time.vector, event.vector) ~ groups.vector, rho=0) #rho=o means log.rank test 
	a.surv.diff.p.value= 1 - pchisq(surv.diff$chisq, length(surv.diff$n) - 1)
	result<-new("SurvDiff", surv.diff.chi.squared=surv.diff$chisq, surv.diff.p.value=a.surv.diff.p.value)
	return (result)
}

survDifftoString <- function(surv.diff.object){
	elements<-c(paste("survDiff G-RHO p-value: ", surv.diff.object@surv.diff.p.value))
	return (elements)
	
}

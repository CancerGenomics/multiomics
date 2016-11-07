###############################################################################
# Class containing the result of evaluating the correlation between expression data and follow up data. It uses coxph test. It doesnt need to do clustering
# Author: Matias
###############################################################################
library("survival")

setClass("CoxphResult",
		representation(gen="character", coxph.coef="numeric", coxph.exp.coef="numeric", coxph.Rsquare="numeric", coxph.concordance="numeric", coxph.log.rank.p.value="numeric", wald.test.score="numeric",  wald.test.p.value="numeric"))

example<-new("CoxphResult", coxph.coef=1, coxph.exp.coef=2, coxph.Rsquare=3, coxph.concordance=4, coxph.log.rank.p.value=5,  wald.test.p.value=6)

#It receives time, event and expression vector and it generates a PDF with the survival and survdiff p-value.
#It uses cut2 for doing the cluster.
multiomics.coxph.for.one.gene <- function(gen="", time.vector, event.vector, expression.vector){
	coxph <- coxph(Surv(time.vector,event.vector) ~ expression.vector,method="breslow")
	return (create.multiomics.coxph(gen, coxph))
}


coxph.writeInFile<- function(multiomics.coxph, output.file){
	createFolderIfDoesntExist(output.file)
	sink(output.file)
	cat(paste("---------------------------------------\n"))
	cat(paste("Coxph between ", multiomics.coxph@gen, "and follow up data"))
	cat(paste("\n---------------------------------------\n"))
	cat("coxph.coef: ", multiomics.coxph@coxph.coef)
	cat("\n")
	cat("coxph.p.value: ", multiomics.coxph@coxph.exp.coef)
	cat("\n")
	cat("\n")
	cat("wald.p.value: ", multiomics.coxph@wald.test.p.value)
	cat("coxph.Rsquare: ", multiomics.coxph@coxph.Rsquare)
	cat("coxph.concordance: ", multiomics.coxph@coxph.concordance)
	cat("coxph.log.rank.p.value: ", multiomics.coxph@coxph.log.rank.p.value)
	cat("\n")
	sink()
	
}

#It receives time, event and expression vector and it generates a PDF with the survival and survdiff p-value.
#It uses cut2 for doing the cluster.
multiomics.coxph.for.two.genes <- function(gen.pair.name="", time.vector, event.vector, expression.of.candidate.gene, expression.of.potential.partner.gene){
	coxph <- coxph(Surv(time.vector,event.vector) ~ expression.of.candidate.gene + expression.of.potential.partner.gene,method="breslow")
	return (create.multiomics.coxph(gen.pair.name, coxph))	
}




create.multiomics.coxph <- function(gen.name.or.gene.pair.name, coxph){
	summary(coxph)$rsq

	#agregar
	coef<-summary(coxph)$coefficients[1]
	expcoef<-summary(coxph)$coefficients[2]
	rsq<-summary(coxph)$rsq[1]
	concordance<-summary(coxph)$concordance[1]
	waldtest_p_value<-summary(coxph)$waldtest[3]
	log.rank.p.value<-summary(coxph)$sctest[3]

	#score<-summary(coxph)$sctest[1]
	#the.coefficients<-coxph$coefficients
	wald.test<-summary(coxph)$waldtest[1]
	
	
	result<-new("CoxphResult", gen=gen.name.or.gene.pair.name, coxph.coef=coef, coxph.exp.coef=expcoef, coxph.Rsquare=rsq, coxph.concordance=concordance, coxph.log.rank.p.value=log.rank.p.value,  wald.test.score=wald.test, wald.test.p.value=waldtest_p_value)
	#result<-new("CoxphResult", gen=gen.name.or.gene.pair.name, coefficients=the.coefficients, coxph.test.score=score, coxph.test.p.value=score_p_value, wald.test.score=wald.test, wald.test.p.value=waldtest_p_value)
	return (result)
}


#It returns each component of coxph as an element in a vector.
coxphToString <- function(coxph.object){
	elements<-c(paste("coxph.coef: ", coxph.object@coxph.coef), paste("coxph.exp.coef: ", coxph.object@coxph.exp.coef), paste("coxph.Rsquare: ", coxph.object@coxph.Rsquare), paste("coxph.concordance: ", coxph.object@coxph.concordance), paste("coxph.log.rank.p.value: ", coxph.object@coxph.log.rank.p.value),  paste("wald.p.value: ", coxph.object@wald.test.p.value))
	return (elements)
}




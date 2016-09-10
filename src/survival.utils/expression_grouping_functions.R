# TODO: Add comment
# 
# Author: Matias
###############################################################################



multiomics.kmeans <-function (expression, number.of.clusters, an.iter.max=1000){
	#The transpose is done because the kmeans function is expecting the genes as columns
	if (class(expression)=='matrix') expression<-t(expression)
	
	fit <- kmeans(expression,centers=number.of.clusters,iter.max = an.iter.max)
	groups <- fit$cluster
	return(groups)
}



multiomics.cut2 <-function (expression, number.of.clusters, an.iter.max=1000){
	library(Hmisc)
	xgene <- cbind(expression)
	quantile(xgene)
	groups <- cut2(xgene, g=number.of.clusters)
	return(groups)
}


#It must recieve always a mtrix with the vectors of both files
multiomics.G1PairUpG2PairDown <-function (expression.matrix, number.of.clusters, an.iter.max=1000){
	if (class(expression.matrix)!='matrix') stop (paste("The function multiomics.G1PairUpG2PairDown is receving a: ", class(expression.matrix) ," and it must receive a matrix wwitht the expression of the pair gene")) 

	#Le sumo esta constante al quantil, para que los valores que coinciden con el cuantil se cosnideen down y no up. De esa forma me garantizo que los que son up es porque son mayor que el 70% de todos los valores.
	constant=0.0000001
	
	candidate.gene.expression<-expression.matrix[1,]
	if (all(candidate.gene.expression==0)) stop("Candidate gene is full of zeros")
	candidate.quantile<-quantile(candidate.gene.expression, probs=c(.50))+constant
	
	
	partner.gene.expression<-expression.matrix[2,]
	if (all(partner.gene.expression==0)) stop("Partner gene is full of zeros")
	partner.gene.quantile<-quantile(partner.gene.expression, probs=c(.50))+constant

	candidate.upOrDown<-factor(findInterval(candidate.gene.expression, c(-Inf,candidate.quantile,Inf)), labels=c("DOWN","UP")) 
	otherGene.upOrDown<-factor(findInterval(partner.gene.expression, c(-Inf,partner.gene.quantile,Inf)), labels=c("DOWN","UP"))
	
	groups<-c(1:length(partner.gene.expression))
	for (i in 1:length(partner.gene.expression)) {
		if ((candidate.upOrDown[i]=="UP") && (otherGene.upOrDown[i]=="UP")) groups[i]<-"up" else groups[i]<-"down" 
	}
	return(groups)
}

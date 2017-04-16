#The methylation cgXXX is in methylation files of TCGA.
#We have also a file which correlates cgXXX with gene. Take into account that a cgXXX can do methilation over multiple genes and that a gene can be methylazed by multiple cgXXX.
#For doing it, we read a gen from the mrnaFile, then we look in the "methylationPlatformFile" all the cgXXX corresponding to this gene, and then we look for those CGXXX in the methilation file. It keeps the pairs (gene, cgXXX) which correlates negatively, considering the mrna file and the methylation file.
#Take into account that the platform file should Have two columns: GEN - CG. If a gen has goto multiple cgs, there will be many rows for this gene

methXMrnas <- function(mrna, meth, meth.platform, output.path="~/", 
                      output.file.name="methylationXMrna.csv",
                      r.minimium=0.7, 
                      pearsons.method = "pearson", 
                      inc.progress = F){
  
  ptm <- proc.time()
  total.rows=nrow(mrna)*15
  print(paste("Running pipeline methylaion_X_mrnas with", r.minimium, 
              "threshold and pearson's method:", pearsons.method, sep=" "))
  
  
  # The result matix is created
  res <- matrix(nrow=total.rows,ncol=5)
  colnames(res)<-(c("Gene","Location", "methylation-id", "Methylation-mRNA correlation", "p-value"))
  
  actual<-1
  actual.n.correlated<-1
  
  ids.in.meth.dataset<-meth[,1]
  
  print("Start process!")
  for (i in 1:nrow(mrna)) {
    actual<-actual+1
    actual.gen<-as.character(mrna[i,1])
    #position.in.cnv.dataset<-which(ids.in.cnv.dataset == actual.gen)
    #In meths i have each cg, which do methialtion over the gene.
    cgs.for.this.gene<-get.cgs.for.this.genes(actual.gen,meth.platform)
    #Se queda con el primer CNV
    if (length(cgs.for.this.gene)>0){
      for (actual.cg.index in 1:length(cgs.for.this.gene)){
            actual.cg<-cgs.for.this.gene[actual.cg.index]
            #position.in.cnv.dataset<-position.in.cnv.dataset[1]
            actual.mrna<-mrna[i,2:ncol(mrna)]
            #position.in.meth.dataset<-which(ids.in.meth.dataset == actual.cg)
            actual.meth<-meth[meth[,1] == actual.cg,2:ncol(meth)]
            if (nrow(actual.meth)>0){
                  resultado.pearson<-cor.test(as.numeric(actual.mrna),
                                              as.numeric(actual.meth), 
                                              method = pearsons.method)
                  if (!is.na(abs(resultado.pearson$estimate))) {
                    if ((abs(resultado.pearson$estimate) > r.minimium) & (resultado.pearson$estimate<0)) {
                      location<-getGeneLocation(actual.gen);
                      newValue<-c(as.character(actual.gen), location, actual.cg,
                                  resultado.pearson$estimate, resultado.pearson$p.value)
                      res[actual.n.correlated,1:5] <- newValue
                      actual.n.correlated<-actual.n.correlated+1
                    }
                  }
      
                  if ((actual)%%500==0)print(paste("analised ", actual, " from ", total.rows))
                  if ((actual)%%1000==0) {
                    elapsedTime <- (proc.time() - ptm)[3]
                    print(paste(
                      "elapsed time: (seconds)", format2Print(elapsedTime), 
                      " - (minutes)", format2Print(elapsedTime/60), 
                      " - (hours)", format2Print(elapsedTime/60/60)
                      
                    ))
                    remainingTime <- ((total.rows*elapsedTime)/actual) - elapsedTime
                    print(paste("estimated remaining time (seconds)", format2Print(remainingTime),
                                " - (minutes)", format2Print(remainingTime/60), 
                                " - (hours)", format2Print(remainingTime/60/60)
                    ))
                  }
                }
            }
    }
    if(inc.progress) {
      incProgress(1/nrow(mrna));
    }	
  }
  
  # deleting useless and unused rows
  res <- res[c(1:actual.n.correlated-1),c(1:5)]
  
  #if (!(folder.exists(output.path))) {dir.create(output.path)}
  file.path<-paste(output.path, output.file.name, sep="")
  write.table(res, file.path, sep="\t", row.names=FALSE, 
              col.names=TRUE, quote=FALSE)
  print(proc.time() - ptm)
  
  return (convertVectorToMatrix(res))
  
}



get.cgs.for.this.genes <- function(actual.gen,meth.platform){
    as.character((meth.platform[meth.platform$gene == actual.gen,])[,2])
} 



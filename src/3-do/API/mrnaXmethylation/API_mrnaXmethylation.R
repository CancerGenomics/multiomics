#The methylation cgXXX is in methylation files of TCGA.
#We have also a file which correlates cgXXX with gene. Take into account that a cgXXX can do methilation over multiple genes and that a gene can be methylazed by multiple cgXXX.
#For doing it, we read a gen from the mrnaFile, then we look in the "methylationPlatformFile" all the cgXXX corresponding to this gene, and then we look for those CGXXX in the methilation file. It keeps the pairs (gene, cgXXX) which correlates negatively, considering the mrna file and the methylation file.
#Take into account that the platform file should Have two columns: GEN - CG. If a gen has goto multiple cgs, there will be many rows for this gene

methXMrnas <- function(mrna, meth, meth.platform, output.path="~/", 
                       output.file.name="methylationXMrna.csv",
                       r.minimium=0.7, 
                       pearsons.method = "pearson", 
                       inc.progress = F,keep.pos.cor=T, keep.neg.cor=F){
  
  ptm <- proc.time()
  total.rows=nrow(mrna)*15
  print(paste("Running pipeline methylaion_X_mrnas with", r.minimium, 
              "threshold and pearson's method:", pearsons.method, sep=" "))
  
  
  ###MDB: 26/2/2018 - P.ADJUST
  #Columns are: "Gene","Location", "methylation-id", "Methylation-mRNA correlation", "p-value", "p_value_fdr_adjustedMeth_Mrna_Correlation", "ID"
  num.of.result.columns<-7
  position.of.adjusted.p.value<-6
  ####
  
  
  # The result matix is created
  res <- matrix(nrow=total.rows,ncol=num.of.result.columns)
  colnames(res)<-(c("Gene","Location", "methylation-id", "Methylation-mRNA correlation", "p-value", "p-value-fdr-adjusted", "ID"))
  
  ###MDB: 26/2/2018 - P.ADJUST - Start on 0
  actual<-0
  actual.n.correlated<-1
  print("Start process!")
  
  ###MDB: 26/2/2018 - P.ADJUST
  p.values.all<-c()
  p.values.positions.of.correlated.pairs<-c()
  ids<-c()
  ###
  
  ids.in.meth.dataset<-meth[,1]
  print("Start process!")
  for (i in 1:nrow(mrna)) {
    
    actual.gen<-as.character(mrna[i,1])
    #position.in.cnv.dataset<-which(ids.in.cnv.dataset == actual.gen)
    #In meths i have each cg, which do methialtion over the gene.
    cgs.for.this.gene<-get.cgs.for.this.genes(actual.gen,meth.platform)
    #Se queda con el primer CNV
    
    if (length(cgs.for.this.gene)>0){
      
      for (actual.cg.index in 1:length(cgs.for.this.gene)){
        actual<-actual+1
        actual.cg<-cgs.for.this.gene[actual.cg.index]
        #position.in.cnv.dataset<-position.in.cnv.dataset[1]
        actual.mrna<-mrna[i,2:ncol(mrna)]
        #position.in.meth.dataset<-which(ids.in.meth.dataset == actual.cg)
        actual.meth<-meth[meth[,1] == actual.cg,2:ncol(meth)]
        if (nrow(actual.meth)>0){
          resultado.pearson<-cor.test(as.numeric(actual.mrna),
                                      as.numeric(actual.meth), 
                                      method = pearsons.method)
          
          ###MDB: 26/2/2018 - P.ADJUST
          p.values.all<-append(p.values.all, resultado.pearson$p.value)
          #id<-paste(actual, actual.gen, actual.mirna, sep="-")
          id<-actual
          ids<-append(ids, id)
   
          
          if (!is.na(abs(resultado.pearson$estimate))) {

            if (abs(resultado.pearson$estimate) > r.minimium) {

              if ((keep.pos.cor==T && resultado.pearson$estimate>0) || ((keep.neg.cor==T && resultado.pearson$estimate<0))){ 
                location<-getGeneLocation(actual.gen);
                newValue<-c(as.character(actual.gen), location, actual.cg,
                            resultado.pearson$estimate, resultado.pearson$p.value, -9999, id)
                res[actual.n.correlated,1:num.of.result.columns] <- newValue
                actual.n.correlated<-actual.n.correlated+1
                ###MDB: 26/2/2018 - P.ADJUST
                p.values.positions.of.correlated.pairs<-append(p.values.positions.of.correlated.pairs, id)
              }
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

  if (length(p.values.all)>0){
    ###MDB: 26/2/2018 - P.ADJUST
    p.values.adjusted.fdr<-p.adjust(p.values.all, method="fdr", n=length(p.values.all))
    names(p.values.adjusted.fdr)<-ids
    
    ###MDB: 26/2/2018 - P.ADJUST
    res[res[,"ID"] %in% p.values.positions.of.correlated.pairs, position.of.adjusted.p.value]<-p.values.adjusted.fdr[as.character(p.values.positions.of.correlated.pairs)]
    ####
  }
  
  # deleting useless and unused rows
  res <- res[c(1:actual.n.correlated-1),c(1:num.of.result.columns)]
  
  #if (!(folder.exists(output.path))) {dir.create(output.path)}
  file.path<-paste(output.path, output.file.name, sep="")
  write.table(res, file.path, sep="\t", row.names=FALSE, 
              col.names=TRUE, quote=FALSE)
  print(proc.time() - ptm)
  
  return (convertVectorToMatrix(res))
  
}

methXMrnasWCGNA <- function(mrna, meth, meth.platform, output.path="~/", 
                            output.file.name="methylationXMrna.csv",
                            r.minimium=0.7,
                            inc.progress = F, keep.pos.cor=F, keep.neg.cor=T){
  
  library("WGCNA")
  library("reshape2")
  library("data.table")
  
  ptm <- proc.time()
  print(paste("Running pipeline methylaion_X_mrnas with", r.minimium, "threshold", sep=" "))
  
  final.data.frame <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(final.data.frame) <- c("x", "y", "correlation", "p.value", "p.value.fdr.adjusted")
  
  # For each gen on MRNA file
  for (i in 1:nrow(mrna)) {
    # Build a dataframe with actual gen row only
    actual.mrna<-mrna[i,2:ncol(mrna)]
    actual.gen<-as.character(mrna[i,1])
    
    # Get a list of all CG associated to the actual gen using the file meth.platform
    cg.genes <- subset(meth.platform, gene == actual.gen)
    
    if (nrow(cg.genes) > 0) {
      # Create dataframe with all files on meth file which key is a value on cg.genes
      current.gen.meth.values.by.cg <- subset(meth, meth %in% cg.genes$cg)
      if (nrow(current.gen.meth.values.by.cg) > 0) {
        # set row names using first row values, then remove first column
        rownames(current.gen.meth.values.by.cg) <- current.gen.meth.values.by.cg[,1]
        current.gen.meth.values.by.cg <- current.gen.meth.values.by.cg[,2:ncol(current.gen.meth.values.by.cg)]
        
        print("Correlation....")
        # calcultate correlation using wcgna
        correlation.result <-correlation.with.wcgna(actual.mrna, current.gen.meth.values.by.cg,r.minimium, keep.pos.cor=keep.pos.cor, keep.neg.cor=keep.neg.cor)
        # colnames(correlation.result)<-(c("Gene","Location", "CNV_mRNA_Correlation", "p-value", "p_value_fdr_adjusted"))
        final.data.frame <- rbind(final.data.frame, correlation.result)
      }
    }
  }
  
  colnames(final.data.frame)<-(c("Gene", "methylation-id", "Methylation_mRNA_correlation", "p-value", "p_value_fdr_adjusted"))
  
  # For each row, calculate methylation-id position
  if (nrow(final.data.frame) > 0) {
    final.data.frame$Location <- apply(final.data.frame[,1], 1, getGeneLocationFromFactor) 
    final.data.frame <- final.data.frame[,c("Gene","Location", "methylation-id", "Methylation_mRNA_correlation", "p-value", "p_value_fdr_adjusted")]
  }
  
  # Write the result to a file
  write.to.file(final.data.frame, output.path, output.file.name)
  
  print(proc.time() - ptm)
  
  return (as.matrix(final.data.frame))
}


get.cgs.for.this.genes <- function(actual.gen,meth.platform){
  as.character((meth.platform[meth.platform$gene == actual.gen,])[,2])
} 


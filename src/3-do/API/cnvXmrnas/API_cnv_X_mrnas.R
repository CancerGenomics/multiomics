#If the CNV has got repeated genes, it will use the first row in the order.
#It will compare all the genes that will be in both files. That is why the fdr is applied over a vector of size = number of elements in common between both files.
CnvXMrnas <- function(mrna, cnv, output.path="~/", 
                                                     output.file.name="cnvXMrna.csv",
                                                     r.minimium=0.7, 
                                                     pearsons.method = "pearson", 
                                                     inc.progress = F){
  
  ptm <- proc.time()
  total.rows=nrow(mrna)
  print(paste("Running pipeline CNv_X_mrnas with", r.minimium, 
              "threshold and pearson's method:", pearsons.method, sep=" "))
  
  ###MDB: 27/2/2018 - P.ADJUST
  #Columns are: "Gene","Location", "CNV-mRNA correlation", "p-value", "p_value_fdr_adjusted_CNV-mRNA_Correlation", "ID"
  num.of.result.columns<-6
  position.of.adjusted.p.value<-5
  ####
  
  ###MDB: 27/2/2018 - P.ADJUST - Start on 0
  actual<-0
  actual.n.correlated<-1
  print("Start process!")
  
  ###MDB: 27/2/2018 - P.ADJUST
  p.values.all<-c()
  p.values.positions.of.correlated.pairs<-c()
  ids<-c()
  ###
  
  # The result matix is created
  res <- matrix(nrow=total.rows,ncol=num.of.result.columns)
  colnames(res)<-(c("Gene","Location", "CNV-mRNA correlation", "p-value", "p_value_fdr_adjusted", "ID"))
  

  ids.in.cnv.dataset<-cnv[,1]
  
  print("Start process!")
  for (i in 1:nrow(mrna)) {
    actual<-actual+1
    actual.gen<-as.character(mrna[i,1])
    position.in.cnv.dataset<-which(ids.in.cnv.dataset == actual.gen)
    #Se queda con el primer CNV
    if (length(position.in.cnv.dataset)>0){
        position.in.cnv.dataset<-position.in.cnv.dataset[1]
        actual.mrna<-mrna[i,2:ncol(mrna)]
        actual.cnv<-cnv[position.in.cnv.dataset,2:ncol(cnv)]
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
        resultado.pearson<-cor.test(as.numeric(actual.mrna),
                                    as.numeric(actual.cnv), 
                                    method = pearsons.method)
        
        ###MDB: 27/2/2018 - P.ADJUST
        p.values.all<-append(p.values.all, resultado.pearson$p.value)
        id<-actual
        ids<-append(ids, id)
        
        if (!is.na(abs(resultado.pearson$estimate))) {
          if (abs(resultado.pearson$estimate) > r.minimium) {
            location<-getGeneLocation(actual.gen);
            newValue<-c(as.character(actual.gen), location,
                        resultado.pearson$estimate, resultado.pearson$p.value, -9999, id)
            ###MDB: 26/2/2018 - P.ADJUST
            res[actual.n.correlated,1:num.of.result.columns] <- newValue
            actual.n.correlated<-actual.n.correlated+1
            
            ###MDB: 26/2/2018 - P.ADJUST
            p.values.positions.of.correlated.pairs<-append(p.values.positions.of.correlated.pairs, id)
          }
        }
    }
    if(inc.progress) {
      incProgress(1/nrow(cnv));
    }	
  }
  
  ###MDB: 26/2/2018 - P.ADJUST
  p.values.adjusted.fdr<-p.adjust(p.values.all, method="fdr", n=length(p.values.all))
  names(p.values.adjusted.fdr)<-ids
  
  ###MDB: 26/2/2018 - P.ADJUST
  res[res[,"ID"] %in% p.values.positions.of.correlated.pairs, position.of.adjusted.p.value]<-p.values.adjusted.fdr[as.character(p.values.positions.of.correlated.pairs)]
  ####
  
  ###MDB: 27/2/2018 - P.ADJUST
  # deleting useless and unused rows
  res <- res[c(1:actual.n.correlated-1),c(1:num.of.result.columns)]
  
  #if (!(folder.exists(output.path))) {dir.create(output.path)}
  file.path<-paste(output.path, output.file.name, sep="")
  write.table(res, file.path, sep="\t", row.names=FALSE, 
              col.names=TRUE, quote=FALSE)
  print(proc.time() - ptm)
  
  return (convertVectorToMatrix(res))
  
}

CnvXMrnasWCGNA <- function(expression, cnv, output.path="~/", 
                      output.file.name="cnvXMrna.csv",
                      r.minimium=0.7, 
                      pearsons.method = "pearson", 
                      inc.progress = F){

  ###MDB: 26/3/2018
  library("WGCNA")
  library("reshape2")
  library("data.table")
  
  ####Number of rows to evaluate (number of mrnas * number of mirnas)
  ptm <- proc.time()
  print(paste("Running pipeline with", r.minimium,"threshold and pearson's method:", pearsons.method, sep=" "))
  
  #Organize the matrix to keep it in the following format
  #Mirna or mrna as columns, samples as rows. mirna names o mrna names as column names, and sample names as row names.
  row.names(expression)<-expression[,1]
  row.names(cnv)<-cnv[,1]
  expression<-expression[,2:ncol(expression)]
  mirna<-cnv[,2:ncol(cnv)]
  
  # calcultate correlation using wcgna
  correlation.result <-correlation.with.wcgna(expression, mirna,r.minimium)
  colnames(correlation.result)<-(c("Gene","Location", "CNV_mRNA_Correlation", "p-value", "p_value_fdr_adjusted"))
  
  # Apply function to set Gene Location on resulting dataframe. This is specific to this Process
  correlation.result[,2] <- apply(correlation.result[,2], 1, getGeneLocationFromFactor)
  
  # Write the result to a file
  write.to.file(correlation.result, output.path, output.file.name)
  
  print(proc.time() - ptm)
  
  return (as.matrix(correlation.result))
  
}





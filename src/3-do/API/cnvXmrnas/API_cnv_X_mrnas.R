#If the CNV has got repeated genes, it will use the first row in the order.

CnvXMrnas <- function(mrna, cnv, output.path="~/", 
                                                     output.file.name="cnvXMrna.csv",
                                                     r.minimium=0.7, 
                                                     pearsons.method = "pearson", 
                                                     inc.progress = F){
  
  ptm <- proc.time()
  total.rows=nrow(mrna)
  print(paste("Running pipeline CNv_X_mrnas with", r.minimium, 
              "threshold and pearson's method:", pearsons.method, sep=" "))
  
  # The result matix is created
  res <- matrix(nrow=total.rows,ncol=3)
  colnames(res)<-(c("Gene symbol","CNV-Mrna correlation", "p_value of correlation"))
  
  actual<-1
  actual.n.correlated<-1
  
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
        if (!is.na(abs(resultado.pearson$estimate))) {
          if (abs(resultado.pearson$estimate) > r.minimium) {
            newValue<-c(as.character(actual.gen),
                        resultado.pearson$estimate, resultado.pearson$p.value)
            res[actual.n.correlated,1:3] <- newValue
            actual.n.correlated<-actual.n.correlated+1
          }
        }
    }
    if(inc.progress) {
      incProgress(1/nrow(cnv));
    }	
  }
  
  # deleting useless and unused rows
  res <- res[c(1:actual.n.correlated-1),c(1:3)]
  
  #if (!(folder.exists(output.path))) {dir.create(output.path)}
  file.path<-paste(output.path, output.file.name, sep="")
  print(file.path)
  write.table(res, file.path, sep="\t", row.names=FALSE, 
              col.names=TRUE, quote=FALSE)
  print(proc.time() - ptm)
  
  return (res)
  
    
}

#Revisar porque tiene bugs. No funciona correctamente.
removeRowsWithLowSD<-function(input, minimumsd, outputpath){
  result_temp<-matrix(nrow = nrow(input), ncol = ncol(input)-1)
  result<-matrix(nrow = nrow(input), ncol = ncol(input))

  genesWithHighSD<-c(1:nrow(input))

  cantgenes<-1

  colnames<-as.character(input[,1])
  
  for (i in 1:nrow(input)) {
    sd<-sd(input[i,2:ncol(input)])
    if (sd>minimumsd) {
        result_temp[cantgenes,] <- as.numeric(input[i,2:ncol(input)])
        genesWithHighSD[cantgenes] <- colnames[i]
        cantgenes=cantgenes+1
      }
    print(i)
  }
  
  result <- cbind(genesWithHighSD[1:cantgenes-1],result_temp[1:cantgenes-1,])
  write.table(file = outputpath, x = result, append = F,quote = F, row.names = F, col.names = F, sep="\t")
}
  

#Revisar porque tiene bugs. No funciona correctamente.
DivideMirnaFromMRNAJustHighSD<-function(input, minimumsd, outputpath){

  result_gene<-matrix(nrow = nrow(input), ncol = ncol(input))
  result_mirna<-matrix(nrow = nrow(input), ncol = ncol(input))
  
  genenames<-c(1:nrow(input))
  mirnanames<-c(1:nrow(input))
  tmp_result_gene<-matrix(nrow = nrow(input), ncol = (ncol(input)-1))
  tmp_result_mirna<-matrix(nrow = nrow(input), ncol = (ncol(input)-1))
  
  
  cantmirnas<-1
  cantgenes<-1
  colnames<-as.character(input$col1)
  for (i in 1:nrow(input)) {
    #row<-unlist(class(input[i,]))
    sd<-sd(input[i,2:ncol(input)])
    if (sd>minimumsd) {
      if (grepl('^mir', tolower(colnames[i]))) {
        #result_mirna[cantmirnas,]<-append(colnames[i], as.numeric(input[i,2:ncol(input)]))
        tmp_result_mirna[cantmirnas,] <- as.numeric(input[i,2:ncol(input)])
        mirnanames[cantmirnas] <- colnames[i]
        cantmirnas=cantmirnas+1
      }
      else {
        #result_gene[cantgenes,]<-append(colnames[i], as.numeric(input[i,2:ncol(input)]))
        tmp_result_gene[cantgenes,] <- as.numeric(input[i,2:ncol(input)])
        genenames[cantgenes] <- colnames[i]
        cantgenes=cantgenes+1
      }
    }
    print(i)
  }
  
  print("tmp_result_mirna")
  print(tmp_result_mirna)
  print("tmp_result_gene")
  print(tmp_result_gene)
  print("mirnanames")
  print(mirnanames)
  print("genenames")
  print(genenames)
  
  
  result_mirna <- append(mirnanames[1:cantmirnas-1],tmp_result_mirna[1:cantmirnas-1,])
  print("result_mirna")
  print(result_mirna)
  
  result_gene <- append(genenames[1:cantgenes-1],tmp_result_gene[1:cantgenes-1,])
  print("result_gene")
  print(result_gene)
  
  #row.names(result)<-row.names(input)
  #colnames(result_mirna)<-colnames(input)
  #colnames(result_gene)<-colnames(input)
  
  #Agregarle la columna de genes
  
  #if (cantgenes-1==1) result_gene<-t(as.matrix((result_gene[1:cantgenes-1,])))
  #else result_gene<-(result_gene[1:cantgenes-1,])
  
  #if (cantmirnas-1==1) result_mirna<-t(as.matrix((result_mirna[1:cantmirnas-1,])))
  #else if (cantmirnas-1>=2) result_mirna<-result_mirna[1:cantmirnas-1,]
  
  write.table(file = paste(outputpath, "outmirna.txt", sep=""), x = result_mirna, append = T,quote = F, row.names = F, col.names = F)
  write.table(file = paste(outputpath, "outgene.txt", sep=""), x = result_gene, append = T,quote = F, row.names = F, col.names = F)
}





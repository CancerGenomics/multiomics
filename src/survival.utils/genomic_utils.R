getGeneLocation <- function(gene.id){
  library(org.Hs.eg.db)
  cols = c("SYMBOL","GENENAME","MAP")

  tryCatch({
    annots <- select(org.Hs.eg.db,keys=gene.id,columns=cols,keytype="SYMBOL")  
    if (nrow(annots)>1) location<-annots[1,]$MAP
    else location<-annots$MAP
    
  },
  error = function(err){
    location<-"notFound"
  })
}


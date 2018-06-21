correlation.with.wcgna <- function(x, y, minimium, keep.pos.cor=T, keep.neg.cor=T) {
  
  ### Enable parallel processing for WCGNA Correlation
  enableWGCNAThreads()
  
  # Calculate correlation between x and y using  WCGNA
  correlation.start <- proc.time()
  # transpose matrix before correlation
  x.transposed <-t(x)
  y.transposed <- t(y)
  
  x.transposed.numeric<-apply(x.transposed, 2, as.numeric)
  y.transposed.numeric<-apply(y.transposed, 2, as.numeric)
  
  print("Running correlation using WCGNA.")
  cor.and.pvalue <- corAndPvalue(x.transposed.numeric, y.transposed.numeric)
  cat("CORRELATION TIME: : ", (proc.time() - correlation.start)["elapsed"], "\n")
  
  merge.start <- proc.time()
  # correlation result into a dataframe
  cor.melt<-melt(cor.and.pvalue$cor)
  colnames(cor.melt) <- c("x","y","correlation")
  # Filter rows with a correlation that does not meet the minimum required value
  cor.melt <- subset(cor.melt, abs(correlation) > minimium)
  
  if (!keep.pos.cor)   cor.melt <- subset(cor.melt, correlation < 0)
  if (!keep.neg.cor)   cor.melt <- subset(cor.melt, correlation > 0)
  
  #pvalue result into a dataframe
  p.melt<-melt(cor.and.pvalue$p)
  colnames(p.melt) <- c("x","y","p.value")
  
  padj.melt<-p.melt
  all.p.values<-p.melt[,3]
  padj.melt[,3] = p.adjust(all.p.values, length(all.p.values), method = "fdr")
  colnames(padj.melt) <- c("x","y","p.value.fdr.adjusted")
  
  # Merge all three previos tables into a single dataframe.
  # Transform data.frame's to data.table's to improve merge performance. Once the merge is done, transform back to data.frame,
  # as the caller code needs a data.frame to work
  temp.table <- merge(as.data.table(cor.melt), as.data.table(p.melt), by=c("x","y"))
  result.table <-merge(temp.table, as.data.table(padj.melt),  by=c("x","y"))
  cat("MERGE TIME: : ", (proc.time() - merge.start)["elapsed"], "\n")
  
  return(result.table)
}


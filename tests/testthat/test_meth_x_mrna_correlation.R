context("Meth x Mrna Correlation Test")

test_that("WCGNA correlation gives same results as previous implementation", {
  
  sourceBaseLocation <- gsub("tests/testthat", "", getwd(), ignore.case  = TRUE)
  source(paste(sourceBaseLocation, "/src/3-do/API/mrnaXmethylation/API_mrnaXmethylation.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_correlation_wcgna.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/survival.utils/methylation.platforms.utils.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/survival.utils/file_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
  
  working.path=paste(sourceBaseLocation, "/test/examples/methylation_X_mrnas/",sep="")
  mrna.dif.expr.path.file="mrna-with-extra-samples.csv"
  meth.file="meth.csv"
  my.predicted.cut.off=10
  the.output.path=tempdir()
  
  mrna.dif.expr.path<-paste(working.path, mrna.dif.expr.path.file, sep="")
  meth.dif.expr.path<-paste(working.path, meth.file, sep="")
  mrna.dif.expr <- readMrnaExpressionFile(mrna.dif.expr.path)
  meth.dif.expr <- readMethylationFile(meth.dif.expr.path)
  
  intersection<-keepSameColumns(mrna.dif.expr,meth.dif.expr)
  mrna.dif.expr<-(intersection[[1]])
  meth.dif.expr<-(intersection[[2]])
  
  path.platform<-paste(sourceBaseLocation, "/resources/methilation.platforms/illuminaMethyl450_hg19_GPL16304.txt", sep="")
  meth.platform.table = getMethylationPlatformTableForPipeline("HumanMethylation450 BeadChip", path.platform)
  
  # Test with a r.minimum greater than 0
  calculated.legacy <- methXMrnas(mrna.dif.expr, meth.dif.expr, meth.platform.table, output.path=the.output.path,
                                       output.file.name="cnvXMrna.csv",
                                       r.minimium=0.2, 
                                       pearsons.method = "pearson",
                                       keep.pos.cor=T, keep.neg.cor=T,
                                       inc.progress = F)
  
  calculated.wcgna <- methXMrnasWCGNA(mrna.dif.expr, meth.dif.expr, meth.platform.table, output.path=the.output.path,
                                       output.file.name="cnvXMrna.csv",
                                       r.minimium=0.2, 
                                       keep.pos.cor=T, keep.neg.cor=T,
                                       inc.progress = F)
  
  # Remove ID column before comparing and transform matrix to data.frame
  calculated.legacy <- calculated.legacy[,1:6]
  
  # Validate data frame sizes are equal
  expect_equal(dim(calculated.legacy),dim(calculated.wcgna))
  
  # Compare First Row
  compare.rows.meth(calculated.wcgna, calculated.legacy, 1)
  # Compare Second Row
  compare.rows.meth(calculated.wcgna, calculated.legacy, 2)
  # Compare Third Row
  compare.rows.meth(calculated.wcgna, calculated.legacy, 3)
  
})

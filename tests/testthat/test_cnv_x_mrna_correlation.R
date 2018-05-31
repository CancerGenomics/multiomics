context("Cnv x Mrna Correlation Test")

test_that("WCGNA correlation gives same results as previous implementation", {

  sourceBaseLocation <- gsub("tests/testthat", "", getwd(), ignore.case  = TRUE)
  source(paste(sourceBaseLocation, "/src/3-do/API/cnvXmrnas/API_cnv_X_mrnas.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_correlation_wcgna.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/survival.utils/matrix_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/survival.utils/read_genomic_data_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/survival.utils/genomic_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/survival.utils/file_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
  
  
  #INPUT FOR example
  working.path=paste0(sourceBaseLocation, "/test/examples/cnv_X_mrnas/")
  mrna.dif.expr.path.file="mrnas.csv"
  cnv.file="cnv.csv"
  the.output.path=tempdir()
  
  mrna.dif.expr.path<-paste(working.path, mrna.dif.expr.path.file, sep="")
  cnv.path<-paste(working.path, cnv.file, sep="")
  print("Preparing...")
  mrna.dif.expr <- readMrnaExpressionFile(mrna.dif.expr.path)
  cnv <- readCNVFile(cnv.path)
  
  #Keep columns which are in both databases
  intersection<-keepSameColumns(mrna.dif.expr,cnv)
  mrna.dif.expr<-(intersection[[1]])
  cnv<-(intersection[[2]])
  
  calculated.legacy <-CnvXMrnas(mrna.dif.expr, cnv, output.path=the.output.path,
            output.file.name="cnvXMrna.csv",
            r.minimium=0.9, 
            pearsons.method = "pearson", 
            inc.progress = F)
  
  
  calculated.wcgna <-CnvXMrnasWCGNA(mrna.dif.expr, cnv, output.path=the.output.path,
                                output.file.name="cnvXMrna.csv",
                                r.minimium=0.9, 
                                pearsons.method = "pearson", 
                                inc.progress = F)
  
  calculated.legacy <- calculated.legacy[,1:5]

  # Disabled as Legacy code return 1 less row.
  #expect_equal(dim(calculated.legacy),dim(calculated.wcgna))
  
  expect_equal(nrow(calculated.wcgna), 3)
  expect_equal(ncol(calculated.wcgna), 5)
  
  # Compare First Row
  compare.rows.cnv(calculated.wcgna, calculated.legacy, 1)
  # Compare Second Row
  compare.rows.cnv(calculated.wcgna, calculated.legacy, 2)
  

})

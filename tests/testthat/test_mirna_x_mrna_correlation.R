context("Mirna x Mrna Correlation Test")

test_that("WCGNA correlation gives same results as previous implementation", {
  
  # Required packages for testing
  #install.packages("testthat")
  #install.packages("roxygen2")
  
  sourceBaseLocation <- gsub("tests/testthat", "", getwd(), ignore.case  = TRUE)
  source(paste(sourceBaseLocation, "/src/3-do/API/mirnasRegulatingMrnas/API_multiomics_for_finding_mirnas_regulating_mrnas.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_correlation_wcgna.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/survival.utils/matrix_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/survival.utils/read_genomic_data_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
  source(paste(sourceBaseLocation, "/src/survival.utils/file_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
  
  
  working.path=paste(sourceBaseLocation, "/test/examples/miRnas_regulating_mRnas/",sep="")
  mrna.dif.expr.path.file="mrnas.csv"
  mirna.dif.expr.path.file="mirnas.csv"
  my.predicted.cut.off=10
  maturemirna.x.mrna.correlation.file.name="inputStep2-matureMirnaXmrna.csv"
  just.betters.maturemirna.X.mrna.considering.mirna.databases="inputStep3_justBettersMirnaXmRNA_ConsideringMirnaDatabases.csv"

  mrna.dif.expr.path<-paste(working.path, mrna.dif.expr.path.file, sep="")
  mirna.dif.expr.path<-paste(working.path, mirna.dif.expr.path.file, sep="")
  mrna.dif.expr <- readMrnaExpressionFile(mrna.dif.expr.path)
  mirna.dif.expr <- readMirnaExpressionFile(mirna.dif.expr.path)
  
  
  # Test with a r.minimum equal to 0
  calculated.wcgna <- CalculateCorrelationsMirnaMrnaUsingWCGNA(mrna.dif.expr,mirna.dif.expr, working.path, 
                                                         output.file.name=maturemirna.x.mrna.correlation.file.name,
                                                         r.minimium=0.0)

  calculated.legacy <- CalculateCorrelationsMirnaMrna(mrna.dif.expr,mirna.dif.expr, working.path, 
                                output.file.name=maturemirna.x.mrna.correlation.file.name,
                               r.minimium=0.0)
  
  # Remove ID column before comparing and transform matrix to data.frame
  calculated.legacy <- calculated.legacy[,1:5]
  calculated.legacy <- as.data.frame(calculated.legacy)
  
  # Validate data frame sizes are equal
  expect_equal(dim(calculated.legacy),dim(calculated.wcgna))
  
  
  # Test with a r.minimum greater than 0
  calculated.wcgna <- CalculateCorrelationsMirnaMrnaUsingWCGNA(mrna.dif.expr,mirna.dif.expr, working.path, 
                                                                  output.file.name=maturemirna.x.mrna.correlation.file.name,
                                                                  r.minimium=0.8)
  
  calculated.legacy <- CalculateCorrelationsMirnaMrna(mrna.dif.expr,mirna.dif.expr, working.path, 
                                                     output.file.name=maturemirna.x.mrna.correlation.file.name,
                                                     r.minimium=0.8)
  
  # Remove ID column before comparing and transform matrix to data.frame
  calculated.legacy <- calculated.legacy[,1:5]

  # Validate data frame sizes are equal
  expect_equal(dim(calculated.legacy),dim(calculated.wcgna))
  
  # Compare First Row
  compare.rows.mirna(calculated.wcgna, calculated.legacy, 1)
  # Compare Second Row
  compare.rows.mirna(calculated.wcgna, calculated.legacy, 2)
  
})

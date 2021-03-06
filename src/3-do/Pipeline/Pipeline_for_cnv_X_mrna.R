# path adonde estan los archivos de entrada, se puede cambiar por un path absoluto
sourceBaseLocation=getwd()

source(paste(sourceBaseLocation, "/src/2-load/load_multiomics.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_private_data_validation.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_private_multimir_interaction.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_correlation_wcgna.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/mirnasRegulatingMrnas/API_multiomics_for_finding_mirnas_regulating_mrnas.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/survivalGeneByGene/API_multiomics_for_survival_gene_by_gene.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/survivalGeneByGene/API_multiomics_for_survival_gene_by_gene.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/cnvXmrnas/API_cnv_X_mrnas.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/matrix_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/read_genomic_data_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/genomic_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/file_utils.R",sep=""), echo=FALSE, encoding="Cp1252")



###########CONFIG#################
#INPUT FOR example
working.path=paste0(sourceBaseLocation, "/test/examples/cnv_X_mrnas/")
mrna.dif.expr.path.file="mrnas-with-extra-samples.csv"
cnv.file="cnv.csv"
working.path=paste(sourceBaseLocation, "/test/examples/cnv_X_mrnas/",sep="")

#YOUR INPUT
#working.path="D:\\matias\\academia\\investigacion\\medicina personalizada\\8-DatosGenomica\\2016-09-13---paper multiomics\\CNV\\"
#mrna.dif.expr.path.file="mRNA_breast_750s_mostvariables.csv"
#cnv.file="CNV_breast_750s.csv"
#the.output.path="D:\\matias\\academia\\investigacion\\medicina personalizada\\8-DatosGenomica\\2016-09-13---paper multiomics\\CNV\\"


#####DO  
mrna.dif.expr.path<-paste(working.path, mrna.dif.expr.path.file, sep="")
cnv.path<-paste(working.path, cnv.file, sep="")
print("Preparing...")

##use this if you want to test with files
#mrna.dif.expr.path<-"C:\\Users\\matia\\Desktop\\temp2\\Pathway_Paradigm_mRNA"
#cnv.path<-"C:\\Users\\matia\\Desktop\\temp2\\Gistic2_CopyNumber_Gistic2_all_data_by_genes_small"

mrna.dif.expr <- readMrnaExpressionFile(mrna.dif.expr.path)
cnv <- readCNVFile(cnv.path)



#Keep columns which are in both databases
intersection<-keepSameColumns(mrna.dif.expr,cnv)
mrna.dif.expr<-(intersection[[1]])
cnv<-(intersection[[2]])


correlation <- CnvXMrnas(mrna.dif.expr, cnv, output.path=working.path,
                      output.file.name="cnvXMrna.csv",
                      r.minimium=0.9, 
                      pearsons.method = "pearson", 
                      inc.progress = F)


correlation2 <- CnvXMrnasWCGNA(mrna.dif.expr, cnv, output.path=working.path,
          output.file.name="cnvXMrna.csv",
          r.minimium=0.9, 
          pearsons.method = "pearson", 
          inc.progress = F)

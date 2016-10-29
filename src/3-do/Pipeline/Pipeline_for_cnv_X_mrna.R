# path adonde estan los archivos de entrada, se puede cambiar por un path absoluto
sourceBaseLocation=getwd()

source(paste(sourceBaseLocation, "/src/2-load/load_multiomics.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_private_data_validation.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_private_multimir_interaction.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/mirnasRegulatingMrnas/API_multiomics_for_finding_mirnas_regulating_mrnas.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/survivalGeneByGene/API_multiomics_for_survival_gene_by_gene.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/survivalGeneByGene/API_multiomics_for_survival_gene_by_gene.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/cnvXmrnas/API_cnv_X_mrnas.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/matrix_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/read_genomic_data_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/genomic_utils.R",sep=""), echo=FALSE, encoding="Cp1252")

###########CONFIG#################
#INPUT FOR example
#working.path=paste(sourceBaseLocation, "/examples/cnv_X_mrnas/",sep="")
#mrna.dif.expr.path.file="mrnas.csv"
#cnv.file="cnv.csv"
#the.output.path="d:\\temp\\"

#YOUR INPUT
working.path="D:\\matias\\academia\\investigacion\\medicina personalizada\\8-DatosGenomica\\2016-09-13---paper multiomics\\CNV\\"
mrna.dif.expr.path.file="mRNA_breast_750s_mostvariables.csv"
cnv.file="CNV_breast_750s.csv"
the.output.path="D:\\matias\\academia\\investigacion\\medicina personalizada\\8-DatosGenomica\\2016-09-13---paper multiomics\\CNV\\"


#####DO  
mrna.dif.expr.path<-paste(working.path, mrna.dif.expr.path.file, sep="")
cnv.path<-paste(working.path, cnv.file, sep="")
print("Preparing...")
mrna.dif.expr <- readMrnaExpressionFile(mrna.dif.expr.path)
cnv <- readCNVFile(cnv.path)
print("Checking if both files has the same samples in the same order...")

CnvXMrnas(mrna.dif.expr, cnv, output.path=the.output.path,
                      output.file.name="cnvXMrna-Spearman-0.5.csv",
                      r.minimium=0.5, 
                      pearsons.method = "spearman", 
                      inc.progress = F)



=======
# path adonde estan los archivos de entrada, se puede cambiar por un path absoluto
sourceBaseLocation=getwd()

source(paste(sourceBaseLocation, "/src/2-load/load_multiomics.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_private_data_validation.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_private_multimir_interaction.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/mirnasRegulatingMrnas/API_multiomics_for_finding_mirnas_regulating_mrnas.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/survivalGeneByGene/API_multiomics_for_survival_gene_by_gene.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/survivalGeneByGene/API_multiomics_for_survival_gene_by_gene.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/cnvXmrnas/API_cnv_X_mrnas.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/matrix_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/read_genomic_data_utils.R",sep=""), echo=FALSE, encoding="Cp1252")

###########CONFIG#################
#INPUT FOR example
working.path=paste(sourceBaseLocation, "/examples/cnv_X_mrnas/",sep="")
mrna.dif.expr.path.file="mrnas.csv"
cnv.file="cnv.csv"
the.output.path="c:\\_tmp\\"

#YOUR INPUT
#working.path="D:\\matias\\academia\\investigacion\\medicina personalizada\\8-DatosGenomica\\2016-09-13---paper multiomics\\CNV\\"
#mrna.dif.expr.path.file="mRNA_breast_750s_mostvariables.csv"
#cnv.file="CNV_breast_750s.csv"
#the.output.path="D:\\matias\\academia\\investigacion\\medicina personalizada\\8-DatosGenomica\\2016-09-13---paper multiomics\\CNV\\"


#####DO  
mrna.dif.expr.path<-paste(working.path, mrna.dif.expr.path.file, sep="")
cnv.path<-paste(working.path, cnv.file, sep="")
print("Preparing...")
mrna.dif.expr <- readMrnaExpressionFile(mrna.dif.expr.path)
cnv <- readCNVFile(cnv.path)

CnvXMrnas(mrna.dif.expr, cnv, output.path=the.output.path,
                      output.file.name="cnvXMrna.csv",
                      r.minimium=0.2, 
                      pearsons.method = "pearson", 
                      inc.progress = F)




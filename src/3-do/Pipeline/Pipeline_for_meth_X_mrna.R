# path adonde estan los archivos de entrada, se puede cambiar por un path absoluto
sourceBaseLocation=getwd()

source(paste(sourceBaseLocation, "/src/2-load/load_multiomics.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_private_data_validation.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/Private/multiomics_private_multimir_interaction.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/mirnasRegulatingMrnas/API_multiomics_for_finding_mirnas_regulating_mrnas.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/survivalGeneByGene/API_multiomics_for_survival_gene_by_gene.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/mrnaXmethylation/API_mrnaXmethylation.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/3-do/API/cnvXmrnas/API_cnv_X_mrnas.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/matrix_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/read_genomic_data_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/genomic_utils.R",sep=""), echo=FALSE, encoding="Cp1252")
source(paste(sourceBaseLocation, "/src/survival.utils/methylation.platforms.utils.R",sep=""), echo=FALSE, encoding="Cp1252")



###########CONFIG#################
#INPUT FOR example
working.path=paste0(sourceBaseLocation, "/test/examples/methylation_X_mrnas/")
mrna.dif.expr.path.file="mrna-with-extra-samples.csv"
meth.file="meth.csv"
the.output.path=tempdir()

#YOUR INPUT
#working.path="D:\\matias\\academia\\investigacion\\medicina personalizada\\8-DatosGenomica\\2016-09-13---paper multiomics\\CNV\\"
#mrna.dif.expr.path.file="mRNA_breast_750s_mostvariables.csv"
#cnv.file="CNV_breast_750s.csv"
#the.output.path="D:\\matias\\academia\\investigacion\\medicina personalizada\\8-DatosGenomica\\2016-09-13---paper multiomics\\CNV\\"


#####DO  
mrna.dif.expr.path<-paste(working.path, mrna.dif.expr.path.file, sep="")
meth.path<-paste(working.path, meth.file, sep="")
print("Preparing...")
mrna.dif.expr <- readMrnaExpressionFile(mrna.dif.expr.path)
meth <- readMethylationFile(meth.path)

#Keep columns which are in both databases
intersection<-keepSameColumns(mrna.dif.expr,meth)
mrna.dif.expr<-(intersection[[1]])
meth<-(intersection[[2]])

path.platform<-paste(sourceBaseLocation, "/resources/methilation.platforms/illuminaMethyl450_hg19_GPL16304.txt", sep="")

res<-methXMrnas(mrna.dif.expr, meth, getMethylationPlatformTableForPipeline("HumanMethylation450 BeadChip", path.platform), output.path=the.output.path,
          output.file.name="cnvXMrna.csv",
          r.minimium=0.2, 
          pearsons.method = "pearson", 
          inc.progress = F)

dim(res)


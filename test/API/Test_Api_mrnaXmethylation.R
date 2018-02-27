source(paste0(getwd(),"/src/3-do/API/mrnaXmethylation/API_mrnaXmethylation.R"))
source(paste0(getwd(),"/src/survival.utils/matrix_utils.R"))
source(paste0(getwd(),"/src/survival.utils/genomic_utils.R"))
source(paste0(getwd(),"/src/survival.utils/read_genomic_data_utils.R"))

do.test <- function(meth.platform){
  mrna<-readMrnaExpressionFile(paste0(getwd(),"/test/examples/methylation_X_mrnas/mrnas.csv"))
  meth<-readMethylationFile(paste0(getwd(),"/test/examples/methylation_X_mrnas/meth.csv"))
  methXMrnas(mrna, meth, meth.platform, output.path=tempdir(), 
             output.file.name="methylationXMrna.csv",
             r.minimium=0.3, 
             pearsons.method = "pearson", 
             inc.progress = F)
  
}

#Ejemplo con el file de la plataforma completo
#El resultado debe ser
    #Gene	Location	methylation-id	Methylation-mRNA correlation	p-value
    #HPSE	4q21.3	cg25428494	-0.711446450679181	1.6598986905973e-06
    #HPSE	4q21.3	cg13508557	-0.802838430583129	6.54656276791991e-09
    #AUTS2	7q11.22	cg17027195	-0.802838430583129	6.54656276791991e-09
meth.platform<-read.table(paste0(getwd(),"/resources/methilation.platforms/illuminaMethyl450_hg19_GPL16304.txt"), header = TRUE, sep="\t")
do.test(meth.platform)



  #Ejemplo con un file de methylation platform preparado para el test.
#El resultado debe ser:
    #Gene	Location	methylation-id	Methylation-mRNA correlation	p-value
    #HPSE	4q21.3	cg25428494	-0.711446450679181	1.6598986905973e-06
    #HPSE	4q21.3	cg13508557	-0.802838430583129	6.54656276791991e-09
    #AUTS2	7q11.22	cg17027195	-0.802838430583129	6.54656276791991e-09
    #AUTS2	7q11.22	cg13508557	-0.802838430583129	6.54656276791991e-09
meth.platform<-read.table(paste0(getwd(),"/test/examples/methylation_X_mrnas/cg-gene.txt"), header = TRUE)
do.test(meth.platform)



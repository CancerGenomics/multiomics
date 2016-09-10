# TODO: Add comment
# 
# Author: Matias
###############################################################################

files<-c()
API_multiomics_for_finding_mirnas_regulating_mrnas.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\API\\mirnasRegulatingMrnas\\API_multiomics_for_finding_mirnas_regulating_mrnas.R"
API_multiomics_for_survival_gene_by_gene_UsageSample.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\API\\survivalGeneByGene\\usageSamples\\API_multiomics_for_survival_gene_by_gene_UsageSample.R"
API_multiomics_for_survival_gene_by_gene.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\API\\survivalGeneByGene\\API_multiomics_for_survival_gene_by_gene.R"
API_multiomics_for_survival_multiple_genes.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\API\\survivalMultipleGenes\\API_multiomics_for_survival_multiple_genes.R"
API_Pam_50.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\API\\API_Pam_50.R"
Pipeline_for_survival_gene_by_gene.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\Pipeline\\Pipeline_for_survival_gene_by_gene.R"
Pipeline_multiomics_for_mirnas_regulating_mrnas.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\Pipeline\\Pipeline_multiomics_for_mirnas_regulating_mrnas.R"
Pipeline_Pam_50.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\Pipeline\\Pipeline_Pam_50.R"

subtypePrediction_distributed.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\Private\\bioclassifier\\subtypePrediction_distributed.R"
subtypePrediction_functions.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\Private\\bioclassifier\\subtypePrediction_functions.R"
multiomics_private_data_validation.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\Private\\multiomics_private_data_validation.R"
multiomics_private_multimir_interaction.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\3-do\\Private\\multiomics_private_multimir_interaction.R"

ConcordanceIndexEntity.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\survival.entities\\ConcordanceIndexEntity.R"
CoxphEntity.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\survival.entities\\CoxphEntity.R"
ExpressionXSurvivalEntity.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\survival.entities\\ExpressionXSurvivalEntity.R"
SurvdiffEntity.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\survival.entities\\SurvdiffEntity.R"
expression_grouping_functions.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\survival.utils\\expression_grouping_functions.R"
survival_utils.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\survival.utils\\survival_utils.R"

Test_expression_grouping_functions.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\test\\survival.utils\\Test_expression_grouping_functions.R"
Config_SurvivalGeneByGene_Test.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\test\\survivalGeneByGene\\Config_SurvivalGeneByGene_Test.R"
Test_SurvivalGeneByGene.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\test\\survivalGeneByGene\\Test_SurvivalGeneByGene.R"
Assertions.R<-"D:\\desarrollo\\workspaces\\R\\multiomics\\testing.utils\\Assertions.R"

files<-c()
files<-c(files, API_multiomics_for_finding_mirnas_regulating_mrnas.R)
files<-c(files, API_multiomics_for_finding_mirnas_regulating_mrnas.R)
#files<-c(files, API_multiomics_for_survival_gene_by_gene_UsageSample.R)
files<-c(files, API_multiomics_for_survival_gene_by_gene.R)
files<-c(files, API_multiomics_for_survival_multiple_genes.R)
files<-c(files, API_Pam_50.R)

#files<-c(files, Pipeline_for_survival_gene_by_gene.R)
#files<-c(files, Pipeline_multiomics_for_mirnas_regulating_mrnas.R)
#files<-c(files, Pipeline_Pam_50.R)

files<-c(files, subtypePrediction_distributed.R)
files<-c(files, subtypePrediction_functions.R)
files<-c(files, multiomics_private_data_validation.R)
files<-c(files, multiomics_private_multimir_interaction.R)

files<-c(files, ConcordanceIndexEntity.R)
files<-c(files, CoxphEntity.R)
files<-c(files, ExpressionXSurvivalEntity.R)
files<-c(files, SurvdiffEntity.R)
files<-c(files, expression_grouping_functions.R)
files<-c(files, survival_utils.R)

#files<-c(files, Test_expression_grouping_functions.R)
#files<-c(files, Config_SurvivalGeneByGene_Test.R)
#files<-c(files, Test_SurvivalGeneByGene.R)
#files<-c(files, Assertions.R)

package.skeleton(name = "multiomics",
		path = "d:\\temp", force = FALSE,
		code_files = files)


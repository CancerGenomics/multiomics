# source('directoryInput.R')
library("shiny")
library("DT")
library("shinyBS")
library("shinyjs")
library("multiMiR")
library("survival")

library("org.Hs.eg.db")
library("XenaR")

source("../survival.utils/matrix_utils.R")
source("../survival.utils/read_genomic_data_utils.R")
source("../survival.utils/genomic_utils.R")
source("../survival.utils/general_utils.R")
source("../survival.utils/expression_grouping_functions.R")
source("../survival.utils/survival_utils.R")

source("../survival.entities/ExpressionXSurvivalEntity.R")
source("../survival.entities/CoxphEntity.R")
source("../survival.entities/ValidationResult.R")
source("../survival.entities/SurvdiffEntity.R")

source("../3-do/Private/multiomics_private_data_validation.R")
source("../3-do/Private/multiomics_private_multimir_interaction.R")
source("../3-do/Private/multiomics_private_tcga.R")

source("../3-do/API/mirnasRegulatingMrnas/API_multiomics_for_finding_mirnas_regulating_mrnas.R")
source("../3-do/API/survivalGeneByGene/API_multiomics_for_survival_gene_by_gene.R")
source("../3-do/API/cnvXmrnas/API_cnv_X_mrnas.R")

source("../survival.utils/methylation.platforms.utils.R")
source("../3-do/API/mrnaXmethylation/API_mrnaXmethylation.R")


## ---------------------------
##
## Script name: all_in_one
##
## Purpose of script:
##
## Author: Diego Salazar
##
## Date Created: 2021-10-04
##
## Copyright (c) Diego Salazar, 2021
## Email: das4019@med.cornell.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------
## --------------------------------- ##
### ... Libraries and functions ... ###
## -------------------------------- ###
library(tmle)
library(caret)
library(dplyr)
library(purrr)
library(SuperLearner)
source("G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/give_O_GBM_and_LGG.R")

## ------------------------------------------- ##
### ... central_path and pathFeatureLabel ... ###
## ------------------------------------------- ##
if (identical(character(0), dir("D:/pathFeatureLabelKernels"))) {
  central_path <- "E:"
} else {
  central_path <- "D:"
}

pathFeatureLabel <- paste0(central_path, "/pathFeatureLabelKernels")

## ----------------------------------------- ##
### ... Load data and pre-Pross W matrix ... ###
## ---------------------------------------- ###
O <- give_O_GBM_LGG(pathFeatureLabel, dimReduc=FALSE, outTMLE3=TRUE, index=5, Y_as_NaN = FALSE)


Y <- O[[1]]; A <- O[[2]]; Delta <- O[[3]]; W <- O[[4]]; clinical <- O[[5]];
A <- A[!Y == 0]; Delta <- as.numeric(Delta[!Y == 0]); W <- W[!Y == 0,]; clinical <- clinical[!Y == 0,]; Y <- Y[!Y == 0];
Y.in.years <- Y * 12 / 365
delta.20.month <- ifelse(Y.in.years < 20 & Delta == 1, 1, 0)

## -------------------------- ##
### ... Propensity score ... ###
## -------------------------- ##
(mtry_seq = floor(sqrt(ncol(clinical)) * c(0.5, 1, 2)))
learners <- create.Learner("SL.ranger", params = list(num.trees = 1000, mtry = mtry_seq))
sl_rf <- CV.SuperLearner(Y = A, X = clinical, family = binomial(), cvControl = list(V = 10),
                     SL.library = c(learners$names,"SL.ranger"))
summary(sl_rf)
g1W <- sl_rf$SL.predict

## ---------------------------- ##
### ... Outcome prediction ... ###
## ---------------------------- ##
malacardGenes <- read.csv(file = paste0(pathFeatureLabel, "/co-mod_bases_datos_importantes/GeneCards-SearchResults_glioma.csv"))
malacardGenes <- paste0("mrna_", malacardGenes$Gene.Symbol[malacardGenes$Relevance.score > quantile(malacardGenes$Relevance.score, .75)])
# malacardGenes <- paste0("mrna_", malacardGenes$Gene.Symbol)

### choose a sample
W.noclinic <- W %>%
  dplyr::select(., -starts_with("clinical_"))

W.cnv <- W.noclinic %>%
  dplyr::select(., starts_with("cnv_"))
# W.cnv <- W.cnv[,sample(ncol(W.cnv), 20)]; dim(W.cnv)
W.mrna <- W.noclinic %>%
  dplyr::select(., starts_with("mrna_"))
W.mrna <- W.mrna[,colnames(W.mrna) %in% malacardGenes]
# W.mrna <- W.mrna[,sample(ncol(W.mrna), 20)]; dim(W.mrna)

W <- clinical %>%
  dplyr::bind_cols(W.cnv, W.mrna, delta.20.month)
W[1:5,1:20];dim(W)
# O <- data.frame(Y,A,delta.20.month,W); dim(O)
# colnames(O)[3] <- "delta"
# write.csv(O, file = paste0(pathFeatureLabel, "/co-md_records/O_Y_A_Delta_W.csv"))

(mtry_seq = floor(sqrt(ncol(W)) * c(0.5, 1, 2)))
learners <- create.Learner("SL.ranger", params = list(num.trees = 1000, mtry = mtry_seq))
sl_rf <- SuperLearner(Y = Y, X = W, family = gaussian(), cvControl = list(V = 10),
                         SL.library = c(learners$names,"SL.ranger"))

Q1W <- predict(sl_rf,newX=data.frame(A=1, W))$pred
Q0W <- predict(sl_rf,newX=data.frame(A=0, W))$pred

Q <- data.frame(Q0W,Q1W)

SL.library <- c("SL.xgboost")

result.Qcgc <- tmle(Y, A, W, Delta=delta.20.month, Q = Q, g1W = g1W,
                    Q.SL.library = SL.library, g.SL.library = SL.library)

### ----------------------- ###
# library(gtsummary)
# trial <- clinical %>% dplyr::bind_cols(data.frame(A,Delta))
# 
# trial %>% tbl_summary(by=A,
#                       statistic = list(all_continuous() ~ "{mean} ({sd})"),
#                       label = list(clinical_primary_diagnosis ~ "Primary diagnosis", 
#                                    clinical_age_at_index ~ "Age",
#                                    clinical_gender ~ "Gender",
#                                    clinical_IDH_codel_subtype ~ "IDH codel subtype",
#                                    clinical_MGMT_promoter_status ~ "MGMT promoter status",
#                                    clinical_ATRX_status ~ "ATRX status",
#                                    clinical_Chr7_gain_and_Chr10_loss ~ "Chr7 gain / Chr10 loss",
#                                    clinical_Karnofsky_Performance_Score ~ "Karnofsky",
#                                    clinical_Mutation_Count ~ "Mutation count",
#                                    Delta ~ "Censoring"
#                                    ),
#                       missing_text = "(Missing)") %>% 
#   bold_labels() %>% 
#   add_overall() %>% 
#   add_n() %>% 
#   add_p() %>% add_q()


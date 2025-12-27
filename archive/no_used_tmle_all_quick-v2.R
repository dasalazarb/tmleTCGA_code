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
## see: https://benkeser.github.io/sllecture/
## ------------------------------------ ##
### ... 0. Libraries and functions ... ###
## ----------------------------------- ###
set.seed(1)
library(tmle)
library(caret)
library(dplyr)
library(purrr)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(SuperLearner)
options(mc.cores = 4)
library(doParallel)
registerDoParallel(cores=4)

## ---------------------- ##
### ... 1. Load data ... ###
## ---------------------- ##
if (!identical(character(0), dir("G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/"))) {
  ruta <- "G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/"
} else {
  ruta <- "G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/"
}

source(paste0(ruta, "func_extractVars.R"))
source(paste0(ruta, "func_learners_v2.R"))
source(paste0(ruta, "func_SuperLearners.R"))
source(paste0(ruta, "tmle_ps_and_missmech.R"))

## ----------------------------------- ##
### ... 5. Estimation of E(Y|A,W) ... ###
## ----------------------------------- ##
## Todas las estrategias
sl__model1 <- CV.SuperLearner(Y = Y[A==1], X = W.new[A==1,], family = binomial(), control = list(saveFitLibrary = TRUE),
                              # For a real analysis we would use V = 10.
                              cvControl = list(V = 10, stratifyCV = TRUE), innerCvControl = list(list(V=10)),
                              parallel = "multicore", verbose = TRUE, method = "method.NNloglik",
                              SL.library = purrr::reduce(list(
                                # clinical info
                                list("SL.clinical.rf"), #list("SL.clinical.mars"),
                                # # early integration
                                list("SL.early.rf"), #list("SL.early.mars"),
                                # # # late integration
                                list("SL.Late.rf"), #list("SL.Late.mars"), # list("SL.LinearComb.logit"), list("SL.LinearComb.rf"),
                                # # # intermediate integration
                                list("SL.Int.rf.encoder"), #list("SL.Int.mars.encoder"), # list("SL.Int.logit.encoder"), list("SL.Int.rf.encoder"),
                                # # # # intermediate integration
                                list("SL.Int.rf.pca"), #list("SL.Int.mars.pca"), # list("SL.Int.logit.pca"), list("SL.Int.rf.pca"),
                                # # # # intermediate integration
                                list("SL.Int.rf.nmf")#, list("SL.Int.mars.nmf") # list("SL.Int.logit.nmf"), list("SL.Int.rf.nmf")
                              ), append))

# saveRDS(object = sl__model1, file = paste0(ruta, "os_lgg_data_var", nvar, "/tmle_sl__model1_var", nvar, ".rds"))
saveRDS(object = sl__model1, file = paste0("C:/Users/da.salazarb/Documents/tmle/tmle_sl__model1_var", nvar, ".rds"))

## Todas las estrategias
sl__model0 <- CV.SuperLearner(Y = Y[A==0], X = W.new[A==0,], family = binomial(), control = list(saveFitLibrary = TRUE),
                              # For a real analysis we would use V = 10.
                              cvControl = list(V = 10, stratifyCV = TRUE), innerCvControl = list(list(V=10)),
                              parallel = "multicore", verbose = TRUE, method = "method.NNloglik",
                              SL.library = purrr::reduce(list(
                                # clinical info
                                list("SL.clinical.rf"), #list("SL.clinical.mars"),
                                # # early integration
                                list("SL.early.rf"), #list("SL.early.mars"),
                                # # # late integration
                                list("SL.Late.rf"), #list("SL.Late.mars"), # list("SL.LinearComb.logit"), list("SL.LinearComb.rf"),
                                # # # intermediate integration
                                list("SL.Int.rf.encoder"), #list("SL.Int.mars.encoder"), # list("SL.Int.logit.encoder"), list("SL.Int.rf.encoder"),
                                # # # # intermediate integration
                                list("SL.Int.rf.pca"), #list("SL.Int.mars.pca"), # list("SL.Int.logit.pca"), list("SL.Int.rf.pca"),
                                # # # # intermediate integration
                                list("SL.Int.rf.nmf")#, list("SL.Int.mars.nmf") # list("SL.Int.logit.nmf"), list("SL.Int.rf.nmf")
                              ), append))

# saveRDS(object = sl__model0, file = paste0(ruta, "os_lgg_data_var", nvar, "/tmle_sl__model0_var", nvar, ".rds"))
saveRDS(object = sl__model0, file = paste0("C:/Users/da.salazarb/Documents/tmle/tmle_sl__model0_var", nvar, ".rds"))

## Quienes tienen nulo
df_null <- data.frame(matrix(0, nrow = length(sl__model1$AllSL), ncol = length(sl__model1$AllSL[[1]]$libraryNames)))
for (i in 1:length(sl__model0$AllSL)) {
  for (j in 1:length(sl__model0$AllSL[[i]]$libraryNames)) {
    df_null[i,j] <- is.null(sl__model0$AllSL[[i]]$fitLibrary[[j]]$object)
  }
}

# libs_final <- names(sl__model0$AllSL[[1]]$fitLibrary)[which(apply(df_null, 2, sum) == 0)]
## para q solo salgan estos, los demas no interesan porque se quiere comparar contra SL.clinical.model
libs_final <- c("SL.clinical.rf_All")

tmle_comparison <- data.frame(matrix(data = 0, nrow = 0, ncol = 32)); dim(tmle_comparison)


for (i in libs_final) {
  ## 
  Q1W <- apply(data.frame(lapply(sl__model1$AllSL, function(x) SuperLearner::predict.SuperLearner(x$fitLibrary[[i]]$object, newdata = W.new, type = "response")$pred)), 1, mean)
  Q0W <- apply(data.frame(lapply(sl__model0$AllSL, function(x) SuperLearner::predict.SuperLearner(x$fitLibrary[[i]]$object, newdata = W.new, type = "response")$pred)), 1, mean)
  
  # Q1W <- SuperLearner::predict.SuperLearner(sl__model1$AllSL$`1`$fitLibrary[[i]]$object, newdata = W.new, type = "response")$pred
  # Q0W <- SuperLearner::predict.SuperLearner(sl__model0$AllSL$`1`$fitLibrary[[i]]$object, newdata = W.new, type = "response")$pred
  Q <- data.frame(Q0W,Q1W)
  head(Q); tail(Q); dim(Q)
  
  ## ------------------------ ##
  ### ... g. tmle update ... ###
  ## ------------------------ ##
  result.Qcgc <- tmle(Y = Y, A = A, W = W.new, verbose = TRUE, Delta = Delta,
                      Q = Q, g1W = g1W, pDelta1 = pDelta1, family = "binomial")
  
  ATE <- data.frame(result.Qcgc$estimates[2]); ATE[is.na(ATE)] <- 0; ATE <- ATE %>% tidyr::pivot_wider(names_from = ATE.CI, values_from = ATE.CI)
  RR <- data.frame(result.Qcgc$estimates[3]); RR[is.na(RR)] <- 0; RR <- RR %>% tidyr::pivot_wider(names_from = RR.CI, values_from = RR.CI)
  OR <- data.frame(result.Qcgc$estimates[4]); OR[is.na(OR)] <- 0; OR <- OR %>% tidyr::pivot_wider(names_from = OR.CI, values_from = OR.CI)
  ATT <- data.frame(result.Qcgc$estimates[6]); ATT[is.na(ATT)] <- 0; ATT <- ATT %>% tidyr::pivot_wider(names_from = ATT.CI, values_from = ATT.CI)
  ATC <- data.frame(result.Qcgc$estimates[7]); ATC[is.na(ATC)] <- 0; ATC <- ATC %>% tidyr::pivot_wider(names_from = ATC.CI, values_from = ATC.CI)
  
  tmle_comparison[i,1] <- i#sl__model1$libraryNames[i]
  tmle_comparison[i,2:6] <- ATE
  tmle_comparison[i,7:12] <- RR
  tmle_comparison[i,13:18] <- OR
  tmle_comparison[i,19:24] <- ATT
  tmle_comparison[i,25:30] <- ATC
  tmle_comparison[i,31] <- apply(data.frame(lapply(sl__model1$AllSL, function(x) x$cvRisk)), 1, mean, na.rm =TRUE)[i]
  tmle_comparison[i,32] <- apply(data.frame(lapply(sl__model0$AllSL, function(x) x$cvRisk)), 1, mean, na.rm =TRUE)[i]
  
  if (i == libs_final[length(libs_final)]) {
    colnames(tmle_comparison) <- c("learner_strategy", "ATE.psi", "ATE.var.psi", "ATE.pvalue", "ATE.CI.LO", "ATE.CI.UP", "RR.psi", "RR.pvalue", "RR.log.psi", "RR.var.log.psi", "RR.CI.LO", "RR.CI.UP",
                                   "OR.psi", "OR.pvalue", "OR.log.psi", "OR.var.log.psi", "OR.CI.LO", "OR.CI.UP","ATT.psi", "ATT.converged", "ATT.var.psi", "ATT.pvalue", "ATT.CI.LO", "ATT.CI.UP",
                                   "ATC.psi", "ATC.converged", "ATC.var.psi", "ATC.pvalue", "ATC.CI.LO", "ATC.CI.UP", "cvRisk_sl_model1", "cvRisk_sl_model0")
  }
  
  # print(tmle_comparison[i,])
}

tmle_comparison$learner_strategy <- gsub("_All", "", tmle_comparison$learner_strategy)
# tmle_comparison$learner_strategy <- gsub("Late", "Int", tmle_comparison$learner_strategy)

## superlearner
df_super1 <- list()
k <- 1
for (i in sl__model1$AllS) {
  for (j in libs_final) {
    df_super1[[j]] <- SuperLearner::predict.SuperLearner(i$fitLibrary[[j]]$object, newdata = W.new, type = "response")$pred * sl__model1$coef[k,j]
  }
  k <- k + 1
}

Q1W <- apply(data.frame(df_super1), 1, mean)

df_super0 <- list()
k <- 1
for (i in sl__model0$AllS) {
  for (j in libs_final) {
    df_super0[[j]] <- SuperLearner::predict.SuperLearner(i$fitLibrary[[j]]$object, newdata = W.new, type = "response")$pred * sl__model0$coef[k,j]
  }
  k <- k + 1
}

Q0W <- apply(data.frame(df_super0), 1, mean)

Q <- data.frame(Q0W,Q1W)
head(Q); tail(Q); dim(Q)

result.Qcgc <- tmle(Y = Y, A = A, W = W.new, verbose = TRUE, Delta = Delta,
                    Q = Q, g1W = g1W, pDelta1 = pDelta1, family = "binomial")

saveRDS(object = result.Qcgc, file = paste0(ruta, "os_lgg_data_var", nvar, "/tmle_result_var", nvar, ".rds"))

ATE <- data.frame(result.Qcgc$estimates[2]); ATE[is.na(ATE)] <- 0; ATE <- ATE %>% tidyr::pivot_wider(names_from = ATE.CI, values_from = ATE.CI)
RR <- data.frame(result.Qcgc$estimates[3]); RR[is.na(RR)] <- 0; RR <- RR %>% tidyr::pivot_wider(names_from = RR.CI, values_from = RR.CI)
OR <- data.frame(result.Qcgc$estimates[4]); OR[is.na(OR)] <- 0; OR <- OR %>% tidyr::pivot_wider(names_from = OR.CI, values_from = OR.CI)
ATT <- data.frame(result.Qcgc$estimates[6]); ATT[is.na(ATT)] <- 0; ATT <- ATT %>% tidyr::pivot_wider(names_from = ATT.CI, values_from = ATT.CI)
ATC <- data.frame(result.Qcgc$estimates[7]); ATC[is.na(ATC)] <- 0; ATC <- ATC %>% tidyr::pivot_wider(names_from = ATC.CI, values_from = ATC.CI)

tmle_comparison["SuperLearner",1] <- "SuperLearner"
tmle_comparison["SuperLearner",2:6] <- ATE
tmle_comparison["SuperLearner",7:12] <- RR
tmle_comparison["SuperLearner",13:18] <- OR
tmle_comparison["SuperLearner",19:24] <- ATT
tmle_comparison["SuperLearner",25:30] <- ATC

### best model in 1 and best model in 0
which_model1 <- names(which.max(apply(sl__model1$coef, 2, mean)))
which_model0 <- names(which.max(apply(sl__model0$coef, 2, mean)))

Q1W <- apply(data.frame(lapply(sl__model1$AllSL, function(x) SuperLearner::predict.SuperLearner(x$fitLibrary[[ which_model1 ]]$object, newdata = W.new, type = "response")$pred)), 1, mean)
Q0W <- apply(data.frame(lapply(sl__model0$AllSL, function(x) SuperLearner::predict.SuperLearner(x$fitLibrary[[ which_model0 ]]$object, newdata = W.new, type = "response")$pred)), 1, mean)


Q <- data.frame(Q0W,Q1W)
head(Q); tail(Q); dim(Q)

## ------------------------ ##
### ... g. tmle update ... ###
## ------------------------ ##
result.Qcgc <- tmle(Y = Y, A = A, W = W.new, verbose = TRUE, Delta = Delta,
                    Q = Q, g1W = g1W, pDelta1 = pDelta1, family = "binomial")

saveRDS(object = result.Qcgc, file = paste0(ruta, "os_lgg_data_var", nvar, "/tmle_result_var", nvar, "_best_SuperLearners.rds"))

ATE <- data.frame(result.Qcgc$estimates[2]); ATE[is.na(ATE)] <- 0; ATE <- ATE %>% tidyr::pivot_wider(names_from = ATE.CI, values_from = ATE.CI)
RR <- data.frame(result.Qcgc$estimates[3]); RR[is.na(RR)] <- 0; RR <- RR %>% tidyr::pivot_wider(names_from = RR.CI, values_from = RR.CI)
OR <- data.frame(result.Qcgc$estimates[4]); OR[is.na(OR)] <- 0; OR <- OR %>% tidyr::pivot_wider(names_from = OR.CI, values_from = OR.CI)
ATT <- data.frame(result.Qcgc$estimates[6]); ATT[is.na(ATT)] <- 0; ATT <- ATT %>% tidyr::pivot_wider(names_from = ATT.CI, values_from = ATT.CI)
ATC <- data.frame(result.Qcgc$estimates[7]); ATC[is.na(ATC)] <- 0; ATC <- ATC %>% tidyr::pivot_wider(names_from = ATC.CI, values_from = ATC.CI)


modelname <- paste0("M1: ", which_model1, " -", " /n M0: ", which_model0 )

tmle_comparison[modelname,1] <- modelname
tmle_comparison[modelname,2:6] <- ATE
tmle_comparison[modelname,7:12] <- RR
tmle_comparison[modelname,13:18] <- OR
tmle_comparison[modelname,19:24] <- ATT
tmle_comparison[modelname,25:30] <- ATC

save(x = tmle_comparison, file = paste0(ruta, "os_lgg_data_var", nvar, "/tmle_comparison_var", nvar, ".rda"))

## ------------------------------------------------------------ ##
### ... Graficos y m?s graficos de tmle y del superlearner ... ###
## ------------------------------------------------------------ ##

### ATE plot
psi.ate <- ggplot(tmle_comparison, aes(x=reorder(learner_strategy, -ATE.psi), y=ATE.psi)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_text_repel(label = round(tmle_comparison$ATE.psi, digits = 4), size=3) + 
  geom_errorbar(aes(ymin=ATE.CI.LO,ymax=ATE.CI.UP), width=.2,
                position=position_dodge(0.05)) + 
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
        text = element_text(size = 13))

### ATT plot
psi.att <- ggplot(tmle_comparison, aes(x=reorder(learner_strategy, -ATE.psi), y=ATT.psi)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_text_repel(label = round(tmle_comparison$ATT.psi, digits = 4), size=3) + 
  geom_errorbar(aes(ymin=ATT.CI.LO,ymax=ATT.CI.UP), width=.2,
                position=position_dodge(0.05)) + 
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
        text = element_text(size = 13))

### ATC plot
psi.atc <- ggplot(tmle_comparison, aes(x=reorder(learner_strategy, -ATE.psi), y=ATC.psi)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_text_repel(label = round(tmle_comparison$ATC.psi, digits = 4), size=3) + 
  geom_errorbar(aes(ymin=ATC.CI.LO,ymax=ATC.CI.UP), width=.2,
                position=position_dodge(0.05)) + 
  theme(axis.text.x = element_text(angle = 0, hjust=.75,vjust=0.95), 
        text = element_text(size = 13)) + xlab("Method")

psi <- ggarrange(psi.ate, psi.att, psi.atc, # + rremove("x.text"),
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)

### RR plot
psi.rr <- ggplot(tmle_comparison, aes(x=reorder(learner_strategy, -ATE.psi), y=RR.psi)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept=1, linetype="dashed", color = "black") + 
  geom_text_repel(label = round(tmle_comparison$RR.psi, digits = 4), size=3) + 
  geom_errorbar(aes(ymin=RR.CI.LO,ymax=RR.CI.UP), width=.2,
                position=position_dodge(0.05)) + 
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
        text = element_text(size = 13))

### OR plot
psi.or <- ggplot(tmle_comparison, aes(x=reorder(learner_strategy, -ATE.psi), y=OR.psi)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept=1, linetype="dashed", color = "black") + 
  geom_text_repel(label = round(tmle_comparison$OR.psi, digits = 4), size=3) + 
  geom_errorbar(aes(ymin=OR.CI.LO,ymax=OR.CI.UP), width=.2,
                position=position_dodge(0.05)) + 
  theme(axis.text.x = element_text(angle = 0, hjust=.75,vjust=0.95), 
        text = element_text(size = 13)) + xlab("Method")

or.rr <- ggarrange(psi.rr, psi.or, # + rremove("x.text"),
          labels = c("D", "E"),
          ncol = 1, nrow = 2)

###first set of plots

jpeg(filename = paste0(ruta, "os_lgg_data_var", nvar, "/fig_ATE_RR_OR_var", nvar, ".jpg"), width = 13, height = 8, units = "in", res=300)
ggarrange(psi, or.rr, ncol = 2)
dev.off()

# tiff(filename = paste0(ruta, "os_lgg_data_var", nvar, "/fig_ATE_RR_OR_var", nvar, ".tiff"), width = 13, height = 8, units = "in", res=300)
# ggarrange(psi, or.rr, ncol = 2)
# dev.off()

pdf(file = paste0(ruta, "os_lgg_data_var", nvar, "/fig_ATE_RR_OR_var", nvar, ".pdf"), width = 13, height = 8)
ggarrange(psi, or.rr, ncol = 2)
dev.off()

###second set of plots

jpeg(filename = paste0(ruta, "os_lgg_data_var", nvar, "/fig_ps_and_cvRisk_var", nvar, ".jpg"), width = 15, height = 10, units = "in", res=300)
ggarrange(ps_plot, ggarrange(plot(sl__model0) + theme_bw() + labs(title="E[Y|A=0,W]"), plot(sl__model1) + theme_bw() + labs(title="E[Y|A=1,W]"), ncol = 1, labels = c("B", "C")), ncol = 2, labels = c("A"))
dev.off()

# tiff(filename = paste0(ruta, "os_lgg_data_var", nvar, "/fig_ps_and_cvRisk_var", nvar, ".tiff"), width = 15, height = 10, units = "in", res=300)
# ggarrange(ps_plot, ggarrange(plot(sl__model0) + theme_bw() + labs(title="E[Y|A=0,W]"), plot(sl__model1) + theme_bw() + labs(title="E[Y|A=1,W]"), ncol = 1, labels = c("B", "C")), ncol = 2, labels = c("A"))
# dev.off()

pdf(file = paste0(ruta, "os_lgg_data_var", nvar, "/fig_ps_and_cvRisk_var", nvar, ".pdf"), width = 15, height = 10)
ggarrange(ps_plot, ggarrange(plot(sl__model0) + theme_bw() + labs(title="E[Y|A=0,W]"), plot(sl__model1) + theme_bw() + labs(title="E[Y|A=1,W]"), ncol = 1, labels = c("B", "C")), ncol = 2, labels = c("A"))
dev.off()

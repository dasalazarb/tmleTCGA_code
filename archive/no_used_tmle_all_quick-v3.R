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
library(tmle)
library(caret)
library(dplyr)
library(purrr)
library(ggplot2)
library(SuperLearner)

## ---------------------- ##
### ... 1. Load data ... ###
## ---------------------- ##
source("C:/Users/da.salazarb/Documents/tmle/tmleTCGA/func_extractVars.R")
source("C:/Users/da.salazarb/Documents/tmle/tmleTCGA/func_learners_v2.R")
source("C:/Users/da.salazarb/Documents/tmle/tmleTCGA/func_SuperLearners.R")

## --------------------------------------------------------------------- ##
### ... 2. Conditional treatment assignment probabilities, P(A=1|W) ... ###
## --------------------------------------------------------------------- ##
clinical <- W %>% 
  dplyr::select(.,starts_with("clinical")) %>%
  # dplyr::select(-clinical_primary_diagnosis) %>%
  dplyr::mutate_if(is.character, as.factor, ordered = TRUE)

# tune = list(num.trees = 5000,
#             max_depth = c(25, 200),
#             mtry_seq = floor(sqrt(ncol(clinical)) * c(0.5, 1, 2)))
# 
# learners <- create.Learner("SL.ranger", tune = tune)
# sl_log <- mcSuperLearner(Y = A, X = clinical, family = binomial(), cvControl = list(V = 10),
#                         SL.library = c(learners$names))

lrn_log <- create.Learner("SL.glm")
# lrn_log <- create.Learner("SL.gam")
sl_log <- mcSuperLearner(Y = A, X = clinical, family = binomial(), cvControl = list(V = 10),
                         SL.library = c(lrn_log$names))

summary(sl_log)
g1W <- sl_log$SL.predict

## ----------------------- ##
### ... Plot P(A=1|W) ... ###
## ----------------------- ##
ps <- data.frame(cbind(g1W,treat=A))
head(ps,20)

# PS plot
ggplot() +
  # Top
  geom_density(data=subset(ps, treat==0), aes(x = V1, y = ..density..), fill="#69b3a2" ) +
  # geom_label( aes(x=0.4, y=250, label="RAD+TMZ"), color="#69b3a2") +
  # Bottom
  geom_density(data=subset(ps, treat==1), aes(x = V1, y = -..density..), fill= "#404080") #+ scale_y_continuous(breaks=round(seq(min(-2.5), max(3.5), by=0.2),2))
# geom_label( aes(x=0.5, y=-100, label="RAD"), color="#404080")

## ---------------------------------------------------------------------------------------------------- ##
### ... 3. Conditional probabilities for missingness mechanism, P(Delta=1|A=0,W), P(Delta=1|A=1,W) ... ###
## ---------------------------------------------------------------------------------------------------- ##
lrn_log <- create.Learner("SL.glm")
sl_rf <- mcSuperLearner(Y = Delta, X = data.frame(A=as.factor(A), clinical), family = binomial(), cvControl = list(V = 10),
                        SL.library = c(lrn_log$names))

# tune = list(num.trees = 1000,
#             max_depth = c(25, 200),
#             mtry_seq = floor(sqrt(ncol(clinical)) * c(0.5, 1, 2)))
# 
# learners <- create.Learner("SL.ranger", tune = tune)
# sl_rf <- mcSuperLearner(Y = Delta, X = data.frame(A=as.factor(A), clinical), family = binomial(), cvControl = list(V = 10),
#                         SL.library = c(learners$names))
summary(sl_rf)

d1W <- predict(sl_rf,newdata=data.frame(clinical, A=factor(1)), type="response")$pred
d0W <- predict(sl_rf,newdata=data.frame(clinical, A=factor(0)), type = "response")$pred

pDelta1 <- data.frame(d0W,d1W)
head(pDelta1); dim(pDelta1)


## ----------------------------------- ##
### ... 5. Estimation of E(Y|A,W) ... ###
## ----------------------------------- ##

### caso 1: random forest -> Resultado: Funciona
tune = list(num.trees = 5000, 
            max_depth = 50)

lrn_rf <- create.Learner("SL.ranger", params = tune)
lrn_rf.clinical <- create.Learner("SL.rf.clinical", params = tune)
lrn_rf.mrna <- create.Learner("SL.rf.mrna", params = tune)
lrn_rf.mrna.pca <- create.Learner("SL.rf.mrna.pca", params = tune)
lrn_rf.mrna.encoder <- create.Learner("SL.rf.mrna.encoder", params = tune)
lrn_rf.mrna.nmf <- create.Learner("SL.rf.mrna.nmf", params = tune)
lrn_rf.cnv <- create.Learner("SL.rf.cnv", params = tune)
lrn_rf.cnv.pca <- create.Learner("SL.rf.cnv.pca", params = tune)
lrn_rf.cnv.encoder <- create.Learner("SL.rf.cnv.encoder", params = tune)
lrn_rf.cnv.nmf <- create.Learner("SL.rf.cnv.nmf", params = tune)


### caso 2: logistic regression 
lrn_log <- create.Learner("SL.glm")
lrn_log.clinical <- create.Learner("SL.glm.clinical", params = tune)
lrn_log.mrna <- create.Learner("SL.glm.mrna", params = tune)
lrn_log.mrna.pca <- create.Learner("SL.glm.mrna.pca", params = tune)
lrn_log.mrna.encoder <- create.Learner("SL.glm.mrna.encoder", params = tune)
lrn_log.mrna.nmf <- create.Learner("SL.glm.mrna.nmf", params = tune)
lrn_log.cnv <- create.Learner("SL.glm.cnv", params = tune)
lrn_log.cnv.pca <- create.Learner("SL.glm.cnv.pca", params = tune)
lrn_log.cnv.encoder <- create.Learner("SL.glm.cnv.encoder", params = tune)
lrn_log.cnv.nmf <- create.Learner("SL.glm.cnv.nmf", params = tune)


### caso 3: xgboost -> Resultado: Funciona
tune = list(ntrees = 5000,
            max_depth = 50, 
            shrinkage = 0.001)

lrn_xgb <- create.Learner("SL.xgboost", params = tune, detailed_names = TRUE, name_prefix = "xgb")
lrn_xgb.clinical <- create.Learner("SL.xgboost.clinical", detailed_names = TRUE, name_prefix = "xgb.clinical")
lrn_xgb.mrna <- create.Learner("SL.xgboost.mrna", detailed_names = TRUE, name_prefix = "xgb.mrna")
lrn_xgb.mrna.pca <- create.Learner("SL.xgboost.mrna.pca", detailed_names = TRUE, name_prefix = "xgb.mrna.pca")
lrn_xgb.mrna.encoder <- create.Learner("SL.xgboost.mrna.encoder", detailed_names = TRUE, name_prefix = "xgb.mrna.encoder")
lrn_xgb.mrna.nmf <- create.Learner("SL.xgboost.mrna.nmf", detailed_names = TRUE, name_prefix = "xgb.mrna.nmf")
lrn_xgb.cnv <- create.Learner("SL.xgboost.cnv", detailed_names = TRUE, name_prefix = "xgb.cnv")
lrn_xgb.cnv.pca <- create.Learner("SL.xgboost.cnv.pca", detailed_names = TRUE, name_prefix = "xgb.cnv.pca")
lrn_xgb.cnv.encoder <- create.Learner("SL.xgboost.cnv.encoder", detailed_names = TRUE, name_prefix = "xgb.cnv.encoder")
lrn_xgb.cnv.nmf <- create.Learner("SL.xgboost.cnv.nmf", detailed_names = TRUE, name_prefix = "xgb.cnv.nmf")


## caso 4: mars -> Resultado: Funciona
tune = list(degree = 2,
            # nprune = 5,
            Use.beta.cache=FALSE, 
            newvar.penalty=0.01, 
            thresh=0.01) # Note: earth's beta cache would require 3.63 GB, so forcing Use.beta.cache=FALSE.

lrn_mars <- create.Learner("SL.earth.FeatureSelection", params = tune, name_prefix = "mars")
lrn_mars.clinical <- create.Learner("SL.earth.clinical", params = tune, name_prefix = "mars.clinical")
lrn_mars.mrna <- create.Learner("SL.earth.mrna", params = tune, name_prefix = "mars.mrna")
lrn_mars.mrna.pca <- create.Learner("SL.earth.mrna.pca", params = tune, name_prefix = "mars.mrna.pca")
lrn_mars.mrna.encoder <- create.Learner("SL.earth.mrna.encoder", params = tune, name_prefix = "mars.mrna.encoder")
lrn_mars.mrna.nmf <- create.Learner("SL.earth.mrna.nmf", params = tune, name_prefix = "mars.mrna.nmf")
lrn_mars.cnv <- create.Learner("SL.earth.cnv", params = tune, name_prefix = "mars.cnv")
lrn_mars.cnv.pca <- create.Learner("SL.earth.cnv.pca", params = tune, name_prefix = "mars.cnv.pca")
lrn_mars.cnv.encoder <- create.Learner("SL.earth.cnv.encoder", params = tune, name_prefix = "mars.cnv.encoder")
lrn_mars.cnv.nmf <- create.Learner("SL.earth.cnv.nmf", params = tune, name_prefix = "mars.cnv.nmf")


## Todas las estrategias
sl__model <- mcSuperLearner(Y = Y, X = data.frame(W.new, clinical_A=as.factor(A)), family = binomial(), cvControl = list(V = 3),
                            SL.library = c(
                              ##### clinical info
                              # lrn_mars.clinical$names, 
                              lrn_rf.clinical$names,
                              # lrn_xgb.clinical$names,
                              # lrn_log.clinical$names,

                              ##### early integration
                              # lrn_mars$names,
                              lrn_rf$names,
                              # lrn_xgb$names,
                              # lrn_log$names,

                              ##### late integration
                              # lrn_mars.mrna$names, lrn_mars.cnv$names,
                              lrn_rf.mrna$names, lrn_rf.cnv$names,
                              # lrn_xgb.mrna$names, lrn_xgb.cnv$names,
                              # lrn_log.mrna$names, lrn_log.cnv$names,

                              ##### intermediate integration (encoder)
                              # lrn_mars.mrna.encoder$names, lrn_mars.cnv.encoder$names,
                              lrn_rf.mrna.encoder$names, lrn_rf.cnv.encoder$names,
                              # lrn_xgb.mrna.encoder$names, lrn_xgb.cnv.encoder$names,
                              # lrn_log.mrna.encoder$names, lrn_log.cnv.encoder$names,

                              ##### intermediate integration (pca)
                              # lrn_mars.mrna.pca$names, lrn_mars.cnv.pca$names,
                              lrn_rf.mrna.pca$names, lrn_rf.cnv.pca$names,
                              # lrn_xgb.mrna.pca$names, lrn_xgb.cnv.pca$names,
                              # lrn_log.mrna.pca$names, lrn_log.cnv.pca$names,

                              ##### intermediate integration (nmf)
                              # lrn_mars.mrna.nmf$names, lrn_mars.cnv.nmf$names
                              lrn_rf.mrna.nmf$names, lrn_rf.cnv.nmf$names
                              # lrn_xgb.mrna.nmf$names, lrn_xgb.cnv.nmf$names
                              # lrn_log.mrna.nmf$names, lrn_log.cnv.nmf$names
                             ))


Q1W <- SuperLearner::predict.SuperLearner(sl__model, newdata = data.frame(W.new, clinical_A=factor(1)), type = "response", onlySL = T)$pred
Q0W <- SuperLearner::predict.SuperLearner(sl__model, newdata = data.frame(W.new, clinical_A=factor(0)), type = "response", onlySL = T)$pred
Q <- data.frame(Q0W,Q1W)
head(Q); tail(Q); dim(Q)

## ------------------------ ##
### ... g. tmle update ... ###
## ------------------------ ##
result.Qcgc <- tmle(Y = Y, A = A, W = W.new, verbose = TRUE, Delta = Delta,
                    Q = Q, g1W = g1W, pDelta1 = pDelta1, family = "binomial")

result.Qcgc

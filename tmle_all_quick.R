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
# SL.gam
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
            max_depth = 50, 
            mtry_seq = floor(sqrt(ncol(W))))

lrn_rf <- create.Learner("SL.ranger", params = tune)


### caso 2: logistic regression 
lrn_log <- create.Learner("SL.glm")


### caso 3: xgboost -> Resultado: Funciona
tune = list(ntrees = 5000,
            max_depth = 50, 
            shrinkage = 0.001)

lrn_xgb <- create.Learner("SL.xgboost", params = tune, detailed_names = TRUE, name_prefix = "xgb")


## caso 4: mars -> Resultado: Funciona
tune = list(degree = 2,
            nprune = 5,
            Use.beta.cache=FALSE) # Note: earth's beta cache would require 3.63 GB, so forcing Use.beta.cache=FALSE.

lrn_mars <- create.Learner("SL.earth", params = tune, name_prefix = "mars")


## Todas las estrategias
sl__model1 <- mcSuperLearner(Y = Y[A==1], X = W.new[A==1,], family = binomial(), cvControl = list(V = 3),
                            SL.library = purrr::reduce(list(
                              # clinical info
                              lapply(lrn_mars$names, function(x) c(x, "func_clinical")), 
                              lapply(lrn_rf$names, function(x) c(x, "func_clinical")), 
                              lapply(lrn_xgb$names, function(x) c(x, "func_clinical")), 
                              lapply(lrn_log$names, function(x) c(x, "func_clinical")), 
                              # # early integration
                              # lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest")),
                              # lapply(lrn_rf$names, function(x) c(x)),
                              # lapply(lrn_xgb$names, function(x) c(x)),
                              # lapply(lrn_log$names, function(x) c(x)),
                              # # late integration
                              # lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna")),
                              # lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv")),
                              # lapply(lrn_rf$names, function(x) c(x, "func_mrna")),
                              # lapply(lrn_rf$names, function(x) c(x, "func_cnv")),
                              # lapply(lrn_xgb$names, function(x) c(x, "func_mrna")),
                              # lapply(lrn_xgb$names, function(x) c(x, "func_cnv")),
                              # lapply(lrn_log$names, function(x) c(x, "func_mrna")),
                              # lapply(lrn_log$names, function(x) c(x, "func_cnv")),
                              # intermediate integration (encoder)
                              lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna.encoder")),
                              lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv.encoder")),
                              # lapply(lrn_rf$names, function(x) c(x, "func_mrna.encoder")),
                              # lapply(lrn_rf$names, function(x) c(x, "func_cnv.encoder")),
                              lapply(lrn_xgb$names, function(x) c(x, "func_mrna.encoder")),
                              lapply(lrn_xgb$names, function(x) c(x, "func_cnv.encoder")),
                              # lapply(lrn_log$names, function(x) c(x, "func_mrna.encoder")),
                              # lapply(lrn_log$names, function(x) c(x, "func_cnv.encoder")),
                              # intermediate integration (pca)
                              lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna.pca")),
                              lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv.pca")),
                              # lapply(lrn_rf$names, function(x) c(x, "func_mrna.pca")),
                              # lapply(lrn_rf$names, function(x) c(x, "func_cnv.pca")),
                              lapply(lrn_xgb$names, function(x) c(x, "func_mrna.pca")),
                              lapply(lrn_xgb$names, function(x) c(x, "func_cnv.pca")),
                              # lapply(lrn_log$names, function(x) c(x, "func_mrna.pca")),
                              # lapply(lrn_log$names, function(x) c(x, "func_cnv.pca")),
                              # intermediate integration (nmf)
                              lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna.nmf")),
                              lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv.nmf")),
                              # lapply(lrn_rf$names, function(x) c(x, "func_mrna.nmf")),
                              # lapply(lrn_rf$names, function(x) c(x, "func_cnv.nmf")),
                              lapply(lrn_xgb$names, function(x) c(x, "func_mrna.nmf")),
                              lapply(lrn_xgb$names, function(x) c(x, "func_cnv.nmf"))#,
                              # lapply(lrn_log$names, function(x) c(x, "func_mrna.nmf")),
                              # lapply(lrn_log$names, function(x) c(x, "func_cnv.nmf"))
                              ), append))

## Todas las estrategias
sl__model0 <- mcSuperLearner(Y = Y[A==0], X = W.new[A==0,], family = binomial(), cvControl = list(V = 3),
                             SL.library = purrr::reduce(list(
                               # clinical info
                               lapply(lrn_mars$names, function(x) c(x, "func_clinical")), 
                               lapply(lrn_rf$names, function(x) c(x, "func_clinical")), 
                               lapply(lrn_xgb$names, function(x) c(x, "func_clinical")), 
                               lapply(lrn_log$names, function(x) c(x, "func_clinical")), 
                               # # early integration
                               # lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest")),
                               # lapply(lrn_rf$names, function(x) c(x)),
                               # lapply(lrn_xgb$names, function(x) c(x)),
                               # lapply(lrn_log$names, function(x) c(x)),
                               # # late integration
                               # lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna")),
                               # lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv")),
                               # lapply(lrn_rf$names, function(x) c(x, "func_mrna")),
                               # lapply(lrn_rf$names, function(x) c(x, "func_cnv")),
                               # lapply(lrn_xgb$names, function(x) c(x, "func_mrna")),
                               # lapply(lrn_xgb$names, function(x) c(x, "func_cnv")),
                               # lapply(lrn_log$names, function(x) c(x, "func_mrna")),
                               # lapply(lrn_log$names, function(x) c(x, "func_cnv")),
                               # intermediate integration (encoder)
                               lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna.encoder")),
                               lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv.encoder")),
                               # lapply(lrn_rf$names, function(x) c(x, "func_mrna.encoder")),
                               # lapply(lrn_rf$names, function(x) c(x, "func_cnv.encoder")),
                               lapply(lrn_xgb$names, function(x) c(x, "func_mrna.encoder")),
                               lapply(lrn_xgb$names, function(x) c(x, "func_cnv.encoder")),
                               # lapply(lrn_log$names, function(x) c(x, "func_mrna.encoder")),
                               # lapply(lrn_log$names, function(x) c(x, "func_cnv.encoder")),
                               # intermediate integration (pca)
                               lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna.pca")),
                               lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv.pca")),
                               # lapply(lrn_rf$names, function(x) c(x, "func_mrna.pca")),
                               # lapply(lrn_rf$names, function(x) c(x, "func_cnv.pca")),
                               lapply(lrn_xgb$names, function(x) c(x, "func_mrna.pca")),
                               lapply(lrn_xgb$names, function(x) c(x, "func_cnv.pca")),
                               # lapply(lrn_log$names, function(x) c(x, "func_mrna.pca")),
                               # lapply(lrn_log$names, function(x) c(x, "func_cnv.pca")),
                               # intermediate integration (nmf)
                               lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna.nmf")),
                               lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv.nmf")),
                               # lapply(lrn_rf$names, function(x) c(x, "func_mrna.nmf")),
                               # lapply(lrn_rf$names, function(x) c(x, "func_cnv.nmf")),
                               lapply(lrn_xgb$names, function(x) c(x, "func_mrna.nmf")),
                               lapply(lrn_xgb$names, function(x) c(x, "func_cnv.nmf"))#,
                               # lapply(lrn_log$names, function(x) c(x, "func_mrna.nmf")),
                               # lapply(lrn_log$names, function(x) c(x, "func_cnv.nmf"))
                             ), append))


Q1W <- SuperLearner::predict.SuperLearner(sl__model1, newdata = W.new, type = "response")$pred
Q0W <- SuperLearner::predict.SuperLearner(sl__model0, newdata = W.new, type = "response")$pred
Q <- data.frame(Q0W,Q1W)
head(Q); tail(Q); dim(Q)

## ------------------------ ##
### ... g. tmle update ... ###
## ------------------------ ##
result.Qcgc <- tmle(Y = Y, A = A, W = W.new, verbose = TRUE, Delta = Delta,
                    Q = Q, g1W = g1W, pDelta1 = pDelta1, family = "binomial")

result.Qcgc

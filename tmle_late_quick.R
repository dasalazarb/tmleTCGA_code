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

func_mrna <- function(X,...){
  returnCols <- rep(FALSE, ncol(X))
  # returnCols[grep("mrna", names(X))] <- TRUE
  returnCols[grepl("A|mrna_|clinical_", names(X))] <- TRUE
  return(returnCols)
}

func_cnv <- function(X,...){
  returnCols <- rep(FALSE, ncol(X))
  # returnCols[grep("cnv", names(X))] <- TRUE
  returnCols[grepl("A|cnv_|clinical_", names(X))] <- TRUE
  return(returnCols)
}

## ---------------------- ##
### ... 1. Load data ... ###
## ---------------------- ##
O <- read.csv(file = "C:/Users/da.salazarb/Documents/tmle/tmleTCGA/O_Y_A_Delta_W.csv") ## <- Cambiar ruta del archivo!
O[1:5,1:5]; dim(O)

Y <- ifelse(O$Y < 365, 1, 0); Delta <- ifelse(O$delta == 0 & O$Y < 365, 0, 1); A <- O$A; W <- O[,!(colnames(O) %in% c("Y", "X", "delta", "...3078"))];

W <- W %>% 
  dplyr::mutate(clinical_primary_diagnosis=as.factor(clinical_primary_diagnosis), clinical_gender = as.factor(clinical_gender), 
                clinical_IDH_codel_subtype=as.factor(clinical_IDH_codel_subtype), clinical_MGMT_promoter_status=as.factor(clinical_MGMT_promoter_status), 
                clinical_ATRX_status=as.factor(clinical_ATRX_status), clinical_Chr7_gain_and_Chr10_loss=as.factor(clinical_Chr7_gain_and_Chr10_loss), 
                A=as.factor(A))

## --------------------------------------------------------------------- ##
### ... 2. Conditional treatment assignment probabilities, P(A=1|W) ... ###
## --------------------------------------------------------------------- ##
clinical <- W %>% 
  dplyr::select(.,starts_with("clinical"))

tune = list(num.trees = 1000,
            max_depth = c(25, 200), # seq(10,200,100), 
            mtry_seq = floor(sqrt(ncol(clinical)) * c(0.5, 1, 2)))

learners <- create.Learner("SL.ranger", tune = tune)
sl_rf <- mcSuperLearner(Y = A, X = clinical, family = binomial(), cvControl = list(V = 10),
                        SL.library = c(learners$names))
summary(sl_rf)
g1W <- sl_rf$SL.predict

## ----------------------------------------------------------------------------------------------------- ##
### ... 3. Conditional probabilities for missingness mechanism, P(Delta=1|A=0,W), P(Delta=1|A=1,W) ... ###
## ----------------------------------------------------------------------------------------------------- ##
tune = list(num.trees = 1000,
            max_depth = c(25, 200), # seq(10,200,100),  
            mtry_seq = floor(sqrt(ncol(clinical)) * c(0.5, 1, 2)))

learners <- create.Learner("SL.ranger", tune = tune)
sl_rf <- mcSuperLearner(Y = Delta, X = clinical, family = binomial(), cvControl = list(V = 10),
                        SL.library = c(learners$names))
summary(sl_rf)

d1W <- predict(sl_rf,newX=data.frame(learners, A=1), type="response")$pred
d0W <- predict(sl_rf,newX=data.frame(learners, A=0), type = "response")$pred

pDelta1 <- data.frame(d0W,d1W)
head(pDelta1)

## -------------------------- ##
### ... 4. plot P(A=1|W) ... ###
## -------------------------- ##
ps <- data.frame(cbind(g1W,treat=A))
head(ps)

# PS plot
ggplot() +
  # Top
  geom_density(data=subset(ps, treat==0), aes(x = g1W, y = ..density..), fill="#69b3a2" ) +
  # geom_label( aes(x=0.4, y=250, label="RAD+TMZ"), color="#69b3a2") +
  # Bottom
  geom_density(data=subset(ps, treat==1), aes(x = g1W, y = -..density..), fill= "#404080") #+ scale_y_continuous(breaks=round(seq(min(-2.5), max(3.5), by=0.2),2))
# geom_label( aes(x=0.5, y=-100, label="RAD"), color="#404080")

## ----------------------------------- ##
### ... 5. Estimation of E(Y|A,W) ... ###
## ----------------------------------- ##
W_ <- W
# names(W_)  <- paste0('x', 1:ncol(W_))
W_[1:5,1:5]; dim(W_)

### caso 1: random forest -> Resultado: Funciona
SL.LinearComb.rf <- function(Y, X, family, ...) {
  # load required packages
  require("ranger")
  
  tune = list(num.trees = 5000, 
              max_depth = c(25,200), 
              mtry_seq = floor(sqrt(ncol(W_)) * c(0.5, 1, 2)))
  
  lrn_rf <- create.Learner("SL.ranger", tune = tune)
  
  sl__model <- mcSuperLearner(Y = Y, X = X, family = family$family, #cvControl = list(V = 10),
                              SL.library = append(lapply(lrn_rf$names, function(x) c(x, "func_mrna")), lapply(lrn_rf$names, function(x) c(x, "func_cnv"))),
                              ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl__model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl__model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

# sl_rf_late <- mcSuperLearner(Y = Y, X = W_, family = binomial(), cvControl = list(V = 10), 
#                         SL.library = append(lapply(lrn_rf$names, function(x) c(x, "func_mrna")), lapply(lrn_rf$names, function(x) c(x, "func_cnv"))))


### caso 2: kernelSVM -> Resultado: Funciona
SL.LinearComb.SVM <- function(Y, X, family, ...) {
  # load required packages
  require("kernlab")
  
  tune <- list(kernel = c("rbfdot", "laplacedot"), 
               kpar = list(list(sigma=2), 
                           list(sigma=3)))
  lrn_rbf <- create.Learner("SL.ksvm", tune = tune, name_prefix = "rbf")
  
  tune <- list(kernel = c("polydot"), 
               kpar = list(list(degree=2), 
                           list(degree=3)))
  lrn_poly <- create.Learner("SL.ksvm", tune = tune, name_prefix = "poly")
  
  tune <- list(kernel = c("tanhdot"), 
               kpar = list(list(scale=2, offset=1),
                           list(scale=3, offset=1)))
  lrn_tanh <- create.Learner("SL.ksvm", tune = tune, name_prefix = "tanh")
  
  tune <- list(kernel = c("anovadot"), 
               kpar = list(list(sigma=2, degree=2),
                           list(sigma=3, degree=2)))
  lrn_anova <- create.Learner("SL.ksvm", tune = tune, name_prefix = "anova")
  
  sl__model <- mcSuperLearner(Y = Y, X = X, family = family$family, #cvControl = list(V = 10),
                              SL.library = purrr::reduce(list(lapply(lrn_rbf$names, function(x) c(x, "func_mrna")), lapply(lrn_rbf$names, function(x) c(x, "func_cnv")),
                                                              lapply(lrn_poly$names, function(x) c(x, "func_mrna")), lapply(lrn_poly$names, function(x) c(x, "func_cnv")),
                                                              lapply(lrn_tanh$names, function(x) c(x, "func_mrna")), lapply(lrn_tanh$names, function(x) c(x, "func_cnv")),
                                                              lapply(lrn_anova$names, function(x) c(x, "func_mrna")), lapply(lrn_anova$names, function(x) c(x, "func_cnv"))), append),
                              ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl__model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl__model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

# sl_ksvm_late <- mcSuperLearner(Y = Y, X = W_, family = binomial(), cvControl = list(V = 10),
# SL.library = purrr::reduce(list(lapply(lrn_rbf$names, function(x) c(x, "func_mrna")), lapply(lrn_rbf$names, function(x) c(x, "func_cnv")),
#                       lapply(lrn_poly$names, function(x) c(x, "func_mrna")), lapply(lrn_poly$names, function(x) c(x, "func_cnv")),
#                       lapply(lrn_tanh$names, function(x) c(x, "func_mrna")), lapply(lrn_tanh$names, function(x) c(x, "func_cnv")),
#                       lapply(lrn_anova$names, function(x) c(x, "func_mrna")), lapply(lrn_anova$names, function(x) c(x, "func_cnv"))), append) )


### caso 3: xgboost -> Resultado: Funciona
SL.LinearComb.xgb <- function(Y, X, family, ...) {
  # load required packages
  require("xgboost")
  
  tune = list(ntrees = 5000,
              max_depth = c(25,100), 
              shrinkage = c(0.001, 0.01))
  
  lrn_xgb <- create.Learner("SL.xgboost", tune = tune, detailed_names = TRUE, name_prefix = "xgb")
  
  sl__model <- mcSuperLearner(Y = Y, X = X, family = family$family, #cvControl = list(V = 10),
                              SL.library = append(lapply(lrn_xgb$names, function(x) c(x, "func_mrna")), lapply(lrn_xgb$names, function(x) c(x, "func_cnv"))),
                              ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl__model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl__model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

# sl_xgb_late <- mcSuperLearner(Y = Y, X = W_, family = binomial(), cvControl = list(V = 10),
#                          SL.library = append(lapply(lrn_xgb$names, function(x) c(x, "func_mrna")), lapply(lrn_xgb$names, function(x) c(x, "func_cnv"))))


## caso 4: mars -> Resultado: Funciona
SL.LinearComb.mars <- function(Y, X, family, ...) {
  # load required packages
  require("earth")
  
  tune = list(degree = c(2,3),
              nprune = c(2,3,5),
              Use.beta.cache=FALSE) # Note: earth's beta cache would require 3.63 GB, so forcing Use.beta.cache=FALSE.
  
  new.screen.randomForest_mrna <- function(...) {
    screen.randomForest(..., nVar = 200) & func_mrna(X, ...)
  }
  
  new.screen.randomForest_cnv <- function(...) {
    screen.randomForest(..., nVar = 200) & func_cnv(X, ...)
  }
  
  lrn_mars <- create.Learner("SL.earth", tune = tune, name_prefix = "mars")
  
  sl__model <- mcSuperLearner(Y = Y, X = X, family = family$family, #cvControl = list(V = 10),
                              SL.library = append(lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest_mrna")), 
                                                  lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest_cnv"))),
                              ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl__model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl__model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

# sl_mars_late <- mcSuperLearner(Y = Y, X = W_, family = binomial(), cvControl = list(V = 10),
#                           SL.library = append(lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest_mrna")), 
#                                               lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest_cnv"))))


sl_model_late <- CV.SuperLearner(Y = Y, X = W_, family = binomial(), 
                                 # For a real analysis we would use V = 10.
                                 cvControl = list(V = 5), innerCvControl = list(list(V=5)),
                                 parallel = "multicore", 
                                 SL.library = list("SL.LinearComb.rf", "SL.LinearComb.SVM", "SL.LinearComb.xgb", "SL.LinearComb.mars") )

Q1W <- SuperLearner::predict.SuperLearner(sl_model_late,newdata=data.frame(A=1, W_ %>% dplyr::select(-"A")), type="response")$pred
Q0W <- SuperLearner::predict.SuperLearner(sl_model_late,newdata=data.frame(A=0, W_ %>% dplyr::select(-"A")), type = "response")$pred

Q <- data.frame(Q0W,Q1W)
head(Q); dim(Q)

## ------------------------ ##
### ... g. tmle update ... ###
## ------------------------ ##
result.Qcgc <- tmle(Y = Y, A = A, W = W_, verbose = TRUE, Delta = Delta,
                    Q = Q, g1W = g1W, pDelta1 = pDelta1, family = "binomial")

result.Qcgc

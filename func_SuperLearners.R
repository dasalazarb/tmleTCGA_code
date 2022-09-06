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
#### ------------------------------------------------ ###
#### ------------------------------------------------ ###
#### ------------------------------------------------ ###
### caso 1: random forest -> Resultado: Funciona
SL.clinical.rf <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.clinical.rf "))
  # load required packages
  require("ranger")
  
  tune = list(num.trees = c(1000, 5000), 
              min.node.size = c(5, 25, 100),
              mtry = floor(sqrt(sum(func_clinical(X))) * c(0.5,0.75,1)))
  
  lrn_rf <- create.Learner("SL.ranger", tune = tune)
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.rf.mrna", "SL.rf.cnv"),
                               SL.library = lapply(lrn_rf$names, function(x) c(x, "func_clinical")),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

### caso 3: xgboost -> Resultado: Funciona
SL.clinical.xgb <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.clinical.xgb "))
  # load required packages
  require("xgboost")
  
  tune = list(ntrees = 2000,
              max_depth = 50, 
              shrinkage = 0.001)
  
  lrn_xgb <- create.Learner("SL.xgboost", params = tune, detailed_names = TRUE, name_prefix = "xgb")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.xgboost.mrna", "SL.xgboost.cnv"),
                               SL.library = lapply(lrn_xgb$names, function(x) c(x, "func_clinical")),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

## caso 4: mars -> Resultado: Funciona
SL.clinical.mars <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.clinical.mars "))
  require("earth")
  
  tune = list(degree = 1:2, pmethod = "backward",
              nprune = seq(10, 50, length.out = 5) %>% floor(),
              Use.beta.cache=FALSE, newvar.penalty=0.02) # Note: earth's beta cache would require 3.63 GB, so forcing Use.beta.cache=FALSE, newvar.penalty=0.02.
  
  lrn_mars <- create.Learner("SL.earth", tune = tune, name_prefix = "clinical.mars")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = lapply(lrn_mars$names, function(x) c(x, "func_clinical")),
                               # SL.library = append(lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna")), 
                               #                     lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

## caso 6: kernelknn -> Resultado: Funciona
SL.clinical.kknn <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.clinical.kknn "))
  require("KernelKnn")
  
  tune = list(k=c(3,10), method=c("euclidian", "manhattan", "mahalanobis"), weights_function="gaussian", transf_categ_cols = TRUE)
  
  lrn_kknn <- create.Learner("SL.kernelKnn", tune = tune, name_prefix = "KKnn")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = lapply(lrn_kknn$names, function(x) c(x, "func_clinical")),
                               # SL.library = append(lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna")), 
                               #                     lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

## caso 7: bayesglm -> Resultado: Funciona
SL.clinical.bayesglm <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.clinical.bayesglm "))
  require("arm")
  
  lrn_bayes <- create.Learner("SL.bayesglm", name_prefix = "bayesglm")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = lapply(lrn_bayes$names, function(x) c(x, "func_clinical")),
                               # SL.library = append(lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna")), 
                               #                     lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

#### ------------------------------------------------ ###
#### ------------------------------------------------ ###
#### ------------------------------------------------ ###
## caso early: mars -> Resultado: Funciona
SL.early.mars <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.early.mars "))
  # load required packages
  require("earth")
  
  tune = list(degree = 1, pmethod = "backward",
              nprune = seq(10, 50, length.out = 5) %>% floor(),
              Use.beta.cache=FALSE, newvar.penalty=0.02) # Note: earth's beta cache would require 3.63 GB, so forcing Use.beta.cache=FALSE, newvar.penalty=0.02.
  
  lrn_mars <- create.Learner("SL.earth", tune = tune, name_prefix = "early.mars")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest")),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso late: kernelSVM -> Resultado: Funciona
SL.early.SVM <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.early.SVM "))
  # load required packages
  require("kernlab")
  
  tune <- list(kernel = c("rbfdot"), 
               kpar = list(list(sigma=2)))
  lrn_rbf <- create.Learner("SL.ksvm", params = tune, name_prefix = "early.rbf")
  
  tune <- list(kernel = c("laplacedot"), 
               kpar = list(list(sigma=2)))
  lrn_laplace <- create.Learner("SL.ksvm", params = tune, name_prefix = "early.laplace")
  
  tune <- list(kernel = c("polydot"), 
               kpar = list(list(degree=2)))
  lrn_poly <- create.Learner("SL.ksvm", params = tune, name_prefix = "early.poly")
  
  tune <- list(kernel = c("tanhdot"), 
               kpar = list(list(scale=2, offset=1)))
  lrn_tanh <- create.Learner("SL.ksvm", params = tune, name_prefix = "early.tanh")
  
  tune <- list(kernel = c("anovadot"), 
               kpar = list(list(sigma=2, degree=2)))
  lrn_anova <- create.Learner("SL.ksvm", params = tune, name_prefix = "early.anova")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = purrr::reduce(list(lapply(lrn_rbf$names, function(x) c(x, "new.screen.randomForest")),
                                                               lapply(lrn_laplace$names, function(x) c(x, "new.screen.randomForest")),
                                                               lapply(lrn_poly$names, function(x) c(x, "new.screen.randomForest")), 
                                                               lapply(lrn_tanh$names, function(x) c(x, "new.screen.randomForest")),
                                                               lapply(lrn_anova$names, function(x) c(x, "new.screen.randomForest"))), append),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso early: xgboost -> Resultado: Funciona
SL.early.xgb <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.early.xgb "))
  require("xgboost")
  
  tune = list(ntrees = 2000,
              max_depth = 50, 
              shrinkage = 0.001)
  
  lrn_xgb <- create.Learner("SL.xgboost", params = tune, detailed_names = TRUE, name_prefix = "xgb")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.xgboost.mrna", "SL.xgboost.cnv"),
                               SL.library = lapply(lrn_xgb$names, function(x) c(x)),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso early: random forest -> Resultado: Funciona
SL.early.rf <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.early.rf "))
  # load required packages
  require("ranger")
  
  tune = list(num.trees = c(1000, 5000), 
              min.node.size = c(5, 25, 100),
              mtry = floor(sqrt(ncol(X)) * c(0.5,0.75,1)))
  
  lrn_rf <- create.Learner("SL.ranger", tune = tune)
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.rf.mrna", "SL.rf.cnv"),
                               SL.library = lapply(lrn_rf$names, function(x) c(x)),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso early: KKnn -> Resultado: Funciona
SL.early.KKnn <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.early.KKnn "))
  require("KernelKnn")
  
  tune = list(k=c(3,10), method=c("euclidian", "manhattan", "mahalanobis"), weights_function="gaussian", transf_categ_cols = TRUE)
  
  lrn_kknn <- create.Learner("SL.kernelKnn", tune = tune, name_prefix = "KKnn")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.rf.mrna", "SL.rf.cnv"),
                               SL.library = lapply(lrn_kknn$names, function(x) c(x)),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso early: bayesglm -> Resultado: Funciona
SL.early.bayesglm <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.early.bayesglm "))
  require("arm")
  
  lrn_bayesglm <- create.Learner("SL.bayesglm", name_prefix = "bayesglm")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.rf.mrna", "SL.rf.cnv"),
                               SL.library = lapply(lrn_bayesglm$names, function(x) c(x)),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

#### ------------------------------------------------ ###
#### ------------------------------------------------ ###
#### ------------------------------------------------ ###
### caso late: random forest -> Resultado: Funciona
SL.Late.rf <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Late.rf "))
  # load required packages
  require("ranger")
  
  tune = list(num.trees = c(1000, 5000), 
              min.node.size = c(5, 25, 100),
              mtry = floor(sqrt(ncol(X)) * c(0.5,0.75,1)))
  
  lrn_rf <- create.Learner("SL.ranger", tune = tune)
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.rf.mrna", "SL.rf.cnv"),
                               SL.library = append(lapply(lrn_rf$names, function(x) c(x, "new.screen.rf.mrna")), lapply(lrn_rf$names, function(x) c(x, "new.screen.rf.cnv"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso late: kernelSVM -> Resultado: Funciona
SL.Late.SVM <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Late.SVM "))
  # load required packages
  require("kernlab")
  
  tune <- list(kernel = c("rbfdot"), 
               kpar = list(list(sigma=2)))
  lrn_rbf <- create.Learner("SL.ksvm", params = tune, name_prefix = "late.rbf")
  
  tune <- list(kernel = c("laplacedot"), 
               kpar = list(list(sigma=2)))
  lrn_laplace <- create.Learner("SL.ksvm", params = tune, name_prefix = "late.laplace")
  
  tune <- list(kernel = c("polydot"), 
               kpar = list(list(degree=2)))
  lrn_poly <- create.Learner("SL.ksvm", params = tune, name_prefix = "late.poly")
  
  tune <- list(kernel = c("tanhdot"), 
               kpar = list(list(scale=2, offset=1)))
  lrn_tanh <- create.Learner("SL.ksvm", params = tune, name_prefix = "late.tanh")
  
  tune <- list(kernel = c("anovadot"), 
               kpar = list(list(sigma=2, degree=2)))
  lrn_anova <- create.Learner("SL.ksvm", params = tune, name_prefix = "late.anova")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = purrr::reduce(list(lapply(lrn_rbf$names, function(x) c(x, "new.screen.rf.mrna")), lapply(lrn_rbf$names, function(x) c(x, "new.screen.rf.cnv")),
                                                               lapply(lrn_laplace$names, function(x) c(x, "new.screen.rf.mrna")), lapply(lrn_laplace$names, function(x) c(x, "new.screen.rf.cnv")),
                                                               lapply(lrn_poly$names, function(x) c(x, "new.screen.rf.mrna")), lapply(lrn_poly$names, function(x) c(x, "new.screen.rf.cnv")),
                                                               lapply(lrn_tanh$names, function(x) c(x, "new.screen.rf.mrna")), lapply(lrn_tanh$names, function(x) c(x, "new.screen.rf.cnv")),
                                                               lapply(lrn_anova$names, function(x) c(x, "new.screen.rf.mrna")), lapply(lrn_anova$names, function(x) c(x, "new.screen.rf.cnv"))), append),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso late: xgboost -> Resultado: Funciona
SL.Late.xgb <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Late.xgb "))
  # load required packages
  require("xgboost")
  
  tune = list(ntrees = 2000,
              max_depth = 50, 
              shrinkage = 0.001)
  
  lrn_xgb <- create.Learner("SL.xgboost", params = tune, detailed_names = TRUE, name_prefix = "xgb")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.xgboost.mrna", "SL.xgboost.cnv"),
                               SL.library = append(lapply(lrn_xgb$names, function(x) c(x, "new.screen.rf.mrna")), lapply(lrn_xgb$names, function(x) c(x, "new.screen.rf.cnv"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


## caso late: mars -> Resultado: Funciona
SL.Late.mars <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Late.mars "))
  # load required packages
  require("earth")
  
  tune = list(degree = 1, pmethod = "backward",
              nprune = seq(10, 50, length.out = 5) %>% floor(),
              Use.beta.cache=FALSE, newvar.penalty=0.02) # Note: earth's beta cache would require 3.63 GB, so forcing Use.beta.cache=FALSE, newvar.penalty=0.02.
  
  lrn_mars <- create.Learner("SL.earth", tune = tune, name_prefix = "late.mars")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.earth.mrna","SL.earth.cnv"),
                               SL.library = append(lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.mrna")),
                                                   lapply(lrn_mars$names, function(x) c(x, "new.screen.rf.cnv"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso late: KKnn -> Resultado: Funciona
SL.Late.KKnn <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Late.KKnn "))
  # load required packages
  require("KernelKnn")
  
  tune = list(k=c(3,10), method=c("euclidian", "manhattan", "mahalanobis"), weights_function="gaussian", transf_categ_cols = TRUE)
  
  lrn_kknn <- create.Learner("SL.kernelKnn", tune = tune, name_prefix = "KKnn")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.rf.mrna", "SL.rf.cnv"),
                               SL.library = append(lapply(lrn_kknn$names, function(x) c(x, "new.screen.rf.mrna")), lapply(lrn_kknn$names, function(x) c(x, "new.screen.rf.cnv"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso late: bayesglm -> Resultado: Funciona
SL.Late.bayesglm <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Late.bayesglm "))
  # load required packages
  require("arm")
  
  lrn_bayesglm <- create.Learner("SL.bayesglm", name_prefix = "bayesglm")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.rf.mrna", "SL.rf.cnv"),
                               SL.library = append(lapply(lrn_bayesglm$names, function(x) c(x, "new.screen.rf.mrna")), lapply(lrn_bayesglm$names, function(x) c(x, "new.screen.rf.cnv"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

#### ------------------------------------------------ ###
#### ------------------------------------------------ ###
#### ------------------------------------------------ ###
#################### encoder
### caso encoder: random forest -> Resultado: Funciona
SL.Int.rf.encoder <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.rf.encoder "))
  # load required packages
  require("ranger")
  
  tune = list(num.trees = c(1000, 5000), 
              min.node.size = c(5, 25, 100),
              mtry = floor(sqrt(sum(func_mrna.encoder(X))) * c(0.5,0.75,1)))
  
  lrn_rf.mrna <- create.Learner("SL.ranger", tune = tune)
  
  tune = list(num.trees = c(1000, 5000), 
              min.node.size = c(5, 25, 100),
              mtry = floor(sqrt(sum(func_cnv.encoder(X))) * c(0.5,0.75,1)))
  
  lrn_rf.cnv <- create.Learner("SL.ranger", tune = tune)
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_rf.mrna$names, function(x) c(x, "func_mrna.encoder")), lapply(lrn_rf.cnv$names, function(x) c(x, "func_cnv.encoder"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso encoder: xgboost -> Resultado: Funciona
SL.Int.xgb.encoder <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.xgb.encoder "))
  require("xgboost")
  
  tune = list(ntrees = 2000,
              max_depth = 50, 
              shrinkage = 0.001)
  
  lrn_xgb <- create.Learner("SL.xgboost", params = tune, detailed_names = TRUE, name_prefix = "xgb")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               # SL.library = c("SL.xgboost.mrna.encoder", "SL.xgboost.cnv.encoder"),
                               SL.library = append(lapply(lrn_xgb$names, function(x) c(x, "func_mrna.encoder")), lapply(lrn_xgb$names, function(x) c(x, "func_cnv.encoder"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso encoder: kernelSVM -> Resultado: Funciona
SL.Int.SVM.encoder <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.SVM.encoder "))
  # load required packages
  require("kernlab")
  
  tune <- list(kernel = c("rbfdot"), 
               kpar = list(list(sigma=2)))
  lrn_rbf <- create.Learner("SL.ksvm", params = tune, name_prefix = "encoder.rbf")
  
  tune <- list(kernel = c("laplacedot"), 
               kpar = list(list(sigma=2)))
  lrn_laplace <- create.Learner("SL.ksvm", params = tune, name_prefix = "encoder.laplace")
  
  tune <- list(kernel = c("polydot"), 
               kpar = list(list(degree=2)))
  lrn_poly <- create.Learner("SL.ksvm", params = tune, name_prefix = "encoder.poly")
  
  tune <- list(kernel = c("tanhdot"), 
               kpar = list(list(scale=2, offset=1)))
  lrn_tanh <- create.Learner("SL.ksvm", params = tune, name_prefix = "encoder.tanh")
  
  tune <- list(kernel = c("anovadot"), 
               kpar = list(list(sigma=2, degree=2)))
  lrn_anova <- create.Learner("SL.ksvm", params = tune, name_prefix = "encoder.anova")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = purrr::reduce(list(lapply(lrn_rbf$names, function(x) c(x, "func_mrna.encoder")), lapply(lrn_rbf$names, function(x) c(x, "func_cnv.encoder")),
                                                               lapply(lrn_laplace$names, function(x) c(x, "func_mrna.encoder")), lapply(lrn_laplace$names, function(x) c(x, "func_cnv.encoder")),
                                                               lapply(lrn_poly$names, function(x) c(x, "func_mrna.encoder")), lapply(lrn_poly$names, function(x) c(x, "func_cnv.encoder")),
                                                               lapply(lrn_tanh$names, function(x) c(x, "func_mrna.encoder")), lapply(lrn_tanh$names, function(x) c(x, "func_cnv.encoder")),
                                                               lapply(lrn_anova$names, function(x) c(x, "func_mrna.encoder")), lapply(lrn_anova$names, function(x) c(x, "func_cnv.encoder"))), append),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


## caso encoder: mars -> Resultado: Funciona
SL.Int.mars.encoder <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.mars.encoder "))
  require("earth")
  
  tune = list(degree = 1, pmethod = "backward",
              nprune = seq(10, 50, length.out = 5) %>% floor(),
              Use.beta.cache=FALSE, newvar.penalty=0.02) # Note: earth's beta cache would require 3.63 GB, so forcing Use.beta.cache=FALSE, newvar.penalty=0.02.
  
  lrn_mars <- create.Learner("SL.earth", tune = tune, name_prefix = "int.mars.encoder")
  
  new.screen.randomForest_mrna <- function(...) {
    screen.randomForest(..., nVar = 250) & func_mrna.encoder(X, ...)
  }
  
  new.screen.randomForest_cnv <- function(...) {
    screen.randomForest(..., nVar = 250) & func_cnv.encoder(X, ...)
  }
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest_mrna")), 
                                                   lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest_cnv"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso encoder: KKnn -> Resultado: Funciona
SL.Int.KKnn.encoder <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.KKnn.encoder "))
  # load required packages
  require("KernelKnn")
  
  tune = list(k=c(3,10), method=c("euclidian", "manhattan", "mahalanobis"), weights_function="gaussian", transf_categ_cols = TRUE)
  
  lrn_kknn.mrna <- create.Learner("SL.kernelKnn", tune = tune, name_prefix = "KKnn.mrna")
  
  tune = list(k=c(3,10), method=c("euclidian", "manhattan", "mahalanobis"), weights_function="gaussian", transf_categ_cols = TRUE)
  
  lrn_kknn.cnv <- create.Learner("SL.kernelKnn", tune = tune, name_prefix = "KKnn.cnv")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_kknn.mrna$names, function(x) c(x, "func_mrna.encoder")), lapply(lrn_kknn.cnv$names, function(x) c(x, "func_cnv.encoder"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso encoder: bayesglm -> Resultado: Funciona
SL.Int.bayesglm.encoder <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.bayesglm.encoder "))
  # load required packages
  require("arm")
  
  lrn_bayesglm.mrna <- create.Learner("SL.bayesglm", name_prefix = "bayesglm.mrna")
  
  lrn_bayesglm.cnv <- create.Learner("SL.bayesglm", name_prefix = "bayesglm.cnv")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_bayesglm.mrna$names, function(x) c(x, "func_mrna.encoder")), lapply(lrn_bayesglm.cnv$names, function(x) c(x, "func_cnv.encoder"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

#################### pca
### caso pca: random forest -> Resultado: Funciona
SL.Int.rf.pca <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.rf.pca "))
  # load required packages
  require("ranger")
  
  tune = list(num.trees = c(1000, 5000), 
              min.node.size = c(5, 25, 100),
              mtry = floor(sqrt(sum(func_mrna.pca(X))) * c(0.5,0.75,1)))
  
  lrn_rf.mrna <- create.Learner("SL.ranger", tune = tune)
  
  tune = list(num.trees = c(1000, 5000), 
              min.node.size = c(5, 25, 100),
              mtry = floor(sqrt(sum(func_cnv.pca(X))) * c(0.5,0.75,1)))
  
  lrn_rf.cnv <- create.Learner("SL.ranger", tune = tune)
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_rf.mrna$names, function(x) c(x, "func_mrna.pca")), lapply(lrn_rf.cnv$names, function(x) c(x, "func_cnv.pca"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso pca: xgboost -> Resultado: Funciona
SL.Int.xgb.pca <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.xgb.pca "))
  # load required packages
  require("xgboost")
  
  tune = list(ntrees = 2000,
              max_depth = 50, 
              shrinkage = 0.001)
  
  lrn_xgb <- create.Learner("SL.xgboost", params = tune, detailed_names = TRUE, name_prefix = "xgb")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_xgb$names, function(x) c(x, "func_mrna.pca")), lapply(lrn_xgb$names, function(x) c(x, "func_cnv.pca"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso pca: kernelSVM -> Resultado: Funciona
SL.Int.SVM.pca <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.SVM.pca "))
  # load required packages
  require("kernlab")
  
  tune <- list(kernel = c("rbfdot"), 
               kpar = list(list(sigma=2)))
  lrn_rbf <- create.Learner("SL.ksvm", params = tune, name_prefix = "pca.rbf")
  
  tune <- list(kernel = c("laplacedot"), 
               kpar = list(list(sigma=2)))
  lrn_laplace <- create.Learner("SL.ksvm", params = tune, name_prefix = "pca.laplace")
  
  tune <- list(kernel = c("polydot"), 
               kpar = list(list(degree=2)))
  lrn_poly <- create.Learner("SL.ksvm", params = tune, name_prefix = "pca.poly")
  
  tune <- list(kernel = c("tanhdot"), 
               kpar = list(list(scale=2, offset=1)))
  lrn_tanh <- create.Learner("SL.ksvm", params = tune, name_prefix = "pca.tanh")
  
  tune <- list(kernel = c("anovadot"), 
               kpar = list(list(sigma=2, degree=2)))
  lrn_anova <- create.Learner("SL.ksvm", params = tune, name_prefix = "pca.anova")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = purrr::reduce(list(lapply(lrn_rbf$names, function(x) c(x, "func_mrna.pca")), lapply(lrn_rbf$names, function(x) c(x, "func_cnv.pca")),
                                                               lapply(lrn_laplace$names, function(x) c(x, "func_mrna.pca")), lapply(lrn_laplace$names, function(x) c(x, "func_cnv.pca")),
                                                               lapply(lrn_poly$names, function(x) c(x, "func_mrna.pca")), lapply(lrn_poly$names, function(x) c(x, "func_cnv.pca")),
                                                               lapply(lrn_tanh$names, function(x) c(x, "func_mrna.pca")), lapply(lrn_tanh$names, function(x) c(x, "func_cnv.pca")),
                                                               lapply(lrn_anova$names, function(x) c(x, "func_mrna.pca")), lapply(lrn_anova$names, function(x) c(x, "func_cnv.pca"))), append),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


## caso pca: mars -> Resultado: Funciona
SL.Int.mars.pca <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.mars.pca "))
  require("earth")
  
  tune = list(degree = 1, pmethod = "backward",
              nprune = seq(10, 50, length.out = 5) %>% floor(),
              Use.beta.cache=FALSE, newvar.penalty=0.02) # Note: earth's beta cache would require 3.63 GB, so forcing Use.beta.cache=FALSE, newvar.penalty=0.02.
  
  lrn_mars <- create.Learner("SL.earth", tune = tune, name_prefix = "int.mars.pca")
  
  new.screen.randomForest_mrna <- function(...) {
    screen.randomForest(..., nVar = 250) & func_mrna.pca(X, ...)
  }
  
  new.screen.randomForest_cnv <- function(...) {
    screen.randomForest(..., nVar = 250) & func_cnv.pca(X, ...)
  }
  
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest_mrna")), 
                                                   lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest_cnv"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso pca: KKnn -> Resultado: Funciona
SL.Int.KKnn.pca <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.KKnn.pca "))
  # load required packages
  require("KernelKnn")
  
  tune = list(k=c(3,10), method=c("euclidian", "manhattan", "mahalanobis"), weights_function="gaussian", transf_categ_cols = TRUE)
  
  lrn_kknn.mrna <- create.Learner("SL.kernelKnn", tune = tune, name_prefix = "KKnn.mrna")
  
  tune = list(k=c(3,10), method=c("euclidian", "manhattan", "mahalanobis"), weights_function="gaussian", transf_categ_cols = TRUE)
  
  lrn_kknn.cnv <- create.Learner("SL.kernelKnn", tune = tune, name_prefix = "KKnn.cnv")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_kknn.mrna$names, function(x) c(x, "func_mrna.pca")), lapply(lrn_kknn.cnv$names, function(x) c(x, "func_cnv.pca"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

### caso pca: bayesglm -> Resultado: Funciona
SL.Int.bayesglm.pca <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.bayesglm.pca "))
  # load required packages
  require("arm")
  
  lrn_bayesglm.mrna <- create.Learner("SL.bayesglm", name_prefix = "bayesglm.mrna")
  
  lrn_bayesglm.cnv <- create.Learner("SL.bayesglm", name_prefix = "bayesglm.cnv")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_bayesglm.mrna$names, function(x) c(x, "func_mrna.pca")), lapply(lrn_bayesglm.cnv$names, function(x) c(x, "func_cnv.pca"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

#################### nmf
### caso nmf: random forest -> Resultado: Funciona
SL.Int.rf.nmf <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.rf.nmf "))
  # load required packages
  require("ranger")
  
  tune = list(num.trees = c(1000, 5000), 
              min.node.size = c(5, 25, 100),
              mtry = floor(sqrt(sum(func_mrna.nmf(X))) * c(0.5,0.75,1)))
  
  lrn_rf.mrna <- create.Learner("SL.ranger", tune = tune)
  
  tune = list(num.trees = c(1000, 5000), 
              min.node.size = c(5, 25, 100),
              mtry = floor(sqrt(sum(func_cnv.nmf(X))) * c(0.5,0.75,1)))
  
  lrn_rf.cnv <- create.Learner("SL.ranger", tune = tune)
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_rf.mrna$names, function(x) c(x, "func_mrna.nmf")), lapply(lrn_rf.cnv$names, function(x) c(x, "func_cnv.nmf"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso nmf: xgboost -> Resultado: Funciona
SL.Int.xgb.nmf <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training thelearner SL.Int.xgb.nmf "))
  # load required packages
  require("xgboost")
  
  tune = list(ntrees = 2000,
              max_depth = 50, 
              shrinkage = 0.001)
  
  lrn_xgb <- create.Learner("SL.xgboost", params = tune, detailed_names = TRUE, name_prefix = "xgb")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_xgb$names, function(x) c(x, "func_mrna.nmf")), lapply(lrn_xgb$names, function(x) c(x, "func_cnv.nmf"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


## caso nmf: mars -> Resultado: Funciona
SL.Int.mars.nmf <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training thelearner SL.Int.mars.nmf "))
  # load required packages
  require("earth")
  
  tune = list(degree = 1, pmethod = "backward",
              nprune = seq(10, 50, length.out = 5) %>% floor(),
              Use.beta.cache=FALSE, newvar.penalty=0.02) # Note: earth's beta cache would require 3.63 GB, so forcing Use.beta.cache=FALSE, newvar.penalty=0.02.
  
  lrn_mars <- create.Learner("SL.earth", tune = tune, name_prefix = "int.mars.nmf")
  
  new.screen.randomForest_mrna <- function(...) {
    screen.randomForest(..., nVar = 250) & func_mrna.nmf(X, ...)
  }
  
  new.screen.randomForest_cnv <- function(...) {
    screen.randomForest(..., nVar = 250) & func_cnv.nmf(X, ...)
  }
  
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest_mrna")), 
                                                   lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest_cnv"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso nmf: kernelSVM -> Resultado: Funciona
SL.Int.SVM.nmf <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.SVM.nmf "))
  # load required packages
  require("kernlab")
  
  tune <- list(kernel = c("rbfdot"), 
               kpar = list(list(sigma=2)))
  lrn_rbf <- create.Learner("SL.ksvm", params = tune, name_prefix = "nmf.rbf")
  
  tune <- list(kernel = c("laplacedot"), 
               kpar = list(list(sigma=2)))
  lrn_laplace <- create.Learner("SL.ksvm", params = tune, name_prefix = "nmf.laplace")
  
  tune <- list(kernel = c("polydot"), 
               kpar = list(list(degree=2)))
  lrn_poly <- create.Learner("SL.ksvm", params = tune, name_prefix = "nmf.poly")
  
  tune <- list(kernel = c("tanhdot"), 
               kpar = list(list(scale=2, offset=1)))
  lrn_tanh <- create.Learner("SL.ksvm", params = tune, name_prefix = "nmf.tanh")
  
  tune <- list(kernel = c("anovadot"), 
               kpar = list(list(sigma=2, degree=2)))
  lrn_anova <- create.Learner("SL.ksvm", params = tune, name_prefix = "nmf.anova")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = purrr::reduce(list(lapply(lrn_rbf$names, function(x) c(x, "func_mrna.nmf")), lapply(lrn_rbf$names, function(x) c(x, "func_cnv.nmf")),
                                                               lapply(lrn_laplace$names, function(x) c(x, "func_mrna.nmf")), lapply(lrn_laplace$names, function(x) c(x, "func_cnv.nmf")),
                                                               lapply(lrn_poly$names, function(x) c(x, "func_mrna.nmf")), lapply(lrn_poly$names, function(x) c(x, "func_cnv.nmf")),
                                                               lapply(lrn_tanh$names, function(x) c(x, "func_mrna.nmf")), lapply(lrn_tanh$names, function(x) c(x, "func_cnv.nmf")),
                                                               lapply(lrn_anova$names, function(x) c(x, "func_mrna.nmf")), lapply(lrn_anova$names, function(x) c(x, "func_cnv.nmf"))), append),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso nmf: KKnn -> Resultado: Funciona
SL.Int.KKnn.nmf <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.KKnn.nmf "))
  require("KernelKnn")
  
  tune = list(k=c(3,10), method=c("euclidian", "manhattan", "mahalanobis"), weights_function="gaussian", transf_categ_cols = TRUE)
  
  lrn_kknn.mrna <- create.Learner("SL.kernelKnn", tune = tune, name_prefix = "KKnn.mrna")
  
  tune = list(k=c(3,10), method=c("euclidian", "manhattan", "mahalanobis"), weights_function="gaussian", transf_categ_cols = TRUE)
  
  lrn_kknn.cnv <- create.Learner("SL.kernelKnn", tune = tune, name_prefix = "KKnn.cnv")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_kknn.mrna$names, function(x) c(x, "func_mrna.nmf")), lapply(lrn_kknn.cnv$names, function(x) c(x, "func_cnv.nmf"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}


### caso nmf: bayesglm -> Resultado: Funciona
SL.Int.bayesglm.nmf <- function(Y, X, family, ...) {
  print(paste0(outcome, "_", cancer, "_data_var", nvar, " ... -> -> -> Training the learner SL.Int.bayesglm.nmf "))
  require("arm")
  
  lrn_bayesglm.mrna <- create.Learner("SL.bayesglm", name_prefix = "bayesglm.mrna")
  
  lrn_bayesglm.cnv <- create.Learner("SL.bayesglm", name_prefix = "bayesglm.cnv")
  
  sl___model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE), #cvControl = list(stratifyCV = TRUE),
                               SL.library = append(lapply(lrn_bayesglm.mrna$names, function(x) c(x, "func_mrna.nmf")), lapply(lrn_bayesglm.cnv$names, function(x) c(x, "func_cnv.nmf"))),
                               method = "method.NNloglik", ...)
  
  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <-  sl___model$SL.predict
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = sl___model)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.template'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}
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
### mars with selection of variables
SL.earth.FeatureSelection <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                                       nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                                       nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) 
{
  X <- X[,new.screen.randomForest(X)]
  require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan, 
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}

### clinical
SL.rf.clinical <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
                        write.forest = TRUE, probability = family$family == "binomial", 
                        min.node.size = ifelse(family$family == "gaussian", 5, 1), 
                        replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = T, ...) 
{
  X <- X[,func_clinical(X)]
  require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, 
                                               X), num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = probability, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.xgboost.clinical <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
                             max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                             nthread = 1, verbose = 0, save_period = NULL, ...) 
{
  X <- X[,func_clinical(X)]
  names(X)  <- paste0('x', 1:ncol(X))
  require("xgboost")
  if (packageVersion("xgboost") < 0.6) 
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    if (packageVersion("xgboost") >= "1.1.1.1") {
      objective <- "reg:squarederror"
    }
    else {
      objective <- "reg:linear"
    }
    model = xgboost::xgboost(data = xgmat, objective = objective, 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period, eval_metric = "logloss")
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

SL.earth.clinical <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                           nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                           nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) 
{
  X <- X[,func_clinical(X)]
  require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan, 
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}

SL.glm.clinical <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  X <- X[,func_clinical(X)]
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  fit.glm <- glm(Y ~ ., data = X, family = family, weights = obsWeights, 
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

### rf.mrna + xgboost.mrna + earth.mrna
SL.rf.mrna <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
                            write.forest = TRUE, probability = family$family == "binomial", 
                            min.node.size = ifelse(family$family == "gaussian", 5, 1), 
                            replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = T, ...) 
{
  X <- X[,func_mrna(X)]
  require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, 
                                               X), num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = probability, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.xgboost.mrna <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
                                 max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                                 nthread = 1, verbose = 0, save_period = NULL, ...) 
{
  X <- X[,func_mrna(X)]
  names(X)  <- paste0('x', 1:ncol(X))
  require("xgboost")
  if (packageVersion("xgboost") < 0.6) 
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    if (packageVersion("xgboost") >= "1.1.1.1") {
      objective <- "reg:squarederror"
    }
    else {
      objective <- "reg:linear"
    }
    model = xgboost::xgboost(data = xgmat, objective = objective, 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period, eval_metric = "logloss")
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

SL.earth.mrna <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                               nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                               nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) 
{
  X <- X[,func_mrna(X)]
  require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan, 
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}

SL.glm.mrna <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  X <- X[,func_mrna(X)]
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  fit.glm <- glm(Y ~ ., data = X, family = family, weights = obsWeights, 
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.rf.mrna.pca <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
                       write.forest = TRUE, probability = family$family == "binomial", 
                       min.node.size = ifelse(family$family == "gaussian", 5, 1), 
                       replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = T, ...) 
{
  X <- X[,func_mrna.pca(X)]
  require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, 
                                               X), num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = probability, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.rf.mrna.encoder <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
                       write.forest = TRUE, probability = family$family == "binomial", 
                       min.node.size = ifelse(family$family == "gaussian", 5, 1), 
                       replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = T, ...) 
{
  X <- X[,func_mrna.encoder(X)]
  require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, 
                                               X), num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = probability, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.rf.mrna.nmf <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
                       write.forest = TRUE, probability = family$family == "binomial", 
                       min.node.size = ifelse(family$family == "gaussian", 5, 1), 
                       replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = T, ...) 
{
  X <- X[,func_mrna.nmf(X)]
  require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, 
                                               X), num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = probability, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.xgboost.mrna.pca <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
          max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
          nthread = 1, verbose = 0, save_period = NULL, ...) 
{
  X <- X[,func_mrna.pca(X)]
  names(X)  <- paste0('x', 1:ncol(X))
  require("xgboost")
  if (packageVersion("xgboost") < 0.6) 
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    if (packageVersion("xgboost") >= "1.1.1.1") {
      objective <- "reg:squarederror"
    }
    else {
      objective <- "reg:linear"
    }
    model = xgboost::xgboost(data = xgmat, objective = objective, 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period, eval_metric = "logloss")
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

SL.xgboost.mrna.encoder <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
                            max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                            nthread = 1, verbose = 0, save_period = NULL, ...) 
{
  X <- X[,func_mrna.encoder(X)]
  names(X)  <- paste0('x', 1:ncol(X))
  require("xgboost")
  if (packageVersion("xgboost") < 0.6) 
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    if (packageVersion("xgboost") >= "1.1.1.1") {
      objective <- "reg:squarederror"
    }
    else {
      objective <- "reg:linear"
    }
    model = xgboost::xgboost(data = xgmat, objective = objective, 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period, eval_metric = "logloss")
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

SL.xgboost.mrna.nmf <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
                                     max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                                     nthread = 1, verbose = 0, save_period = NULL, ...) 
{
  X <- X[,func_mrna.nmf(X)]
  names(X)  <- paste0('x', 1:ncol(X))
  require("xgboost")
  if (packageVersion("xgboost") < 0.6) 
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    if (packageVersion("xgboost") >= "1.1.1.1") {
      objective <- "reg:squarederror"
    }
    else {
      objective <- "reg:linear"
    }
    model = xgboost::xgboost(data = xgmat, objective = objective, 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period, eval_metric = "logloss")
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

SL.earth.mrna.pca <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
          nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
          nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) 
{
  X <- X[,func_mrna.pca(X)]
  require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan, 
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}

SL.earth.mrna.encoder <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                               nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                               nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) 
{
  X <- X[,func_mrna.encoder(X)]
  require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan, 
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}

SL.earth.mrna.nmf <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                                   nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                                   nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) 
{
  X <- X[,func_mrna.nmf(X)]
  require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan, 
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}

SL.glm.mrna.pca <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  X <- X[,func_mrna.pca(X)]
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  fit.glm <- glm(Y ~ ., data = X, family = family, weights = obsWeights, 
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glm.mrna.encoder <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  X <- X[,func_mrna.encoder(X)]
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  fit.glm <- glm(Y ~ ., data = X, family = family, weights = obsWeights, 
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glm.mrna.nmf <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  X <- X[,func_mrna.nmf(X)]
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  fit.glm <- glm(Y ~ ., data = X, family = family, weights = obsWeights, 
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

### rf.cnv + xgboost.cnv + earth.cnv
SL.rf.cnv <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
                       write.forest = TRUE, probability = family$family == "binomial", 
                       min.node.size = ifelse(family$family == "gaussian", 5, 1), 
                       replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = T, ...) 
{
  X <- X[,func_cnv(X)]
  require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, 
                                               X), num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = probability, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.xgboost.cnv <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
                            max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                            nthread = 1, verbose = 0, save_period = NULL, ...) 
{
  X <- X[,func_cnv(X)]
  names(X)  <- paste0('x', 1:ncol(X))
  require("xgboost")
  if (packageVersion("xgboost") < 0.6) 
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    if (packageVersion("xgboost") >= "1.1.1.1") {
      objective <- "reg:squarederror"
    }
    else {
      objective <- "reg:linear"
    }
    model = xgboost::xgboost(data = xgmat, objective = objective, 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period, eval_metric = "logloss")
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

SL.earth.cnv <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                          nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                          nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) 
{
  X <- X[,func_cnv(X)]
  require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan, 
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}

SL.glm.cnv <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  X <- X[,func_cnv(X)]
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  fit.glm <- glm(Y ~ ., data = X, family = family, weights = obsWeights, 
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}


SL.rf.cnv.pca <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
                           write.forest = TRUE, probability = family$family == "binomial", 
                           min.node.size = ifelse(family$family == "gaussian", 5, 1), 
                           replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = T, ...) 
{
  X <- X[,func_cnv.pca(X)]
  require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, 
                                               X), num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = probability, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.rf.cnv.encoder <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
                               write.forest = TRUE, probability = family$family == "binomial", 
                               min.node.size = ifelse(family$family == "gaussian", 5, 1), 
                               replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = T, ...) 
{
  X <- X[,func_cnv.encoder(X)]
  require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, 
                                               X), num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = probability, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.rf.cnv.nmf <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
                           write.forest = TRUE, probability = family$family == "binomial", 
                           min.node.size = ifelse(family$family == "gaussian", 5, 1), 
                           replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = T, ...) 
{
  X <- X[,func_cnv.nmf(X)]
  require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, 
                                               X), num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = probability, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.xgboost.cnv.pca <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
                                max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                                nthread = 1, verbose = 0, save_period = NULL, ...) 
{
  X <- X[,func_cnv.pca(X)]
  names(X)  <- paste0('x', 1:ncol(X))
  require("xgboost")
  if (packageVersion("xgboost") < 0.6) 
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    if (packageVersion("xgboost") >= "1.1.1.1") {
      objective <- "reg:squarederror"
    }
    else {
      objective <- "reg:linear"
    }
    model = xgboost::xgboost(data = xgmat, objective = objective, 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period, eval_metric = "logloss")
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

SL.xgboost.cnv.encoder <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
                                    max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                                    nthread = 1, verbose = 0, save_period = NULL, ...) 
{
  X <- X[,func_cnv.encoder(X)]
  names(X)  <- paste0('x', 1:ncol(X))
  require("xgboost")
  if (packageVersion("xgboost") < 0.6) 
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    if (packageVersion("xgboost") >= "1.1.1.1") {
      objective <- "reg:squarederror"
    }
    else {
      objective <- "reg:linear"
    }
    model = xgboost::xgboost(data = xgmat, objective = objective, 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period, eval_metric = "logloss")
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

SL.xgboost.cnv.nmf <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
                                max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                                nthread = 1, verbose = 0, save_period = NULL, ...) 
{
  X <- X[,func_cnv.nmf(X)]
  names(X)  <- paste0('x', 1:ncol(X))
  require("xgboost")
  if (packageVersion("xgboost") < 0.6) 
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    if (packageVersion("xgboost") >= "1.1.1.1") {
      objective <- "reg:squarederror"
    }
    else {
      objective <- "reg:linear"
    }
    model = xgboost::xgboost(data = xgmat, objective = objective, 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, nthread = nthread, 
                             params = params, save_period = save_period, eval_metric = "logloss")
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

SL.earth.cnv.pca <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                              nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                              nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) 
{
  X <- X[,func_cnv.pca(X)]
  require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan, 
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}

SL.earth.cnv.encoder <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                                  nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                                  nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) 
{
  X <- X[,func_cnv.encoder(X)]
  require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan, 
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}

SL.earth.cnv.nmf <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                              nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                              nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) 
{
  X <- X[,func_cnv.nmf(X)]
  require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan, 
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}

SL.glm.cnv.pca <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  X <- X[,func_cnv.pca(X)]
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  fit.glm <- glm(Y ~ ., data = X, family = family, weights = obsWeights, 
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glm.cnv.encoder <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  X <- X[,func_cnv.encoder(X)]
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  fit.glm <- glm(Y ~ ., data = X, family = family, weights = obsWeights, 
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glm.cnv.nmf <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  X <- X[,func_cnv.nmf(X)]
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  fit.glm <- glm(Y ~ ., data = X, family = family, weights = obsWeights, 
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

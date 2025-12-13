# SuperLearner wrappers with sensible defaults for different feature subsets.
# Each function constrains the design matrix before fitting to keep pipelines clear.

#' Multivariate adaptive regression splines with RF-based feature selection
SL.earth.FeatureSelection <- function(Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3,
                                      nk = max(21, 2 * ncol(X) + 1), pmethod = "backward",
                                      nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) {
  selected_mask <- SuperLearner::screen.randomForest(Y = Y, X = X, family = family, nVar = 300)
  selected <- X[, selected_mask, drop = FALSE]
  if (!requireNamespace("earth", quietly = TRUE)) stop("earth package is required")

  fit <- if (family$family == "gaussian") {
    earth::earth(x = selected, y = Y, degree = degree, nk = nk, penalty = penalty,
                 pmethod = pmethod, nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan)
  } else {
    earth::earth(x = selected, y = Y, degree = degree, nk = nk, penalty = penalty,
                 pmethod = pmethod, nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan,
                 glm = list(family = binomial))
  }
  list(pred = predict(fit, newdata = newX, type = "response"), fit = structure(list(object = fit), class = "SL.earth"))
}

#' Ranger random forest restricted to clinical variables
SL.rf.clinical <- function(Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))),
                           write.forest = TRUE, probability = family$family == "binomial",
                           min.node.size = ifelse(family$family == "gaussian", 5, 1),
                           replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = FALSE, ...) {
  subset_X <- X[, func_clinical(X), drop = FALSE]
  if (!requireNamespace("ranger", quietly = TRUE)) stop("ranger package is required")
  if (family$family == "binomial") Y <- as.factor(Y)
  if (is.matrix(subset_X)) subset_X <- data.frame(subset_X)

  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, subset_X), num.trees = num.trees, mtry = mtry,
                        min.node.size = min.node.size, replace = replace, sample.fraction = sample.fraction,
                        case.weights = obsWeights, write.forest = write.forest, probability = probability,
                        num.threads = num.threads, verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") pred <- pred[, "1"]
  list(pred = pred, fit = structure(list(object = fit, verbose = verbose), class = "SL.ranger"))
}

#' XGBoost restricted to clinical variables
SL.xgboost.clinical <- function(Y, X, newX, family, obsWeights, id, ntrees = 1000, max_depth = 4,
                                shrinkage = 0.1, minobspernode = 10, params = list(), nthread = 1,
                                verbose = 0, save_period = NULL, ...) {
  subset_X <- X[, func_clinical(X), drop = FALSE]
  names(subset_X) <- paste0("x", seq_len(ncol(subset_X)))
  if (!requireNamespace("xgboost", quietly = TRUE)) stop("xgboost package is required")
  if (!is.matrix(subset_X)) subset_X <- model.matrix(~ . - 1, subset_X)

  xgmat <- xgboost::xgb.DMatrix(data = subset_X, label = Y, weight = obsWeights)
  objective <- switch(family$family, gaussian = if (packageVersion("xgboost") >= "1.1.1.1") "reg:squarederror" else "reg:linear",
                      binomial = "binary:logistic", multinomial = "multi:softmax", stop("Unsupported family"))
  model <- xgboost::xgboost(data = xgmat, objective = objective, nrounds = ntrees, max_depth = max_depth,
                            min_child_weight = minobspernode, eta = shrinkage, verbose = verbose, nthread = nthread,
                            params = params, save_period = save_period, eval_metric = if (family$family == "binomial") "logloss" else NULL,
                            num_class = if (family$family == "multinomial") length(unique(Y)) else NULL)

  if (!is.matrix(newX)) newX <- model.matrix(~ . - 1, newX)
  pred <- predict(model, newdata = newX)
  list(pred = pred, fit = structure(list(object = model), class = "SL.xgboost"))
}

#' Earth model on clinical variables
SL.earth.clinical <- function(Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3,
                              nk = max(21, 2 * ncol(X) + 1), pmethod = "backward",
                              nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) {
  subset_X <- X[, func_clinical(X), drop = FALSE]
  if (!requireNamespace("earth", quietly = TRUE)) stop("earth package is required")
  fit <- if (family$family == "gaussian") {
    earth::earth(x = subset_X, y = Y, degree = degree, nk = nk, penalty = penalty,
                 pmethod = pmethod, nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan)
  } else {
    earth::earth(x = subset_X, y = Y, degree = degree, nk = nk, penalty = penalty,
                 pmethod = pmethod, nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan,
                 glm = list(family = binomial))
  }
  list(pred = predict(fit, newdata = newX, type = "response"), fit = structure(list(object = fit), class = "SL.earth"))
}

#' GLM restricted to clinical variables
SL.glm.clinical <- function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
  subset_X <- X[, func_clinical(X), drop = FALSE]
  if (is.matrix(subset_X)) subset_X <- as.data.frame(subset_X)
  fit <- stats::glm(Y ~ ., data = subset_X, family = family, weights = obsWeights, model = model)
  if (is.matrix(newX)) newX <- as.data.frame(newX)
  list(pred = predict(fit, newdata = newX, type = "response"), fit = structure(list(object = fit), class = "SL.glm"))
}

#' Ranger random forest restricted to mRNA variables
SL.rf.mrna <- function(Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))),
                       write.forest = TRUE, probability = family$family == "binomial",
                       min.node.size = ifelse(family$family == "gaussian", 5, 1),
                       replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = FALSE, ...) {
  subset_X <- X[, func_mrna(X), drop = FALSE]
  if (!requireNamespace("ranger", quietly = TRUE)) stop("ranger package is required")
  if (family$family == "binomial") Y <- as.factor(Y)
  if (is.matrix(subset_X)) subset_X <- data.frame(subset_X)

  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, subset_X), num.trees = num.trees, mtry = mtry,
                        min.node.size = min.node.size, replace = replace, sample.fraction = sample.fraction,
                        case.weights = obsWeights, write.forest = write.forest, probability = probability,
                        num.threads = num.threads, verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") pred <- pred[, "1"]
  list(pred = pred, fit = structure(list(object = fit, verbose = verbose), class = "SL.ranger"))
}

#' XGBoost restricted to mRNA variables
SL.xgboost.mrna <- function(Y, X, newX, family, obsWeights, id, ntrees = 1000, max_depth = 4,
                            shrinkage = 0.1, minobspernode = 10, params = list(), nthread = 1,
                            verbose = 0, save_period = NULL, ...) {
  subset_X <- X[, func_mrna(X), drop = FALSE]
  names(subset_X) <- paste0("x", seq_len(ncol(subset_X)))
  if (!requireNamespace("xgboost", quietly = TRUE)) stop("xgboost package is required")
  if (!is.matrix(subset_X)) subset_X <- model.matrix(~ . - 1, subset_X)

  xgmat <- xgboost::xgb.DMatrix(data = subset_X, label = Y, weight = obsWeights)
  objective <- switch(family$family, gaussian = if (packageVersion("xgboost") >= "1.1.1.1") "reg:squarederror" else "reg:linear",
                      binomial = "binary:logistic", multinomial = "multi:softmax", stop("Unsupported family"))
  model <- xgboost::xgboost(data = xgmat, objective = objective, nrounds = ntrees, max_depth = max_depth,
                            min_child_weight = minobspernode, eta = shrinkage, verbose = verbose, nthread = nthread,
                            params = params, save_period = save_period, eval_metric = if (family$family == "binomial") "logloss" else NULL,
                            num_class = if (family$family == "multinomial") length(unique(Y)) else NULL)

  if (!is.matrix(newX)) newX <- model.matrix(~ . - 1, newX)
  pred <- predict(model, newdata = newX)
  list(pred = pred, fit = structure(list(object = model), class = "SL.xgboost"))
}

#' Earth model on mRNA variables
SL.earth.mrna <- function(Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3,
                          nk = max(21, 2 * ncol(X) + 1), pmethod = "backward",
                          nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...) {
  subset_X <- X[, func_mrna(X), drop = FALSE]
  if (!requireNamespace("earth", quietly = TRUE)) stop("earth package is required")
  fit <- if (family$family == "gaussian") {
    earth::earth(x = subset_X, y = Y, degree = degree, nk = nk, penalty = penalty,
                 pmethod = pmethod, nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan)
  } else {
    earth::earth(x = subset_X, y = Y, degree = degree, nk = nk, penalty = penalty,
                 pmethod = pmethod, nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan,
                 glm = list(family = binomial))
  }
  list(pred = predict(fit, newdata = newX, type = "response"), fit = structure(list(object = fit), class = "SL.earth"))
}

#' GLM restricted to mRNA variables
SL.glm.mrna <- function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
  subset_X <- X[, func_mrna(X), drop = FALSE]
  if (is.matrix(subset_X)) subset_X <- as.data.frame(subset_X)
  fit <- stats::glm(Y ~ ., data = subset_X, family = family, weights = obsWeights, model = model)
  if (is.matrix(newX)) newX <- as.data.frame(newX)
  list(pred = predict(fit, newdata = newX, type = "response"), fit = structure(list(object = fit), class = "SL.glm"))
}

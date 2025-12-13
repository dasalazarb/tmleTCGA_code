# High-level wrappers to train common SuperLearner configurations without
# duplicating boilerplate code. Each function returns a list with
# predictions and the fitted mcSuperLearner object.

train_superlearner <- function(Y, X, family, learner_defs, screener = NULL, msg = NULL, ...) {
  if (!is.null(msg)) message(msg)
  lib <- if (is.null(screener)) {
    learner_defs$names
  } else {
    lapply(learner_defs$names, function(name) c(name, screener))
  }
  model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE),
                          SL.library = lib, method = "method.NNloglik", ...)
  structure(list(pred = model$SL.predict, fit = structure(list(object = model), class = "SL.template")),
            class = "SL.wrapper")
}

# Clinical learners -----------------------------------------------------
SL.clinical.rf <- function(Y, X, family, ...) {
  tune <- list(num.trees = c(1000, 5000), min.node.size = c(5, 25, 100),
               mtry = floor(sqrt(sum(func_clinical(X))) * c(0.5, 0.75, 1)))
  learner <- SuperLearner::create.Learner("SL.ranger", tune = tune)
  train_superlearner(Y, X, family, learner, screener = "func_clinical",
                    msg = "Training SL.clinical.rf", ...)
}

SL.clinical.xgb <- function(Y, X, family, ...) {
  tune <- list(ntrees = 2000, max_depth = 50, shrinkage = 0.001)
  learner <- SuperLearner::create.Learner("SL.xgboost", params = tune, detailed_names = TRUE, name_prefix = "xgb")
  train_superlearner(Y, X, family, learner, screener = "func_clinical",
                    msg = "Training SL.clinical.xgb", ...)
}

SL.clinical.mars <- function(Y, X, family, ...) {
  tune <- list(degree = 1:2, pmethod = "backward", nprune = floor(seq(10, 50, length.out = 5)),
               Use.beta.cache = FALSE, newvar.penalty = 0.02)
  learner <- SuperLearner::create.Learner("SL.earth", tune = tune, name_prefix = "clinical.mars")
  train_superlearner(Y, X, family, learner, screener = "func_clinical",
                    msg = "Training SL.clinical.mars", ...)
}

SL.clinical.kknn <- function(Y, X, family, ...) {
  tune <- list(k = c(3, 10), method = c("euclidian", "manhattan", "mahalanobis"),
               weights_function = "gaussian", transf_categ_cols = TRUE)
  learner <- SuperLearner::create.Learner("SL.kernelKnn", tune = tune, name_prefix = "KKnn")
  train_superlearner(Y, X, family, learner, screener = "func_clinical",
                    msg = "Training SL.clinical.kknn", ...)
}

SL.clinical.bayesglm <- function(Y, X, family, ...) {
  learner <- SuperLearner::create.Learner("SL.bayesglm", name_prefix = "bayesglm")
  train_superlearner(Y, X, family, learner, screener = "func_clinical",
                    msg = "Training SL.clinical.bayesglm", ...)
}

# Early-stage learners using global screening ---------------------------------
SL.early.mars <- function(Y, X, family, ...) {
  tune <- list(degree = 1, pmethod = "backward", nprune = floor(seq(10, 50, length.out = 5)),
               Use.beta.cache = FALSE, newvar.penalty = 0.02)
  learner <- SuperLearner::create.Learner("SL.earth", tune = tune, name_prefix = "early.mars")
  train_superlearner(Y, X, family, learner, screener = "new.screen.randomForest",
                    msg = "Training SL.early.mars", ...)
}

SL.early.SVM <- function(Y, X, family, ...) {
  tune <- list(kernel = c("rbfdot", "laplacedot", "polydot", "tanhdot", "anovadot"),
               kpar = list(list(sigma = 2), list(sigma = 2), list(degree = 2), list(scale = 2, offset = 1),
                           list(sigma = 2, degree = 2)))
  learner <- SuperLearner::create.Learner("SL.ksvm", params = tune, name_prefix = "early.ksvm")
  train_superlearner(Y, X, family, learner, screener = "new.screen.randomForest",
                    msg = "Training SL.early.SVM", ...)
}

SL.early.xgb <- function(Y, X, family, ...) {
  tune <- list(ntrees = 2000, max_depth = 50, shrinkage = 0.001)
  learner <- SuperLearner::create.Learner("SL.xgboost", params = tune, detailed_names = TRUE, name_prefix = "xgb")
  train_superlearner(Y, X, family, learner, msg = "Training SL.early.xgb", ...)
}

SL.early.rf <- function(Y, X, family, ...) {
  tune <- list(num.trees = c(1000, 5000), min.node.size = c(5, 25, 100), mtry = floor(sqrt(ncol(X)) * c(0.5, 0.75, 1)))
  learner <- SuperLearner::create.Learner("SL.ranger", tune = tune)
  train_superlearner(Y, X, family, learner, msg = "Training SL.early.rf", ...)
}

SL.early.KKnn <- function(Y, X, family, ...) {
  tune <- list(k = c(3, 10), method = c("euclidian", "manhattan", "mahalanobis"),
               weights_function = "gaussian", transf_categ_cols = TRUE)
  learner <- SuperLearner::create.Learner("SL.kernelKnn", tune = tune, name_prefix = "KKnn")
  train_superlearner(Y, X, family, learner, msg = "Training SL.early.KKnn", ...)
}

SL.early.bayesglm <- function(Y, X, family, ...) {
  learner <- SuperLearner::create.Learner("SL.bayesglm", name_prefix = "bayesglm")
  train_superlearner(Y, X, family, learner, msg = "Training SL.early.bayesglm", ...)
}

# Late-stage learners using modality-specific screening -----------------------
SL.Late.rf <- function(Y, X, family, ...) {
  tune <- list(num.trees = c(1000, 5000), min.node.size = c(5, 25, 100), mtry = floor(sqrt(ncol(X)) * c(0.5, 0.75, 1)))
  learner <- SuperLearner::create.Learner("SL.ranger", tune = tune)
  screens <- c("new.screen.rf.mrna", "new.screen.rf.cnv")
  lib <- unlist(lapply(learner$names, function(name) lapply(screens, function(scr) c(name, scr))), recursive = FALSE)
  model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE),
                          SL.library = lib, method = "method.NNloglik", ...)
  structure(list(pred = model$SL.predict, fit = structure(list(object = model), class = "SL.template")),
            class = "SL.wrapper")
}

SL.Late.SVM <- function(Y, X, family, ...) {
  tune <- list(kernel = c("rbfdot", "laplacedot", "polydot", "tanhdot", "anovadot"),
               kpar = list(list(sigma = 2), list(sigma = 2), list(degree = 2), list(scale = 2, offset = 1),
                           list(sigma = 2, degree = 2)))
  learner <- SuperLearner::create.Learner("SL.ksvm", params = tune, name_prefix = "late.ksvm")
  screens <- c("new.screen.rf.mrna", "new.screen.rf.cnv")
  lib <- unlist(lapply(learner$names, function(name) lapply(screens, function(scr) c(name, scr))), recursive = FALSE)
  model <- mcSuperLearner(Y = Y, X = X, family = family$family, control = list(saveFitLibrary = TRUE),
                          SL.library = lib, method = "method.NNloglik", ...)
  structure(list(pred = model$SL.predict, fit = structure(list(object = model), class = "SL.template")),
            class = "SL.wrapper")
}

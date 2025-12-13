# Convenience wrappers to create tuned SuperLearner models.

#' Define a ranger random forest grid using mtry scaled by feature count
#' @param W_ Feature matrix used only to size the mtry grid.
#' @return List returned by SuperLearner::create.Learner
def_rf <- function(W_) {
  mtry_seq <- floor(sqrt(ncol(W_)) * c(0.5, 1, 2))
  SuperLearner::create.Learner("SL.ranger", params = list(num.trees = 1000, mtry = mtry_seq))
}

#' Define kernel SVM learners across multiple kernels
#' @return List of learner definitions for SL.ksvm
def_ksvm <- function() {
  kernel <- c("rbfdot", "polydot", "vanilladot", "tanhdot", "laplacedot", "besseldot", "anovadot", "splinedot")
  kpar <- list(sigma = 0.1, degree = 2)
  SuperLearner::create.Learner("SL.ksvm", params = list(kernel = kernel, kpar = kpar))
}

#' Lightweight xgboost defaults for quick experiments
#' @return List of learner definitions for SL.xgboost
def0_xgboost <- function() {
  params <- list(ntrees = c(10, 20), max_depth = 1:2, shrinkage = c(0.001, 0.01))
  SuperLearner::create.Learner("SL.xgboost", params = params)
}


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
library(caret)
library(SuperLearner)

def_rf <- function(W_) {
  (mtry_seq = floor(sqrt(ncol(W_)) * c(0.5, 1, 2)))
  learners <- create.Learner("SL.ranger", params = list(num.trees = 1000, mtry = mtry_seq))
  return(learners)
}

def_ksvm <- function() {
  (kernel <- c("rbfdot ", "polydot", "vanilladot", "tanhdot", "laplacedot", "besseldot", 
               "anovadot", "splinedot"))
  (kpar <- list(sigma=0.1, degree=2))
  learners <- create.Learner("SL.ksvm", params = list(kernel = kernel, kpar = kpar))
  return(learners)
}

def0_xgboost <- function() {
  (params = list(ntrees = c(10, 20),
        max_depth = 1:2,
        shrinkage = c(0.001, 0.01)))
  learners <- create.Learner("SL.xgboost", params = params)
  return(learners)
}
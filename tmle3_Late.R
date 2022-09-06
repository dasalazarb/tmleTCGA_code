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
# http://benkeser.github.io/sllecture/
# https://rdrr.io/cran/tmle/src/R/tmle.R#sym-.verifyArgs
# https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html#xgboost-hyperparameter-exploration
source("G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/give_O_GBM_and_LGG.R")
source("G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/func_lateIntegration.R")

library(tmle3)
library(sl3)
library(kableExtra)
library(knitr)
# library(hal9001)
library(tidyverse)
library(data.table)
library(SuperLearner)
library(origami)

# ----------------------------- #
# ... load and process data ... #
# ----------------------------- #
#' 1. Set central path
if (identical(character(0), dir("D:/pathFeatureLabelKernels"))) {
  central_path <- "E:"
} else {
  central_path <- "D:"
}

pathFeatureLabel <- paste0(central_path, "/pathFeatureLabelKernels")

index <- 1:5
i <- 1
for (i in index)  {#' 2. Load O data: (Y, A, delta, W, clinical)
  O <- give_O_GBM_LGG(pathFeatureLabel,dimReduc=FALSE,outTMLE3=TRUE, index=i, Y_as_NaN=FALSE)
  Y <- O[[1]]; A <- O[[2]]; Delta <- O[[3]]; W <- O[[4]]; clinical <- O[[5]];
  
  #' #' (optional) Sampling W for quick results
  #' W.noclinic <- W %>%
  #'   dplyr::select(., -starts_with("clinical_"))
  #' 
  #' W.cnv <- W.noclinic %>%
  #'   dplyr::select(., starts_with("cnv_")); dim(W.cnv)
  #' W.mrna <- W.noclinic %>%
  #'   dplyr::select(., starts_with("mrna_")); dim(W.mrna)
  #' 
  #' W <- clinical %>%
  #'   dplyr::bind_cols(W.cnv[,sample(1:ncol(W.cnv), 200)], W.mrna[,sample(1:ncol(W.mrna), 200)])
  #' W[1:5,1:20]; dim(W)
  
  #' 3. Generate O_data
  O_data <- data.table(cbind(as.data.frame(Y), as.data.frame(as.factor(Delta)),
                  as.data.frame(as.factor(A)),W)); colnames(O_data)[2:3] <- c("Delta", "A")
  O_data[1:5,1:10]; dim(O_data)
  
  #' 4. Explore O_data: A vs Y_(with missing) and A vs Y
  table(is.na(O_data$A), is.na(O_data$Y))
  table(as.factor(O_data$A), is.na(O_data$Y))
  
  # ------------------------- #
  # ... tmle3 Components  ... #
  # ------------------------- #
  #' 1. Define the variable roles: W -> A -> Y
  node_list <- list(
    W = colnames(O_data)[!colnames(O_data) %in% c("Y","A","Delta")], 
    A = "A", 
    Y = "Y"
  )
  
  #' #' 1.1. Missingness
  processed <- process_missing(O_data, node_list)
  O_data <- processed$data
  node_list <- processed$node_list
  
  #' 2. Create a "Spec" Object
  ate_spec <- tmle_ATE(
    treatment_level = 0,
    control_level = 1
  )
  
  #' 3. Define the Relevant Super Learners
  # sl3_list_properties()
  # sl3_list_learners(c("binomial"))
  # sl3_list_learners(c("continuous"))
  
  #' Choose base learners and metalearner
  screen_mrna <- Lrnr_pkg_SuperLearner_screener$new("screen.mrna")
  screen_cnv <- Lrnr_pkg_SuperLearner_screener$new("screen.cnv")
  
  # ## grid for RF
  # grid_params <- list(
  #   mtry = seq(500, 1000, length.out = 3) %>% floor(),
  #   max_nodes = c(25,500),
  #   ntree = seq(500, 1000, length.out = 3) %>% floor()
  # )
  # 
  # grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  # learners <- apply(grid, MARGIN = 1, function(tuning_params) {
  #   do.call(Lrnr_randomForest$new, as.list(tuning_params))
  # })
  # 
  # stack_rf <- make_learner(Stack, learners)
  
  # lrnr <- make_learner(Lrnr_randomForest, mtry=100, max_nodes=25, ntree=500)
  
  grid_params <- list(
    max_depth = c(10,20),
    eta = c(0.001, 0.1, 0.3),
    nrounds = c(100,200)
  )
  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
  # grid
  
  lrnr <- apply(grid, MARGIN = 1, function(tuning_params) {
    do.call(Lrnr_xgboost$new, as.list(tuning_params))
  })
  
  # lrnr <- make_learner(Lrnr_xgboost, nrounds=100, eta=0.01, gamma=0.1, max_depth=15)
  lrnr <- make_learner(Stack, lrnr)
  nombreLearner <- lrnr$name
  
  ## separate models for each omic profile
  pipe_rf_mrna <- make_learner(Pipeline, screen_mrna, lrnr)
  pipe_rf_cnv <- make_learner(Pipeline, screen_cnv, lrnr)
  pipe_all <- make_learner(Stack, c(pipe_rf_mrna, pipe_rf_cnv))
  
  ## learners for each part of TMLE
  sl_Y <- Lrnr_sl$new(
    learners = pipe_all # list(lrnr)
    # learners = unlist(list(lrnr_lasso), recursive = TRUE)
  )
  sl_Delta <- Lrnr_sl$new(
    learners = pipe_all # list(lrnr)
    # learners = list(lrnr_mean, lrnr_glm, lrnr_lasso, lrnr_ridge)
  )
  sl_A <- Lrnr_sl$new(
    learners = pipe_all # list(lrnr)
    # learners = list(lrnr_mean, lrnr_glm)
  )
  
  learner_list <- list(A = sl_A, delta_Y = sl_Delta, Y = sl_Y)
  
  # ## 4. Fit the TMLE
  # tmle_fit <- tmle3(ate_spec, O_data, node_list, learner_list)
  
  #' 4. Initial Likelihood
  tmle_task <- ate_spec$make_tmle_task(data = O_data, node_list = node_list)
  # print(tmle_task$npsem)
  initial_likelihood <- ate_spec$make_initial_likelihood(tmle_task,learner_list)
  print(initial_likelihood)
  
  initial_likelihood$get_likelihoods(tmle_task)
  
  write.csv(x = as.matrix(initial_likelihood$get_likelihoods(tmle_task)), file = paste0(pathFeatureLabel, "/co-md_records/tmle_initLikelihoods_lgg_gbm_",i,"_late_", nombreLearner, ".csv"))
  
  #' 5. Targeted Likelihood (updater) - CV-TMLE
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
  tmle_params <- ate_spec$make_params(tmle_task, targeted_likelihood)
  print(tmle_params)
  
  #' 6. Putting it all together
  tmle_fit_manual <- fit_tmle3(
    tmle_task, targeted_likelihood, tmle_params,
    targeted_likelihood$updater
  )
  print(tmle_fit_manual)
  cat("\n")
  
  cat(paste0("Results for imputation No. ", i, "\n"), file = paste0(pathFeatureLabel, "/co-md_records/tmle_finalEstimation_lgg_gbm__late_", nombreLearner, ".csv"), append = TRUE)
  capture.output(tmle_fit_manual, file = paste0(pathFeatureLabel, "/co-md_records/tmle_finalEstimation_lgg_gbm__late_", nombreLearner, ".csv"), append = TRUE)
  cat(paste0("Results for imputation No. ", i, "\n"), file = paste0(pathFeatureLabel, "/co-md_records/tmle_metalearner_lgg_gbm__late_", nombreLearner, ".csv"), append = TRUE)
  capture.output(tmle_fit_manual$likelihood$factor_list[["A"]]$learner, file = paste0(pathFeatureLabel, "/co-md_records/tmle_metalearner_lgg_gbm__late_", nombreLearner, ".csv"), append = TRUE)
  
  }

## imputed_1
# [1] "Cross-validated risk (MSE, squared error loss):"
# learner coefficients      risk         se    fold_sd fold_min_risk fold_max_risk
# 1: Lrnr_randomForest_500_TRUE_5_100_25            1 0.2200385 0.01644308 0.03440609     0.1655623     0.2837483

# 
# 
# ## ------------- ##
# ps <- read.csv(paste0(pathFeatureLabel,"/co-md_records/initial_likelihoods_superlearners_lgg_gbm_2.csv"),row.names = 1)
# ps <- cbind(ps,treat=A,missingness=Delta)
# head(ps)
# 
# # Libraries
# library(ggplot2)
# 
# # Chart
# ggplot() +
#   # Top
#   geom_density(data=subset(ps, treat==0), aes(x = A, y = ..density..), fill="#69b3a2" ) +
#   # geom_label( aes(x=0.4, y=250, label="RAD+TMZ"), color="#69b3a2") +
#   # Bottom
#   geom_density(data=subset(ps, treat==1), aes(x = A, y = -..density..), fill= "#404080") #+ scale_y_continuous(breaks=round(seq(min(-2.5), max(3.5), by=0.2),2))
#   # geom_label( aes(x=0.5, y=-100, label="RAD"), color="#404080")
# 

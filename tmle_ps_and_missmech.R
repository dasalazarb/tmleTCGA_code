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
library(caret)
library(dplyr)
library(purrr)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(SuperLearner)

if (file.exists(paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_sl_ps.rds")) & file.exists(paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_sl_missingness_mechanism.rds")) ) {
  print("Reading files from tmle_ps_miss_mech_coxModel folder")
  load(paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_g1W.rda"))
  load(paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_pDelta1.rda"))
  
  ## ----------------------- ##
  ### ... Plot P(A=1|W) ... ###
  ## ----------------------- ##
  ps <- data.frame(cbind(g1W,treat=A))
  head(ps,20)
  
  # PS plot
  ps_plot <- ggplot() +
    # Top
    geom_density(data=subset(ps, treat==0), aes(x = g1W, y = ..density..), fill="#69b3a2" ) +
    geom_label( aes(x=0.5, y=2, label="RAD"), color="#69b3a2") +
    # Bottom
    geom_density(data=subset(ps, treat==1), aes(x = g1W, y = -..density..), fill= "#404080") +
    geom_label( aes(x=0.5, y=-2, label="RAD+TMZ"), color="#404080") +
    labs(x = "Propensity score estimated", y = "Density")
  
} else {
  ## ----------------------------------- ##
  ### ... 1.1. Cox regression model ... ###
  ## ----------------------------------- ##
  # install.packages(c("survival", "survminer"))
  clinical <- W %>%
    dplyr::select(.,starts_with("clinical")) %>%
    # dplyr::select(-clinical_primary_diagnosis) %>%
    dplyr::mutate_if(is.character, as.factor)

  if (cancer == "lgg") {
    library("survival")
    library("survminer")

    attach(clinical)

    res.cox <- coxph(Surv(time = O$Y, event = Y, type = "right") ~ factor(A) + factor(clinical_primary_diagnosis) + clinical_age_at_index + factor(clinical_gender) +
                       factor(clinical_IDH_codel_subtype) + factor(clinical_MGMT_promoter_status) + factor(clinical_ATRX_status) +
                       factor(clinical_Chr7_gain_and_Chr10_loss) + clinical_Karnofsky_Performance_Score +
                       clinical_Mutation_Count)
    summary(res.cox)
    saveRDS(res.cox, file = paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_cox_model.rds"))
  }


  ## --------------------------------------------------------------------- ##
  ### ... 2. Conditional treatment assignment probabilities, P(A=1|W) ... ###
  ## --------------------------------------------------------------------- ##
  tune = list(num.trees = c(1000, 5000),
              min.node.size = c(5, 25, 100),
              mtry = floor(sqrt(ncol(clinical)) * c(0.5,0.75,1)))

  lrn_rf <- create.Learner("SL.ranger", tune = tune, detailed_names = TRUE)

  tune = list(degree = 1:3, pmethod = "backward",
              nprune = seq(20, 150, length.out = 5) %>% floor(),
              Use.beta.cache=FALSE) # Note: earth's beta cache would require 3.63 GB, so forcing Use.beta.cache=FALSE.

  lrn_mars <- create.Learner("SL.earth", tune = tune, name_prefix = "SL.mars", detailed_names = TRUE)

  lrn_log <- create.Learner("SL.glm", detailed_names = TRUE)

  sl_ps <- CV.SuperLearner(Y = A, X = clinical, family = binomial(),
                           # For a real analysis we would use V = 10.
                           cvControl = list(V = 10, stratifyCV = TRUE), innerCvControl = list(list(V=10, stratifyCV = TRUE)),
                           parallel = "multicore", verbose = TRUE,
                           SL.library = c(lrn_log$names, lrn_rf$names, lrn_mars$names) )

  saveRDS(sl_ps, file = paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_sl_ps.rds"))

  # sl_log <- mcSuperLearner(Y = A, X = clinical, family = binomial(), cvControl = list(V = 10),
  #                          SL.library = c(lrn_log$names))

  plot(sl_ps) + theme_bw()

  summary(sl_ps)
  g1W <- sl_ps$SL.predict
  save(x = g1W, file = paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_g1W.rda"))

  ## ----------------------- ##
  ### ... Plot P(A=1|W) ... ###
  ## ----------------------- ##
  ps <- data.frame(cbind(g1W,treat=A))
  head(ps,20)

  # PS plot
  ps_plot <- ggplot() +
    # Top
    geom_density(data=subset(ps, treat==0), aes(x = g1W, y = ..density..), fill="#69b3a2" ) +
    geom_label( aes(x=0.5, y=2.5, label="RAD"), color="#69b3a2") +
    # Bottom
    geom_density(data=subset(ps, treat==1), aes(x = g1W, y = -..density..), fill= "#404080") +
    geom_label( aes(x=0.5, y=-2.5, label="RAD+TMZ"), color="#404080") +
    labs(x = "Propensity score estimated", y = "Density")

  ## ---------------------------------------------------------------------------------------------------- ##
  ### ... 3. Conditional probabilities for missingness mechanism, P(Delta=1|A=0,W), P(Delta=1|A=1,W) ... ###
  ## ---------------------------------------------------------------------------------------------------- ##
  lrn_log <- create.Learner("SL.glm")
  sl_rf <- mcSuperLearner(Y = Delta, X = data.frame(A=as.factor(A), clinical), family = binomial(), cvControl = list(V = 10),
                          SL.library = c(lrn_log$names))


  summary(sl_rf)

  d1W <- predict(sl_rf,newdata=data.frame(clinical, A=factor(1)), type="response")$pred
  d0W <- predict(sl_rf,newdata=data.frame(clinical, A=factor(0)), type = "response")$pred

  pDelta1 <- data.frame(d0W,d1W)
  head(pDelta1); dim(pDelta1)

  saveRDS(sl_rf, file = paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_sl_missingness_mechanism.rds"))
  save(x = pDelta1, file = paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_pDelta1.rda"))
}

# if (file.exists(paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_new_screen_randomForest.rds")) ) {
#   new.screen.randomForest <- function(...) {
#     load(file = paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_screen_randomForest.rda"))
#     return(rf.screen.vars)
#   }
# } else {
#   
#   rf.screen.vars <- screen.randomForest(Y, W, nVar = 300, family = binomial())
#   
#   save(x = rf.screen.vars, file = paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_screen_randomForest.rda"))
#   
#   new.screen.randomForest <- function(...) {
#     load(file = paste0(ruta, "tmle_ps_miss_mech_coxModel/tmle_screen_randomForest.rda"))
#     return(rf.screen.vars)
#   }
# }


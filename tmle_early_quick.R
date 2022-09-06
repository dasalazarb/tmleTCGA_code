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
## ------------------------------------ ##
### ... a. Libraries and functions ... ###
## ----------------------------------- ###
library(tmle)
library(caret)
library(dplyr)
library(purrr)
library(ggplot2)
library(RhpcBLASctl)
library(SuperLearner)

## ----------------------------------- ##
### ... b. parallelizing options ... ###
## ---------------------------------- ###
num_cores = RhpcBLASctl::get_num_cores()
options(mc.cores = num_cores)
getOption("mc.cores")

## ---------------------- ##
### ... c. Load data ... ###
## ---------------------- ##
O <- read.csv(file = "C:/Users/da.salazarb/Documents/tmle/tmleTCGA/O_Y_A_Delta_W.csv") ## <- Cambiar ruta del archivo!
O[1:5,1:5]; dim(O)

Y <- ifelse(O$Y < 365, 1, 0); Delta <- ifelse(O$delta == 0 & O$Y < 365, 0, 1); A <- O$A; W <- O[,!(colnames(O) %in% c("Y", "X", "delta", "...3078"))];

W <- W %>% 
  dplyr::mutate(clinical_primary_diagnosis=as.factor(clinical_primary_diagnosis), clinical_gender = as.factor(clinical_gender), 
                clinical_IDH_codel_subtype=as.factor(clinical_IDH_codel_subtype), clinical_MGMT_promoter_status=as.factor(clinical_MGMT_promoter_status), 
                clinical_ATRX_status=as.factor(clinical_ATRX_status), clinical_Chr7_gain_and_Chr10_loss=as.factor(clinical_Chr7_gain_and_Chr10_loss), 
                A=as.factor(A))

# ### choose a sample
# clinical <- W %>%
#   dplyr::select(.,starts_with("clinical"))
# 
# W.noclinic <- W %>%
#   dplyr::select(., -starts_with("clinical_"))
# 
# W.cnv <- W.noclinic %>%
#   dplyr::select(., starts_with("cnv_"))
# W.cnv <- W.cnv[,sample(ncol(W.cnv), 50)]; dim(W.cnv)
# W.mrna <- W.noclinic %>%
#   dplyr::select(., starts_with("mrna_"))
# W.mrna <- W.mrna[,sample(ncol(W.mrna), 50)]; dim(W.mrna)
# 
# W <- clinical %>%
#   dplyr::bind_cols(data.frame(A=A),W.cnv, W.mrna) %>% 
#   dplyr::mutate(A=as.factor(A))
# W[1:5,1:20];dim(W)

## --------------------------------------------------------------------- ##
### ... d. Conditional treatment assignment probabilities, P(A=1|W) ... ###
## --------------------------------------------------------------------- ##
clinical <- W %>% 
  dplyr::select(.,starts_with("clinical"))

tune = list(num.trees = 5000,
            max_depth = c(25, 200), # seq(10,200,100), 
            mtry_seq = floor(sqrt(ncol(clinical)) * c(0.5, 1, 2)))

learners <- create.Learner("SL.ranger", tune = tune)
sl_rf <- mcSuperLearner(Y = A, X = clinical, family = binomial(), cvControl = list(V = 10),
                        SL.library = c(learners$names))
summary(sl_rf)
g1W <- sl_rf$SL.predict

## ----------------------- ##
### ... Plot P(A=1|W) ... ###
## ----------------------- ##
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

## ----------------------------------------------------------------------------------------------------- ##
### ... e. Conditional probabilities for missingness mechanism, P(Delta=1|A=0,W), P(Delta=1|A=1,W) ... ###
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

## ----------------------------------- ##
### ... f. Estimation of E(Y|A,W) ... ###
## ----------------------------------- ##
W_ <- W
W_[1:5,1:5]; dim(W_)

### screen with glmnet
new.screen.glmnet <- function(...) {
  screen.glmnet(..., minscreen = 20, alpha = 0.1)
}

### caso 1: random forest
tune = list(num.trees = 5000, 
            max_depth = c(25,200), 
            mtry_seq = floor(sqrt(ncol(W_)) * c(0.5, 1, 2)))

lrn_rf <- create.Learner("SL.ranger", tune = tune)

sl_rf_early <- mcSuperLearner(Y = Y, X = W_, family = binomial(), cvControl = list(V = 10),
                              SL.library = lapply(lrn_rf$names, function(x) c(x, "new.screen.glmnet")))

### caso 2: logistic regression 
lrn_log <- create.Learner("SL.glm")

# sl_log_early <- mcSuperLearner(Y = Y, X = W_, family = binomial(), cvControl = list(V = 10),
#                               SL.library = c(lrn_log$names))


### caso 3: xgboost
tune = list(ntrees = 5000,
            max_depth = seq(10,200,70), 
            shrinkage = c(0.001, 0.01))

lrn_xgb <- create.Learner("SL.xgboost", tune = tune, detailed_names = TRUE, name_prefix = "xgb")

# sl_xgb_early <- mcSuperLearner(Y = Y, X = W_, family = binomial(), cvControl = list(V = 10),
#                       SL.library = c(lrn_xgb$names))

## caso 4: mars -> Resultado: Funciona
tune = list(degree = c(2,3),
            nprune = c(2,3,5),
            Use.beta.cache=FALSE) # Note: earth's beta cache would require 3.63 GB, so forcing Use.beta.cache=FALSE.

new.screen.randomForest <- function(...) {
  screen.randomForest(..., nVar = 300)
}

lrn_mars <- create.Learner("SL.earth", tune = tune, name_prefix = "mars")

# sl_mars_early <- mcSuperLearner(Y = Y, X = W_, family = binomial(), cvControl = list(V = 10),
#                          SL.library = lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest")))

## CV.SuperLearner
sl_model_early <- CV.SuperLearner(Y = Y, X = W_, family = binomial(), 
                                  # For a real analysis we would use V = 10.
                                  cvControl = list(V = 5), innerCvControl = list(list(V=5)),
                                  parallel = "multicore", 
                                  SL.library = purrr::reduce(list(lapply(lrn_mars$names, function(x) c(x, "new.screen.randomForest")),
                                                                  lapply(lrn_rf$names, function(x) c(x)),
                                                                  lapply(lrn_xgb$names, function(x) c(x)),
                                                                  list("SL.LinearComb.SVM")), append) )

Q1W <- SuperLearner::predict.SuperLearner(sl_model_early,newX=data.frame(A=1, W_ %>% dplyr::select(-"A")), type="response")$pred
Q0W <- SuperLearner::predict.SuperLearner(sl_model_early,newX=data.frame(A=0, W_ %>% dplyr::select(-"A")), type = "response")$pred

Q <- data.frame(Q0W,Q1W)
head(Q); dim(Q)

## ------------------------ ##
### ... g. tmle update ... ###
## ------------------------ ##
result.Qcgc <- tmle(Y = Y, A = A, W = W_, verbose = TRUE, Delta = Delta,
                    Q = Q, g1W = g1W, pDelta1 = pDelta1, family = "binomial")

result.Qcgc

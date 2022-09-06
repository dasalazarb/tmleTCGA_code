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
## ---------------------- ##
### ... 1. Load data ... ###
## ---------------------- ##
# O <- read.csv(file = "C:/Users/da.salazarb/Documents/tmle/tmleTCGA/imputation_1.csv") ## <- Cambiar ruta del archivo!
# O <- read.csv(file = "C:/Users/da.salazarb/Documents/tmle/tmleTCGA/O_Y_A_Delta_W.csv") ## <- Cambiar ruta del archivo!

nvar <- 50
time.t <- 730
cancer <- "lgg"
outcome <- "os"

if (!identical(character(0), dir("G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/"))) {
  init_path <- "G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/"
} else {
  init_path <- "G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/"
}


if (cancer == "lgg") {
  path_to_data <- paste0(init_path, outcome, "_", cancer, "_data_var", nvar)
  print(path_to_data)
  
  O <- read.csv(file = paste0(path_to_data, "/imputation_1.csv")) ## <- Cambiar ruta del archivo!
  O[1:5,1:5]; dim(O)
  
  Y <- ifelse(O$Y < time.t & O$delta == 1, 1, 0); Delta <- ifelse(O$delta == 0 & O$Y < time.t, 0, 1); 
  A <- O$A; W <- O %>% dplyr::select(starts_with("clinical_"), starts_with("mrna_"), starts_with("cnv_"))
  
  ## Delta: durante el periodo de observacion, si la persona perdio el seguimiento sin observar el evento entonces Delta==0
  ## Delta -> siguio en seguimiento durante todo el periodo de observacion?
  
  W <- W %>% 
    dplyr::mutate(clinical_primary_diagnosis=as.factor(clinical_primary_diagnosis), clinical_gender=as.factor(clinical_gender), 
                  clinical_IDH_codel_subtype=as.factor(clinical_IDH_codel_subtype), clinical_MGMT_promoter_status=as.factor(clinical_MGMT_promoter_status), 
                  clinical_ATRX_status=as.factor(clinical_ATRX_status), clinical_Chr7_gain_and_Chr10_loss=as.factor(clinical_Chr7_gain_and_Chr10_loss))
} else {
  path_to_data <- paste0(init_path, outcome, "_", cancer, "_data_var", nvar)
  print(path_to_data)
  
  O <- read.csv(file = paste0(path_to_data, "/imputation_1.csv")) ## <- Cambiar ruta del archivo!
  O[1:5,1:5]; dim(O)
  
  Y <- ifelse(O$Y < time.t & O$delta == 1, 1, 0); Delta <- ifelse(O$delta == 0 & O$Y < time.t, 0, 1); 
  A <- O$A; W <- O %>% dplyr::select(starts_with("clinical_"), starts_with("mirna_"), starts_with("cnv_"))
  
  W <- W %>% 
    dplyr::mutate(clinical_VITALSTATUS=as.factor(clinical_VITALSTATUS), clinical_TUMORSTAGE=as.factor(clinical_TUMORSTAGE), 
                  clinical_TUMORGRADE=as.factor(clinical_TUMORGRADE), clinical_TUMORRESIDUALDISEASE=as.factor(clinical_TUMORRESIDUALDISEASE), 
                  clinical_PRIMARYTHERAPYOUTCOMESUCCESS=as.factor(clinical_PRIMARYTHERAPYOUTCOMESUCCESS), clinical_PERSONNEOPLASMCANCERSTATUS=as.factor(clinical_PERSONNEOPLASMCANCERSTATUS), 
                  clinical_ProgressionFreeStatus=as.factor(clinical_ProgressionFreeStatus), clinical_PlatinumStatus=as.factor(clinical_PlatinumStatus)) %>% 
    dplyr::select(-clinical_TUMORSTAGE, -clinical_PRIMARYTHERAPYOUTCOMESUCCESS, -clinical_TUMORRESIDUALDISEASE)
  colnames(W) <- gsub("mirna_", "mrna_", colnames(W))
}
W[10:15,10:15]; dim(W)

dimen.method <- "encoder"
cnv.dimen <- read.csv(paste0(path_to_data, "/encoder_W_cnv.csv"))
mrna.dimen <- read.csv(paste0(path_to_data, "/encoder_W_mrna.csv"))
colnames(cnv.dimen) <- paste0("cnv_", dimen.method,"_", colnames(cnv.dimen))
colnames(mrna.dimen) <- paste0("mrna_", dimen.method, "_", colnames(mrna.dimen))
W_dimen.encoder <- data.frame(cnv.dimen, mrna.dimen) %>% 
  dplyr::select(-cnv_encoder_X,-mrna_encoder_X)

dimen.method <- "pca"
cnv.dimen <- read.csv(paste0(path_to_data, "/pca_W_cnv.csv"))
mrna.dimen <- read.csv(paste0(path_to_data, "/pca_W_mrna.csv"))
colnames(cnv.dimen) <- paste0("cnv_", dimen.method,"_", colnames(cnv.dimen))
colnames(mrna.dimen) <- paste0("mrna_", dimen.method, "_", colnames(mrna.dimen))
W_dimen.pca <- data.frame(cnv.dimen, mrna.dimen) %>% 
  dplyr::select(-cnv_pca_X,-mrna_pca_X)

dimen.method <- "nmf"
cnv.dimen <- read.csv(paste0(path_to_data, "/nmf_W_cnv.csv"))
mrna.dimen <- read.csv(paste0(path_to_data, "/nmf_W_mrna.csv"))
colnames(cnv.dimen) <- paste0("cnv_", dimen.method,"_", colnames(cnv.dimen))
colnames(mrna.dimen) <- paste0("mrna_", dimen.method, "_", colnames(mrna.dimen))
W_dimen.nmf <- data.frame(cnv.dimen, mrna.dimen) %>% 
  dplyr::select(-cnv_nmf_X,-mrna_nmf_X)

W.new <- data.frame(W, W_dimen.encoder, W_dimen.pca, W_dimen.nmf)
W.new[10:15,10:15]; dim(W.new)

## ---------------------------------------------- ##
### ... funciones de extraccion de variables ... ###
## ---------------------------------------------- ##
func_mrna <- function(X,...){
  returnCols <- rep(FALSE, ncol(X))
  # returnCols[grep("mrna", names(X))] <- TRUE
  returnCols[grepl("mrna_|clinical_", names(X)) & !grepl("encoder_|pca_|nmf_|cnv_", names(X))] <- TRUE
  return(returnCols)
}

func_cnv <- function(X,...){
  returnCols <- rep(FALSE, ncol(X))
  # returnCols[grep("cnv", names(X))] <- TRUE
  returnCols[grepl("cnv_|clinical_", names(X)) & !grepl("encoder_|pca_|nmf_|mrna_", names(X))] <- TRUE
  return(returnCols)
}

func_clinical <- function(X,...){
  returnCols <- rep(FALSE, ncol(X))
  # returnCols[grep("cnv", names(X))] <- TRUE
  returnCols[grepl("clinical_", names(X))] <- TRUE
  return(returnCols)
}

func_mrna.encoder <- function(X,...){
  returnCols <- rep(FALSE, ncol(X))
  # returnCols[grep("mrna", names(X))] <- TRUE
  returnCols[grepl("mrna_encoder_|clinical_", names(X))] <- TRUE
  return(returnCols)
}

func_cnv.encoder <- function(X,...){
  returnCols <- rep(FALSE, ncol(X))
  # returnCols[grep("cnv", names(X))] <- TRUE
  returnCols[grepl("cnv_encoder_|clinical_", names(X))] <- TRUE
  return(returnCols)
}

func_mrna.pca <- function(X,...){
  returnCols <- rep(FALSE, ncol(X))
  # returnCols[grep("mrna", names(X))] <- TRUE
  returnCols[grepl("mrna_pca_|clinical_", names(X))] <- TRUE
  return(returnCols)
}

func_cnv.pca <- function(X,...){
  returnCols <- rep(FALSE, ncol(X))
  # returnCols[grep("cnv", names(X))] <- TRUE
  returnCols[grepl("cnv_pca_|clinical_", names(X))] <- TRUE
  return(returnCols)
}

func_mrna.nmf <- function(X,...){
  returnCols <- rep(FALSE, ncol(X))
  # returnCols[grep("mrna", names(X))] <- TRUE
  returnCols[grepl("mrna_nmf_|clinical_", names(X))] <- TRUE
  return(returnCols)
}

func_cnv.nmf <- function(X,...){
  returnCols <- rep(FALSE, ncol(X))
  # returnCols[grep("cnv", names(X))] <- TRUE
  returnCols[grepl("cnv_nmf_|clinical_", names(X))] <- TRUE
  return(returnCols)
}

new.screen.randomForest <- function(...) {
  screen.randomForest(..., nVar = 300)
}


new.screen.rf.mrna <- function (Y, X, family, nVar = 300, ntree = 1000, 
                                mtry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), 
                                nodesize = ifelse(family$family == "gaussian", 5, 3), maxnodes = NULL, ...) {
  
  new.X <- X[,func_mrna(X)]
  
  require("randomForest")
  if (family$family == "gaussian") {
    rank.rf.fit <- randomForest::randomForest(Y ~ ., data = new.X, 
                                              ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  if (family$family == "binomial") {
    rank.rf.fit <- randomForest::randomForest(as.factor(Y) ~ 
                                                ., data = new.X, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  whichVariable <- (rank(-rank.rf.fit$importance) <= nVar)
  
  whichVariable <- colnames(X) %in% colnames(new.X[,whichVariable])
  
  return(whichVariable)
}

new.screen.rf.cnv <- function (Y, X, family, nVar = 300, ntree = 1000, 
                               mtry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), 
                               nodesize = ifelse(family$family == "gaussian", 5, 3), maxnodes = NULL, ...) {
  
  new.X <- X[,func_cnv(X)]
  
  require("randomForest")
  if (family$family == "gaussian") {
    rank.rf.fit <- randomForest::randomForest(Y ~ ., data = new.X, 
                                              ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  if (family$family == "binomial") {
    rank.rf.fit <- randomForest::randomForest(as.factor(Y) ~ 
                                                ., data = new.X, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  whichVariable <- (rank(-rank.rf.fit$importance) <= nVar)
  
  whichVariable <- colnames(X) %in% colnames(new.X[,whichVariable])
  
  return(whichVariable)
}

## ---------------------------------------------- ##
### ... funciones de extraccion de variables ... ###
## ---------------------------------------------- ##
new.screen.rf.mrna.encoder <- function (Y, X, family, nVar = 300, ntree = 1000, 
                                mtry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), 
                                nodesize = ifelse(family$family == "gaussian", 5, 3), maxnodes = NULL, ...) {
  
  new.X <- X[,func_mrna.encoder(X)]
  
  require("randomForest")
  if (family$family == "gaussian") {
    rank.rf.fit <- randomForest::randomForest(Y ~ ., data = new.X, 
                                              ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  if (family$family == "binomial") {
    rank.rf.fit <- randomForest::randomForest(as.factor(Y) ~ 
                                                ., data = new.X, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  whichVariable <- (rank(-rank.rf.fit$importance) <= nVar)
  
  whichVariable <- colnames(X) %in% colnames(new.X[,whichVariable])
  
  return(whichVariable)
}

new.screen.rf.cnv.encoder <- function (Y, X, family, nVar = 300, ntree = 1000, 
                               mtry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), 
                               nodesize = ifelse(family$family == "gaussian", 5, 3), maxnodes = NULL, ...) {
  
  new.X <- X[,func_cnv.encoder(X)]
  
  require("randomForest")
  if (family$family == "gaussian") {
    rank.rf.fit <- randomForest::randomForest(Y ~ ., data = new.X, 
                                              ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  if (family$family == "binomial") {
    rank.rf.fit <- randomForest::randomForest(as.factor(Y) ~ 
                                                ., data = new.X, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  whichVariable <- (rank(-rank.rf.fit$importance) <= nVar)
  
  whichVariable <- colnames(X) %in% colnames(new.X[,whichVariable])
  
  return(whichVariable)
}

new.screen.rf.mrna.pca <- function (Y, X, family, nVar = 300, ntree = 1000, 
                                      mtry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), 
                                      nodesize = ifelse(family$family == "gaussian", 5, 3), maxnodes = NULL, ...) {
  
  new.X <- X[,func_mrna.pca(X)]
  
  require("randomForest")
  if (family$family == "gaussian") {
    rank.rf.fit <- randomForest::randomForest(Y ~ ., data = new.X, 
                                              ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  if (family$family == "binomial") {
    rank.rf.fit <- randomForest::randomForest(as.factor(Y) ~ 
                                                ., data = new.X, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  whichVariable <- (rank(-rank.rf.fit$importance) <= nVar)
  
  whichVariable <- colnames(X) %in% colnames(new.X[,whichVariable])
  
  return(whichVariable)
}

new.screen.rf.cnv.pca <- function (Y, X, family, nVar = 300, ntree = 1000, 
                                     mtry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), 
                                     nodesize = ifelse(family$family == "gaussian", 5, 3), maxnodes = NULL, ...) {
  
  new.X <- X[,func_cnv.pca(X)]
  
  require("randomForest")
  if (family$family == "gaussian") {
    rank.rf.fit <- randomForest::randomForest(Y ~ ., data = new.X, 
                                              ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  if (family$family == "binomial") {
    rank.rf.fit <- randomForest::randomForest(as.factor(Y) ~ 
                                                ., data = new.X, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  whichVariable <- (rank(-rank.rf.fit$importance) <= nVar)
  
  whichVariable <- colnames(X) %in% colnames(new.X[,whichVariable])
  
  return(whichVariable)
}

new.screen.rf.mrna.nmf <- function (Y, X, family, nVar = 300, ntree = 1000, 
                                      mtry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), 
                                      nodesize = ifelse(family$family == "gaussian", 5, 3), maxnodes = NULL, ...) {
  
  new.X <- X[,func_mrna.nmf(X)]
  
  require("randomForest")
  if (family$family == "gaussian") {
    rank.rf.fit <- randomForest::randomForest(Y ~ ., data = new.X, 
                                              ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  if (family$family == "binomial") {
    rank.rf.fit <- randomForest::randomForest(as.factor(Y) ~ 
                                                ., data = new.X, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  whichVariable <- (rank(-rank.rf.fit$importance) <= nVar)
  
  whichVariable <- colnames(X) %in% colnames(new.X[,whichVariable])
  
  return(whichVariable)
}

new.screen.rf.cnv.nmf <- function (Y, X, family, nVar = 300, ntree = 1000, 
                                     mtry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), 
                                     nodesize = ifelse(family$family == "gaussian", 5, 3), maxnodes = NULL, ...) {
  
  new.X <- X[,func_cnv.nmf(X)]
  
  require("randomForest")
  if (family$family == "gaussian") {
    rank.rf.fit <- randomForest::randomForest(Y ~ ., data = new.X, 
                                              ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  if (family$family == "binomial") {
    rank.rf.fit <- randomForest::randomForest(as.factor(Y) ~ 
                                                ., data = new.X, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  whichVariable <- (rank(-rank.rf.fit$importance) <= nVar)
  
  whichVariable <- colnames(X) %in% colnames(new.X[,whichVariable])
  
  return(whichVariable)
}



## ------------------------------------------------ ##

# dimen.method <- "encoder"
# if (dimen.method == "pca") {
#   cnv.dimen <- read.csv("C:/Users/da.salazarb/Documents/tmle/tmleTCGA/pca_W_cnv.csv")
#   mrna.dimen <- read.csv("C:/Users/da.salazarb/Documents/tmle/tmleTCGA/pca_W_mrna.csv")
# } else if (dimen.method == "nmf") {
#   cnv.dimen <- read.csv("C:/Users/da.salazarb/Documents/tmle/tmleTCGA/nmf_W_cnv.csv")
#   mrna.dimen <- read.csv("C:/Users/da.salazarb/Documents/tmle/tmleTCGA/nmf_W_mrna.csv")
# } else {
#   cnv.dimen <- read.csv("C:/Users/da.salazarb/Documents/tmle/tmleTCGA/encoder_W_cnv.csv")
#   mrna.dimen <- read.csv("C:/Users/da.salazarb/Documents/tmle/tmleTCGA/encoder_W_mrna.csv")
# }
# colnames(cnv.dimen) <- paste0("cnv_", dimen.method,"_", colnames(cnv.dimen))
# colnames(mrna.dimen) <- paste0("mrna_", dimen.method, "_", colnames(mrna.dimen))

# W_dimen.reduc <- data.frame(cnv.dimen, mrna.dimen)
# W_dimen.reduc <- W %>% 
# dplyr::select(., starts_with("clinical")) %>% 
# dplyr::bind_cols(mrna.dimen, cnv.dimen) #%>% 
# dplyr::select(-cnv_X, -mrna_X)  %>% 
# dplyr::bind_cols(data.frame(A)) %>% 
# dplyr::mutate(A=as.factor(A))
# W_dimen.reduc[10:15,10:15]; dim(W_dimen.reduc)

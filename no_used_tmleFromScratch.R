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
#----------------Init-------------------
#' Load functions, libraries and others
# .................................... #
# http://benkeser.github.io/sllecture/
# https://rdrr.io/cran/tmle/src/R/tmle.R#sym-.verifyArgs
# https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html#xgboost-hyperparameter-exploration
## EnsembleModels: http://homes.dsi.unimi.it/valentini/papers/ens.review.revised.pdf
## Tomado de: ? https://migariane.github.io/TMLE.nb.html
## check: https://www.khstats.com/blog/tmle/tutorial/
source("G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/func_DimenReduc.R")
# source("G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/func_lateIntegration.R")
source('G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/func_othersFromtmlePackage.R')

library(tmle)
library(caret)
library(dplyr)
library(purrr)
library(glmnet)
library(SuperLearner)
library(caretEnsemble)

#' #----------------Path-------------------
#' #' Set path of files -> pathFeatureLabelKernels
#' # .................................... #
#' if (identical(character(0), dir("D:/pathFeatureLabelKernels"))) {
#'   central_path <- "E:"
#' } else {
#'   central_path <- "D:"
#' }
#' 
#' pathFeatureLabel <- paste0(central_path, "/pathFeatureLabelKernels")
#' 
#' #----------------Data-------------------
#' #' Load the data -> O ~ (Y,A,W_I,Delta)
#' # .................................... #
#' # O <- give_O_GBM_LGG(pathFeatureLabel,dimReduc=TRUE)
#' O <- give_O_LGG(pathFeatureLabel,dimReduc=FALSE)
#' 
#' Y <- O[[1]]; A <- O[[2]]; Delta <- O[[3]]; W_ <- O[[4]]; clinical <- O[[5]];
#' W <- as.matrix(W_[,sample(1:length(colnames(W_)), 10)]);dim(W); class(W)
#' 
#' # Summary O data
#' resumeData(O)
#' 
#' ## ---------------------- ##
### ... 1. Load data ... ###
## ---------------------- ##
O <- read.csv(file = "D:/pathFeatureLabelKernels/co-md_records/O_Y_A_Delta_W.csv") ## <- Cambiar ruta del archivo!
O[1:5,1:5]; dim(O)

Y <- ifelse(O$Y > 365,  1, 0); Delta <- ifelse(O$Y > 365 & O$delta == 1, 1, 0); A <- O$A; W <- O[,!(colnames(O) %in% c("Y", "A", "X", "Delta"))];

## procesar Y
stage1 <- .initStage1(Y,A,Delta,alpha=0.005,maptoYstar=TRUE,family="binomial")
Y <- stage1$Ystar
Qbounds <- stage1$Qbounds

V <- 3 ## numero de k-fold cv
id <- 1:length(Y) ## ids de cada una de las observaciones
n <- length(Y)
folds <- vector(mode = "list", length = V) ## particiones de datos en folds se guarda en forma de lista
uid <- unique(id) # ids unicos
n.id <- length(uid) # longitud de ids unicos
id.split <- split(sample(1:n.id), rep(1:V, length = n.id)) ## crear los grupos dependiendo de V
family <- "binomial" ## familia de Y

#-------------estimateQ------------------
#' Estimating initial regression of Y on A and W
# ..................................... #
cat("\tEstimating initial regression of Y on A and W\n")
cat("\t using SuperLearner\n")
X <- data.frame(A, W) ## datos originales
X00 <- data.frame(A = 0, W) ## como si todos fueran tratados con A=0
X01 <- data.frame(A = 1, W) ## como si todos fueran tratados con A=1
newX <- rbind(X, X00, X01) ## concatenacion de datos

## Learners for Q and G
learners.svm.radial <- create.Learner("SL.svm.radial", tune = list(C = c(0.01, 0.1, 1, 10, 100, 1000), 
                                                  #kpar = c(list(sigma=0.25, degree=1),list(sigma=0.25, degree=2)), 
                                                  sigma = c(0.001, 0.01, 0.1, 1), 
                                                  # kernel = c("rbfdot", "polydot", "anovadot"), 
                                                  degree = c(1, 2, 3)))

learners.svm.poly <- create.Learner("SL.svm.poly", tune = list(C = c(0.01, 0.1, 1, 10, 100, 1000), 
                                                                   #kpar = c(list(sigma=0.25, degree=1),list(sigma=0.25, degree=2)), 
                                                                   # sigma = c(0.001, 0.01, 0.1, 1), 
                                                                   # kernel = c("rbfdot", "polydot", "anovadot"), 
                                                                   degree = c(1, 2, 3)))

Super.Q.SL.library <- c("SL.svm.linear", learners.svm.poly$names, learners.svm.radial$names) ## libreria para entrenar
Super.Q.SL.library <- c("SL.ranger") ## libreria para entrenar
SL.library <- c("SL.ranger") ## libreria para entrenar
# Super.Q.SL.library <- c("SL.xgboost") ## libreria para entrenar

# Q.SL.library <- list(c(Super.Q.SL.library, "func_protein"), 
#                      c(Super.Q.SL.library, "func_mrna"), 
#                      c(Super.Q.SL.library, "func_cnv"))

arglist <- list(Y = Y[Delta == 1], X = X[Delta == 1, ], newX = newX, SL.library = Super.Q.SL.library, 
                cvControl = list(V = V), family = family, control = list(saveFitLibrary = FALSE), 
                id = id[Delta == 1])

suppressWarnings(m <- try(do.call(SuperLearner, arglist))) ## Ajustar el SuperLearner

coef <- m$coef
keepAlg <- which(m$coef > 0)
SL.coef <- m$coef[m$coef > 0]
type <- "SuperLearner, ensemble"
type <- paste0("cv-", type)
SL.library.keep <- .expandLib(Super.Q.SL.library)[keepAlg]
predictions <- rep(NA, nrow(newX))

for (v in 1:V) {
  folds[[v]] <- which(id %in% uid[id.split[[v]]])
  TRAIN <- rep(TRUE, length(Y))
  TRAIN[folds[[v]]] <- FALSE
  arglist <- list(Y = Y[TRAIN & Delta == 1], X = X[TRAIN & Delta == 1, ], newX = newX[!TRAIN, ], SL.library = SL.library.keep, 
                  cvControl = list(V = V), family = family, control = list(saveFitLibrary = FALSE), 
                  id = id[TRAIN & Delta == 1])
  
  suppressWarnings(m <- try(do.call(SuperLearner, arglist)))
  predictions[!TRAIN] <- as.vector(as.matrix(m$library.predict) %*% SL.coef)
}

Q <- matrix(predictions, nrow = n, byrow = FALSE) ## matriz de n x 3
colnames(Q) <- c("QAW", "Q0W", "Q1W")[1:ncol(Q)] ## nombre de variables
Q <- .bound(Q, Qbounds)
head(Q)
tail(Q)

Qfamily <- family
Q <- list(Q = Q, family = Qfamily, coef = coef, type = type, 
          SL.library = SL.library)

#-------------estimateG------------------
#' Estimating initial regression of A on W
# ..................................... #
gbound= 5/sqrt(length(Y))/log(length(Y))
gbound.ATT <- c(gbound, 1-gbound)
gbound <- c(min(gbound), 1)

# Screen covariates for association with Y or residual 
prescreenW.g <- TRUE
min.retain <- 2
RESID <- FALSE
Q.offset <-  plogis(Q$Q[,"QAW"])
retain.W <- .prescreenW.g(stage1$Ystar, A, as.matrix(W), Delta , QAW = Q.offset,  family = family, min.retain, RESID=RESID)

## estimateG
# g <- suppressWarnings(estimateG(d=data.frame(A,W[,retain.W]), g1W, gform, g.SL.library, id=id, V = V, verbose, "treatment mechanism", outcome="A",  discreteSL = g.discreteSL)) 
m <- NULL
coef <- NA
type <- NULL
d <- data.frame(A,W[,retain.W])
g1W <- rep(1,nrow(d)) ## dimension de g1W
message <- "treatment mechanism"
# newdata <- d
type <- paste("No", strsplit(message, " ")[[1]][1])
arglist <- list(Y=d[,1], X=d[,-1, drop=FALSE], newX=d[,-1, drop=FALSE], family="gaussian", SL.library=SL.library, cvControl=list(V=V), id=id)
suppressWarnings(m <- try(do.call(SuperLearner,arglist)))
g1W <- as.vector(m$SL.predict)
g <- list(g1W=g1W, coef=coef, type=type)
g$SL.library=SL.library
g$coef <- m$coef

### paso final
g.z <- NULL
g.z$type="No intermediate variable"
g.z$coef=NA
pDelta1 <- NULL
verbose <- TRUE
g.Deltaform <- NULL
g.Delta.SL.library <- SL.library
g.Delta.discreteSL <- FALSE
# if prescreenW.g, keep all the variables in the g.Deltaform, or keep a minimum of 2 variables. Always keep A.
g.Delta <- suppressWarnings(estimateG(d=data.frame(Delta, Z=1, A, W[,retain.W]), pDelta1, g.Deltaform,# newdata = data.frame(Delta, Z=1, A, W[,retain.W]),
                                      SL.library = g.Delta.SL.library, id=id, V = V, verbose = verbose, "missingness mechanism", outcome="D",  discreteSL= g.Delta.discreteSL)) 
g$bound.ATT <- gbound.ATT
g$bound <- gbound
g1W.total <- .bound(g$g1W*g.Delta$g1W[,"Z0A1"], g$bound)
g0W.total <- .bound((1-g$g1W)*g.Delta$g1W[,"Z0A0"], g$bound)
if(all(g1W.total==0)){g1W.total <- rep(10^-9, length(g1W.total))}
if(all(g0W.total==0)){g0W.total <- rep(10^-9, length(g0W.total))}

target.gwt <- TRUE

if(target.gwt){
  wt <- A/g1W.total + (1-A)/g0W.total
  H1W <- A
  H0W <- 1-A
} else{
  wt <- rep(1, length(A))
  H1W <- A/g1W.total
  H0W <- (1-A)/g0W.total
}

suppressWarnings(epsilon <- coef(glm(stage1$Ystar~-1 + offset(Q$Q[,"QAW"]) + H0W + H1W, family=Q$family, weights = wt, subset=Delta==1)))
epsilon[is.na(epsilon)] <- 0  # needed for EY1 calculation
if (target.gwt){
  Qstar <- Q$Q + c((epsilon[1]*H0W + epsilon[2]*H1W), rep(epsilon[1], length(Y)), rep(epsilon[2], length(Y)))
} else {
  Qstar <- Q$Q + c((epsilon[1]*H0W + epsilon[2]*H1W), epsilon[1]/g0W.total, epsilon[2]/g1W.total)
}
colnames(Qstar) <- c("QAW", "Q0W", "Q1W")
Ystar <- stage1$Ystar
maptoYstar <- TRUE
if (maptoYstar) {
  Qstar <- plogis(Qstar)*diff(stage1$ab)+stage1$ab[1] 
  Q$Q <- plogis(Q$Q)*diff(stage1$ab)+stage1$ab[1]
  Ystar <- Ystar*diff(stage1$ab)+stage1$ab[1]
} else if (family == "poisson"){  	
  Q$Q <- exp(Q$Q)				  
  Qstar <- exp(Qstar)
}
colnames(Q$Q) <- c("QAW", "Q0W", "Q1W")
#Q$Q <- Q$Q[,-1]
res <- calcParameters(Ystar, A, I.Z=rep(1, length(Ystar)), Delta, g1W.total, g0W.total, Qstar, 
                      mu1=mean(Qstar[,"Q1W"]), mu0=mean(Qstar[,"Q0W"]), id, family)
#ATT  & ATC - additive effect: psi, CI, pvalue 
ATT <- ATC <- IC.ATT <- IC.ATE <- NULL
if(length(unique(A)) > 1){
  n.id <- length(unique(id))
  depsilon <-  0.001
  # browser()
  Q.ATT <- (Q$Q- stage1$ab[1])/diff(stage1$ab)  # changed to should be Q$Q
  res.ATT <- try(oneStepATT(Y = stage1$Ystar, A = A, Delta, Q = Q.ATT, g1W = g$g1W, 
                            pDelta1 = g.Delta$g1W[,c("Z0A0", "Z0A1")],
                            depsilon = depsilon, max_iter = max(1000, 2/depsilon), gbounds = gbound.ATT, Qbounds = stage1$Qbound))
  if(!(identical(class(res.ATT), "try-error"))){
    ATT$psi <- res.ATT$psi * diff(stage1$ab) 
    ATT$converged <- res.ATT$conv
    if(n.id < length(id)){
      IC.ATT <- as.vector(by(res.ATT$IC, id, mean)) * diff(stage1$ab) + stage1$ab[1]
    } else {
      IC.ATT <- res.ATT$IC * diff(stage1$ab) + stage1$ab[1]
    }
    ATT$var.psi <- var(IC.ATT)/n.id
    ATT$CI <- c(ATT$psi -1.96 *sqrt(ATT$var.psi), ATT$psi +1.96 *sqrt(ATT$var.psi))
    ATT$pvalue <- 2*pnorm(-abs(ATT$psi/sqrt(ATT$var.psi)))	
  }
  
  res.ATC <- try(oneStepATT(Y = stage1$Ystar, A = 1-A, Delta, Q = cbind(QAW = Q.ATT[,1], Q0W = Q.ATT[,"Q1W"], Q1W = Q.ATT[,"Q0W"]), 
                            g1W = 1-g$g1W, pDelta1 = g.Delta$g1W[,c("Z0A1", "Z0A0")],
                            depsilon = depsilon, max_iter = max(1000, 2/depsilon), gbounds = gbound.ATT, Qbounds = stage1$Qbound))
  
  if(!(identical(class(res.ATC), "try-error"))){
    ATC$psi <- -res.ATC$psi * diff(stage1$ab) 
    ATC$converged <- res.ATC$conv
    if(n.id < length(id)){
      IC.ATC <- as.vector(by(res.ATC$IC, id, mean)) * diff(stage1$ab) + stage1$ab[1]
    } else {
      IC.ATC <- res.ATC$IC * diff(stage1$ab) + stage1$ab[1]
    }
    ATC$var.psi <- var(IC.ATC)/n.id
    ATC$CI <- c(ATC$psi -1.96 *sqrt(ATC$var.psi), ATC$psi +1.96 *sqrt(ATC$var.psi))
    ATC$pvalue <- 2*pnorm(-abs(ATC$psi/sqrt(ATC$var.psi)))	
  }
  res$ATT <- ATT
  res$ATC <- ATC
  res$IC$IC.ATT <- IC.ATT
  res$IC$IC.ATC <- IC.ATC
}

# calculate Rsq - complete case
#  browser()
m.rsq <- glm(Y~ 1, family = family)
Yhat <- predict(m.rsq, newdata = data.frame(A),  type = "response")
Q$Rsq  <-  NULL

if(family == "binomial"){
  Q$Rsq <- 1 - sum(Delta * ( Y * log(Q$Q[,"QAW"]) + (1-Y)*log(1-Q$Q[,"QAW"]))) /
    sum(Delta * ( Y * log(Yhat) + (1-Y)*log(1-Yhat)))		
  Q$Rsq.type <- "pseudo R squared"
} else {
  Q$Rsq <- 1 - sum((Y[Delta == 1]-Q$Q[Delta == 1,"QAW"])^2)/ sum((Y[Delta == 1]-mean(Y[Delta==1]))^2)
  Q$Rsq.type <- "R squared"
}
cvQinit <- TRUE
if (cvQinit){
  Q$Rsq.type <- paste("Cross-validated", Q$Rsq.type)
} else {
  Q$Rsq.type <- paste("Empirical", Q$Rsq.type)
}
Q$Q <- Q$Q[,-1]
returnVal <- list(estimates=res, Qinit=Q, g=g, g.Z=g.z, g.Delta=g.Delta, Qstar=Qstar[,-1], epsilon=epsilon, gbound = gbound, gbound.ATT = gbound.ATT, W.retained = colnames(W)[retain.W]) 

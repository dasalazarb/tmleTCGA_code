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
.initStage1 <- function(Y,A, Delta, alpha, maptoYstar, family){
# .initStage1 <- function(Y,A, Q, Q.Z1=NULL, Delta, Qbounds, alpha, maptoYstar, family){
  Q <- NULL
  Qbounds <- NULL
  if(family=="binomial") {Qbounds <- c(0,1)}
  if(is.null(Qbounds)) {
    if(maptoYstar){ 
      Qbounds <- range(Y[Delta==1])
      Qbounds <- Qbounds + .01*c(-abs(Qbounds[1]),abs(Qbounds[2]))
    } else {
      Qbounds <- c(-Inf, Inf)
    }
  }
  if(!is.null(Q)){
    QAW <- (1-A)*Q[,1] + A*Q[,2]
    Q <- cbind(QAW, Q0W=Q[,1], Q1W=Q[,2])
  }
  
  # if(!is.null(Q.Z1)){
  #   Q <- cbind(Q, Q0W.Z1=Q.Z1[,1], Q1W.Z1=Q.Z1[,2])
  # }
  ab <- c(0,1)
  Ystar <- Y
  if(maptoYstar){ 
    Ystar <- .bound(Y, Qbounds)
    if(!is.null(Q)){
      Q <- .bound(Q, Qbounds)
    }
    if(0 >= alpha | 1 <= alpha){
      alpha <- .995
      warning(paste("\n\talpha must be between 0 and 1, alpha reset to",alpha,"\n"),
              immediate. = TRUE)
    }
    ab <- range(Ystar, na.rm=TRUE)
    Ystar[is.na(Ystar)] <- 0  
    Ystar <- (Ystar-ab[1])/diff(ab)	
    if(!is.null(Q)){Q <- (Q-ab[1])/diff(ab)}
    Qbounds <- c(alpha, 1-alpha)
  }
  return(list(Ystar=Ystar, Q=Q, Qbounds=Qbounds, ab=ab))
} 

.bound <- function(x, bounds){
  x[x>max(bounds)] <- max(bounds)
  x[x<min(bounds)] <- min(bounds)
  return(x)
}

.expandLib <- function(SL.lib) {
  if (is.list(SL.lib)) {
    counts <- sapply(SL.lib, length)
    numExtra <- sum(counts > 2)
    m <- matrix("", nrow = length(SL.lib) + numExtra, 
                ncol = 2)
    rowIndex <- 1
    for (i in 1:length(counts)) {
      if (counts[i] == 1) {
        m[rowIndex, ] <- c(SL.lib[[i]], "All")
      }
      else {
        m[rowIndex, ] <- SL.lib[[i]][1:2]
      }
      rowIndex <- rowIndex + 1
      if (counts[i] > 2) {
        for (j in 3:counts[i]) {
          m[rowIndex, ] <- SL.lib[[i]][c(1, j)]
          rowIndex <- rowIndex + 1
        }
      }
    }
    return(split(m, 1:nrow(m)))
  }
  else {
    return(SL.lib)
  }
}

#---------------.prescreenW.g----------------------
# Screen covariates for association with Y or residual 
# depending on value of RESID flag
# after initial regression to use in estimating g
# selects lasso covariates with non-zero coefficients
# using initial fit as offset, so it's a function of association with residuals
# 	Y - outcome on same scale as Q
#  W - eligible covariates for screening
#  Delta - indicator the outcome is observed
#  QAW - initial estimate of QAW on appropriate scale for offset
# family = "binomial", "gaussian", or "poisson", set by estimateQ
# min.retain = minimum number of of variables to retain
#------------------------------------------------------
.prescreenW.g <- function(Y, A, W, Delta, QAW, family, min.retain, RESID=FALSE){  
  if(NCOL(W) < min.retain){
    min.retain <- NCOL(W)
  }
  #require(glmnet)
  if(RESID){
    offset <- QAW[Delta == 1]
  } else {
    offset <- NULL
  }
  m.lasso <- cv.glmnet(W[Delta == 1,], Y[Delta == 1], family =  family, offset = offset)
  beta <- coef(m.lasso, s = m.lasso$lambda.min)[-1] # ignore the intercept
  retain <- which(abs(beta) > 0)
  if (length(retain) < min.retain ){
    if (length(unique(A)) == 1){
      retain <-  unique( c(retain, order(abs(cor(Delta, W)))[1:min.retain]))[1:min.retain]
    } else {
      retain <- unique( c(retain, order(abs(cor(A, W)))[1:min.retain]))[1:min.retain]
    }
  } 
  return(retain)
}
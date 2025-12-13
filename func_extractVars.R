# Column-selection helpers for clinical, mRNA and CNV features.
# Functions return logical vectors suitable for SuperLearner screening.

#' Select mRNA columns (raw features)
func_mrna <- function(X, ...) {
  grepl("mrna_|clinical_", names(X)) & !grepl("encoder_|pca_|nmf_|cnv_", names(X))
}

#' Select CNV columns (raw features)
func_cnv <- function(X, ...) {
  grepl("cnv_|clinical_", names(X)) & !grepl("encoder_|pca_|nmf_|mrna_", names(X))
}

#' Select clinical variables only
func_clinical <- function(X, ...) {
  grepl("clinical_", names(X))
}

#' Select mRNA embeddings from autoencoders
func_mrna.encoder <- function(X, ...) {
  grepl("mrna_encoder_|clinical_", names(X))
}

#' Select CNV embeddings from autoencoders
func_cnv.encoder <- function(X, ...) {
  grepl("cnv_encoder_|clinical_", names(X))
}

#' Select mRNA PCA components
func_mrna.pca <- function(X, ...) {
  grepl("mrna_pca_|clinical_", names(X))
}

#' Select CNV PCA components
func_cnv.pca <- function(X, ...) {
  grepl("cnv_pca_|clinical_", names(X))
}

#' Select mRNA NMF components
func_mrna.nmf <- function(X, ...) {
  grepl("mrna_nmf_|clinical_", names(X))
}

#' Select CNV NMF components
func_cnv.nmf <- function(X, ...) {
  grepl("cnv_nmf_|clinical_", names(X))
}

#' Wrapper around screen.randomForest with a higher variable budget
new.screen.randomForest <- function(...) {
  SuperLearner::screen.randomForest(..., nVar = 300)
}

#' RandomForest prescreening for mRNA variables only
new.screen.rf.mrna <- function(Y, X, family, nVar = 300, ntree = 1000,
                               mtry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)),
                               nodesize = ifelse(family$family == "gaussian", 5, 3), maxnodes = NULL, ...) {
  selected <- X[, func_mrna(X), drop = FALSE]
  if (!requireNamespace("randomForest", quietly = TRUE)) stop("randomForest package is required")
  if (family$family == "gaussian") {
    rf_fit <- randomForest::randomForest(Y ~ ., data = selected, ntree = ntree, mtry = mtry,
                                         nodesize = nodesize, keep.forest = FALSE, maxnodes = maxnodes)
  } else {
    rf_fit <- randomForest::randomForest(as.factor(Y) ~ ., data = selected, ntree = ntree,
                                         mtry = mtry, nodesize = nodesize, keep.forest = FALSE, maxnodes = maxnodes)
  }
  important <- rank(-rf_fit$importance) <= nVar
  colnames(X) %in% colnames(selected)[important]
}

#' RandomForest prescreening for CNV variables only
new.screen.rf.cnv <- function(Y, X, family, nVar = 300, ntree = 1000,
                              mtry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)),
                              nodesize = ifelse(family$family == "gaussian", 5, 3), maxnodes = NULL, ...) {
  selected <- X[, func_cnv(X), drop = FALSE]
  if (!requireNamespace("randomForest", quietly = TRUE)) stop("randomForest package is required")
  if (family$family == "gaussian") {
    rf_fit <- randomForest::randomForest(Y ~ ., data = selected, ntree = ntree, mtry = mtry,
                                         nodesize = nodesize, keep.forest = FALSE, maxnodes = maxnodes)
  } else {
    rf_fit <- randomForest::randomForest(as.factor(Y) ~ ., data = selected, ntree = ntree,
                                         mtry = mtry, nodesize = nodesize, keep.forest = FALSE, maxnodes = maxnodes)
  }
  important <- rank(-rf_fit$importance) <= nVar
  colnames(X) %in% colnames(selected)[important]
}

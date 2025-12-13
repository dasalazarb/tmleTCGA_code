# Utility functions for common dimensionality-reduction workflows.
# Each function receives a feature matrix and returns a transformed
# representation that can be used in downstream models.

#' Kernel PCA with an RBF kernel
#' @param W Numeric matrix with samples in rows.
#' @return Matrix of transformed principal components.
func_kpca <- function(W) {
  stopifnot(is.matrix(W))
  kpc <- kernlab::kpca(~., data = W, kernel = "rbfdot", kpar = list(sigma = 2), features = 100)
  W.kernel <- kernlab::pcv(kpc)
  rownames(W.kernel) <- rownames(W)
  colnames(W.kernel) <- paste0("PC.", seq_len(ncol(W.kernel)))
  W.kernel
}

#' Standard PCA with automatic rank selection (up to 100 PCs)
#' @param W Numeric matrix with samples in rows.
#' @return PCA scores matrix.
func_pca <- function(W) {
  stopifnot(is.matrix(W))
  rank_limit <- if (nrow(W) > ncol(W)) 100 else NULL
  pca <- stats::prcomp(W, scale. = TRUE, rank. = rank_limit)
  pca$x
}

#' Non-negative matrix factorization
#' @param W Non-negative matrix.
#' @param rango Number of components to extract.
#' @return Basis matrix from the fitted NMF model.
func_nmf <- function(W, rango) {
  stopifnot(is.matrix(W), rango > 0)
  nmf_fit <- NMF::nmf(W, rango, "lee", "nndsvd")
  as.matrix(NMF::basis(nmf_fit))
}

#' Sparse PCA using PMA::SPC
#' @param W Numeric matrix with samples in rows.
#' @return Sparse principal components matrix.
func_sPCA <- function(W) {
  stopifnot(is.matrix(W))
  cv.out <- PMA::SPC.cv(as.matrix(W))
  spc <- PMA::SPC(as.matrix(W), sumabsv = cv.out$bestsumabsv, K = 50)
  spc$u
}

#' Simple dense autoencoder bottleneck representation (2-D)
#' @param W Numeric matrix with samples in rows.
#' @return Matrix with two-dimensional bottleneck embeddings.
func_encoder <- function(W) {
  stopifnot(is.matrix(W))

  model <- keras::keras_model_sequential() %>%
    keras::layer_dense(units = 6, activation = "tanh", input_shape = ncol(W)) %>%
    keras::layer_dense(units = 2, activation = "tanh", name = "bottleneck") %>%
    keras::layer_dense(units = 6, activation = "tanh") %>%
    keras::layer_dense(units = ncol(W))

  model %>% keras::compile(loss = "mean_squared_error", optimizer = "adam")
  model %>% keras::fit(x = W, y = W, epochs = 200, verbose = 0)

  bottleneck <- keras::keras_model(inputs = model$input, outputs = keras::get_layer(model, "bottleneck")$output)
  keras::predict(bottleneck, W)
}

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
library(NMF)
library(ANN2)
library(PMA) # https://rdrr.io/cran/PMA/man/SPC.html
# library(h2o) # https://bradleyboehmke.github.io/HOML/autoencoders.html
library(keras)
library(kernlab) # https://www.rdocumentation.org/packages/kernlab/versions/0.9-29/topics/kpca
library(nsprcomp) # https://www.rdocumentation.org/packages/nsprcomp/versions/0.5.1-2/topics/nsprcomp
library(pcaMethods) # https://www.bioconductor.org/packages/release/bioc/vignettes/pcaMethods/inst/doc/pcaMethods.pdf
# h2o.init( port = 54321, startH2O = TRUE)

func_kpca <- function(W) {
  kpc <- kernlab::kpca(~.,data=W,kernel="rbfdot",
              kpar=list(sigma=2),features=100)
  W.kernel <- kernlab::pcv(kpc)
  row.names(W.kernel) <- row.names(W)
  colnames(W.kernel) <- paste0("PC.", rep(1:dim(W.kernel)[2]))
  return(W.kernel)
}

func_pca <- function(W) {
  if (dim(W)[1] > dim(W)[2]) {
    pca <- prcomp(W, scale = TRUE, rank. = 100)
  } else {
    pca <- prcomp(W, scale = TRUE)
  }
  W.kernel <- pca$x
  return(W.kernel)
}

func_nmf <- function(W,rango) {
  nmf_ <- NMF::nmf(W, rango, "lee", "nndsvd")
  W.kernel <- NMF::basis(nmf_)
  W.kernel <- as.matrix(W.kernel)
  return(W.kernel)
}

func_sPCA <- function(W) {
  # Performs sparse principal components analysis by applying PMD 
  # to a data matrix with lasso ($L_1$) penalty on the columns and no penalty on the rows.
  cv.out <- PMA::SPC.cv(as.matrix(W))
  W.kernel <- PMA::SPC(as.matrix(W), sumabsv=cv.out$bestsumabsv, K=50)
  W.kernel <- W.kernel$u
  return(W.kernel)
}

func_encoder_2ANN <- function() {
  AE <- autoencoder(X = W_,
                    hidden.layers = c(10,3,10),
                    loss.type = 'pseudo-huber',
                    optim.type = 'adam',
                    n.epochs = 5000)
  AE$Rcpp_ANN
}

func_encoder <- function(W) {
  W <- as.matrix(W)
  
  model <- keras::keras_model_sequential()
  model %>%
    layer_dense(units = 6, activation = "tanh", input_shape = ncol(W)) %>%
    layer_dense(units = 2, activation = "tanh", name = "bottleneck") %>%
    layer_dense(units = 6, activation = "tanh") %>%
    layer_dense(units = ncol(W))
  
  # view model layers
  summary(model)
  
  # compile model
  model %>% compile(
    loss = "mean_squared_error", 
    optimizer = "adam"
  )
  
  # fit model
  model %>% keras::fit(
    x = W, 
    y = W, 
    epochs = 2000,
    verbose = 0
  )
  # evaluate the performance of the model
  mse.ae2 <- evaluate(model, W, W)
  mse.ae2
  
  # extract the bottleneck layer
  intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
  intermediate_output <- predict(intermediate_layer_model, W)
  
  ggplot(data.frame(PC1 = intermediate_output[,1], PC2 = intermediate_output[,2]), aes(x = PC1, y = PC2)) + geom_point()
}

# func_nsPCA <- function(W) {
#   W.kernel <- pcaMethods::pca(as.matrix(W), method = "nipals", center = TRUE, nPcs = 50)
# }

# func_encoder <- function(W){
#   features <- as.h2o(W)
#   
#   # Hyperparameter search grid
#   hyper_grid <- list(hidden = list(
#     c(100, 75 ,50, 75, 100),
#     c(75, 50, 25, 50, 75), 
#     c(250, 100, 250),
#     c(50, 25, 10, 25 ,50),
#     c(500, 250, 100, 50, 100, 250, 500)
#   ))
#   
#   # Execute grid search
#   ae_grid <- h2o.grid(
#     algorithm = 'deeplearning',
#     x = seq_along(features),
#     training_frame = features,
#     grid_id = 'autoencoder_grid',
#     autoencoder = TRUE,
#     activation = 'Tanh',
#     hyper_params = hyper_grid,
#     sparse = TRUE,
#     ignore_const_cols = FALSE,
#     seed = 123
#   )
#   
#   # Print grid details
#   a <- h2o.getGrid('autoencoder_grid', sort_by = 'mse', decreasing = FALSE)
#   b <- a@summary_table
#   best_model_id <- b$model_ids[1]
#   best_model <- h2o.getModel(best_model_id)
#   hidden_layer <- ceiling(length(strsplit(b$hidden[1], ",")[[1]])/2)
#   
#   W.kernel <- h2o.deepfeatures(best_model, features, layer = hidden_layer)
#   return(W.kernel)
# }
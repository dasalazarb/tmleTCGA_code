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
library(dplyr)
library(purrr)
library(doParallel)
library(missForest)

read_profile <- function(path_feature_label, filename) {
  read.csv(
    file = file.path(path_feature_label, "co-mod_ccle_tcga_depurados", filename),
    header = TRUE,
    sep = ",",
    row.names = 1
  )
}

give_O_data <- function(pathFeatureLabel) {
  ## --------------------------------- ##
  ## data TCGA
  message("1. loading profiles")
  W_protein <- read_profile(pathFeatureLabel, "tcga_protein_tmle.csv")
  W_mirna <- read_profile(pathFeatureLabel, "tcga_mirna_tmle.csv")
  W_cnv <- read_profile(pathFeatureLabel, "tcga_cnv_tmle.csv")
  W_mrna <- read_profile(pathFeatureLabel, "tcga_mrna_tmle.csv")

  genes_DEA <- read.csv(file.path(pathFeatureLabel, "LGG", "DEA_LGG.csv"), header = TRUE)
  unicos_genes <- intersect(as.vector(genes_DEA$X), colnames(W_mrna))
  W_mrna <- W_mrna[, unicos_genes, drop = FALSE]

  log_dimensions <- function(name, data) {
    message(sprintf("    * init %s dimensions: %s x %s", name, nrow(data), ncol(data)))
  }

  purrr::walk2(
    c("W_protein", "W_mirna", "W_cnv", "W_mrna"),
    list(W_protein, W_mirna, W_cnv, W_mrna),
    log_dimensions
  )

  ### SuperLearner as a early integration method
  W_ <- cbind(W_protein, W_mirna, W_cnv, W_mrna)
  message(sprintf("    ** init W_ dimensions: %s x %s", nrow(W_), ncol(W_)))
  message("---------------------------------")

  ## --------------------------------- ##
  message("2. loading clinical info")
  clinical <- read.csv(
    file.path(pathFeatureLabel, "co-mod_best_results_last_run", "clinical_data_TCGA.csv"),
    header = TRUE,
    row.names = 1
  )
  message(sprintf("    init clinical dimensions: %s x %s", nrow(clinical), ncol(clinical)))
  message("---------------------------------")

  ## --------------------------------- ##
  ## outcome  + right-censored
  message("3. loading Y and delta info")
  Y <- read.csv(
    file.path(pathFeatureLabel, "co-mod_bases_datos_importantes", "survivalTCGA.csv"),
    header = TRUE,
    sep = ","
  )
  Y_delta <- Y %>% ## delta hace referencia a los right-censored obs
    dplyr::select(bcr_patient_barcode, type, PFI, PFI.time) %>%
    dplyr::filter(type == "LGG")

  row.names(Y_delta) <- Y_delta$bcr_patient_barcode
  message(sprintf("    * init Y_delta dimensions: %s x %s", nrow(Y_delta), ncol(Y_delta)))
  message("---------------------------------")

  ## --------------------------------- ##
  ## treatment -> archivo LGG_clinic_data.R ->
  ## radiation+temozolomide (125), radiation (49), temozolomide (28)
  message("4. loading A info")
  A <- read.csv(
    file.path(pathFeatureLabel, "co-mod_best_results_last_run", "radiation_treatment_cycles_TCGA_LGG.csv"),
    header = TRUE,
    sep = ","
  )
  A <- A %>%
    dplyr::filter(treat %in% c("radiation+temozolomide", "radiation")) %>% # , "temozolomide"
    dplyr::mutate(treat = ifelse(treat == "radiation+temozolomide", 1, 0))
  row.names(A) <- A$bcr_patient_barcode
  message(sprintf("    * init A dimensions: %s x %s", nrow(A), ncol(A)))
  message("---------------------------------")

  ## --------------------------------- ##
  ## pacientes unicos
  message("5. common patients")
  pacientes_comunes <- Reduce(
    intersect,
    list(row.names(W_), row.names(A), row.names(Y_delta), row.names(clinical))
  )

  A <- A[pacientes_comunes, ]
  A <- A %>%
    dplyr::select(treat)

  A <- as.vector(A$treat)

  Y_delta <- Y_delta[pacientes_comunes, ]
  Y_delta <- Y_delta %>%
    dplyr::select(PFI, PFI.time)
  Y <- as.numeric(Y_delta$PFI.time)
  delta <- as.numeric(as.vector(Y_delta$PFI))

  clinical <- clinical[pacientes_comunes, ]
  colnames(clinical) <- paste0("clinical_", colnames(clinical))
  W <- W_[pacientes_comunes, ]
  message("---------------------------------")

  ## --------------------------------- ##
  ## imputar faltantes en clinical data
  message("6. to impute missing values in clinical data")
  registerDoParallel(cores = 3)
  W_ <- cbind(W, clinical)
  clinical <- missForest(W_, maxiter = 50, variablewise = TRUE, parallelize = "variables", mtry = 100)
  clinical <- clinical$ximp %>%
    dplyr::select(contains("clinical_"))
  clinical <- fastDummies::dummy_cols(clinical, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
  detach("package:doParallel", unload = TRUE)
  message("---------------------------------")

  ## --------------------------------- ##
  message("#. Final results:")
  message(sprintf("    *** final W.init.concat dimension: %s x %s", nrow(W_), ncol(W_)))
  message(sprintf("    *** final clinical dimension: %s x %s", nrow(clinical), ncol(clinical)))
  message(sprintf(
    "    *** final A dimension: %s - Summary (1: radio+TMZ, 0: radio): %s x %s",
    length(A), summary(as.factor(A))[[2]], summary(as.factor(A))[[1]]
  ))
  message(sprintf("    *** final Y dimension: %s - Summary:", length(Y)))
  message(capture.output(summary(Y)))
  message(sprintf(
    "    *** final delta dimension: %s - Summary (1: obsOutcome, 0: right-censored): %s x %s",
    length(delta), summary(as.factor(delta))[[2]], summary(as.factor(delta))[[1]]
  ))

  return(list(Y, A, delta, W, clinical))
}

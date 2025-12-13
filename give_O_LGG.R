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
library(purrr)
library(dplyr)
library(SuperLearner)
library(TCGAbiolinks)

if (identical(character(0), dir("D:/pathFeatureLabelKernels"))) {
  pathFeatureLabel <- "G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/data"
} else {
  pathFeatureLabel <- "D:"
}

read_profile <- function(path_feature_label, profile) {
  read.csv(file.path(path_feature_label, paste0("tcga_LGG_", profile, "_depurado_tmle.csv")))
}

clean_profile <- function(df) {
  df <- df[match(unique(df$X), df$X), ]
  row.names(df) <- df$X
  df$X <- NULL
  df
}

variance_filter <- function(df, lower = 0.1, upper = 0.9) {
  vars <- apply(df, 2, var)
  bounds <- quantile(vars, c(lower, upper), na.rm = TRUE)
  df[, vars > bounds[1] & vars < bounds[2], drop = FALSE]
}

normalize_and_filter <- function(W) {
  normalized <- apply(W, 2, function(x) {
    range_val <- max(x) - min(x)
    if (range_val == 0) {
      return(rep(0, length(x)))
    }
    (x - min(x)) / range_val
  })
  normalized <- data.frame(normalized)
  zeros <- apply(normalized, 2, function(x) sum(x == 0))
  normalized[, zeros < round(nrow(normalized) * 0.40), drop = FALSE]
}

log_dimensions <- function(name, df) {
  message(sprintf("    * %s dimensions: %s x %s", name, nrow(df), ncol(df)))
}

give_O_GBM_LGG <- function(pathFeatureLabel, dimReduc, outTMLE3, index, Y_as_NaN) {

  if (!file.exists(file.path(pathFeatureLabel, "co-mod_best_results_last_run", "all_TCGA_LGG_GBM.csv"))) {
    message("... 1. Load cnv and mrna profiles + merge these profiles (rbind) ...")
    profiles <- c("mrna", "cnv")
    g <- map(profiles, ~ read_profile(pathFeatureLabel, .x) %>% clean_profile())
    g[[2]] <- t(na.omit(t(g[[2]])))

    g <- map(g, variance_filter, lower = 0.1, upper = 0.9)
    walk2(profiles, g, ~ log_dimensions(.x, .y))

    message("... 2. Load lgg info and build glioma_info ...")
    lgg_info <- read.csv(file.path(pathFeatureLabel, "Y_A_clinical_LGG.csv"), row.names = 1)

    glioma_info <- lgg_info %>%
      dplyr::select(-OS.time, -days_to_therapy_start, -treatment, -Telomere_Maintenance, -Grade) %>%
      dplyr::rename(delta = OS, Y = new.OS.time, A = treatment_0_1)
    colnames(glioma_info)[1:11] <- paste0("clinical_", colnames(glioma_info)[1:11])
    row.names(glioma_info) <- glioma_info$clinical_bcr_patient_barcode

    message("... 3. Column bind of cnv and mrna profiles for lgg ...")
    W_ <- cbind(g[[1]], g[[2]])
    W_$X <- substr(row.names(W_), 1, 12)
    W_ <- W_[match(unique(W_$X), W_$X), ]
    row.names(W_) <- W_$X
    W_$X <- NULL

    message(" ... 4. Common patients for clinical and omic profiles for lgg ...")
    pacientes_comunes <- Reduce(intersect, list(row.names(W_), row.names(glioma_info)))
    write(pacientes_comunes, file.path(pathFeatureLabel, "barcodes_gbm_lgg_to_tmle.csv"))

    W <- W_[pacientes_comunes, ]
    W <- normalize_and_filter(W)
    message("---------------------------------")

    all <- cbind(glioma_info[pacientes_comunes, ], W)
    row.names(all) <- all$clinical_bcr_patient_barcode

    write(all$bcr_patient_barcode, file = file.path(pathFeatureLabel, "barcodes_TCGA_LGG.csv"))

    clinicalData <- all %>%
      dplyr::select(-A, -Y, -delta, -clinical_TERT_promoter_status, -clinical_bcr_patient_barcode) %>%
      dplyr::mutate_if(is.character, as.factor)

    library(missForest)
    library(doParallel)
    registerDoParallel(cores = 7)

    clinicalData.imp <- missForest(clinicalData[all$Y != 0, ], parallelize = "variables", ntree = 200, maxnodes = 25, mtry = 2000, verbose = TRUE)

    for (i in c(50, 100, 500, 1000)) {
      mrna.variables <- colnames(all %>% dplyr::filter(Y != 0) %>% dplyr::select(., starts_with("mrna_")))[
        screen.randomForest(
          Y = ifelse(all$delta[all$Y != 0] != 0, 1, 0),
          X = all %>% dplyr::filter(Y != 0) %>% dplyr::select(., starts_with("mrna_")),
          family = binomial(),
          nVar = i
        )
      ]
      all_no.mrna <- colnames(clinicalData.imp$ximp %>% dplyr::select(., -starts_with("mrna_")))

      all_ <- clinicalData.imp$ximp %>%
        dplyr::select(dplyr::one_of(all_no.mrna), dplyr::one_of(mrna.variables)) %>%
        dplyr::bind_cols(all %>% dplyr::filter(Y != 0) %>% dplyr::select(Y, A, delta)) %>%
        dplyr::relocate(Y, A, delta)

      write.csv(x = all_, file = file.path(pathFeatureLabel, paste0("imputation_", i, ".csv")))
    }

    A <- all$A
    Y <- all$Y
    delta <- all$delta
    W <- all %>% dplyr::select(., starts_with("mrna_"), starts_with("cnv_"))
    clinical <- all %>%
      dplyr::select(., starts_with("clinical_")) %>%
      dplyr::mutate_if(is.character, factor)

    message("\tSummary of O ~ (Y,A,W_I,Delta):")
    log_dimensions("final W.init.concat", W)
    log_dimensions("final clinical", clinical)
    message(sprintf(
      "    *** final A dimension: %s - Summary (0: RAD+TMZ, 1: RAD): %s x %s",
      length(A), summary(as.factor(A))[[2]], summary(as.factor(A))[[1]]
    ))
    message(sprintf("    *** final Y dimension: %s - Summary:", length(Y)))
    message(capture.output(summary(Y)))
    message(sprintf(
      "    *** final delta dimension: %s - Summary (1: obsOutcome, 0: right-censored): %s x %s",
      length(delta), summary(as.factor(delta))[[2]], summary(as.factor(delta))[[1]]
    ))

    return(list(Y, A, delta, W, clinical))

  } else if (dimReduc == TRUE & outTMLE3 == FALSE) {
    message("\t * Loading W with reduction of dimensionality - NMF (func_NMF.py) - ")
    all <- read.csv(file.path(pathFeatureLabel, "co-mod_best_results_last_run", "Y_A_delta_clinical_TCGA_LGG_GBM.csv"), header = TRUE, sep = ",")
    mrna.W <- read.csv(file.path(pathFeatureLabel, "co-mod_best_results_last_run", "NMF_W_mrna_LGG.csv"), header = FALSE, sep = ",")
    colnames(mrna.W) <- paste0("mrna_", colnames(mrna.W))
    protein.W <- read.csv(file.path(pathFeatureLabel, "co-mod_best_results_last_run", "NMF_W_protein_LGG.csv"), header = FALSE, sep = ",")
    colnames(protein.W) <- paste0("protein_", colnames(protein.W))
    cnv.W <- read.csv(file.path(pathFeatureLabel, "co-mod_best_results_last_run", "NMF_W_cnv_LGG.csv"), header = FALSE, sep = ",")
    colnames(cnv.W) <- paste0("cnv_", colnames(cnv.W))
    W <- cbind(mrna.W, protein.W, cnv.W)
    A <- all$A
    Y <- all$Y
    delta <- all$delta
    clinical <- all %>%
      dplyr::select(contains("clinical_"))

    nombres_comunes_grupos <- intersect(clinical$clinical_submitter_id, grupos$short_index)
    grupos <- grupos[grupos$short_index %in% nombres_comunes_grupos, ]

    message("\t * Summary of O ~ (Y,A,W_I,Delta):")
    log_dimensions("final W.init.concat", W)
    log_dimensions("final clinical", clinical)
    message(sprintf(
      "    *** final A dimension: %s - Summary (1: radio+TMZ, 0: radio+TMZ+beva): %s x %s",
      length(A), summary(as.factor(A))[[2]], summary(as.factor(A))[[1]]
    ))
    message("    *** Numero de pacientes por clusters:")
    message(capture.output(summary(grupos$cluster)))
    message(sprintf("    *** final Y dimension: %s - Summary:", length(Y)))
    message(capture.output(summary(Y)))
    message(sprintf(
      "    *** final delta dimension: %s - Summary (1: obsOutcome, 0: right-censored): %s x %s",
      length(delta), summary(as.factor(delta))[[2]], summary(as.factor(delta))[[1]]
    ))

    return(list(Y, A, delta, W, clinical, grupos))
  }
  else if (dimReduc == FALSE & outTMLE3 == FALSE) {
    message("W complete")
    all <- read.csv(file.path(pathFeatureLabel, "co-mod_best_results_last_run", "all_TCGA_LGG_GBM.csv"), header = TRUE, sep = ",")
    W <- all %>%
      dplyr::select(contains("mrna_"), contains("cnv_"))
    A <- all$A
    Y <- all$Y
    delta <- all$delta
    clinical <- all %>%
      dplyr::select(contains("clinical_")) %>%
      dplyr::mutate_if(is.character, as.factor)
    library(fastDummies)
    clinical <- data.matrix(fastDummies::dummy_cols(clinical, remove_first_dummy = TRUE, remove_selected_columns = TRUE))
    W <- data.matrix(W)

    message("#. Final results:")
    log_dimensions("final W.init.concat", W)
    log_dimensions("final clinical", clinical)
    message(sprintf(
      "    *** final A dimension: %s - Summary (0: RAD+TMZ, 1: RAD): %s x %s",
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
  else if (dimReduc == FALSE & outTMLE3 == TRUE) {
    message("O_data input structure for tmle3 -> O = (Y_{w/ NAs}, A, Delta, [W + clinical])")
    all <- read.csv(file.path(pathFeatureLabel, "co-mod_ccle_tcga_depurados", paste0("imputation_", index, ".csv")), header = TRUE, sep = ",")
    W <- all %>%
      dplyr::select(contains("mrna_"), contains("cnv_"))
    A <- all$A
    Y <- all$Y
    delta <- all$delta
    clinical <- all %>%
      dplyr::select(contains("clinical_")) %>%
      dplyr::mutate_if(is.character, as.factor)
    W <- cbind(clinical, W)
    if (Y_as_NaN == TRUE) {
      Y[delta == 0] <- NA
    }

    message("#. Final results:")
    log_dimensions("final (W + clinical)", W)
    log_dimensions("final clinical", clinical)
    message(sprintf(
      "    *** final A dimension: %s - Summary (0: RAD+TMZ, 1: RAD): %s x %s",
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
}

resumeData <- function(O) {
  message("\tSummary of O ~ (Y,A,W_I,Delta):")
  Y <- O[[1]]
  A <- O[[2]]
  Delta <- O[[3]]
  W <- O[[4]]
  clinical <- O[[5]]
  log_dimensions("final W.init.concat", W)
  log_dimensions("final clinical", clinical)
  message(sprintf(
    "    *** final A dimension: %s - Summary (0: RAD+TMZ, 1: RAD): %s x %s",
    length(A), summary(as.factor(A))[[2]], summary(as.factor(A))[[1]]
  ))
  message(sprintf("    *** final Y dimension: %s - Summary:", length(Y)))
  message(capture.output(summary(Y)))
  message(sprintf(
    "\n    *** final delta dimension: %s - Summary (1: obsOutcome, 0: right-censored): %s x %s",
    length(Delta), summary(as.factor(Delta))[[2]], summary(as.factor(Delta))[[1]]
  ))
}

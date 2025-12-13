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
library(TCGAbiolinks)

if (identical(character(0), dir("D:/pathFeatureLabelKernels"))) {
  central_path <- "E:"
} else {
  central_path <- "D:"
}

pathFeatureLabel <- paste0(central_path, "/pathFeatureLabelKernels")

read_project_profile <- function(path_feature_label, project, profile) {
  read.csv(
    file.path(
      path_feature_label,
      "co-mod_ccle_tcga_depurados",
      paste0("tcga_", project, "_", profile, "_depurado_tmle.csv")
    )
  )
}

clean_profile <- function(df) {
  df <- df[match(unique(df$X), df$X), ]
  row.names(df) <- df$X
  df$X <- NULL
  df
}

joint_profiles <- function(lgg_df, gbm_df) {
  common_cols <- intersect(colnames(lgg_df), colnames(gbm_df))
  combined <- rbind(lgg_df[, common_cols], gbm_df[, common_cols])
  clean_profile(combined)
}

variance_filter <- function(df, lower = 0.05, upper = 0.99) {
  df <- t(na.omit(t(df)))
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

prepare_common_patients <- function(W, clinical) {
  common <- Reduce(intersect, list(row.names(W), row.names(clinical)))
  list(W = W[common, , drop = FALSE], clinical = clinical[common, , drop = FALSE], patients = common)
}

prepare_final_outputs <- function(all_df, include_clinical = TRUE) {
  A <- all_df$A
  Y <- all_df$Y
  delta <- all_df$delta
  W <- all_df %>% dplyr::select(starts_with("mrna_"), starts_with("cnv_"))
  clinical <- all_df %>%
    dplyr::select(starts_with("clinical_")) %>%
    dplyr::mutate_if(is.character, factor)

  list(Y = Y, A = A, delta = delta, W = W, clinical = clinical)
}

give_O_GBM_LGG <- function(pathFeatureLabel, dimReduc, outTMLE3, index, Y_as_NaN) {

  if (!file.exists(file.path(pathFeatureLabel, "co-mod_best_results_last_run", "all_TCGA_LGG_GBM.csv"))) {
    message("... 1. Load cnv and mrna profiles + merge these profiles (rbind) ...")
    profiles <- c("mrna", "cnv")
    lgg_profiles <- map(profiles, ~ read_project_profile(pathFeatureLabel, "LGG", .x) %>% clean_profile())
    gbm_profiles <- map(profiles, ~ read_project_profile(pathFeatureLabel, "GBM", .x) %>% clean_profile())

    merged_profiles <- map2(lgg_profiles, gbm_profiles, joint_profiles) %>%
      map(~ variance_filter(.x, lower = 0.05, upper = 0.99))

    walk2(profiles, merged_profiles, ~ log_dimensions(.x, .y))

    message("... 2. Load gbm and lgg info and merge into glioma_info ...")
    lgg_info <- read.csv(file.path(pathFeatureLabel, "co-mod_best_results_last_run", "Y_A_clinical_LGG.csv"), row.names = 1)
    gbm_info <- read.csv(file.path(pathFeatureLabel, "co-mod_best_results_last_run", "Y_A_clinical_GBM.csv"), row.names = 1)

    glioma_info <- rbind(lgg_info, gbm_info) %>%
      dplyr::select(-PFI.time, -days_to_therapy_start, -treatment, -Telomere_Maintenance, -Grade) %>%
      dplyr::rename(delta = PFI, Y = new.PFI.time, A = treatment_0_1)
    colnames(glioma_info)[1:11] <- paste0("clinical_", colnames(glioma_info)[1:11])
    row.names(glioma_info) <- glioma_info$clinical_bcr_patient_barcode

    message("... 3. Column bind of cnv and mrna profiles for both projects (gbm and lgg) ...")
    W_ <- cbind(merged_profiles[[1]], merged_profiles[[2]])
    W_$X <- substr(row.names(W_), 1, 12)
    W_ <- W_[match(unique(W_$X), W_$X), ]
    row.names(W_) <- W_$X
    W_$X <- NULL

    message(" ... 4. Common patients for clinical and omic profiles for gbm and lgg ...")
    patients <- Reduce(intersect, list(row.names(W_), row.names(glioma_info)))
    write(patients, file.path(pathFeatureLabel, "co-mod_best_results_last_run", "barcodes_gbm_lgg_to_tmle.csv"))

    W <- W_[patients, ]
    W <- normalize_and_filter(W)
    message("---------------------------------")

    all <- cbind(glioma_info[patients, ], W)
    row.names(all) <- all$clinical_bcr_patient_barcode

    write.csv(all, file.path(pathFeatureLabel, "co-mod_best_results_last_run", "all_TCGA_LGG_GBM.csv"))
    write(all$bcr_patient_barcode, file = file.path(pathFeatureLabel, "co-mod_best_results_last_run", "barcodes_TCGA_LGG_GBM.csv"))

    clinicalData <- all %>%
      dplyr::select(-A, -Y, -delta, -clinical_TERT_promoter_status) %>%
      dplyr::mutate_if(is.character, as.factor)
    clinicalData[is.na(clinicalData)] <- NA

    for (i in 1:5) {
      clinicalData.imp <- read.csv(file.path(pathFeatureLabel, "co-mod_ccle_tcga_depurados", paste0("imputation_", i, ".csv")))
      write.csv(
        x = data.frame(A = all$A, Y = all$Y, delta = all$delta, clinicalData.imp),
        file = file.path(pathFeatureLabel, "co-mod_ccle_tcga_depurados", paste0("imputation_", i, ".csv"))
      )
    }

    outputs <- prepare_final_outputs(all)

    message("\tSummary of O ~ (Y,A,W_I,Delta):")
    log_dimensions("final W.init.concat", outputs$W)
    log_dimensions("final clinical", outputs$clinical)
    message(sprintf(
      "    *** final A dimension: %s - Summary (0: RAD+TMZ, 1: RAD): %s x %s",
      length(outputs$A), summary(as.factor(outputs$A))[[2]], summary(as.factor(outputs$A))[[1]]
    ))
    message(sprintf("    *** final Y dimension: %s - Summary:", length(outputs$Y)))
    message(capture.output(summary(outputs$Y)))
    message(sprintf(
      "    *** final delta dimension: %s - Summary (1: obsOutcome, 0: right-censored): %s x %s",
      length(outputs$delta), summary(as.factor(outputs$delta))[[2]], summary(as.factor(outputs$delta))[[1]]
    ))

    return(list(outputs$Y, outputs$A, outputs$delta, outputs$W, outputs$clinical))

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
    clinical <- all %>% dplyr::select(contains("clinical_"))

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
    W <- all %>% dplyr::select(contains("mrna_"), contains("cnv_"))
    A <- all$A
    Y <- all$Y
    delta <- all$delta
    clinical <- all %>%
      dplyr::select(contains("clinical_")) %>%
      dplyr::mutate_if(is.character, as.factor)
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
    W <- all %>% dplyr::select(contains("mrna_"), contains("cnv_"))
    A <- all$A
    Y <- all$Y
    delta <- all$delta
    clinical <- all %>% dplyr::select(contains("clinical_")) %>% dplyr::mutate_if(is.character, as.factor)
    W <- cbind(clinical, W)
    if (Y_as_NaN == TRUE) {
      Y[delta == 0 & Y > max(Y[delta == 1])] <- max(Y[delta == 1])
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

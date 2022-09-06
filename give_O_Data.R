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
# install.packages("tmle")
# library(tmle)
# library(caret)
library(dplyr)
library(purrr)
library(doParallel)
library(missForest)
# library(SuperLearner)
# source("C:/Users/da.salazarb/Google Drive/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/func_DimenReduc.R")
# 
# if (identical(character(0), dir("D:/pathFeatureLabelKernels"))) {
#   central_path <- "E"
# } else {
#   central_path <- "D"
# }

# pathFeatureLabel <- paste0(central_path, ":/pathFeatureLabelKernels")

give_O_data <- function(pathFeatureLabel) {
  ## --------------------------------- ##
  ## data TCGA
  print("1. loading profiles")
  W_protein <- read.csv(file = paste0(pathFeatureLabel,"/co-mod_ccle_tcga_depurados/tcga_protein_tmle.csv"), header = TRUE, sep = ",", row.names = 1)
  W_mirna <- read.csv(file = paste0(pathFeatureLabel,"/co-mod_ccle_tcga_depurados/tcga_mirna_tmle.csv"), header = TRUE, sep = ",", row.names = 1)
  W_cnv <- read.csv(file = paste0(pathFeatureLabel,"/co-mod_ccle_tcga_depurados/tcga_cnv_tmle.csv"), header = TRUE, sep = ",", row.names = 1)
  W_mrna <- read.csv(file = paste0(pathFeatureLabel,"/co-mod_ccle_tcga_depurados/tcga_mrna_tmle.csv"), header = TRUE, sep = ",", row.names = 1)
  
  # most_var <- apply(W_mrna, 2, var)
  # unicos_genes <- names(most_var[most_var > quantile(most_var, .5)])
  
  genes_DEA <- read.csv(paste0(central_path, ":/LGG/DEA_LGG.csv"), header = TRUE)
  unicos_genes <- intersect(as.vector(genes_DEA$X), colnames(W_mrna))
  W_mrna <- W_mrna[, unicos_genes]
  
  print(paste0("    * init W_protein dimensions: ", dim(W_protein)[1], " x ", dim(W_protein)[2]))
  print(paste0("    * init W_mirna dimensions: ", dim(W_mirna)[1], " x ", dim(W_mirna)[2]))
  print(paste0("    * init W_cnv dimensions: ", dim(W_cnv)[1], " x ", dim(W_cnv)[2]))
  print(paste0("    * init W_mrna dimensions: ", dim(W_mrna)[1], " x ", dim(W_mrna)[2]))
  
  ### SuperLearner as a early integration method
  W_ <- cbind(W_protein, W_mirna, W_cnv, W_mrna)
  print(paste0("    ** init W_ dimensions: ", dim(W_)[1], " x ", dim(W_)[2]))
  # W_ <- W_[,-as.numeric(which(apply(W_, 2, var) == 0))]
  print("---------------------------------")
  
  ## --------------------------------- ##
  print("2. loading clinical info")
  clinical <- read.csv(paste0(pathFeatureLabel, "/co-mod_best_results_last_run/clinical_data_TCGA.csv"), header = TRUE, row.names = 1)
  print(paste0("    init clinical dimensions: ", dim(clinical)[1], " x ", dim(clinical)[2]))
  print("---------------------------------")
  
  ## --------------------------------- ##
  ## outcome  + right-censored
  print("3. loading Y and delta info")
  Y <- read.csv(paste0(pathFeatureLabel,"/co-mod_bases_datos_importantes/survivalTCGA.csv"), header = TRUE, sep = ",")
  Y_delta <- Y %>% ## delta hace referencia a los right-censored obs
    dplyr::select(bcr_patient_barcode, type, PFI, PFI.time) %>% 
    dplyr::filter(type == "LGG")
  
  row.names(Y_delta) <- Y_delta$bcr_patient_barcode
  print(paste0("    * init Y_delta dimensions: ", dim(Y_delta)[1], " x ", dim(Y_delta)[2]))
  print("---------------------------------")
  
  ## --------------------------------- ##
  ## treatment -> archivo LGG_clinic_data.R -> 
  ## radiation+temozolomide (125), radiation (49), temozolomide (28)
  print("4. loading A info")
  A <- read.csv(paste0(pathFeatureLabel,"/co-mod_best_results_last_run/radiation_treatment_cycles_TCGA_LGG.csv"), header = TRUE, sep = ",")
  A <- A %>% 
    dplyr::filter(treat %in% c("radiation+temozolomide", "radiation")) %>% # , "temozolomide"
    dplyr::mutate(treat=ifelse(treat=="radiation+temozolomide", 1, 0)) 
  row.names(A) <- A$bcr_patient_barcode
  print(paste0("    * init A dimensions: ", dim(A)[1], " x ", dim(A)[2]))
  print("---------------------------------")
  
  ## --------------------------------- ##
  ## pacientes unicos
  print("5. common patients")
  pacientes_comunes <- purrr::reduce(list(row.names(W_), row.names(A), row.names(Y_delta), row.names(clinical)), intersect)
  
  A <- A[pacientes_comunes, ]
  A <- A %>% 
    dplyr::select(treat)
  
  A <- as.vector(A$treat)
  
  Y_delta <- Y_delta[pacientes_comunes, ]
  Y_delta <- Y_delta %>% 
    dplyr::select(PFI, PFI.time)
  Y <- as.numeric(Y_delta[2]$PFI.time)
  delta <- as.numeric(as.vector(Y_delta[1]$PFI))
  
  clinical <- clinical[pacientes_comunes, ]
  colnames(clinical) <- paste0("clinical_", colnames(clinical))
  W <- W_[pacientes_comunes, ]
  print("---------------------------------")
  
  ## --------------------------------- ##
  ## imputar faltantes en clinical data
  print("6. to impute missing values in clinical data")
  registerDoParallel(cores=3)
  W_ <- cbind(W, clinical)
  clinical <- missForest(W_, maxiter = 50, variablewise = TRUE, parallelize = "variables", mtry = 100)
  clinical <- clinical$ximp %>% 
    dplyr::select(contains("clinical_"))
  clinical <- fastDummies::dummy_cols(clinical, remove_selected_columns = TRUE, remove_first_dummy=TRUE)
  detach("package:doParallel", unload=TRUE)
  print("---------------------------------")
  
  ## --------------------------------- ##
  print("#. Final results:")
  print(paste0("    *** final W.init.concat dimension: ", dim(W_)[1], " x ", dim(W_)[2]))
  print(paste0("    *** final clinical dimension: ", dim(clinical)[1], " x ", dim(clinical)[2]))
  print(paste0("    *** final A dimension: ", length(A), " - Summary (1: radio+TMZ, 0: radio): ", summary(as.factor(A))[[2]], " x ", summary(as.factor(A)))[[1]] )
  print(paste0("    *** final Y dimension: ", length(Y), " - Summary: "))
  print(summary(Y))
  print(paste0("    *** final delta dimension: ", length(delta), " - Summary (1: obsOutcome, 0: right-censored): ", summary(as.factor(delta))[[2]], " x ", summary(as.factor(delta))[[1]]) )
  
  return(list(Y, A, delta, W, clinical))
}

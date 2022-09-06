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

# _____________________________________________ #
#### ........... Cargar archivos ........... ####
# _____________________________________________ #
if (identical(character(0), dir("D:/pathFeatureLabelKernels"))) {
  pathFeatureLabel <- "G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/data"
} else {
  pathFeatureLabel <- "D:"
}

# pathFeatureLabel <- central_path

give_O_GBM_LGG <- function(pathFeatureLabel,dimReduc,outTMLE3,index,Y_as_NaN) {
  
  if (identical(FALSE, file.exists(paste0(central_path, "/pathFeatureLabelKernels","/co-mod_best_results_last_run/all_TCGA_LGG.csv")))) {
    # ________________________________________________________________________________________________ #
    #### ........... Cargar archivos de GBM y LGG. Unirlos y procesarlos por varianzas. ........... ####
    # ________________________________________________________________________________________________ #
    cat("... 1. Load cnv and mrna profiles + merge these profiles (rbind) ... /n")
    ## cargar perfiles omicos depurados de gbm y lgg
    profiles <- c("mrna", "cnv") ## mirna de gbm solo tiene 5 variables, por eso no se usa / "protein",
    g <- purrr::map(profiles, function(x) read.csv(file = paste0(pathFeatureLabel,"/tcga_LGG_",x,"_depurado_tmle.csv")))
    g <- purrr::map(g, function(x) x[match(unique(x$X), x$X),])
    
    giveme_dim <- function(profiles, purrrList) {
      return(purrr::pmap(list(profiles, purrrList), function(x,y) cat(paste0("Dimension of ",x," is: ", dim(y)[1], " x ", dim(y)[2], "/n")) ))
    }
    
    giveme_dim(profiles, g)
    
    cleanProfiles <- function(xy) { # funcion para fusionar gbm y lgg
      row.names(xy) <- xy$X # dejar variable X como nombres de filas
      xy$X <- NULL ## eliminar X
      return(xy)
    }
    
    g <- purrr::map(g,cleanProfiles) # unir perfiles
    # g[[1]] <- t(na.omit(t(g[[1]]))) ## eliminar variables con NAN de genes
    g[[2]] <- t(na.omit(t(g[[2]]))) ## eliminar variables con NAN de CNV
    
    g.var <- purrr::map(g, function(x) apply(x,2,var)) # revisar varianza
    # purrr::map(h, hist)
    q.inf <- purrr::map(g.var, function(x) quantile(x,.1,na.rm = TRUE)) ## eliminar pocos variables
    g <- purrr::pmap(list(g,g.var,q.inf),function(x,y,z) x[,which(y > z)])
    # purrr::map(h, dim)
    
    g.var <- purrr::map(g, function(x) apply(x,2,var))
    q.sup <- purrr::map(g.var, function(x) quantile(x,.9,na.rm = TRUE)) ## eliminar muy variables
    g <- purrr::pmap(list(g,g.var,q.sup),function(x,y,z) x[,which(y < z)])
    
    for (i in 1:length(profiles)) {
      print(paste0("    * init profile ", profiles[i], " dimensions: ", dim(g[[i]])[1], " x ", dim(g[[i]])[2]))
    }
    
    # _____________________________________________________________________________________________________________________________________ #
    #### ........... Load information from GBM and LGG: barcodes of patients come from intersection of mrna and cnv profiles ........... ####
    # _____________________________________________________________________________________________________________________________________ #
    cat("... 2. Load gbm and lgg info and merge into glioma_info ... /n")
    lgg_info <- read.csv(paste0(pathFeatureLabel, "/Y_A_clinical_LGG.csv"),row.names = 1); head(lgg_info); dim(lgg_info)
    
    glioma_info <- lgg_info %>% 
      dplyr::select(-OS.time, -days_to_therapy_start, -treatment, -Telomere_Maintenance, -Grade) %>% 
      dplyr::rename(delta=OS, Y=new.OS.time, A=treatment_0_1)
    colnames(glioma_info)[1:11] <- paste0("clinical_", colnames(glioma_info)[1:11])
    row.names(glioma_info) <- glioma_info$clinical_bcr_patient_barcode
    head(glioma_info); dim(glioma_info)
    
    # __________________________________________________________________________________________________________________ #
    #### ........... Fusionar informacion adicional: clinica, Tratamiento (A), varRta (Y), censor (delta) ........... ####
    # __________________________________________________________________________________________________________________ #
    cat("... 3. Column bind of cnv and mrna profiles for both projects (gbm and lgg) ... /n")
    W_ <- cbind(g[[1]],g[[2]])
    W_$X <- sapply(row.names(W_), function(x) substr(x, 1, 12))
    W_ <- W_[match(unique(W_$X), W_$X),]
    row.names(W_) <- W_$X
    W_$X <- NULL
    W_[1:10,1:5]; dim(W_)
    
    # _____________________________________________________________________________________________ #
    #### ........... Common patients for clinical and omic profiles for gbm and lgg ........... ####
    # _____________________________________________________________________________________________ #
    cat(" ... 4. Common patients for clinical and omic profiles for gbm and lgg ... /n")
    cat(" ... 5. Save clean data in '/co-mod_best_results_last_run/all_TCGA_LGG_GBM.csv' ... /n")
    pacientes_comunes <- purrr::reduce(list(row.names(W_), row.names(glioma_info)), intersect)
    write(pacientes_comunes,paste0(pathFeatureLabel,"/barcodes_gbm_lgg_to_tmle.csv"))
    
    W <- W_[pacientes_comunes, ]
    W <- apply(W, 2, function(x){(x-min(x))/(max(x)-min(x))})
    W <- data.frame(W)
    ceros <- apply(W, 2, function(x) sum(x == 0))
    W <- W[, which(ceros < round(dim(W)[1] * .40)) ]
    cat("--------------------------------- /n")
    
    all <- cbind(glioma_info, W)
    row.names(all) <- all$clinical_bcr_patient_barcode
    
    all[1:10,1:15]; dim(all)
    
    write(all$bcr_patient_barcode, file = paste0(pathFeatureLabel,"/barcodes_TCGA_LGG.csv"))
    
    # ________________________________________ #
    #### ........... Imputation ........... ####
    # ________________________________________ #
    # malaCard <- read.csv(file = paste0(pathFeatureLabel,"/co-mod_bases_datos_importantes/GeneCards-SearchResults_glioma.csv"))
    # malaCard <- paste("mrna_", malaCard$Gene.Symbol, sep = "")
    # all <- read.csv(paste0(pathFeatureLabel,"/co-mod_best_results_last_run/all_TCGA_LGG_GBM.csv"))
    # all[1:10,1:15]; dim(all)
    
    clinicalData <- all %>% 
      dplyr::select(-A,-Y,-delta, -clinical_TERT_promoter_status, -clinical_bcr_patient_barcode) %>% 
      # dplyr::select(.,starts_with("clinical_"), one_of(malaCard)) %>% 
      dplyr::mutate_if(is.character, as.factor)
    clinicalData[1:5,1:5]; str(clinicalData[1:5,1:5]); dim(clinicalData)
    
    # library(VIM)
    # library(mice)
    # VIM::aggr(clinicalData %>% dplyr::select(.,starts_with("clinical")), col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(clinicalData), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
    # VIM::matrixplot(clinicalData %>% dplyr::select(.,starts_with("clinical")), sortby = c("clinical_Karnofsky_Performance_Score"), cex.axis=.7)
    
    ## imputation using missForest
    library(missForest)
    library(doParallel)
    registerDoParallel(cores=7)
    
    clinicalData.imp <- missForest(clinicalData[all$Y!=0,], parallelize="variables", ntree = 200, maxnodes = 25, mtry = 2000, verbose = TRUE)
    
    ## Seleccionar X variables de mrna relacionadas con Y aprate de las clinicas y del cnv
    for (i in c(50, 100, 500, 1000)) {
      mrna.variables <- colnames(all %>% dplyr::filter(Y != 0) %>% dplyr::select(., starts_with("mrna_")))[screen.randomForest(Y = ifelse(all$delta[all$Y != 0] != 0, 1, 0), 
                                                                                                     X = all %>% dplyr::filter(Y != 0) %>% dplyr::select(., starts_with("mrna_")), 
                                                                                                     family = binomial(), 
                                                                                                     nVar = i)]
      all_no.mrna <- colnames(clinicalData.imp$ximp %>% dplyr::select(., -starts_with("mrna_")))
      
      all_ <- clinicalData.imp$ximp %>% 
        dplyr::select(dplyr::one_of(all_no.mrna),dplyr::one_of(mrna.variables)) %>% 
        dplyr::bind_cols(all %>% dplyr::filter(Y != 0) %>% dplyr::select(Y, A, delta)) %>% 
        dplyr::relocate(Y,A,delta)
      print(dim(all_))
      
      write.csv(x = all_, file = paste0(pathFeatureLabel, "/imputation_", i, ".csv"))
      # write.csv(x = data.frame(A=all$A, Y=all$Y, delta=all$delta, clinicalData.imp$ximp), file = paste0(pathFeatureLabel, "/imputation_", i, ".csv"))
    }
    
    
    # for (i in 1:5) {
    #   clinicalData.imp <- missForest(clinicalData, parallelize="variables", ntree = 200, maxnodes = 25, mtry = 2000, verbose = TRUE)
    #   write.csv(x = data.frame(A=all$A, Y=all$Y, delta=all$delta, clinicalData.imp$ximp), file = paste0(pathFeatureLabel, "/imputation_", i, ".csv"))
    # }
    
    # for (i in 1:5) {
    #   clinicalData.imp <- read.csv(file = paste0(pathFeatureLabel, "/imputation_", i, ".csv"))
    #   write.csv(x = data.frame(A=all$A, Y=all$Y, delta=all$delta, clinicalData.imp), file = paste0(pathFeatureLabel, "/imputation_", i, ".csv"))
    # }
    
    
    # mean(c(0.2016183, 0.1923476, 0.1949457, 0.1927466 , 0.1997466,  0.1972337 ,  0.2025307 , 0.1885439, 0.1951557 , 0.2016287 , 
    # 0.1918307 ,  0.1958268 , 0.2013338 , 0.1924255 , 0.1981284, 0.2034309 , 0.2044298,  0.2005327 , 0.2028219, 0.1921532, 
    # 0.1975249 , 0.1962466 , 0.1989299 ,  0.196337 , 0.1961408 , 0.1921587,  0.1961409 ))
    # ## Imputation with miselect and mice
    # library(mice)
    # library(miselect)
    # 
    # #creates a list of predictors of missing data with a mininum correlation of 0.1
    # #and at least 50% of useful data
    # predictor.selection <- quickpred(clinicalData, mincor=0.1, minpuc=0.5,method='pearson')
    # mids <- mice(clinicalData, m=5, predictorMatrix = predictor.selection)
    # 
    # # Generate list of completed data.frames
    # dfs <- lapply(1:5, function(i) complete(mids, action = i))
    
    
    
    # ____________________________________________ #
    #### ........... print results ........... ####
    # ____________________________________________ #
    A <- all$A
    Y <- all$Y
    delta <- all$delta
    W <- all %>% 
      dplyr::select(.,starts_with("mrna_"),starts_with("cnv_"))
    # ____________________________________________________________ #
    #### ........... clinical data for gbm and lgg ........... ####
    # ____________________________________________________________ #
    clinical <- all %>% dplyr::na_if("NA") %>% ## old clinicalData
      dplyr::select(.,starts_with("clinical_")) %>% 
      # dplyr::mutate(A=factor(A)) %>%  #, age_at_index=factor(age_at_index)
      dplyr::mutate_if(is.character, factor)
    # dplyr::select(Telomere_Maintenance, MGMT_promoter_status, ATRX_status, IDH_codel_subtype, A)
    # head(clinical)
    
    ## --------------------------------- ##
    cat("/tSummary of O ~ (Y,A,W_I,Delta):/n")
    cat(paste0("    *** final W.init.concat dimension: ", dim(W)[1], " x ", dim(W)[2], " /n"))
    cat(paste0("    *** final clinical dimension: ", dim(clinical)[1], " x ", dim(clinical)[2], " /n"))
    cat(paste0("    *** final A dimension: ", length(A), " - Summary (0: RAD+TMZ, 1: RAD): ", summary(as.factor(A))[[2]], " x ", summary(as.factor(A)))[[1]], " /n" )
    cat(paste0("    *** final Y dimension: ", length(Y), " - Summary: /n"))
    cat(summary(Y))
    cat(" /n")
    cat(paste0("    *** final delta dimension: ", length(delta), " - Summary (1: obsOutcome, 0: right-censored): ", summary(as.factor(delta))[[2]], " x ", summary(as.factor(delta))[[1]]) , " /n")
    
    return(list(Y, A, delta, W, clinical))
    
  } else if (dimReduc == TRUE & outTMLE3 == FALSE) {
    cat("/t * Loading W with reduction of dimensionality - NMF (func_NMF.py) - /n")
    # cat("/t 2.  is used for dimensionality reduction ... try other methods! /n")
    all <- read.csv(paste0(pathFeatureLabel,"/co-mod_best_results_last_run/Y_A_delta_clinical_TCGA_LGG_GBM.csv"), header = TRUE, sep = ",")
    ### los siguientes archivos provienen de func_NMF.py
    mrna.W <- read.csv(paste0(pathFeatureLabel,"/co-mod_best_results_last_run/NMF_W_mrna_LGG.csv"), header = FALSE, sep = ",")
    colnames(mrna.W) <- paste0("mrna_",colnames(mrna.W))
    protein.W <- read.csv(paste0(pathFeatureLabel,"/co-mod_best_results_last_run/NMF_W_protein_LGG.csv"), header = FALSE, sep = ",")
    colnames(protein.W) <- paste0("protein_",colnames(protein.W))
    cnv.W <- read.csv(paste0(pathFeatureLabel,"/co-mod_best_results_last_run/NMF_W_cnv_LGG.csv"), header = FALSE, sep = ",")
    colnames(cnv.W) <- paste0("cnv_",colnames(cnv.W))
    W <- cbind(mrna.W,protein.W,cnv.W)
    A <- all$A
    Y <- all$Y
    delta <- all$delta
    clinical <- all %>% 
      dplyr::select(contains("clinical_"))
    
    ## --------------------------------- ##
    ## Dejando solo los pacientes con grupos
    nombres_comunes_grupos <- intersect(clinical$clinical_submitter_id, grupos$short_index)
    grupos <- grupos[grupos$short_index %in% nombres_comunes_grupos,]
    
    ## --------------------------------- ##
    cat("/t * Summary of O ~ (Y,A,W_I,Delta):/n")
    cat(paste0("    *** final W.init.concat dimension: ", dim(W)[1], " x ", dim(W)[2], "/n"))
    cat(paste0("    *** final clinical dimension: ", dim(clinical)[1], " x ", dim(clinical)[2], "/n"))
    cat(paste0("    *** final A dimension: ", length(A), " - Summary (1: radio+TMZ, 0: radio+TMZ+beva): ", summary(as.factor(A))[[2]], " x ", summary(as.factor(A))[[1]], "/n") )
    cat("    *** Numero de pacientes por clusters: /n")
    print(summary(grupos$cluster))
    cat(paste0("    *** final Y dimension: ", length(Y), " - Summary: /n"))
    cat(paste0("    *** final clinical dimension: ", dim(clinical)[1], " x ", dim(clinical)[2], "/n"))
    print(summary(Y))
    cat(paste0("    *** final delta dimension: ", length(delta), " - Summary (1: obsOutcome, 0: right-censored): ", summary(as.factor(delta))[[2]], " x ", summary(as.factor(delta))[[1]], "/n") )
    
    return(list(Y, A, delta, W, clinical, grupos))
  } 
  else if (dimReduc == FALSE & outTMLE3 == FALSE) {
    cat("W complete /n")
    # all.1 <- read.csv(paste0(pathFeatureLabel,"/co-mod_best_results_last_run/all_TCGA_LGG_GBM.csv"), header = TRUE, sep = ",")
    # all <- read.csv(paste0(pathFeatureLabel,"/co-mod_ccle_tcga_depurados/imputation_1.csv"), header = TRUE, sep = ",")
    # all.2 <- data.frame(Y=all.1$Y, A=all.1$A, delta=all.1$delta, all)
    # all <- all.2
    W <- all %>% 
      dplyr::select(contains("mrna_"), contains("cnv_"))
    A <- all$A
    Y <- all$Y
    delta <- all$delta
    clinical <- all %>% 
      dplyr::select(contains("clinical_")) %>% 
      dplyr::mutate_if(is.character,as.factor)
    library(fastDummies)
    clinical <- data.matrix(fastDummies::dummy_cols(clinical, remove_first_dummy = TRUE,remove_selected_columns = TRUE))
    W <- data.matrix(W)
    
    ## --------------------------------- ##
    cat("#. Final results: /n")
    cat(paste0("    *** final W.init.concat dimension: ", dim(W)[1], " x ", dim(W)[2], " /n"))
    cat(paste0("    *** final clinical dimension: ", dim(clinical)[1], " x ", dim(clinical)[2], " /n"))
    cat(paste0("    *** final A dimension: ", length(A), " - Summary (0: RAD+TMZ, 1: RAD): ", summary(as.factor(A))[[2]], " x ", summary(as.factor(A))[[1]], " /n") )
    cat(paste0("    *** final Y dimension: ", length(Y), " - Summary:  /n"))
    print(summary(Y))
    cat(paste0("    *** final delta dimension: ", length(delta), " - Summary (1: obsOutcome, 0: right-censored): ", summary(as.factor(delta))[[2]], " x ", summary(as.factor(delta))[[1]]) , " /n")
    
    return(list(Y, A, delta, W, clinical))
  }
  else if (dimReduc == FALSE & outTMLE3 == TRUE) {
    cat("O_data input structure for tmle3 -> O = (Y_{w/ NAs}, A, Delta, [W + clinical]) /n")
    # all.1 <- read.csv(paste0(pathFeatureLabel,"/co-mod_best_results_last_run/all_TCGA_LGG_GBM.csv"), header = TRUE, sep = ",")
    all <- read.csv(paste0(pathFeatureLabel,"/co-mod_ccle_tcga_depurados/imputation_",index,".csv"), header = TRUE, sep = ",")
    # all.2 <- data.frame(Y=all.1$Y, A=all.1$A, delta=all.1$delta, all)
    # all <- all.2
    W <- all %>% 
      dplyr::select(contains("mrna_"), contains("cnv_"))
    A <- all$A
    Y <- all$Y
    delta <- all$delta
    # Y[delta == 0][Y[delta == 0] > max(Y[delta == 1])] <- max(Y[delta == 1])
    clinical <- all %>% 
      dplyr::select(contains("clinical_")) %>% 
      dplyr::mutate_if(is.character,as.factor)
    # library(fastDummies)
    # clinical <- data.matrix(fastDummies::dummy_cols(clinical, remove_first_dummy = TRUE,remove_selected_columns = TRUE))
    # W <- data.matrix(W)
    W <- cbind(clinical,W)#; dim(W); class(W)
    if (Y_as_NaN == TRUE) {
      Y[which(delta == 0)] <- NA
    }
    
    ## --------------------------------- ##
    cat("#. Final results: /n")
    cat(paste0("    *** final (W + clinical) dimension: ", dim(W)[1], " x ", dim(W)[2], " /n"))
    cat(paste0("    *** final clinical dimension: ", dim(clinical)[1], " x ", dim(clinical)[2], " /n"))
    cat(paste0("    *** final A dimension: ", length(A), " - Summary (0: RAD+TMZ, 1: RAD): ", summary(as.factor(A))[[2]], " x ", summary(as.factor(A))[[1]], " /n") )
    cat(paste0("    *** final Y dimension: ", length(Y), " - Summary:  /n"))
    print(summary(Y))
    cat(paste0("    *** final delta dimension: ", length(delta), " - Summary (1: obsOutcome, 0: right-censored): ", summary(as.factor(delta))[[2]], " x ", summary(as.factor(delta))[[1]]) , " /n")
    
    return(list(Y, A, delta, W, clinical))
  }
}

resumeData <- function(O) {
  cat("/tSummary of O ~ (Y,A,W_I,Delta):/n")
  Y <- O[[1]]; A <- O[[2]]; Delta <- O[[3]]; W <- O[[4]]; clinical <- O[[5]];
  cat(paste0("    *** final W.init.concat dimension: ", dim(W)[1], " x ", dim(W)[2], "/n"))
  cat(paste0("    *** final clinical dimension: ", dim(clinical)[1], " x ", dim(clinical)[2], "/n"))
  cat(paste0("    *** final A dimension: ", length(A), " - Summary (0: RAD+TMZ, 1: RAD): ", summary(as.factor(A))[[2]], " x ", summary(as.factor(A))[[1]] , "/n"))
  cat(paste0("    *** final Y dimension: ", length(Y), " - Summary: /n"))
  print(summary(Y))
  cat(paste0("/n    *** final delta dimension: ", length(Delta), " - Summary (1: obsOutcome, 0: right-censored): ", summary(as.factor(Delta))[[2]], " x ", summary(as.factor(Delta))[[1]], "/n") )  
}

#### para revisar varianzas de tiempos de nuevo eventos cuando se incluye muerte por tumor vs cuando se excluye muerte por tumor
# surv <- read.csv(paste0(pathFeatureLabel,"/co-mod_bases_datos_importantes/survivalTCGA.csv"), header = TRUE, sep = ",") %>% ## delta hace referencia a los right-censored obs
#   dplyr::select(bcr_patient_barcode, type, PFI, PFI.time) %>% 
#   # dplyr::rename(PFI=PFI, PFI.time=PFI.time) %>% 
#   dplyr::filter(type %in% c("LGG")) %>% 
#   dplyr::select(-type) %>% 
#   dplyr::mutate(PFI.time = as.numeric(PFI.time))
# head(surv)
# 
# clinical.follow.up_temp <- clinical.follow.up %>% 
#   dplyr::select(bcr_patient_barcode, days_to_death)
# 
# 
# datos <- clinical.follow.up_temp %>% 
#   dplyr::inner_join(surv, by="bcr_patient_barcode")
# 
# sd(datos[!is.na(datos$days_to_death),])
# sd(datos[!is.na(datos$days_to_death) & datos$days_to_death != datos$PFI.time,]$PFI.time)

# Varianza de variable respuesta (incluyendo muerte por tumor): 512060.2 -> 715.5838
# Varianza de variable respuesta (excluyendo muerte por tumor): 527293.6 -> 726.1498


#### PScore
# ## https://datascienceplus.com/imputing-missing-data-with-r-mice-package/
# library(VIM)
# library(mice)
# VIM::aggr(clinicalData %>% dplyr::select(.,starts_with("clinical")), col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(clinicalData), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
# VIM::matrixplot(clinicalData %>% dplyr::select(.,starts_with("clinical")), sortby = c("clinical_TERT_promoter_status"), cex.axis=.7)
# 
# ## MCAR vs MAR
# library(gtsummary)
# trial <- clinicalData %>% select(.,starts_with("clinical_"))
# 
# trial %>% tbl_summary(by=clinical_primary_diagnosis,
#                       statistic = list(all_continuous() ~ "{mean} ({sd})"),
#                       # label = list(gender ~ "Gender",
#                       #   IDH_codel_subtype ~ "IDH codel subtype",
#                       #   MGMT_promoter_status ~ "MGMT promoter status",
#                       #   ATRX_status ~ "ATRX status",
#                       #   OS ~ "Censoring",
#                       #   OS.time ~ "Overall survival"
#                       #              ),
#                       missing_text = "(Missing)") %>%
#   bold_labels() %>%
#   add_overall() %>%
#   add_n() %>%
#   add_p() %>% add_q()
# 
# library(finalfit)
# clinical.colnames <- colnames(clinicalData %>% 
#                           dplyr::select(., dplyr::starts_with("clinical_")))
# clinical.colnames <- gsub("clinical_","",clinical.colnames)
# clinicalData %>% summary_factorlist(dependent = clinical.colnames[1],explanatory = clinical.colnames[2:length(clinical.colnames)],p = TRUE, na_include = TRUE)
# library(GGally)
# a <- clinicalData
# colnames(a) <- gsub("clinical_", "", colnames(a))
# a %>% finalfit::missing_pairs(dependent = clinical.colnames[1],explanatory = clinical.colnames[2:length(clinical.colnames)],position = "fill")
# clinicalData %>% finalfit::missing_compare(dependent = clinical.colnames[1],explanatory = clinical.colnames[4]) %>% 
#   knitr::kable(row.names=FALSE, align = c("l", "l", "r", "r", "r")) # Omit when you run
# 
# ## http://www.practicalpropensityscore.com/uploads/9/4/5/3/94532409/chapter5_part1__propensity_score_estimation.r
# # https://us.sagepub.com/en-us/nam/practical-propensity-score-methods-using-r/book241054#preview
# # https://osf.io/nygb5/?view_only=a955dc0fe51a4ecfa9ad1d22b0791086
# #examine the number of missing cases
# # missing.indicator <- data.frame(is.na(clinicalData))
# # propMissing <- apply(missing.indicator,2,mean)
# # 
# # #create dummy missing value indicators
# # names(missing.indicator)[propMissing>0] <- paste(names(clinicalData)[propMissing>0],"NA",sep="")
# # #convert dummy missing indicators from logical to numeric variables
# # for (var in 1:ncol(missing.indicator)) {
# #   missing.indicator[,var] <- as.numeric(missing.indicator[,var]) }
# # 
# # #merge covariate names with missing indicator names
# # ELS.data <- cbind(clinicalData,missing.indicator[,propMissing>0])
# # 
# # #show percentage missing
# # print(round(propMissing,3))
# 
# #======================
# #impute separately treated and untreated groups
# 
# #create storage
# long.imputation <- c()
# 
# #loop through the treatment groups
# for (group in 0:1) {
#   
#   
#   #creates a list of predictors of missing data with a mininum correlation of 0.1
#   #and at least 50% of useful data
#   predictor.selection <- quickpred(subset(clinicalData,A==group), mincor=0.1, minpuc=0.5,method='pearson',
#                                    exclude=c("STU_ID"))
#   
#   
#   #impute variables by from least missing to most missing
#   #Using multiple imputation by chained equations
#   #with predictive mean matching as the univariate imputation method
#   imputation <- mice(subset(clinicalData,A==group), m=10, method="pmm", visitSequence="monotone",
#                      predictorMatrix = predictor.selection)
#   
#   
#   #extract stacked data files
#   long.imputation <- rbind(long.imputation,complete(imputation, action="long"))
#   
# } #finish loop
# # 
# # #extract a single imputation dataset
# # imputation1 <- subset(long.imputation, subset=.imp==1)
# # 
# # #create a list of all imputed datasets
# library(mitools) #load package to combine imputed datasets
# # #this object was specifically designed to be analyzed with the survey package
# allImputations <- imputationList(list(
#   subset(long.imputation, subset=.imp==1),
#   subset(long.imputation, subset=.imp==2),
#   subset(long.imputation, subset=.imp==3),
#   subset(long.imputation, subset=.imp==4),
#   subset(long.imputation, subset=.imp==5),
#   subset(long.imputation, subset=.imp==6),
#   subset(long.imputation, subset=.imp==7),
#   subset(long.imputation, subset=.imp==8),
#   subset(long.imputation, subset=.imp==9),
#   subset(long.imputation, subset=.imp==10)))
# 
# ##clinical names
# # clinical.names <- colnames(glioma_info %>% dplyr::select(-clinical_bcr_patient_barcode, -delta, -A, -Y))
# # paste(clinical.names, collapse = " + ")
# 
# # formula.ps <- formula(paste("A~", paste(clinical.names, collapse = "+"),sep = ""))
# 
# # modelFit1 <- with(data=imputation1, expr=glm(formula.ps, family = binomial))
# 
# psModelAll <- with(data=allImputations, expr=glm(A ~ clinical_primary_diagnosis + clinical_age_at_index + clinical_gender + 
#                                                   clinical_IDH_codel_subtype + clinical_MGMT_promoter_status + clinical_ATRX_status + 
#                                                   clinical_Chr7_gain_and_Chr10_loss + clinical_TERT_promoter_status + 
#                                                   clinical_Karnofsky_Performance_Score + clinical_Mutation_Count,
#                                                  family = binomial))
# summary(pool(psModelAll))
# pScoresAll <- sapply(psModelAll, fitted)
# pScoresMean <- apply(pScoresAll,1,mean)
# allImputations <- update(allImputations, pScores = pScoresMean)
# 
# ### one single model
# # imputation1 <- subset(long.imputation, subset=.imp==1)
# # ps.colnames <- colnames(imputation1 %>% 
# #                           dplyr::select(-.imp, -.id, -A))
# 
# # formula(paste("A~", paste(ps.colnames, collapse = "+"), sep = ""))
# # psModel1 <- glm(formula(paste("A~", paste(ps.colnames, collapse = "+"), sep = "")), family = binomial, data = imputation1)
# 
# psModel1 <- glm(A ~ clinical_primary_diagnosis + clinical_age_at_index + clinical_gender + 
#                   clinical_IDH_codel_subtype + clinical_MGMT_promoter_status + clinical_ATRX_status + 
#                   clinical_Chr7_gain_and_Chr10_loss + clinical_TERT_promoter_status + 
#                   clinical_Karnofsky_Performance_Score + clinical_Mutation_Count, 
#                 family = binomial, data = imputation1)
# pScores <-  fitted(psModel1)
# imputation1$pScores <-  pScores
# 
# ggplot() +
#   # Top
#   geom_density(data=subset(imputation1, A==0), aes(x = pScores, y = ..density..), fill="#69b3a2" ) +
#   # geom_label( aes(x=0.4, y=250, label="RAD+TMZ"), color="#69b3a2") +
#   # Bottom
#   geom_density(data=subset(imputation1, A==1), aes(x = pScores, y = -..density..), fill= "#404080") + scale_y_continuous(breaks=round(seq(min(-2.5), max(3.5), by=0.2),2))
#   # geom_label( aes(x=0.5, y=-100, label="RAD"), color="#404080")
# 
# ## with trees
# library(party)
# psModelAll <- with(data=allImputations, expr=ctree(A ~ clinical_primary_diagnosis + clinical_age_at_index + clinical_gender + 
#                                                    clinical_IDH_codel_subtype + clinical_MGMT_promoter_status + clinical_ATRX_status + 
#                                                    clinical_Chr7_gain_and_Chr10_loss + clinical_TERT_promoter_status + 
#                                                    clinical_Karnofsky_Performance_Score + clinical_Mutation_Count + 
#                                                    clinical_Grade))
# 
# summary(pool(psModel1))

### ---- ###

# ### Propensity Score ###
# library(tableone)
# colnames(all)
# baselinevars <- c("primary_diagnosis", "age_at_index", "gender", "IDH_codel_subtype", "MGMT_promoter_status", "ATRX_status", "Telomere_Maintenance","Chr7_gain_and_Chr10_loss")
# all <- all %>% dplyr::mutate_if(is.character, as.factor)
# tabla <- CreateTableOne(vars=baselinevars, 
#                data=all, strata = "A", 
#                includeNA = TRUE, 
#                test = FALSE, smd=TRUE)
# 
# print(tabla, smd = TRUE)
# 
# ## Step 1: PS estimation
# ps.formula <- as.formula(paste("I(A == 1)", "~", 
#                                paste(baselinevars, collapse = "+")))
# ps.formula
# 
# ps.fit <- glm(ps.formula, family = "binomial", 
#               data=all,na.action = na.omit)
# 
# all$ps <- predict(ps.fit, newdata = all, type="response",na.action = na.omit)

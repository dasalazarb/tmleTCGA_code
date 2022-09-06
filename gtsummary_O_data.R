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
library(gt)
library(gtsummary)
O <- read.csv("G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/os_lgg_data_var50/imputation_1.csv"); dim(O)

O$delta <- ifelse(O$delta == 0 & O$Y < 730, 0, 1); O$Y <- ifelse(O$Y < 730 & O$delta == 1, 1, 0); 

as_kable_extra(O %>% 
  select(Y,A,delta, starts_with("clinical_")) %>% 
  mutate(A = ifelse(A==1, "RAD+TMZ", "RAD"), 
         clinical_gender = ifelse(clinical_gender == "female", "Female", "Male"), 
         clinical_ATRX_status = ifelse(clinical_ATRX_status == "WT", "Wild-type", "Mutant")) %>% 
  tbl_summary(by = A, 
              label = list(Y ~ "Overall survival", 
                           delta ~ "Censoring", 
                           clinical_primary_diagnosis ~ "Primary Diagnosis", 
                           clinical_age_at_index ~ "Age at index", 
                           clinical_gender ~ "Gender", 
                           clinical_IDH_codel_subtype ~ "IDH/codelation subtype", 
                           clinical_MGMT_promoter_status ~ "MGMT promoter status", 
                           clinical_ATRX_status ~ "ATRX status", 
                           clinical_Chr7_gain_and_Chr10_loss ~ "Chr7 gain and Chr10 loss", 
                           clinical_Karnofsky_Performance_Score	~ "Karnofsky performance score", 
                           clinical_Mutation_Count ~ "Mutation count")) %>% 
  add_p() %>% 
  modify_header(label ~ "**Variable**") %>% 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Treatment Received**") %>%
  bold_labels(), format = "latex") #%>% huxtable::to_latex()#gt::as_latex()

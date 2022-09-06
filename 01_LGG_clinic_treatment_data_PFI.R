# https://bioc.ism.ac.jp/packages/3.2/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html#tcgaquery_clinic-tcgaquery_clinicfilt-working-with-clinical-data
# https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html
# Avastin FDA Approval History: https://www.drugs.com/history/avastin.html
library(TCGAbiolinks)
library(dplyr)
library(gtsummary)
library(DT)

if (identical(character(0), dir("D:/pathFeatureLabelKernels"))) {
  central_path <- "E:"
} else {
  central_path <- "D:"
}

pathFeatureLabel <- paste0(central_path, "/pathFeatureLabelKernels")

### ---------------------- ###
#### ..... LGG info ..... ####
### ---------------------- ###
barcode_LGG <- read.table(paste0(pathFeatureLabel, "/co-mod_best_results_last_run/barcodes_TCGA_LGG.txt"))
barcode_LGG <- sapply(barcode_LGG$V1, function(x) substr(x, 1, 12))
head(barcode_LGG); length(barcode_LGG)

query <- GDCquery(project = "TCGA-LGG",
                  data.category = "Clinical",
                  file.type = "xml",
                  barcode = barcode_LGG)
## info from TCGABiolinks ##
GDCdownload(query)
clinical.patient <- GDCprepare_clinic(query,"patient"); names(clinical.patient); head(clinical.patient)
clinical.drug <- GDCprepare_clinic(query,"drug"); head(clinical.drug); dim(clinical.drug)
clinical.radiation <- GDCprepare_clinic(query,"radiation"); names(clinical.radiation); head(clinical.radiation)
clinical.follow.up <- GDCprepare_clinic(query,"follow_up"); names(clinical.follow.up); head(clinical.follow.up)

### ----------------------------------- ###
#### ..... Cleaning data for LGG ..... ####
### ----------------------------------- ###
### ... 1. drugs ... ###
drug_correction <- read.csv(file = paste0(pathFeatureLabel, "/co-mod_bases_datos_importantes/DrugCorrection1.csv")); head(drug_correction)
drug_correction <- drug_correction %>% 
  dplyr::mutate(OldName=gsub(" ", "", tolower(OldName)), Correction=gsub(" ", "", tolower(Correction)))

drug_lgg <- clinical.drug %>% 
  dplyr::mutate(drug_name=gsub(" ", "", tolower(drug_name))) %>% 
  dplyr::mutate(drug_name=unlist(lapply(drug_name, function(x) ifelse(x %in% drug_correction$OldName, unique(drug_correction[drug_correction$OldName %in% x,]$Correction), x)) ) ) %>% 
  dplyr::filter(!grepl("nos", drug_name)) %>% 
  dplyr::filter(!grepl("^.{0}$", drug_name)) %>% 
  dplyr::rename(pharmaceutical_therapy_drug_name=drug_name)
head(drug_lgg);dim(drug_lgg); sort(unique(drug_lgg$pharmaceutical_therapy_drug_name))

drug_lgg <- as.data.frame(drug_lgg %>% 
                            dplyr::group_by(bcr_patient_barcode, pharmaceutical_therapy_drug_name) %>% 
                            dplyr::filter(days_to_drug_therapy_start == min(days_to_drug_therapy_start)))
head(drug_lgg);dim(drug_lgg)

write.csv(drug_lgg, file=paste0(pathFeatureLabel, "/co-mod_best_results_last_run/drug_treatment_TCGA_LGG.csv"), row.names = FALSE)

### ... 2. Radiation ... ###
radiation_lgg <- clinical.radiation %>% 
  dplyr::select(bcr_patient_barcode, days_to_radiation_therapy_start, days_to_radiation_therapy_end, radiation_type, regimen_indication) %>% 
  dplyr::mutate(radiation_type=tolower(radiation_type)) %>% 
  dplyr::select(-regimen_indication)

radiation_lgg$radiation_type[radiation_lgg$radiation_type == "external beam"] <- "radiation"
radiation_lgg$radiation_type[radiation_lgg$radiation_type == ""] <- "radiation"
radiation_lgg$radiation_type[radiation_lgg$radiation_type == "radioisotope"] <- "radiation"
radiation_lgg$radiation_type[radiation_lgg$radiation_type == "external"] <- "radiation"
radiation_lgg$radiation_type[radiation_lgg$radiation_type == "other"] <- "radiation"

radiation_lgg <- unique(radiation_lgg)
colnames(radiation_lgg)[colnames(radiation_lgg) %in% "radiation_type"] <- "treat"

head(radiation_lgg);dim(radiation_lgg)

## ... 3. treatment: RAD + DRUGs ... ##
treatment_lgg <- rbind(drug_lgg %>% 
                         dplyr::select(bcr_patient_barcode, days_to_drug_therapy_start, pharmaceutical_therapy_drug_name, -days_to_drug_therapy_end) %>% 
                         dplyr::rename(days_to_therapy_start=days_to_drug_therapy_start, treat=pharmaceutical_therapy_drug_name), 
                       radiation_lgg %>% 
                         dplyr::select(-days_to_radiation_therapy_end) %>% 
                         dplyr::rename(days_to_therapy_start=days_to_radiation_therapy_start)) %>% 
  dplyr::inner_join(read.csv(paste0(pathFeatureLabel,"/co-mod_bases_datos_importantes/survivalTCGA.csv"), header = TRUE, sep = ",") %>% ## delta hace referencia a los right-censored obs
                      # dplyr::select(bcr_patient_barcode, type, OS, OS.time) %>% 
                      dplyr::select(bcr_patient_barcode, type, PFI, PFI.time) %>% 
                      dplyr::rename(PFI=PFI, PFI.time=PFI.time) %>% 
                      dplyr::filter(type %in% c("LGG")) %>% 
                      dplyr::select(-type) %>% 
                      dplyr::mutate(OS.time = as.numeric(PFI.time)), by="bcr_patient_barcode") %>% 
  dplyr::mutate(days_to_therapy_start = ifelse(days_to_therapy_start <= 0, 9999, days_to_therapy_start)) %>% 
  dplyr::filter(!is.na(days_to_therapy_start)) %>% 
  dplyr::arrange(bcr_patient_barcode, days_to_therapy_start) %>% 
  dplyr::mutate(treat = ifelse(is.na(days_to_therapy_start) | is.na(PFI.time),treat,
                               ifelse(days_to_therapy_start > PFI.time, paste0("**",treat, "**"), treat))
  ) %>% 
  dplyr::filter(!grepl("\\*\\*", treat)) %>% 
  dplyr::group_by(bcr_patient_barcode) %>%
  dplyr::mutate(treatment=paste(unique(sort(treat)), collapse="+")) %>%
  dplyr::ungroup() %>%
  # dplyr::filter(
  #   grepl("^radiation\\+temozolomide$", treatment) |
  #     # grepl("^bevacizumab\\+radiation\\+temozolomide$",treatment) |
  #     # grepl("bevacizumab.*temozolomide.*radiation|bevacizumab.*radiation*.temozolomide|temozolomide.*bevacizumab.*radiation|temozolomide.*radiation.*bevacizumab|radiation.*temozolomide.*bevacizumab|radiation.*bevacizumab.*temozolomide",treatment) |
  #     # grepl("^bevacizumab\\+temozolomide$",treatment) | 
  #   grepl("^radiation$",treatment) |
  #   grepl("^temozolomide$",treatment) 
  #   # (grepl("^bevacizumab$",treatment) & grepl("^radiation$",treatment) & grepl("^temozolomide$",treatment)) |
  #   # grepl("^bevacizumab\\+radiation$",treatment)
  # ) %>% 
  dplyr::group_by(bcr_patient_barcode) %>% 
  dplyr::mutate(days_to_therapy_start = max(days_to_therapy_start,na.rm = TRUE)) %>%
  dplyr::ungroup() %>% 
  dplyr::select(-treat) %>% 
  dplyr::distinct()

head(treatment_lgg, 20); dim(treatment_lgg); table(treatment_lgg$treatment, treatment_lgg$PFI)

# tratamientos <- treatment_lgg %>% 
#   dplyr::select(bcr_patient_barcode, treatment) %>% 
#   dplyr::distinct() %>% 
#   dplyr::group_by(treatment) %>% 
#   dplyr::summarise(n()) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::select(treatment, `n()`) %>% 
#   dplyr::arrange(desc(`n()`)) %>% 
#   dplyr::distinct()

# write.csv(tratamientos, file=paste0(pathFeatureLabel, "/co-mod_best_results_last_run/tratamientos_lgg.csv"), row.names = FALSE)
write.csv(treatment_lgg, file=paste0(pathFeatureLabel, "/co-mod_best_results_last_run/radiation_treatment_cycles_TCGA_LGG.csv"), row.names = FALSE)

## ... 4. Informacion clinica ... ##
# Usar linea de abajo por si sale error de paquete TCGABiolinks
# clinical <- read.csv(paste0(pathFeatureLabel, "/co-mod_best_results_last_run/radiation_treatment_cycles_TCGA_LGG.csv"))
clinical <- GDCquery_clinic(project = "TCGA-LGG", type = "clinical")
clinical <- clinical %>% 
  dplyr::filter(submitter_id %in% barcode_LGG) %>% 
  dplyr::mutate(tissue_or_organ_of_origin=ifelse(is.na(tissue_or_organ_of_origin), "Brain, NOS", tissue_or_organ_of_origin)) %>% 
  dplyr::select(submitter_id, primary_diagnosis, gender, age_at_index, tissue_or_organ_of_origin) %>% 
    dplyr::inner_join(readxl::read_excel(paste0(pathFeatureLabel, "/co-mod_bases_datos_importantes/Ceccarelli_data_subtypes_lgg_TCGA.xlsx")) %>% 
                        dplyr::select(Case, Histology, IDH_codel_subtype, MGMT_promoter_status, TERT_promoter_status, ATRX_status, Telomere_Maintenance,`Chr_7_gain/Chr_10_loss`, 
                                      Karnofsky_Performance_Score, Mutation_Count, Grade) %>% 
                        dplyr::rename(Chr7_gain_and_Chr10_loss = `Chr_7_gain/Chr_10_loss`) %>% 
                        dplyr::rename(submitter_id=Case), by = "submitter_id"
                      ) %>% 
  dplyr::select(submitter_id, primary_diagnosis, age_at_index, gender, IDH_codel_subtype, MGMT_promoter_status, ATRX_status, Telomere_Maintenance, Chr7_gain_and_Chr10_loss, 
                TERT_promoter_status, Karnofsky_Performance_Score, Mutation_Count, Grade) %>% 
  dplyr::mutate_if(is.character, as.factor) %>% 
  dplyr::rename(bcr_patient_barcode=submitter_id)

write.csv(clinical, file=paste0(pathFeatureLabel, "/co-mod_best_results_last_run/clinical_data_TCGA_LGG.csv"), row.names = TRUE)

Y_A_clinical <- dplyr::inner_join(clinical, treatment_lgg, by = "bcr_patient_barcode") %>%
  dplyr::filter(
    grepl("^radiation\\+temozolomide$", treatment) |
      # grepl("^bevacizumab\\+radiation\\+temozolomide$",treatment) |
      # grepl("bevacizumab.*temozolomide.*radiation|bevacizumab.*radiation*.temozolomide|temozolomide.*bevacizumab.*radiation|temozolomide.*radiation.*bevacizumab|radiation.*temozolomide.*bevacizumab|radiation.*bevacizumab.*temozolomide",treatment) |
      # grepl("^bevacizumab\\+temozolomide$",treatment) |
      grepl("^radiation$",treatment)
    # grepl("^temozolomide$",treatment) |
    # (grepl("^bevacizumab$",treatment) & grepl("^radiation$",treatment) & grepl("^temozolomide$",treatment)) |
    # grepl("^bevacizumab\\+radiation$",treatment)
  ) %>% 
  dplyr::mutate(treatment_0_1 = ifelse(grepl("^radiation\\+temozolomide$",treatment), 0,1)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(new.OS.time = PFI.time - days_to_therapy_start)
head(Y_A_clinical); dim(Y_A_clinical); table(Y_A_clinical$treatment, Y_A_clinical$PFI)

write.csv(Y_A_clinical, file=paste0(pathFeatureLabel, "/co-mod_best_results_last_run/Y_A_clinical_LGG.csv"), row.names = TRUE)

# Y_A_clinical <- dplyr::inner_join(clinical, treatment_lgg, by = "bcr_patient_barcode") %>% 
#   dplyr::mutate(treatment_0_1 = ifelse(grepl("^radiation\\+temozolomide$",treatment), 0,1)) %>% 
#   dplyr::distinct() %>% 
#   dplyr::mutate(new.PFI.time = PFI.time - days_to_therapy_start)
# 
# Y_A_clinical[grepl("^radiation\\+temozolomide$", Y_A_clinical$treatment),"treatment"] <- "RAD+TMZ"
# Y_A_clinical[grepl("^radiation$", Y_A_clinical$treatment),"treatment"] <- "RAD"
# 
# head(Y_A_clinical); dim(Y_A_clinical); table(Y_A_clinical$treatment, Y_A_clinical$PFI)
# 
# trial <- A_clinical %>% select(primary_diagnosis, treatment, gender, 
#                                IDH_codel_subtype, MGMT_promoter_status, ATRX_status, Telomere_Maintenance, OS, OS.time) %>%
#   dplyr::mutate_all(as.vector) %>% 
#   dplyr::mutate(OS=as.factor(OS))
# 
# trial %>% tbl_summary(by=treatment,
#                       statistic = list(all_continuous() ~ "{mean} ({sd})"),
#                       label = list(gender ~ "Gender",
#                         IDH_codel_subtype ~ "IDH codel subtype", 
#                         MGMT_promoter_status ~ "MGMT promoter status", 
#                         ATRX_status ~ "ATRX status", 
#                         OS ~ "Censoring",
#                         OS.time ~ "Overall survival"
#                                    ), 
#                       missing_text = "(Missing)") %>% 
#   bold_labels() %>% 
#   add_overall() %>% 
#   add_n() %>% 
#   add_p() %>% add_q()
# 
# # write(names(barcode_LGG[barcode_LGG %in% clinical.follow.up_$bcr_patient_barcode]), paste0(pathFeatureLabel,"/co-mod_best_results_last_run/barcodes_LGG_debugged.txt"))
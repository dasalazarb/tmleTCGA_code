# https://bioc.ism.ac.jp/packages/3.2/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html#tcgaquery_clinic-tcgaquery_clinicfilt-working-with-clinical-data
# https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html
# Avastin FDA Approval History: https://www.drugs.com/history/avastin.html
library(TCGAbiolinks)
library(dplyr)
library(gtsummary)
library(DT)

detect_central_path <- function(primary_path, fallback_path) {
  if (identical(character(0), dir(primary_path))) fallback_path else primary_path
}

clean_drug_names <- function(clinical_drug, drug_correction) {
  correction_lookup <- setNames(drug_correction$Correction, drug_correction$OldName)

  clinical_drug %>%
    mutate(drug_name = gsub(" ", "", tolower(drug_name)),
           drug_name = coalesce(correction_lookup[drug_name], drug_name)) %>%
    filter(!grepl("nos", drug_name), !grepl("^.{0}$", drug_name)) %>%
    rename(pharmaceutical_therapy_drug_name = drug_name) %>%
    group_by(bcr_patient_barcode, pharmaceutical_therapy_drug_name) %>%
    filter(days_to_drug_therapy_start == min(days_to_drug_therapy_start)) %>%
    ungroup()
}

standardize_radiation <- function(clinical_radiation) {
  clinical_radiation %>%
    select(bcr_patient_barcode, days_to_radiation_therapy_start, days_to_radiation_therapy_end, radiation_type, regimen_indication) %>%
    mutate(radiation_type = case_when(
      is.na(radiation_type) ~ "radiation",
      radiation_type %in% c("external beam", "radioisotope", "external", "other", "") ~ "radiation",
      TRUE ~ tolower(radiation_type)
    )) %>%
    select(-regimen_indication, -days_to_radiation_therapy_end) %>%
    distinct() %>%
    rename(treat = radiation_type, days_to_therapy_start = days_to_radiation_therapy_start)
}

build_treatment <- function(drug_df, radiation_df, survival_df, time_col) {
  time_col <- rlang::ensym(time_col)

  survival_df <- survival_df %>%
    filter(type %in% c("LGG")) %>%
    select(-type) %>%
    mutate(survival_time = as.numeric(.data[[rlang::as_string(time_col)]]))

  bind_rows(
    drug_df %>%
      select(bcr_patient_barcode, days_to_drug_therapy_start, pharmaceutical_therapy_drug_name) %>%
      rename(days_to_therapy_start = days_to_drug_therapy_start, treat = pharmaceutical_therapy_drug_name),
    radiation_df
  ) %>%
    inner_join(survival_df, by = "bcr_patient_barcode") %>%
    mutate(days_to_therapy_start = ifelse(days_to_therapy_start <= 0, 9999, days_to_therapy_start)) %>%
    filter(!is.na(days_to_therapy_start)) %>%
    arrange(bcr_patient_barcode, days_to_therapy_start) %>%
    mutate(treat = ifelse(is.na(days_to_therapy_start) | is.na(survival_time), treat,
                          ifelse(days_to_therapy_start > survival_time, paste0("**", treat, "**"), treat))) %>%
    filter(!grepl("\\*\\*", treat)) %>%
    group_by(bcr_patient_barcode) %>%
    mutate(
      treatment = paste(unique(sort(treat)), collapse = "+"),
      days_to_therapy_start = max(days_to_therapy_start, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(-treat) %>%
    distinct()
}

central_path <- detect_central_path(
  "D:/pathFeatureLabelKernels",
  "G:/.shortcut-targets-by-id/1xYVRN5UwMFXmj3OdNbFJmArfObc9kWSH/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/data"
)

pathFeatureLabel <- central_path

### ---------------------- ###
#### ..... LGG info ..... ####
### ---------------------- ###
barcode_LGG <- read.table(paste0(pathFeatureLabel, "/barcodes_TCGA_LGG.txt"))
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
drug_correction <- read.csv(file = paste0(pathFeatureLabel, "/DrugCorrection1.csv")) %>%
  mutate(OldName = gsub(" ", "", tolower(OldName)),
         Correction = gsub(" ", "", tolower(Correction)))

drug_lgg <- clean_drug_names(clinical.drug, drug_correction)
head(drug_lgg);dim(drug_lgg); sort(unique(drug_lgg$pharmaceutical_therapy_drug_name))

write.csv(drug_lgg, file=paste0(pathFeatureLabel, "/drug_treatment_TCGA_LGG.csv"), row.names = FALSE)

### ... 2. Radiation ... ###
radiation_lgg <- standardize_radiation(clinical.radiation)

head(radiation_lgg);dim(radiation_lgg)

## ... 3. treatment: RAD + DRUGs ... ##
treatment_lgg <- build_treatment(
  drug_lgg,
  radiation_lgg,
  read.csv(paste0(pathFeatureLabel, "/survivalTCGA.csv"), header = TRUE, sep = ",") %>%
    select(bcr_patient_barcode, type, OS, OS.time),
  time_col = OS.time
)

head(treatment_lgg, 20); dim(treatment_lgg); table(treatment_lgg$treatment, treatment_lgg$OS)

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
write.csv(treatment_lgg, file=paste0(pathFeatureLabel, "/radiation_treatment_cycles_TCGA_LGG.csv"), row.names = FALSE)

## ... 4. Informacion clinica ... ##
# Usar linea de abajo por si sale error de paquete TCGABiolinks
# clinical <- read.csv(paste0(pathFeatureLabel, "/co-mod_best_results_last_run/radiation_treatment_cycles_TCGA_LGG.csv"))
clinical <- GDCquery_clinic(project = "TCGA-LGG", type = "clinical")
clinical <- clinical %>% 
  dplyr::filter(submitter_id %in% barcode_LGG) %>% 
  dplyr::mutate(tissue_or_organ_of_origin=ifelse(is.na(tissue_or_organ_of_origin), "Brain, NOS", tissue_or_organ_of_origin)) %>% 
  dplyr::select(submitter_id, primary_diagnosis, gender, age_at_index, tissue_or_organ_of_origin) %>% 
    dplyr::inner_join(readxl::read_excel(paste0(pathFeatureLabel, "/Ceccarelli_data_subtypes_lgg_TCGA.xlsx")) %>% 
                        dplyr::select(Case, Histology, IDH_codel_subtype, MGMT_promoter_status, TERT_promoter_status, ATRX_status, Telomere_Maintenance,`Chr_7_gain/Chr_10_loss`, 
                                      Karnofsky_Performance_Score, Mutation_Count, Grade) %>% 
                        dplyr::rename(Chr7_gain_and_Chr10_loss = `Chr_7_gain/Chr_10_loss`) %>% 
                        dplyr::rename(submitter_id=Case), by = "submitter_id"
                      ) %>% 
  dplyr::select(submitter_id, primary_diagnosis, age_at_index, gender, IDH_codel_subtype, MGMT_promoter_status, ATRX_status, Telomere_Maintenance, Chr7_gain_and_Chr10_loss, 
                TERT_promoter_status, Karnofsky_Performance_Score, Mutation_Count, Grade) %>% 
  dplyr::mutate_if(is.character, as.factor) %>% 
  dplyr::rename(bcr_patient_barcode=submitter_id)

write.csv(clinical, file=paste0(pathFeatureLabel, "/clinical_data_TCGA_LGG.csv"), row.names = TRUE)

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
  dplyr::mutate(new.OS.time = OS.time - days_to_therapy_start)
head(Y_A_clinical); dim(Y_A_clinical); table(Y_A_clinical$treatment, Y_A_clinical$OS)

write.csv(Y_A_clinical, file=paste0(pathFeatureLabel, "/Y_A_clinical_LGG.csv"), row.names = TRUE)


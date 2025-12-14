## ---------------------------
##
## Script name: synthetic_integration_examples
##
## Purpose of script: generar datos sintéticos reproducibles para ejemplos de
## integración temprana, intermedia y tardía en pipelines TMLE.
##
## Author: ChatGPT
##
## Date Created: 2024-09-14
## Last Updated: 2025-02-06
##
## ---------------------------
##
## Notas:
## - No depende de datos externos; todo se genera con distribuciones normales y
##   categóricas simples.
## - Guarda un archivo .rds con las tres vistas de integración si se ejecuta
##   directamente el script.
## ---------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(stats)
})

scale_numeric <- function(x) {
  as.numeric(scale(x))
}

#' Genera un conjunto de datos sintéticos con múltiples modalidades
#' @param n Observaciones a simular
#' @param seed Semilla para reproducibilidad
#' @param p_mrna Número de variables para la matriz de expresión génica
#' @param p_cnv Número de variables para la matriz de CNV
#' @return Lista con Y, A, Delta y vistas de integración temprana, intermedia y tardía
#' @export
#' @examples
#' data_syn <- generate_synthetic_integration_data(n = 300, seed = 20240201)
#' names(data_syn)
#' glimpse(data_syn$W_intermediate)
generate_synthetic_integration_data <- function(n = 400, seed = 20240915, p_mrna = 50, p_cnv = 30) {
  set.seed(seed)

  clinical <- tibble(
    age = rnorm(n, mean = 60, sd = 8),
    gender = sample(c("Female", "Male"), n, replace = TRUE, prob = c(0.45, 0.55)),
    subtype = sample(c("IDHmut", "IDHwt", "NotReported"), n, replace = TRUE, prob = c(0.35, 0.55, 0.10))
  )

  # Probabilidad de tratamiento dependiente de clínicas
  lin_trt <- -0.2 + 0.02 * (clinical$age - 60) +
    0.6 * (clinical$gender == "Female") + 0.4 * (clinical$subtype == "IDHmut") +
    rnorm(n, 0, 0.35)
  g1W <- plogis(lin_trt)
  A <- rbinom(n, 1, g1W)

  # Matrices ómicas sintéticas (sin correlación complicada para rapidez)
  mrna <- matrix(rnorm(n * p_mrna, mean = 0.3 * A, sd = 1), nrow = n, ncol = p_mrna)
  colnames(mrna) <- sprintf("mrna_%02d", seq_len(p_mrna))

  cnv <- matrix(rnorm(n * p_cnv, mean = 0.2 * A, sd = 1.2), nrow = n, ncol = p_cnv)
  colnames(cnv) <- sprintf("cnv_%02d", seq_len(p_cnv))

  # Riesgo verdadero con combinación clínica + ómica
  mrna_signal <- rowMeans(mrna[, 1:5, drop = FALSE])
  cnv_signal <- rowMeans(cnv[, 1:3, drop = FALSE])

  risk_score <- -1 + 0.5 * A + 0.4 * scale_numeric(mrna_signal) +
    0.3 * scale_numeric(cnv_signal) + 0.25 * (clinical$gender == "Female") -
    0.02 * (clinical$age - 60)

  event_prob <- plogis(risk_score)
  Y <- rbinom(n, 1, event_prob)

  # Mecanismo de censura ligero
  censor_lp <- -0.2 + 0.3 * A - 0.25 * (clinical$subtype == "NotReported") + 0.01 * (clinical$age - 60)
  Delta <- rbinom(n, 1, plogis(censor_lp))

  clinical_fac <- clinical %>%
    mutate(across(where(is.character), as.factor))

  W_early <- bind_cols(clinical_fac, as_tibble(mrna), as_tibble(cnv))

  # Integración intermedia: PCA por modalidad y unión con clínicas
  mrna_pca <- prcomp(mrna, center = TRUE, scale. = TRUE, rank. = 5)
  cnv_pca <- prcomp(cnv, center = TRUE, scale. = TRUE, rank. = 3)

  W_intermediate <- bind_cols(
    clinical_fac,
    as_tibble(mrna_pca$x, .name_repair = ~sprintf("mrna_PC%02d", seq_along(.x))),
    as_tibble(cnv_pca$x, .name_repair = ~sprintf("cnv_PC%02d", seq_along(.x)))
  )

  # Integración tardía: listas separadas para ensamblado posterior
  W_late <- list(
    clinical = clinical_fac,
    mrna = as_tibble(mrna),
    cnv = as_tibble(cnv)
  )

  list(
    Y = Y,
    A = A,
    Delta = Delta,
    W_early = W_early,
    W_intermediate = W_intermediate,
    W_late = W_late
  )
}

if (sys.nframe() == 0) {
  synthetic_data <- generate_synthetic_integration_data()
  dir.create("data", showWarnings = FALSE, recursive = TRUE)
  saveRDS(synthetic_data, file = file.path("data", "synthetic_integration_example.rds"))
  message("Datos sintéticos guardados en data/synthetic_integration_example.rds")
}

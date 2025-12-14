#!/usr/bin/env Rscript
# Ejemplo autocontenido para lanzar un flujo de integración temprana.
# - Usa datos sintéticos de synthetic_integration_examples.R
# - Reutiliza los envoltorios de SuperLearner definidos en func_SuperLearners.R
# - Calcula los componentes mínimos (g1W, pDelta1, Q) y ejecuta tmle.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(SuperLearner)
  library(tmle)
})

# Cargar las utilidades del repositorio
source("synthetic_integration_examples.R")
source("func_extractVars.R")
source("func_SuperLearners.R")

# 1) Datos ---------------------------------------------------------------
data_syn <- generate_synthetic_integration_data(n = 250, seed = 20240201)

Y <- data_syn$Y
A <- data_syn$A
Delta <- data_syn$Delta
W <- data_syn$W_early %>%
  mutate(A = A)

# Renombrar columnas clínicas para compatibilidad con los screeners
W <- W %>%
  rename_with(~ paste0("clinical_", .x), .cols = c("age", "gender", "subtype")) %>%
  mutate(across(c("clinical_gender", "clinical_subtype", "A"), as.factor))

# 2) Probabilidad de tratamiento g(A=1|W) -------------------------------
clinical <- W %>% select(starts_with("clinical"))

sl_g <- SL.early.rf(Y = A, X = clinical, family = binomial())
g1W <- sl_g$pred

# 3) Mecanismo de observación p(Delta=1|A,W) ----------------------------
sl_delta <- SL.early.rf(Y = Delta, X = clinical, family = binomial())

# Predecir para A=0 y A=1
pDelta1 <- tibble(
  d0W = predict(sl_delta$fit$object, newdata = mutate(clinical, A = factor(0)), type = "response")$pred,
  d1W = predict(sl_delta$fit$object, newdata = mutate(clinical, A = factor(1)), type = "response")$pred
)

# 4) Modelo de resultado Q(W,A) ----------------------------------------
sl_Q <- SL.early.rf(Y = Y, X = W, family = binomial())

Q0W <- predict(sl_Q$fit$object, newdata = mutate(W %>% select(-A), A = 0), type = "response")$pred
Q1W <- predict(sl_Q$fit$object, newdata = mutate(W %>% select(-A), A = 1), type = "response")$pred

Q_list <- list(Q0W = Q0W, Q1W = Q1W)

# 5) TMLE ---------------------------------------------------------------
result <- tmle(
  Y = Y,
  A = A,
  W = W,
  Q = Q_list,
  g1W = g1W,
  Delta = Delta,
  pDelta1 = pDelta1,
  family = "binomial",
  verbose = TRUE
)

print(result)

#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: run_tmle_pipeline
##
## Purpose: Ejecutar un pipeline TMLE de punta a punta. Si no se dispone
## de datos propios, genera datos sintéticos y ejecuta el flujo completo
## con bibliotecas sencillas de SuperLearner.
##
## Uso:
##   Rscript pipelines/run_tmle_pipeline.R                # usa TMLE_TCGA_DATA si existe, o datos sintéticos
##   Rscript pipelines/run_tmle_pipeline.R --data path.csv
##   Rscript pipelines/run_tmle_pipeline.R --simulate --n 500 --seed 123
##
## Salida:
##   - Imprime el ATE y el IC 95% en consola.
##   - Guarda un resumen en results/tmle_summary.txt
##
## Dependencias: tmle, SuperLearner, dplyr, readr
## ---------------------------

suppressPackageStartupMessages({
  library(tmle)
  library(SuperLearner)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tools)
})

source(file.path("synthetic_integration_examples.R"))

## ---------------------------
## Utilidades CLI
## ---------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  list(
    data = get_option_value(args, "--data"),
    simulate = "--simulate" %in% args,
    n = as.integer(get_option_value(args, "--n", default = NA)),
    seed = as.integer(get_option_value(args, "--seed", default = NA))
  )
}

get_option_value <- function(args, option, default = NULL) {
  if (option %in% args) {
    pos <- match(option, args)
    if (pos < length(args)) {
      return(args[pos + 1])
    }
  }
  default
}

## ---------------------------
## Carga o simulación de datos
## ---------------------------
load_or_simulate_data <- function(data_path = NULL, n = 300L, seed = 20240201L) {
  if (!is.null(data_path)) {
    if (!file.exists(data_path)) {
      message("No se encontró el archivo en ", data_path, "; se usarán datos sintéticos.")
    } else {
      ext <- tolower(file_ext(data_path))
      if (ext == "rds") {
        obj <- readRDS(data_path)
        if (is.list(obj) && all(c("Y", "A", "Delta", "W_early") %in% names(obj))) {
          return(list(Y = obj$Y, A = obj$A, Delta = obj$Delta %||% rep(1, length(obj$Y)), W = obj$W_early))
        } else if (is.data.frame(obj) && all(c("Y", "A") %in% names(obj))) {
          return(standardize_frame(obj))
        }
      }

      if (ext %in% c("csv", "txt")) {
        df <- read_csv(data_path, show_col_types = FALSE)
        if (all(c("Y", "A") %in% names(df))) {
          return(standardize_frame(df))
        }
      }

      warning("El archivo no tiene el formato esperado; se generan datos sintéticos.")
    }
  }

  synthetic <- generate_synthetic_integration_data(n = n, seed = seed)
  list(Y = synthetic$Y, A = synthetic$A, Delta = synthetic$Delta, W = synthetic$W_early)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

standardize_frame <- function(df) {
  df <- as.data.frame(df)
  if (!"Delta" %in% names(df)) df$Delta <- 1
  W <- df[, !(names(df) %in% c("Y", "A", "Delta")), drop = FALSE]
  W <- W %>% mutate(across(where(is.character), as.factor))
  list(Y = df$Y, A = df$A, Delta = df$Delta, W = W)
}

## ---------------------------
## Ejecución TMLE
## ---------------------------
run_tmle <- function(data_bundle) {
  set.seed(123)
  tmle(
    Y = data_bundle$Y,
    A = data_bundle$A,
    W = data_bundle$W,
    Delta = data_bundle$Delta,
    family = "binomial",
    Q.SL.library = c("SL.glm", "SL.mean"),
    g.SL.library = c("SL.glm", "SL.mean"),
    verbose = FALSE
  )
}

## ---------------------------
## Main
## ---------------------------
main <- function() {
  args <- parse_args()
  data_path <- args$data %||% Sys.getenv("TMLE_TCGA_DATA", unset = NA)
  using_synthetic <- isTRUE(args$simulate)

  if (!is.na(data_path) && !using_synthetic) {
    bundle <- load_or_simulate_data(data_path)
  } else {
    bundle <- load_or_simulate_data(n = ifelse(is.na(args$n), 300L, args$n), seed = ifelse(is.na(args$seed), 20240201L, args$seed))
    using_synthetic <- TRUE
  }

  message(ifelse(using_synthetic, "Ejecutando con datos sintéticos.", paste0("Ejecutando con archivo: ", data_path)))
  fit <- run_tmle(bundle)

  results <- tibble(
    psi = fit$estimates$ATE$psi,
    ci_lower = fit$estimates$ATE$CI[1],
    ci_upper = fit$estimates$ATE$CI[2],
    used_synthetic = using_synthetic
  )

  dir.create("results", showWarnings = FALSE)
  write_lines(capture.output(print(results)), file.path("results", "tmle_summary.txt"))
  print(results)
}

if (identical(sys.nframe(), 0L)) {
  main()
}

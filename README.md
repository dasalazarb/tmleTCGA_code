# tmleTCGA

Repository with scripts to experiment with Targeted Maximum Likelihood Estimation (TMLE) and machine-learning models applied to TCGA cancer data (mainly GBM and LGG). The code is written in both R and Python and covers everything from building feature matrices to evaluating treatment effects with different omics-integration strategies.

## Requirements

### R
- R >= 4.0.
- Main packages: `tmle`, `SuperLearner`, `caret`, `dplyr`, `ggplot2`, `ranger`, `xgboost`, `earth`, `NMF`, `PMA`, `keras`, `ANN2`, `kernlab`, `nsprcomp`, `pcaMethods`, `gtsummary`.
- Some scripts load Windows absolute paths (for example, `C:/Users/...`). Adjust them to your environment before running.

### Python (optional)
- Python 3.8+ with `optuna`, `pandas`, `scikit-learn`, and `tensorflow` for the dimensionality-reduction experiments in `reducDimenMethods_Optuna.py`.

## Repository structure

- **pipelines/run_tmle_pipeline.R**: end-to-end CLI pipeline that consumes a user-provided dataset (CSV/RDS) or automatically generates synthetic data to estimate the ATE with TMLE and SuperLearner.
- **quick_tmle_early_integration.R / quick_tmle_late_integration.R / quick_tmle_intermediate_integration.R / tmle_all_quick-v4.R**: TMLE pipelines for early, late, and intermediate integration of clinical and omics data.
- **tmle_ps_and_missmech.R**: helper functions for TMLE and missingness mechanisms with different configurations.
- **func_learners.R**, **func_learners_v2.R**, **func_SuperLearners.R**: definitions of `SuperLearner` wrappers (random forest, xgboost, MARS, GLM, etc.) and preconfigured ensembles.
- **func_DimenReduc.R**: utilities for dimensionality reduction (PCA, KPCA, NMF, autoencoders) on expression and CNV matrices.
- **synthetic_integration_examples.R**: generates synthetic data for quick examples of early, intermediate, and late integration without external files.
- **func_extractVars.R**, **func_othersFromtmlePackage.R**: helper functions to prepare input variables and maintain compatibility with the `tmle` package.
- **give_O_*.R**, **01_LGG_clinic_treatment_data_*.R**: scripts to prepare clinical/omics data and build input matrices (O = {Y, A, Δ, W}).
- **gtsummary_O_data.R**, **plot_comparison_var50_to_1000.R**: result summaries and visualizations.
- **reducDimenMethods_Optuna.py**: hyperparameter search with Optuna for dimensionality-reduction models in Python.
- **archive/**: deprecated or exploratory scripts preserved for reference.
- **LICENSE**: MIT license.

## Suggested workflow

1. **Prepare clinical and treatment data**: run `01_LGG_clinic_treatment_data_OS.R` or `01_LGG_clinic_treatment_data_PFI.R` to build survival and treatment tables.
2. **Build feature matrices**: use the `give_O_*.R` scripts to combine omics results (mRNA, CNV) with clinical variables to obtain `Y`, `A`, `Δ`, and `W`.
3. **Reduce dimensionality (optional)**: apply `func_DimenReduc.R` or `reducDimenMethods_Optuna.py` to create compressed representations (PCA, NMF, autoencoders) before training models.
4. **Define learners and ensembles**: load `func_learners_v2.R` and `func_SuperLearners.R` to access wrappers tailored to each data modality (clinical, mRNA, CNV) and integration combinations.
5. **Run TMLE**: choose your integration strategy:
   - Early integration: `quick_tmle_early_integration.R` or variants.
   - Intermediate integration (e.g., autoencoders or PCA): `quick_tmle_intermediate_integration.R` and the `tmle_all_quick-v*` versions that combine multiple modalities.
   - Late integration: `quick_tmle_late_integration.R` or `tmle3_Late.R`.
6. **Explore results**: review plots (for example, propensity-score and risk densities) produced in the TMLE scripts or `plot_comparison_var50_to_1000.R`, and summary tables from `gtsummary_O_data.R`.

## Quick R run

### End-to-end CLI pipeline (recomendada)

```bash
# Usa un archivo CSV/RDS con columnas Y, A, Delta y el resto de covariables
Rscript pipelines/run_tmle_pipeline.R --data /ruta/a/datos.csv

# O bien genera datos sintéticos si no tienes archivos locales
Rscript pipelines/run_tmle_pipeline.R --simulate --n 400 --seed 20240301
```

El script guarda un resumen en `results/tmle_summary.txt` y acepta la variable de entorno `TMLE_TCGA_DATA` como fuente por defecto.

### Flujos de notebook/analítico

```r
# 1) Load functions and data
source("func_extractVars.R")
source("func_learners_v2.R")
source("func_SuperLearners.R")
# Adjust the data source paths inside each script to match your environment.
# quick_tmle_early_integration.R uses the environment variable TMLE_TCGA_DATA; if
# the file does not exist, it will automatically generate a synthetic dataset to
# run the early-integration example end-to-end.

# 2) Run a TMLE pipeline (example: early integration)
source("quick_tmle_early_integration.R")

# 3) Review printed outputs or plots produced by the script
```

To use the full ensembles with multiple integration strategies, run `tmle_all_quick-v4.R`, which builds `SuperLearner` libraries with clinical, omics, and dimensionality-reduction learners.

### Self-contained example (no external dependencies)

If your environment does not have R or packages installed, you can try a simplified early-integration flow with a pure Python script:

```bash
python early_integration_example.py
```

The script generates synthetic data, estimates treatment probabilities, missingness mechanisms, and potential outcomes using logistic regression implemented only with the standard library. It prints the sample size, the mean propensity scores, and an ATE estimate via g-computation to validate the full run.

### Self-contained R example for early integration

If you want to see how the functions in the R scripts are chained to estimate `g1W`, `pDelta1`, and `Q` and to call `tmle()` without depending on external files, run:

```r
Rscript early_integration_quickstart.R
```

This flow performs the following:

- `synthetic_integration_examples.R`: generates `Y`, `A`, `Delta`, and the `W_early` view.
- `func_extractVars.R`: provides the screeners `func_clinical`/`new.screen.randomForest` used by the learners.
- `func_SuperLearners.R`: uses `SL.early.rf` to fit the treatment probability, missingness, and outcome models.
- `tmle`: combines `g1W`, `pDelta1`, and `Q` to deliver the TMLE estimate of the treatment effect.

## Minimal example with synthetic data

If you want to test the structure of the three integration strategies without preparing real data, run:

```r
source("synthetic_integration_examples.R")
# generates a data/synthetic_integration_example.rds file with:
# - W_early: clinical + omics concatenated
# - W_intermediate: clinical + principal components by modality
# - W_late: list separated by modality
```

You can then connect those views to your learners in `quick_tmle_early_integration.R`, `quick_tmle_intermediate_integration.R`, or `quick_tmle_late_integration.R` to validate the full flow.

## Data

The scripts expect processed matrices for gene expression (mRNA), copy-number variation (CNV), and clinical variables derived from TCGA. The repository does not include data because of size/privacy constraints; add the input files at the paths referenced in the scripts (`C:/Users/...` in the examples) or adapt the sources to your environment.

## License

This project is distributed under the MIT license. See the `LICENSE` file for details.

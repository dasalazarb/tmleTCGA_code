# tmleTCGA

Repositorio con scripts para experimentar con Targeted Maximum Likelihood Estimation (TMLE) y modelos de *machine learning* aplicados a datos de cáncer de la colección TCGA (principalmente GBM y LGG). El código está organizado en R y Python y cubre desde la construcción de matrices de características hasta la evaluación de efectos de tratamiento con diferentes estrategias de integración de datos ómicos.

## Requisitos

### R
- R \>= 4.0.
- Paquetes principales: `tmle`, `SuperLearner`, `caret`, `dplyr`, `ggplot2`, `ranger`, `xgboost`, `earth`, `NMF`, `PMA`, `keras`, `ANN2`, `kernlab`, `nsprcomp`, `pcaMethods`, `gtsummary`.
- Algunos scripts cargan rutas absolutas de Windows (por ejemplo, `C:/Users/...`). Ajústalas a tu entorno antes de ejecutar.

### Python (opcional)
- Python 3.8+ con `optuna`, `pandas`, `scikit-learn`, `tensorflow` para los experimentos de reducción de dimensionalidad (`reducDimenMethods_Optuna.py`).

## Estructura del repositorio

- **quick_tmle_early_integration.R / quick_tmle_late_integration.R / quick_tmle_intermediate_integration.R / tmle_all_quick-v4.R**: pipelines TMLE con estrategias de integración temprana, tardía e intermedia de datos clínicos y ómicos.
- **tmle_ps_and_missmech.R**: funciones auxiliares para TMLE y mecanismos de falta con diferentes configuraciones.
- **no_used_*.R**: versiones antiguas o dependientes de rutas absolutas que se conservan solo como referencia histórica.
- **func_learners.R**, **func_learners_v2.R**, **func_SuperLearners.R**: definiciones de *wrappers* para `SuperLearner` (random forest, xgboost, MARS, GLM, etc.) y ensamblados preconfigurados.
- **func_DimenReduc.R**: utilidades para reducción de dimensionalidad (PCA, KPCA, NMF, autoencoders) sobre matrices de expresión y CNV.
- **synthetic_integration_examples.R**: genera datos sintéticos listos para ejemplos rápidos de integración temprana, intermedia y tardía sin depender de archivos externos.
- **func_extractVars.R**, **func_othersFromtmlePackage.R**: funciones auxiliares para preparar variables de entrada y compatibilidad con el paquete `tmle`.
- **give_O_*.R**, **01_LGG_clinic_treatment_data_*.R**: scripts de preparación de datos clínicos/ómicos y generación de matrices de entrada (*O* = {Y, A, Δ, W}).
- **gtsummary_O_data.R**, **plot_comparison_var50_to_1000.R**: resúmenes y visualizaciones de resultados.
- **reducDimenMethods_Optuna.py**: búsqueda de hiperparámetros con Optuna para modelos de reducción de dimensionalidad en Python.
- **LICENSE**: licencia MIT.

## Flujo de trabajo sugerido

1. **Preparar datos clínicos y de tratamiento**: ejecutar `01_LGG_clinic_treatment_data_OS.R` o `01_LGG_clinic_treatment_data_PFI.R` para construir tablas de supervivencia y tratamientos.
2. **Construir matrices de características**: usar los scripts `give_O_*.R` para combinar resultados de ómicas (mRNA, CNV) con las variables clínicas y obtener `Y`, `A`, `Δ` y `W`.
3. **Reducir dimensionalidad (opcional)**: aplicar `func_DimenReduc.R` o `reducDimenMethods_Optuna.py` para crear representaciones comprimidas (PCA, NMF, autoencoders) antes de entrenar modelos.
4. **Definir *learners* y ensamblados**: cargar `func_learners_v2.R` y `func_SuperLearners.R` para disponer de *wrappers* específicos por modalidad de datos (clínico, mRNA, CNV) y combinaciones de integración.
5. **Ejecutar TMLE**: elegir la estrategia de integración deseada:
  - Integración temprana: `quick_tmle_early_integration.R` o variantes.
  - Integración intermedia (p. ej. autoencoders o PCA): `quick_tmle_intermediate_integration.R` y versiones `tmle_all_quick-v*` que combinan varias modalidades.
  - Integración tardía: `quick_tmle_late_integration.R` o `tmle3_Late.R`.
6. **Explorar resultados**: revisar gráficos (por ejemplo, densidades de puntaje de propensión y riesgos) generados en los scripts TMLE o `plot_comparison_var50_to_1000.R`, y tablas resumen de `gtsummary_O_data.R`.

## Ejecución rápida en R

```r
# 1) Cargar funciones y datos
source("func_extractVars.R")
source("func_learners_v2.R")
source("func_SuperLearners.R")
# Ajusta las rutas de origen de datos dentro de cada script según tu entorno.
# quick_tmle_early_integration.R usa la variable de entorno TMLE_TCGA_DATA; si
# el archivo no existe, generará automáticamente un dataset sintético para
# ejecutar el ejemplo de integración temprana de principio a fin.

# 2) Ejecutar un pipeline TMLE (ejemplo: integración temprana)
source("quick_tmle_early_integration.R")

# 3) Visualizar las salidas impresas o gráficos producidos por el script
```

Para usar los ensamblados completos con múltiples estrategias de integración, ejecuta `tmle_all_quick-v4.R`, que arma bibliotecas `SuperLearner` con learners clínicos, ómicos y reducciones de dimensionalidad.

### Ejemplo autocontenido (sin dependencias externas)

Si tu entorno no tiene R ni paquetes instalados, puedes probar un flujo simplificado de integración temprana con un script en Python puro:

```bash
python early_integration_example.py
```

El script genera datos sintéticos, estima probabilidades de tratamiento, mecanismo de falta y resultados potenciales usando regresión logística implementada solo con la biblioteca estándar. Al final imprime el tamaño de la muestra, la media de los puntajes de propensión y una estimación del ATE por g-computation para validar el recorrido de principio a fin.

## Ejemplo mínimo con datos sintéticos

Si quieres probar la estructura de las tres estrategias de integración sin preparar datos reales, ejecuta:

```r
source("synthetic_integration_examples.R")
# genera un archivo data/synthetic_integration_example.rds con:
# - W_early: clínicas + ómicas concatenadas
# - W_intermediate: clínicas + componentes principales por modalidad
# - W_late: lista separada por modalidad
```

Luego puedes conectar esas vistas a tus learners en `quick_tmle_early_integration.R`, `quick_tmle_intermediate_integration.R` o `quick_tmle_late_integration.R` para validar el flujo completo.

## Datos

Los scripts esperan matrices procesadas de expresión génica (mRNA), variaciones en el número de copias (CNV) y variables clínicas derivadas de TCGA. El repositorio no incluye datos por restricciones de tamaño/privacidad; agrega los archivos de entrada en las rutas referenciadas en los scripts (`C:/Users/...` en los ejemplos) o adapta las fuentes a tu entorno.

## Licencia

Este proyecto se distribuye bajo la licencia MIT. Consulta el archivo `LICENSE` para más detalles.

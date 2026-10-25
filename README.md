# GRADE de HIPEC en la primera recaída de cáncer de ovario sensible a platino

> **Aviso: prueba técnica.** Este repositorio se preparó como ejercicio de un proceso de selección. Los datos se extrajeron manualmente de tres ensayos publicados y de curvas de Kaplan-Meier digitalizadas. Los resultados ilustran el método y **no constituyen un análisis formal, una revisión sistemática ni una recomendación clínica**.

Pipeline reproducible en R para evaluar con GRADE la evidencia de tres ensayos clínicos aleatorizados sobre quimioterapia hipertérmica intraperitoneal (HIPEC) con platino, añadida a la citorreducción secundaria.

## Pregunta

En mujeres con primera recaída de cáncer epitelial de ovario sensible a platino y citorreducción secundaria completa o casi completa, ¿añadir HIPEC con platino, frente a no añadirla, mejora la supervivencia y a qué costo en toxicidad?

| Ensayo | Diseño | n | HIPEC |
| --- | --- | --- | --- |
| [CHIPOR](https://doi.org/10.1016/S1470-2045(24)00531-X), Lancet Oncol 2024 | Fase III | 415 | Cisplatino; quimioterapia antes de la cirugía |
| [MSK Team Ovary](https://doi.org/10.1200/JCO.21.00605), J Clin Oncol 2021 | Fase II aleatorizado | 98 | Carboplatino 800 mg/m² |
| [HORSE/MITO-18](https://doi.org/10.1200/JCO.24.00686), J Clin Oncol 2025 | Fase III | 167 | Cisplatino |

## Qué hace el pipeline

- **Enfoque GRADE parcialmente contextualizado** con umbrales múltiples por 1.000 pacientes (Zeng et al. 2021–2022), anclados a ESMO-MCBS, y RoB 2. La imprecisión se califica contando los umbrales que cruza el IC fuera del target.
- **Metaanálisis frecuentista**: REML (principal), DerSimonian-Laird, efecto fijo y HKSJ, más omitir un estudio y diferencia de riesgo de Mantel-Haenszel.
- **Metaanálisis bayesiano** (`bayesmeta`) con prior half-normal sobre τ. Da la probabilidad posterior de cada zona de decisión y la sensibilidad a cuatro priors.
- **Grilla de sensibilidad de umbrales** (factor 0,5–2,0), comparación con MID único y tamaño óptimo de información.
- **Reconstrucción de IPD** desde las curvas de Kaplan-Meier (algoritmo de Guyot, `IPDfromKM`), con:
  - validación contra lo publicado;
  - Cox de un paso y prueba de Schoenfeld;
  - HR por tramos;
  - RMST a 60 meses;
  - modelos de Royston-Parmar.
- **Semáforo de riesgo de sesgo RoB 2.**

## Estructura

```
hipec_grade/
├── R/
│   ├── 00_configuracion.R              # paquetes, rutas, umbrales, riesgos basales, juicios GRADE
│   ├── 01_importacion_curado.R         # lectura del Excel, curado y 6 chequeos de validación
│   ├── 02_metaanalisis_frecuentista.R  # REML/DL/FE/HKSJ, omitir uno, RD-MH, forest plots
│   ├── 03_metaanalisis_bayesiano.R     # bayesmeta, probabilidades por zona, priors
│   ├── 04_umbrales_certeza_grade.R     # motor GRADE, MID único, grilla de umbrales, OIS
│   ├── 05a_preparacion_km_digitizado.R # limpieza de las curvas digitalizadas
│   ├── 05_reconstruccion_ipd_km.R      # Guyot, validación, Cox, tramos, RMST, Royston-Parmar
│   ├── 06_riesgo_sesgo_rob2.R          # semáforo y resumen RoB 2
│   └── 99_ejecutar_todo.R              # corre todo en orden y guarda sessionInfo()
├── data/
│   ├── raw/                            # extracción de datos (Excel) y curvas digitalizadas
│   └── derived/                        # tablas curadas (se regeneran)
└── output/
    ├── tablas/                         # CSV de resultados
    └── figuras/                        # PNG
```

## Cómo correrlo

Requiere R 4.2 o posterior (los scripts usan UTF-8). Instale los paquetes una vez:

```r
install.packages(c("here", "readxl", "dplyr", "tidyr", "stringr", "purrr", "readr", "tibble",
                   "forcats", "metafor", "bayesmeta", "survival", "survRM2", "IPDfromKM",
                   "flexsurv", "ggplot2", "patchwork", "scales"))
```

Abra `hipec_grade.Rproj` y ejecute:

```r
source(here::here("R", "99_ejecutar_todo.R"))
```

El script 01 se detiene si los datos curados no coinciden con los totales publicados. Los valores esperados de cada análisis están en la cabecera de los scripts 02 y 04.

## Resultados principales (corrida del 25 de octubre de 2026)

| Desenlace | Efecto relativo (IC 95%) | Diferencia por 1.000 | Certeza |
| --- | --- | --- | --- |
| Mortalidad a 5 años, cisplatino | HR 0,76 (0,60–0,96) | 99 menos (178 a 15 menos) | Moderada |
| Mortalidad a 5 años, todos los platinos | HR 0,86 (0,63–1,17) | 55 menos (163 menos a 58 más) | Baja |
| EA grado ≥3 | RR 1,74 (1,37–2,22) | 157 más (78 a 256 más) | Moderada |
| Falla renal grado ≥3 | RR 6,70 (2,02–22,2) | 82 más (15 a 306 más) | Baja |
| Mortalidad postoperatoria | 0/338 frente a 3/342 | 8 menos (9 menos a 15 más) | Baja |
| Progresión o muerte a 2 años | HR 1,04 (0,72–1,50) | 13 más (122 menos a 132 más) | Muy baja |

## Limitaciones conocidas

- Los tres ensayos no salen de una búsqueda sistemática, y los umbrales se fijaron sin panel de expertos.
- La reconstrucción de IPD pasa la validación en 2 de 6 curvas. Las de CHIPOR tienen un error de calibración del eje Y detectado al compararlas con la tabla de vida publicada: deben redigitalizarse antes de interpretar el HR por tramos o la SLP a nivel IPD.

## Fuentes de datos y derechos

El repositorio incluye solo la extracción de datos propia y las coordenadas digitalizadas, no los artículos. Consulte los originales en los enlaces de la tabla de ensayos.

# ==================================================================================================
# 00_configuracion.R
# HIPEC en cáncer de ovario recurrente sensible a platino: pipeline de síntesis de evidencia GRADE
# Configuración general: paquetes, rutas, umbrales de decisión, riesgos basales, juicios GRADE
# ==================================================================================================
#
# AVISO: PRUEBA TÉCNICA. Este pipeline se escribió como ejercicio de un proceso de selección. Los
# datos se extrajeron manualmente de tres artículos publicados y de curvas de Kaplan-Meier
# digitalizadas. Los resultados ilustran el método; NO constituyen un análisis formal, una
# revisión sistemática ni una recomendación clínica.
#
# __INDEX___
#   1. Paquetes
#   2. Rutas (here::here; raíz del proyecto = carpeta con hipec_grade.Rproj o .here)
#   3. Umbrales de decisión (enfoque parcialmente contextualizado, umbrales múltiples)
#   4. Riesgos basales (grupo control) para convertir efectos relativos en absolutos
#   5. Juicios GRADE distintos de la imprecisión (fijados por los revisores; ver protocolo)
#   6. Tema gráfico, paleta y aviso de prueba técnica
#   7. Resumen
# ==================================================================================================


# __INDEX___ 1. Paquetes ---------------------------------------------------------------------------
# 
install.packages(c("here", "readxl", "dplyr", "tidyr", "stringr", "purrr", "readr", "tibble",
                    "metafor", "bayesmeta", "survival", "survRM2", "IPDfromKM", "flexsurv",
                    "ggplot2", "patchwork", "scales", "forcats"))

hipec_paquetes_vector <- c(
  "here", "readxl", "dplyr", "tidyr", "stringr", "purrr", "readr", "tibble", "forcats",
  "metafor", "bayesmeta", "survival", "survRM2", "IPDfromKM", "flexsurv",
  "ggplot2", "patchwork", "scales"
)

hipec_paquetes_faltantes_vector <- hipec_paquetes_vector[
  !vapply(hipec_paquetes_vector, requireNamespace, logical(1), quietly = TRUE)
]

if (length(hipec_paquetes_faltantes_vector) > 0) {
  stop(
    "Faltan paquetes: ", paste(hipec_paquetes_faltantes_vector, collapse = ", "),
    ". Instálelos con install.packages() antes de correr el pipeline."
  )
}

  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)


set.seed(20261024)
message("[00] Paquetes cargados.")


# __INDEX___ 2. Rutas ------------------------------------------------------------------------------

hipec_rutas_lista <- list(
  excel_crudo    = here::here("data", "raw", "Ensayos_HIPEC_Ovario_Recurrente_Datos_GRADE.xlsx"),
  km_digitizado  = here::here("data", "raw", "km_digitizado"),
  derivados      = here::here("data", "derived"),
  tablas         = here::here("output", "tablas"),
  figuras        = here::here("output", "figuras"),
  modelos        = here::here("output", "modelos")
)

purrr::walk(
  hipec_rutas_lista[c("derivados", "tablas", "figuras", "modelos", "km_digitizado")],
  dir.create, recursive = TRUE, showWarnings = FALSE
)

if (!file.exists(hipec_rutas_lista$excel_crudo)) {
  stop("No se encontró el libro de Excel en: ", hipec_rutas_lista$excel_crudo)
}

message("[00] Rutas resueltas desde la raíz del proyecto: ", here::here())


# __INDEX___ 3. Umbrales de decisión ---------------------------------------------------------------
# Enfoque GRADE parcialmente contextualizado (Zeng et al., J Clin Epidemiol 2021;
# guías GRADE 34-35).
# Ningún panel definió valores, así que los umbrales se anclan a referencias externas y su
# influencia se explora con una sensibilidad (script 04). Unidades: eventos por 1000.
# Zonas (valor absoluto de la diferencia de riesgo): trivial < t_pequeno <= pequeño < t_moderado
#                                                    <= moderado < t_grande <= grande
# Anclajes:
#   - Desenlaces de supervivencia: ESMO-MCBS v1.1 (formulario 2a) da grado 3 a una ganancia de
#     5-10 puntos porcentuales en supervivencia a un tiempo fijo y grado 4 a >= 10 puntos; se toman
#     2 puntos como techo de lo trivial.
#   - Eventos adversos graves: daños transitorios y en su mayoría reversibles -> bandas x2.
#   - Falla renal grave: daño potencialmente permanente -> mismas bandas que supervivencia.
#   - Muerte postoperatoria: desenlace catastrófico y raro -> bandas 4 veces más estrechas.

hipec_umbrales_tbl <- tribble(
  ~desenlace_id, ~desenlace_etiqueta, ~horizonte, ~t_pequeno, ~t_moderado, ~t_grande,
  "os",          "Mortalidad (supervivencia global)",  "5 años",    20,         50,          100,
  "pfs",         "Progresión o muerte (SLP)",          "2 años",    20,         50,          100,
  "ea_g3",       "Eventos adversos grado >=3",         "30-60 días", 50,        100,         200,
  "renal_g3",    "Falla renal grado >=3",              "60 días",   20,         50,          100,
  "mort_postop", "Mortalidad postoperatoria",          "30-60 días", 5,         10,          20
)

# Sensibilidad del script 04: escalamiento multiplicativo de todas las bandas
hipec_umbrales_escala_vector <- seq(0.5, 2.0, by = 0.25)

message("[00] Umbrales de decisión definidos para ", nrow(hipec_umbrales_tbl), " desenlaces.")


# __INDEX___ 4. Riesgos basales --------------------------------------------------------------------
# Los efectos absolutos de los desenlaces de tiempo al evento usan S1 = S0^HR (riesgos
# proporcionales), con S0 la supervivencia del brazo control a un tiempo fijo. Basal principal =
# ensayo con más peso; el alternativo se reporta como sensibilidad. Valores de las figuras.

hipec_riesgo_basal_tbl <- tribble(
  ~desenlace_id, ~escenario, ~fuente, ~supervivencia_control, ~ic_inf, ~ic_sup,
  "os", "principal", "Control de CHIPOR, SG a 5 años", 0.380, 0.31, 0.45,
  "os", "alternativo", "Control de HORSE, SG a 5 años", 0.616, 0.508, 0.724,
  "pfs", "principal", "Control de MSK, SLP a 24 meses", 0.283, 0.165, 0.413,
  "pfs", "alternativo", "Control de HORSE, SLP a 2 años", 0.473, 0.365, 0.581
)
# Desenlaces binarios: el riesgo combinado del brazo control se calcula en el script 01.

message("[00] Riesgos basales fijados (desenlaces de tiempo al evento).")


# __INDEX___ 5. Juicios GRADE distintos de la imprecisión ------------------------------------------
# Las bajadas por riesgo de sesgo, inconsistencia, evidencia indirecta y sesgo de publicación son
# juicios generados por mi en la revisión (justificados en el protocolo). La imprecisión la calcula el script 04 a
# partir del IC y de los umbrales anteriores. `ic_imprecision` indica qué intervalo se usa:
# "fijo" cuando la heterogeneidad ya se penalizó en inconsistencia (evita contarla dos veces).

hipec_juicios_grade_tbl <- tribble(
  ~analisis_id, ~desenlace_id, ~rob, ~inconsistencia, ~indirecta, ~publicacion, ~ic_imprecision,
  "os_todos",          "os",          0,    -1,              0,          0,            "fijo",
  "os_cisplatino",     "os",          0,     0,              0,          0,            "aleatorio",
  "pfs_todos",         "pfs",        -1,    -1,              0,          0,            "fijo",
  "pfs_cisplatino",    "pfs",        -1,     0,              0,          0,            "aleatorio",
  "ea_g3_todos",       "ea_g3",       0,     0,             -1,          0,            "aleatorio",
  "renal_g3_chipor",   "renal_g3",    0,     0,             -1,          0,            "aleatorio",
  "mort_postop_todos", "mort_postop", 0,     0,              0,          0,            "aleatorio"
)

message("[00] Juicios GRADE cargados para ", nrow(hipec_juicios_grade_tbl), " análisis.")


# __INDEX___ 6. Tema gráfico, paleta y aviso de prueba técnica -------------------------------------

hipec_paleta_vector <- c(
  "CHIPOR"   = "#1F5A8A",
  "MSK"      = "#C0504D",
  "HORSE"    = "#4E8F5A",
  "Combinado" = "#222222"
)

hipec_paleta_zonas_vector <- c(
  "Beneficio grande"   = "#1B5E20",
  "Beneficio moderado" = "#4CAF50",
  "Beneficio pequeño"  = "#A5D6A7",
  "Trivial"            = "#E0E0E0",
  "Daño pequeño"       = "#FFCC80",
  "Daño moderado"      = "#FB8C00",
  "Daño grande"        = "#B71C1C"
)

hipec_tema_ggplot <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold"),
    plot.title.position = "plot",
    plot.caption      = element_text(colour = "grey40", hjust = 0),
    legend.position   = "bottom",
    plot.margin       = margin(8, 24, 8, 8),
    strip.text        = element_text(face = "bold")
  )

theme_set(hipec_tema_ggplot)

# Pie de figura común: todas las figuras se rotulan como resultado de una prueba técnica
hipec_aviso_prueba_texto <- paste(
  "Prueba técnica con datos extraídos de publicaciones; no constituye un análisis formal",
  "ni una recomendación clínica."
)


# __INDEX___ 7. Resumen ----------------------------------------------------------------------------
message("[00] AVISO: ", hipec_aviso_prueba_texto)
message(
  "[00] Configuración lista | desenlaces: ", nrow(hipec_umbrales_tbl),
  " | análisis con juicios GRADE: ", nrow(hipec_juicios_grade_tbl),
  " | valores en la grilla de escalamiento: ", length(hipec_umbrales_escala_vector), "."
)

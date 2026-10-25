# ==================================================================================================
# 99_ejecutar_todo.R
# Corre el pipeline completo en orden, cada script en un entorno limpio, y registra tiempo y estado.
# Los scripts 05a, 05 y 06 son opcionales: los dos primeros necesitan las curvas de Kaplan-Meier
# digitalizadas y el 06 dibuja el riesgo de sesgo RoB 2.
# ==================================================================================================
#
# AVISO: PRUEBA TÉCNICA. Datos extraídos manualmente de tres artículos publicados; los resultados
# ilustran el método y NO constituyen un análisis formal ni una recomendación clínica.
#
# __INDEX___
#   1. Definición del pipeline
#   2. Ejecución secuencial con captura de errores
#   3. Información de la sesión (reproducibilidad)
#   4. Resumen
# ==================================================================================================


# __INDEX___ 1. Definición del pipeline ------------------------------------------------------------
pipeline_scripts_tbl <- tibble::tribble(
  ~script,                              ~obligatorio,
  "00_configuracion.R",                 TRUE,
  "01_importacion_curado.R",            TRUE,
  "02_metaanalisis_frecuentista.R",     TRUE,
  "03_metaanalisis_bayesiano.R",        TRUE,
  "04_umbrales_certeza_grade.R",        TRUE,
  "05a_preparacion_km_digitizado.R",    FALSE,
  "05_reconstruccion_ipd_km.R",         FALSE,
  "06_riesgo_sesgo_rob2.R",             FALSE
)


# __INDEX___ 2. Ejecución secuencial ---------------------------------------------------------------
pipeline_estado_tbl <- purrr::pmap_dfr(pipeline_scripts_tbl, function(script, obligatorio) {
  message("\n========== Ejecutando ", script, " ==========")
  inicio_tiempo <- Sys.time()
  resultado_ejecucion <- tryCatch(
    {
      source(here::here("R", script), local = new.env(parent = globalenv()))
      "ok"
    },
    error = function(error_obj) {
      message("[99] ", script, " falló: ", conditionMessage(error_obj))
      if (obligatorio) stop("[99] Falló un script obligatorio; el pipeline se detuvo.",
                            call. = FALSE)
      "omitido"
    }
  )
  tibble::tibble(
    script   = script,
    estado   = resultado_ejecucion,
    segundos = round(as.numeric(difftime(Sys.time(), inicio_tiempo, units = "secs")), 1)
  )
})


# __INDEX___ 3. Información de la sesión -----------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  here::here("output", "sessionInfo.txt")
)
readr::write_csv(pipeline_estado_tbl, here::here("output", "pipeline_estado.csv"))


# __INDEX___ 4. Resumen ----------------------------------------------------------------------------
print(pipeline_estado_tbl)
message(
  "\n[99] Pipeline terminado | ok: ", sum(pipeline_estado_tbl$estado == "ok"),
  " | omitidos: ", sum(pipeline_estado_tbl$estado == "omitido"),
  " | tiempo total: ", sum(pipeline_estado_tbl$segundos), " s",
  "\n[99] Recordatorio: prueba técnica; los resultados no constituyen un análisis formal."
)

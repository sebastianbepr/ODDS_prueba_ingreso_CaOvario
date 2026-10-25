# ==================================================================================================
# 05a_preparacion_km_digitizado.R
# Convierte las exportaciones crudas de WebPlotDigitizer a formato (tiempo_meses, supervivencia)

# ==================================================================================================
#
# AVISO: PRUEBA TÉCNICA. Curvas digitalizadas a mano de figuras publicadas; los resultados
# ilustran el método y NO constituyen un análisis formal ni una recomendación clínica.
#
# __INDEX___
#   1. Configuración y mapa de archivos
#   2. Función auxiliar: leer y limpiar un brazo digitalizado
#   3. Limpieza de todos los brazos y escritura de los insumos del script 05
#   4. Tabla de control de calidad (medianas publicadas frente a digitalizadas)
#   5. Resumen
# ==================================================================================================
#
# Reglas de limpieza (documentadas en el protocolo, sección 7):
#   a) tiempo = x * máximo del eje (los ejes de CHIPOR y MSK terminan en 96 meses; HORSE en 120);
#   b) los tiempos < 0 se llevan a 0 y la supervivencia se recorta a [0, 1] (ruido de calibración);
#   c) los puntos se ordenan por tiempo y se obliga a que la supervivencia no suba (cummin): elimina
#      los clics fuera de orden en las caídas verticales;
#   d) la curva empieza en (0, 1);
#   e) si quedan pacientes en riesgo en el último tiempo publicado pero los clics terminan antes, la
#      curva se extiende plana hasta ese tiempo (el escalón sigue abierto en la figura).


# __INDEX___ 1. Configuración y mapa de archivos ---------------------------------------------------
source(here::here("R", "00_configuracion.R"))

km_originales_ruta <- file.path(hipec_rutas_lista$km_digitizado, "wpd_originales")
hipec_numeros_riesgo_tbl <- readRDS(file.path(hipec_rutas_lista$derivados, "numeros_riesgo.rds"))
hazard_curado_tbl        <- readRDS(file.path(hipec_rutas_lista$derivados, "hazard_curado.rds"))

curvas_mapa_tbl <- tribble(
  ~estudio, ~desenlace_id, ~brazo,    ~eje_maximo, ~archivo_original,
  "CHIPOR", "os",          "HIPEC",   96,  "Supervivencia global HIPEC.csv",
  "CHIPOR", "os",          "Control", 96,  "Supervivencia global control.csv",
  "CHIPOR", "pfs",         "HIPEC",   96,  "Progresi_n libre supervivencia_ HIPEC.csv",
  "CHIPOR", "pfs",         "Control", 96,  "Progresion libre supervivencia_ control.csv",
  "MSK",    "os",          "HIPEC",   96,  "MSK_OS_HIPEC.csv",
  "MSK",    "os",          "Control", 96,  "MSK_OS_CONTROL.csv",
  "MSK",    "pfs",         "HIPEC",   96,  "MSK_PFS_HIPEC.csv",
  "MSK",    "pfs",         "Control", 96,  "MSK_PFS_CONTROL.csv",
  "HORSE",  "os",          "HIPEC",   120, "POST_RECURRENCIA_SUPER_EN_INTEN_TO_TREAT_HIPEC.csv",
  "HORSE",  "os",          "Control", 120, "POST_RECURRENCIA_SUPER_EN_INTEN_TO_TREAT_CONTROL.csv",
  "HORSE",  "pfs",         "HIPEC",   120, "Progres_free_superv_HIPEC.csv",
  "HORSE",  "pfs",         "Control", 120, "Progres_free_control.csv"
) %>%
  mutate(archivo_original = file.path(km_originales_ruta, archivo_original))

stopifnot(all(file.exists(curvas_mapa_tbl$archivo_original)))
message("[05a] Brazos digitalizados encontrados: ", nrow(curvas_mapa_tbl))


# __INDEX___ 2. Función auxiliar: leer y limpiar un brazo digitalizado -----------------------------
curva_digitizada_limpiar <- function(estudio, desenlace_id, brazo, eje_maximo, archivo_original) {
  curva_cruda_tbl <- read_delim(
    archivo_original, delim = ";", col_names = c("x_fraccion", "supervivencia"),
    locale = locale(decimal_mark = ","), show_col_types = FALSE
  ) %>%
    filter(!is.na(x_fraccion), !is.na(supervivencia))

  riesgo_ultimo_tbl <- hipec_numeros_riesgo_tbl %>%
    filter(estudio == !!estudio, desenlace_id == !!desenlace_id, brazo == !!brazo) %>%
    filter(n_riesgo > 0) %>%
    slice_max(tiempo_meses, n = 1)

  curva_limpia_tbl <- curva_cruda_tbl %>%
    mutate(
      tiempo_meses  = pmin(pmax(x_fraccion * eje_maximo, 0), eje_maximo),
      supervivencia = pmin(pmax(supervivencia, 0), 1)
    ) %>%
    arrange(tiempo_meses, desc(supervivencia)) %>%
    mutate(
      supervivencia_cruda = supervivencia,
      supervivencia       = cummin(supervivencia)
    ) %>%
    filter(tiempo_meses > 0) %>%
    select(tiempo_meses, supervivencia, supervivencia_cruda)

  curva_limpia_tbl <- bind_rows(
    tibble(tiempo_meses = 0, supervivencia = 1, supervivencia_cruda = 1),
    curva_limpia_tbl
  )

  extension_aplicada <- riesgo_ultimo_tbl$tiempo_meses > max(curva_limpia_tbl$tiempo_meses)
  if (extension_aplicada) {
    curva_limpia_tbl <- bind_rows(
      curva_limpia_tbl,
      tibble(tiempo_meses = riesgo_ultimo_tbl$tiempo_meses,
             supervivencia = min(curva_limpia_tbl$supervivencia),
             supervivencia_cruda = NA_real_)
    )
  }

  # Mediana de Kaplan-Meier leída de la curva limpia (primer tiempo con supervivencia <= 0,5)
  mediana_digitizada <- curva_limpia_tbl %>%
    filter(supervivencia <= 0.5) %>%
    summarise(mediana = if (n() == 0) NA_real_ else min(tiempo_meses)) %>%
    pull(mediana)

  write_csv(
    curva_limpia_tbl %>% select(tiempo_meses, supervivencia),
    file.path(hipec_rutas_lista$km_digitizado,
              sprintf("%s_%s_%s.csv", estudio, desenlace_id, brazo))
  )

  tibble(
    estudio, desenlace_id, brazo, eje_maximo,
    puntos_crudos         = nrow(curva_cruda_tbl),
    puntos_limpios        = nrow(curva_limpia_tbl),
    tiempo_max_digitizado = max(curva_cruda_tbl$x_fraccion) * eje_maximo,
    tiempo_negativo_n     = sum(curva_cruda_tbl$x_fraccion < 0),
    supervivencia_mayor_1 = sum(curva_cruda_tbl$supervivencia > 1),
    subidas_corregidas_n  = sum(curva_limpia_tbl$supervivencia_cruda >
                                  curva_limpia_tbl$supervivencia + 1e-9, na.rm = TRUE),
    subida_maxima         = max(curva_limpia_tbl$supervivencia_cruda -
                                  curva_limpia_tbl$supervivencia, na.rm = TRUE),
    extension_hasta       = if (extension_aplicada) riesgo_ultimo_tbl$tiempo_meses else NA_real_,
    supervivencia_final   = min(curva_limpia_tbl$supervivencia),
    mediana_digitizada
  )
}


# __INDEX___ 3. Limpieza de todos los brazos e insumos del script 05 -------------------------------
curvas_control_calidad_tbl <- pmap_dfr(curvas_mapa_tbl, curva_digitizada_limpiar)


# __INDEX___ 4. Tabla de control de calidad --------------------------------------------------------
curvas_control_calidad_tbl <- curvas_control_calidad_tbl %>%
  left_join(
    hazard_curado_tbl %>%
      select(estudio, desenlace_id, mediana_hipec, mediana_control) %>%
      pivot_longer(c(mediana_hipec, mediana_control), names_to = "brazo",
                   values_to = "mediana_publicada") %>%
      mutate(brazo = if_else(brazo == "mediana_hipec", "HIPEC", "Control")),
    by = c("estudio", "desenlace_id", "brazo")
  ) %>%
  mutate(
    mediana_diferencia = mediana_digitizada - mediana_publicada,
    mediana_aceptable  = abs(mediana_diferencia) <= 2
  )

write_csv(curvas_control_calidad_tbl,
          file.path(hipec_rutas_lista$tablas, "05a_control_calidad_km.csv"))
print(curvas_control_calidad_tbl %>%
        select(estudio, desenlace_id, brazo, mediana_publicada, mediana_digitizada,
               subidas_corregidas_n, extension_hasta))


# __INDEX___ 5. Resumen ----------------------------------------------------------------------------
message(
  "[05a] Listo | brazos escritos: ", nrow(curvas_control_calidad_tbl),
  " | puntos fuera de orden corregidos: ", sum(curvas_control_calidad_tbl$subidas_corregidas_n),
  " | brazos extendidos al último tiempo de riesgo: ",
  sum(!is.na(curvas_control_calidad_tbl$extension_hasta)),
  " | medianas a <= 2 meses: ", sum(curvas_control_calidad_tbl$mediana_aceptable), "/",
  nrow(curvas_control_calidad_tbl)
)

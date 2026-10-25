# ==================================================================================================
# 01_importacion_curado.R
# ==================================================================================================
#
# AVISO: PRUEBA TÉCNICA. Datos extraídos manualmente de tres artículos publicados; los resultados
# ilustran el método y NO constituyen un análisis formal ni una recomendación clínica.
#
# __INDEX___
#   1. Configuración
#   2. Funciones auxiliares de lectura
#   3. Lectura de la hoja 3 (matriz de resultados crudos) como texto
#   4. Resultados de tiempo al evento (Tablas 3A y 3B) -> hazard ratios en escala log
#   5. Resultados binarios de seguridad (Tabla 3C) -> conteos 2x2
#   6. Extracción complementaria que no está en el libro (tomada de los PDF)
#   7. Números en riesgo de las figuras de Kaplan-Meier publicadas (insumo del script 05)
#   8. Chequeos de validación
#   9. Escritura de tablas derivadas
#  10. Resumen
# ==================================================================================================

# __INDEX___ 1. Configuración ----------------------------------------------------------------------
source(here::here("R", "00_configuracion.R"))
message("[01] Inicio de la importación y curado de: ", basename(hipec_rutas_lista$excel_crudo))


# __INDEX___ 2. Funciones auxiliares de lectura ----------------------------------------------------

# Extrae el primer par "inferior - superior" de un intervalo en texto como "0.56 - 0.96 (MV: ...)"
intervalo_limites_extraer <- function(intervalo_texto_vector) {
  intervalo_match_matriz <- str_match(
    intervalo_texto_vector,
    "^\\s*([0-9]+\\.?[0-9]*)\\s*-\\s*([0-9]+\\.?[0-9]*|NE)"
  )
  tibble(
    ic_inf = suppressWarnings(as.numeric(intervalo_match_matriz[, 2])),
    ic_sup = suppressWarnings(as.numeric(intervalo_match_matriz[, 3]))
  )
}

# Extrae los conteos de un texto como "20 / 207 (9.7%)"
conteo_fraccion_extraer <- function(fraccion_texto_vector) {
  fraccion_match_matriz <- str_match(fraccion_texto_vector, "^\\s*(\\d+)\\s*/\\s*(\\d+)")
  tibble(
    eventos = as.integer(fraccion_match_matriz[, 2]),
    total   = as.integer(fraccion_match_matriz[, 3])
  )
}

# Unifica las etiquetas de estudio entre hojas
estudio_etiqueta_normalizar <- function(estudio_texto_vector) {
  case_when(
    str_detect(estudio_texto_vector, "CHIPOR") ~ "CHIPOR",
    str_detect(estudio_texto_vector, "MSK")    ~ "MSK",
    str_detect(estudio_texto_vector, "HORSE")  ~ "HORSE",
    TRUE                                       ~ NA_character_
  )
}

# Devuelve las filas de datos de un bloque que empieza en `titulo_patron` y termina antes de "Total"
bloque_filas_extraer <- function(hoja_texto_tbl, titulo_patron) {
  bloque_inicio_indice <- which(str_detect(hoja_texto_tbl$col_01, titulo_patron))
  if (length(bloque_inicio_indice) != 1) {
    stop("Bloque no encontrado o duplicado: ", titulo_patron)
  }
  hoja_restante_tbl <- hoja_texto_tbl %>% slice((bloque_inicio_indice + 2):n())
  bloque_fin_indice <- which(str_detect(hoja_restante_tbl$col_01, "^Total"))[1]
  hoja_restante_tbl %>% slice(seq_len(bloque_fin_indice - 1))
}


# __INDEX___ 3. Lectura de la hoja 3 como texto ----------------------------------------------------
resultados_crudos_texto_tbl <- readxl::read_excel(
  hipec_rutas_lista$excel_crudo,
  sheet        = "3_Datos_Crudos_Resultados",
  range        = "A1:N30",
  col_names    = FALSE,
  col_types    = "text",
  .name_repair = ~ sprintf("col_%02d", seq_along(.x))
)

message("[01] Hoja 3 leída: ", nrow(resultados_crudos_texto_tbl), " filas x ",
        ncol(resultados_crudos_texto_tbl), " columnas (todas como texto).")


# __INDEX___ 4. Resultados de tiempo al evento -----------------------------------------------------
# Columnas: 01 estudio | 02 eventos I | 03 N I | 05 eventos C | 06 N C | 08 mediana I |
#           09 IC mediana I | 10 mediana C | 11 IC mediana C | 12 HR | 13 IC HR | 14 p

supervivencia_bloque_curar <- function(bloque_tbl, desenlace_id) {
  bloque_tbl %>%
    transmute(
      estudio            = estudio_etiqueta_normalizar(col_01),
      desenlace_id       = desenlace_id,
      eventos_hipec      = as.integer(col_02),
      n_hipec            = as.integer(col_03),
      eventos_control    = as.integer(col_05),
      n_control          = as.integer(col_06),
      mediana_hipec      = as.numeric(col_08),
      mediana_control    = as.numeric(col_10),
      hr                 = as.numeric(col_12),
      hr_ic_texto        = col_13,
      p_texto            = col_14
    ) %>%
    bind_cols(intervalo_limites_extraer(.$hr_ic_texto)) %>%
    rename(hr_ic_inf = ic_inf, hr_ic_sup = ic_sup) %>%
    mutate(
      log_hr            = log(hr),
      se_log_hr         = (log(hr_ic_sup) - log(hr_ic_inf)) / (2 * qnorm(0.975)),
      asimetria_ic_log  = abs(log_hr - (log(hr_ic_inf) + log(hr_ic_sup)) / 2),
      proporcion_hipec  = eventos_hipec / n_hipec,
      proporcion_control = eventos_control / n_control
    )
}

supervivencia_global_tbl <- resultados_crudos_texto_tbl %>%
  bloque_filas_extraer("^Tabla 3A") %>%
  supervivencia_bloque_curar("os")

supervivencia_libre_progresion_tbl <- resultados_crudos_texto_tbl %>%
  bloque_filas_extraer("^Tabla 3B") %>%
  supervivencia_bloque_curar("pfs")

# La SLP peritoneal de CHIPOR está debajo de la fila "Total" de 3B; se conserva como secundaria
supervivencia_peritoneal_tbl <- resultados_crudos_texto_tbl %>%
  filter(str_detect(col_01, "PFS Peritoneal")) %>%
  supervivencia_bloque_curar("pfs_peritoneal")

hipec_diseno_estudios_tbl <- tribble(
  ~estudio, ~farmaco_hipec,  ~quimio_secuencia,       ~aleatorizacion_momento,          ~fase,
  "CHIPOR", "Cisplatino",    "Antes de la cirugía",   "Tras la citorreducción (CC0/CC1)", "III",
  "MSK",    "Carboplatino",  "Después de la cirugía", "Antes de completar la citorreducción", "II",
  "HORSE",  "Cisplatino",    "Después de la cirugía", "Tras la citorreducción (CC0/CC1)", "III"
)

hazard_curado_tbl <- bind_rows(
  supervivencia_global_tbl,
  supervivencia_libre_progresion_tbl,
  supervivencia_peritoneal_tbl
) %>%
  left_join(hipec_diseno_estudios_tbl, by = "estudio") %>%
  mutate(estudio = factor(estudio, levels = c("CHIPOR", "MSK", "HORSE")))

message("[01] Filas de tiempo al evento curadas: ", nrow(hazard_curado_tbl))


# __INDEX___ 5. Resultados binarios de seguridad ---------------------------------------------------
seguridad_bloque_tbl <- resultados_crudos_texto_tbl %>% bloque_filas_extraer("^Tabla 3C")

seguridad_ea_g3_tbl <- seguridad_bloque_tbl %>%
  transmute(
    estudio         = estudio_etiqueta_normalizar(col_01),
    desenlace_id    = "ea_g3",
    eventos_hipec   = as.integer(col_02),
    n_hipec         = as.integer(col_03),
    eventos_control = as.integer(col_05),
    n_control       = as.integer(col_06),
    ventana_escala  = col_12,
    dato_inferido   = FALSE
  )

seguridad_fraccion_curar <- function(bloque_tbl, columna_hipec, columna_control, desenlace_id) {
  hipec_conteo_tbl   <- conteo_fraccion_extraer(bloque_tbl[[columna_hipec]])
  control_conteo_tbl <- conteo_fraccion_extraer(bloque_tbl[[columna_control]])
  tibble(
    estudio         = estudio_etiqueta_normalizar(bloque_tbl$col_01),
    desenlace_id    = desenlace_id,
    eventos_hipec   = hipec_conteo_tbl$eventos,
    n_hipec         = hipec_conteo_tbl$total,
    eventos_control = control_conteo_tbl$eventos,
    n_control       = control_conteo_tbl$total,
    ventana_escala  = bloque_tbl$col_12
  )
}

seguridad_renal_g3_tbl <- seguridad_fraccion_curar(seguridad_bloque_tbl, "col_08", "col_09",
                                                   "renal_g3") %>%
  # Nota (c): el cero de MSK se infiere de la ausencia de eventos renales en su Tabla 2
  mutate(dato_inferido = estudio == "MSK")

seguridad_mortalidad_tbl <- seguridad_fraccion_curar(seguridad_bloque_tbl, "col_10", "col_11",
                                                     "mort_postop") %>%
  mutate(dato_inferido = FALSE)

binario_curado_tbl <- bind_rows(seguridad_ea_g3_tbl, seguridad_renal_g3_tbl,
                                seguridad_mortalidad_tbl) %>%
  left_join(hipec_diseno_estudios_tbl, by = "estudio") %>%
  mutate(
    estudio            = factor(estudio, levels = c("CHIPOR", "MSK", "HORSE")),
    proporcion_hipec   = eventos_hipec / n_hipec,
    proporcion_control = eventos_control / n_control
  )

message("[01] Filas binarias de seguridad curadas: ", nrow(binario_curado_tbl))


# __INDEX___ 6. Tasas a tiempo fijo, conteos post hoc con tiosulfato y subgrupos preespecificados.
# Página = página de la revista.

hipec_extraccion_complementaria_tbl <- tribble(
  ~estudio, ~variable,                         ~valor_hipec,       ~valor_control, ~fuente,
  "CHIPOR", "SG a 5 años (%)",                 "46",               "38",           "Fig 2A, p1556",
  "CHIPOR", "Renal G>=3 sin tiosulfato (n/N)", "19/156",           "2/154",        "p1558",
  "CHIPOR", "Renal G>=3 con tiosulfato (n/N)", "1/51",             "1/54",         "p1558",
  "CHIPOR", "HR SG, BRCA mutado",              "1.36 (0.59-3.14)", NA,             "Fig 3, p1558",
  "CHIPOR", "HR SG, BRCA no mutado",           "0.71 (0.48-1.03)", NA,             "Fig 3, p1558",
  "CHIPOR", "HR SLP extraperitoneal",          "0.83 (0.59-1.17)", NA,             "Fig 2D, p1557",
  "CHIPOR", "QLQ-C30 completo, 6/12 m (%)",    "34 / 19",          NA,             "p1558",
  "MSK",    "SG a 24 meses (%)",               "89.6",             "95.8",         "Tabla 4",
  "MSK",    "SLP a 24 meses (%, KM)",          "16.3",             "28.3",         "Tabla 3",
  "MSK",    "Resección completa (%)",          "82",               "94",           "Tabla 1",
  "HORSE",  "SG a 5 años (%)",                 "75.9",             "61.6",         "p855",
  "HORSE",  "SLP a 2 años (%)",                "51.1",             "47.3",         "p855",
  "HORSE",  "HR SLP, histología no serosa",    "0.29 (0.11-0.75)", NA,             "Fig 3, p856",
  "HORSE",  "HR SLP, BRCA positivo",           "1.19 (0.64-2.17)", NA,             "Fig 3, p856"
)


# __INDEX___ 7. Números en riesgo de las figuras de Kaplan-Meier publicadas ------------------------
# Transcritos de: CHIPOR Fig 2A-B; MSK Fig 3A-B; HORSE Fig 2A-B. Tiempo en meses.
# Los usa IPDfromKM::preprocess() junto con las coordenadas digitalizadas (script 05).

numeros_riesgo_crear <- function(estudio, desenlace_id, brazo, tiempos_vector, riesgo_vector) {
  tibble(estudio, desenlace_id, brazo, tiempo_meses = tiempos_vector, n_riesgo = riesgo_vector)
}

hipec_numeros_riesgo_tbl <- bind_rows(
  numeros_riesgo_crear("CHIPOR", "os",  "HIPEC",   seq(0, 96, 12),
                       c(207, 193, 158, 113, 87, 65, 44, 25, 15)),
  numeros_riesgo_crear("CHIPOR", "os",  "Control", seq(0, 96, 12),
                       c(208, 194, 150, 106, 75, 49, 28, 13, 7)),
  numeros_riesgo_crear("CHIPOR", "pfs", "HIPEC",   seq(0, 96, 12),
                       c(207, 87, 45, 26, 19, 14, 11, 6, 5)),
  numeros_riesgo_crear("CHIPOR", "pfs", "Control", seq(0, 96, 12),
                       c(208, 86, 32, 22, 15, 8, 5, 4, 0)),
  numeros_riesgo_crear("MSK",    "os",  "HIPEC",   seq(0, 84, 12),
                       c(49, 48, 40, 26, 13, 6, 2, 1)),
  numeros_riesgo_crear("MSK",    "os",  "Control", seq(0, 84, 12),
                       c(49, 48, 44, 29, 17, 6, 5, 3)),
  numeros_riesgo_crear("MSK",    "pfs", "HIPEC",   seq(0, 60, 12),
                       c(49, 30, 8, 2, 2, 2)),
  numeros_riesgo_crear("MSK",    "pfs", "Control", seq(0, 84, 12),
                       c(49, 33, 12, 7, 6, 4, 4, 3)),
  numeros_riesgo_crear("HORSE",  "os",  "HIPEC",   seq(0, 120, 12),
                       c(82, 81, 77, 67, 61, 52, 35, 24, 12, 7, 6)),
  numeros_riesgo_crear("HORSE",  "os",  "Control", seq(0, 120, 12),
                       c(85, 82, 76, 65, 52, 45, 31, 21, 13, 9, 2)),
  numeros_riesgo_crear("HORSE",  "pfs", "HIPEC",   seq(0, 120, 12),
                       c(82, 69, 41, 31, 23, 17, 10, 8, 4, 1, 1)),
  numeros_riesgo_crear("HORSE",  "pfs", "Control", seq(0, 120, 12),
                       c(85, 66, 41, 25, 21, 19, 15, 11, 5, 5, 1))
)


# __INDEX___ 8. Chequeos de validación y control de calidad---------------------------------------
# Cada chequeo compara los datos curados con los totales que reportan los artículos.

validacion_registro_tbl <- tribble(
  ~chequeo, ~esperado, ~observado,
  "Muertes en SG, 3 ensayos",
  377,
  sum(supervivencia_global_tbl$eventos_hipec + supervivencia_global_tbl$eventos_control),
  "Pacientes en SG, 3 ensayos",
  680,
  sum(supervivencia_global_tbl$n_hipec + supervivencia_global_tbl$n_control),
  "Eventos de SLP, 3 ensayos",
  589,
  sum(supervivencia_libre_progresion_tbl$eventos_hipec +
        supervivencia_libre_progresion_tbl$eventos_control),
  "EA grado >=3, brazo HIPEC",
  124,
  sum(seguridad_ea_g3_tbl$eventos_hipec),
  "EA grado >=3, brazo control",
  72,
  sum(seguridad_ea_g3_tbl$eventos_control),
  "Números en riesgo en t = 0 iguales al tamaño de cada brazo",
  1,
  as.numeric(all(
    hipec_numeros_riesgo_tbl %>% filter(tiempo_meses == 0) %>% pull(n_riesgo) %in%
      c(207, 208, 49, 82, 85)
  ))
) %>%
  mutate(aprobado = esperado == observado)

# Los IC publicados del HR son simétricos en escala log; una asimetría > 0,03 sugiere un error
asimetria_alerta_tbl <- hazard_curado_tbl %>%
  filter(asimetria_ic_log > 0.03) %>%
  select(estudio, desenlace_id, hr, hr_ic_inf, hr_ic_sup, asimetria_ic_log)

if (nrow(asimetria_alerta_tbl) > 0) {
  message("[01] ADVERTENCIA: intervalos de HR asimétricos (revisar el redondeo):")
  print(asimetria_alerta_tbl)
}

if (!all(validacion_registro_tbl$aprobado)) {
  print(validacion_registro_tbl)
  stop("[01] Falló la validación: los datos curados no coinciden con los totales publicados.")
}

message("[01] Los ", nrow(validacion_registro_tbl), " chequeos de validación pasaron.")


# __INDEX___ 9. Escritura de tablas derivadas ------------------------------------------------------
saveRDS(hazard_curado_tbl, file.path(hipec_rutas_lista$derivados, "hazard_curado.rds"))
saveRDS(binario_curado_tbl, file.path(hipec_rutas_lista$derivados, "binario_curado.rds"))
saveRDS(hipec_numeros_riesgo_tbl, file.path(hipec_rutas_lista$derivados, "numeros_riesgo.rds"))
saveRDS(hipec_extraccion_complementaria_tbl,
        file.path(hipec_rutas_lista$derivados, "extraccion_complementaria.rds"))

write_csv(hazard_curado_tbl,       file.path(hipec_rutas_lista$tablas, "01_hazard_curado.csv"))
write_csv(binario_curado_tbl,      file.path(hipec_rutas_lista$tablas, "01_binario_curado.csv"))
write_csv(validacion_registro_tbl, file.path(hipec_rutas_lista$tablas, "01_validacion.csv"))


# __INDEX___ 10. Resumen ---------------------------------------------------------------------------
message(
  "[01] Listo | filas de tiempo al evento: ", nrow(hazard_curado_tbl),
  " | filas binarias: ", nrow(binario_curado_tbl),
  " | filas de números en riesgo: ", nrow(hipec_numeros_riesgo_tbl),
  " | datos inferidos marcados: ", sum(binario_curado_tbl$dato_inferido)
)

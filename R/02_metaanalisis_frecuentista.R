# ==================================================================================================
# 02_metaanalisis_frecuentista.R
# Combinación frecuentista de los tres ensayos: efectos aleatorios por varianza inversa sobre
# log(HR) y log(RR), con análisis de sensibilidad de efecto fijo, DerSimonian-Laird y
# Hartung-Knapp-Sidik-Jonkman, influencia de omitir un estudio, diferencia de riesgo de
# Mantel-Haenszel para desenlaces escasos y forest plots.
# ==================================================================================================
#
# AVISO: PRUEBA TÉCNICA. Datos extraídos manualmente de tres artículos publicados; los resultados
# ilustran el método y NO constituyen un análisis formal ni una recomendación clínica.
#
# __INDEX___
#   1. Configuración e insumos curados
#   2. Mapa de análisis (qué ensayos entran en cada análisis)
#   3. Función auxiliar: ajustar y ordenar una familia de modelos para un análisis
#   4. Análisis de tiempo al evento (SG, SLP; todos los ensayos y solo cisplatino)
#   5. Análisis binarios (EA grado >=3, falla renal, mortalidad postoperatoria)
#   6. Influencia de omitir un estudio
#   7. Forest plots
#   8. Escritura de resultados
#   9. Resumen
# ==================================================================================================

# __INDEX___ 1. Configuración e insumos curados ----------------------------------------------------
source(here::here("R", "00_configuracion.R"))

hazard_curado_tbl  <- readRDS(file.path(hipec_rutas_lista$derivados, "hazard_curado.rds"))
binario_curado_tbl <- readRDS(file.path(hipec_rutas_lista$derivados, "binario_curado.rds"))

message("[02] Insumos curados cargados.")


# __INDEX___ 2. Mapa de análisis -------------------------------------------------------------------
hipec_mapa_analisis_tbl <- tribble(
  ~analisis_id,        ~desenlace_id, ~tipo,     ~estudios,
  "os_todos",          "os",          "hazard",  list(c("CHIPOR", "MSK", "HORSE")),
  "os_cisplatino",     "os",          "hazard",  list(c("CHIPOR", "HORSE")),
  "pfs_todos",         "pfs",         "hazard",  list(c("CHIPOR", "MSK", "HORSE")),
  "pfs_cisplatino",    "pfs",         "hazard",  list(c("CHIPOR", "HORSE")),
  "ea_g3_todos",       "ea_g3",       "binario", list(c("CHIPOR", "MSK", "HORSE")),
  "renal_g3_chipor",   "renal_g3",    "binario", list(c("CHIPOR")),
  "mort_postop_todos", "mort_postop", "binario", list(c("CHIPOR", "MSK", "HORSE"))
) %>%
  mutate(estudios = map(estudios, 1))


# __INDEX___ 3. Función auxiliar: ajustar y ordenar una familia de modelos -------------------------

# Ordena un modelo de metafor en la escala de razón
modelo_resumen_extraer <- function(modelo_rma, analisis_id, modelo_etiqueta) {
  prediccion_lista <- predict(modelo_rma, transf = exp)
  tibble(
    analisis_id   = analisis_id,
    modelo        = modelo_etiqueta,
    k             = modelo_rma$k,
    estimador     = exp(as.numeric(modelo_rma$b)),
    ic_inf        = exp(modelo_rma$ci.lb),
    ic_sup        = exp(modelo_rma$ci.ub),
    pi_inf        = if (!is.null(prediccion_lista$pi.lb)) prediccion_lista$pi.lb else NA_real_,
    pi_sup        = if (!is.null(prediccion_lista$pi.ub)) prediccion_lista$pi.ub else NA_real_,
    tau2          = if (!is.null(modelo_rma$tau2)) modelo_rma$tau2 else NA_real_,
    i2            = if (!is.null(modelo_rma$I2)) modelo_rma$I2 else NA_real_,
    q             = if (!is.null(modelo_rma$QE)) modelo_rma$QE else NA_real_,
    q_p           = if (!is.null(modelo_rma$QEp)) modelo_rma$QEp else NA_real_
  )
}

# Ajusta REML (principal), DL, efecto fijo y HKSJ a efectos en escala log con sus errores estándar
familia_modelos_ajustar <- function(log_efecto_vector, se_vector, estudio_vector, analisis_id) {
  if (length(log_efecto_vector) == 1) {
    modelo_unico <- metafor::rma(yi = log_efecto_vector, sei = se_vector, method = "FE",
                                 slab = estudio_vector)
    return(list(
      modelos = list(fijo = modelo_unico),
      resumen = modelo_resumen_extraer(modelo_unico, analisis_id, "ensayo único")
    ))
  }
  modelos_lista <- list(
    reml = metafor::rma(yi = log_efecto_vector, sei = se_vector, method = "REML",
                        slab = estudio_vector),
    dl   = metafor::rma(yi = log_efecto_vector, sei = se_vector, method = "DL",
                        slab = estudio_vector),
    fijo = metafor::rma(yi = log_efecto_vector, sei = se_vector, method = "FE",
                        slab = estudio_vector),
    hksj = metafor::rma(yi = log_efecto_vector, sei = se_vector, method = "REML",
                        test = "knha", slab = estudio_vector)
  )
  modelos_etiquetas_vector <- c(reml = "REML (principal)", dl = "DerSimonian-Laird",
                                fijo = "Efecto fijo", hksj = "REML + HKSJ")
  list(
    modelos = modelos_lista,
    resumen = imap_dfr(modelos_lista, ~ modelo_resumen_extraer(
      .x, analisis_id, modelos_etiquetas_vector[[.y]]
    ))
  )
}


# __INDEX___ 4. Análisis de tiempo al evento -------------------------------------------------------
hazard_ajustes_lista <- hipec_mapa_analisis_tbl %>%
  filter(tipo == "hazard") %>%
  { set_names(pmap(., function(analisis_id, desenlace_id, tipo, estudios) {
      analisis_datos_tbl <- hazard_curado_tbl %>%
        filter(desenlace_id == !!desenlace_id, estudio %in% estudios)
      familia_modelos_ajustar(analisis_datos_tbl$log_hr, analisis_datos_tbl$se_log_hr,
                              as.character(analisis_datos_tbl$estudio), analisis_id)
    }), .$analisis_id) }

message("[02] Análisis de tiempo al evento ajustados: ", length(hazard_ajustes_lista))


# __INDEX___ 5. Análisis binarios ------------------------------------------------------------------
# Riesgo relativo en escala log con corrección de continuidad de 0,5.

binario_ajustes_lista <- hipec_mapa_analisis_tbl %>%
  filter(tipo == "binario") %>%
  { set_names(pmap(., function(analisis_id, desenlace_id, tipo, estudios) {
      analisis_datos_tbl <- binario_curado_tbl %>%
        filter(desenlace_id == !!desenlace_id, estudio %in% estudios)
      riesgo_relativo_tbl <- metafor::escalc(
        measure = "RR",
        ai = eventos_hipec, n1i = n_hipec, ci = eventos_control, n2i = n_control,
        data = analisis_datos_tbl, drop00 = TRUE, slab = as.character(estudio)
      ) %>%
        as_tibble() %>%
        filter(!is.na(yi))
      ajuste_lista <- familia_modelos_ajustar(
        riesgo_relativo_tbl$yi, sqrt(riesgo_relativo_tbl$vi),
        as.character(riesgo_relativo_tbl$estudio), analisis_id
      )
      # Diferencia de riesgo de Mantel-Haenszel con todos los ensayos (incluye dobles ceros)
      ajuste_lista$diferencia_riesgo_mh <- metafor::rma.mh(
        measure = "RD",
        ai = eventos_hipec, n1i = n_hipec, ci = eventos_control, n2i = n_control,
        data = analisis_datos_tbl
      )
      ajuste_lista$riesgo_control_agrupado <-
        sum(analisis_datos_tbl$eventos_control) / sum(analisis_datos_tbl$n_control)
      ajuste_lista
    }), .$analisis_id) }

diferencia_riesgo_mh_tbl <- imap_dfr(binario_ajustes_lista, ~ tibble(
  analisis_id      = .y,
  dr_mh_por_mil    = 1000 * as.numeric(.x$diferencia_riesgo_mh$b),
  dr_mh_ic_inf     = 1000 * .x$diferencia_riesgo_mh$ci.lb,
  dr_mh_ic_sup     = 1000 * .x$diferencia_riesgo_mh$ci.ub,
  riesgo_control   = .x$riesgo_control_agrupado
))

message("[02] Análisis binarios ajustados: ", length(binario_ajustes_lista))


# __INDEX___ 6. Influencia de omitir un estudio ----------------------------------------------------
influencia_omitir_uno_tbl <- c(hazard_ajustes_lista, binario_ajustes_lista) %>%
  keep(~ "reml" %in% names(.x$modelos) && .x$modelos$reml$k >= 3) %>%
  imap_dfr(function(ajuste_lista, analisis_id) {
    omitir_uno_resultado <- metafor::leave1out(ajuste_lista$modelos$reml, transf = exp)
    tibble(
      analisis_id     = analisis_id,
      estudio_omitido = omitir_uno_resultado$slab,
      estimador       = omitir_uno_resultado$estimate,
      ic_inf          = omitir_uno_resultado$ci.lb,
      ic_sup          = omitir_uno_resultado$ci.ub,
      i2              = omitir_uno_resultado$I2
    )
  })


# __INDEX___ 7. Forest plots -----------------------------------------------------------------------
# Forest en ggplot: cuadros por ensayo según el peso REML, diamante combinado, intervalo de
# predicción.

forest_datos_preparar <- function(ajuste_lista, etiqueta_efecto) {
  modelo_principal <- ajuste_lista$modelos[[if ("reml" %in% names(ajuste_lista$modelos)) "reml"
                                            else "fijo"]]
  estudios_tbl <- tibble(
    fila      = modelo_principal$slab,
    estimador = exp(modelo_principal$yi),
    ic_inf    = exp(modelo_principal$yi - qnorm(0.975) * sqrt(modelo_principal$vi)),
    ic_sup    = exp(modelo_principal$yi + qnorm(0.975) * sqrt(modelo_principal$vi)),
    peso      = weights(modelo_principal),
    tipo_fila = "Ensayo"
  )
  agrupado_tbl <- ajuste_lista$resumen %>%
    filter(modelo %in% c("REML (principal)", "ensayo único")) %>%
    transmute(fila = "Combinado (REML)", estimador, ic_inf, ic_sup, peso = NA_real_,
              tipo_fila = "Combinado", pi_inf, pi_sup)
  bind_rows(estudios_tbl, agrupado_tbl) %>%
    mutate(
      fila = factor(fila, levels = rev(c(modelo_principal$slab, "Combinado (REML)"))),
      etiqueta_texto = sprintf("%.2f (%.2f-%.2f)", estimador, ic_inf, ic_sup),
      etiqueta_efecto = etiqueta_efecto
    )
}

forest_graficar <- function(forest_tbl, titulo_texto) {
  # Columna de etiquetas a la derecha de todos los intervalos (escala log)
  etiqueta_x <- max(c(forest_tbl$ic_sup, forest_tbl$pi_sup), na.rm = TRUE) * 1.15
  limite_inf <- min(c(forest_tbl$ic_inf, forest_tbl$pi_inf), na.rm = TRUE) * 0.9
  ggplot(forest_tbl, aes(x = estimador, y = fila)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = ic_inf, xmax = ic_sup), width = 0.15, orientation = "y") +
    geom_errorbar(data = ~ filter(.x, tipo_fila == "Combinado", !is.na(pi_inf)),
                  aes(xmin = pi_inf, xmax = pi_sup), width = 0, linewidth = 2.5,
                  alpha = 0.15, orientation = "y") +
    geom_point(data = ~ filter(.x, tipo_fila == "Ensayo"),
               aes(size = peso, colour = fila), shape = 15, show.legend = FALSE) +
    geom_point(data = ~ filter(.x, tipo_fila == "Combinado"), shape = 18, size = 6) +
    geom_text(aes(x = etiqueta_x, label = etiqueta_texto), hjust = 0, size = 3.4) +
    scale_x_log10(breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3),
                  limits = c(limite_inf, etiqueta_x * 2), expand = c(0, 0)) +
    scale_colour_manual(values = hipec_paleta_vector) +
    scale_size_continuous(range = c(3, 8)) +
    labs(
      title = titulo_texto,
      x = paste0(unique(forest_tbl$etiqueta_efecto), " (escala log); < 1 favorece a HIPEC"),
      y = NULL,
      caption = paste0("Banda sombreada: intervalo de predicción del 95% (REML). ",
                       "Tamaño del cuadro: peso REML.\n", hipec_aviso_prueba_texto)
    )
}

forest_titulos_vector <- c(
  os_todos          = "Supervivencia global - HIPEC con cualquier platino (3 ECA)",
  os_cisplatino     = "Supervivencia global - HIPEC con cisplatino (2 ECA)",
  pfs_todos         = "Supervivencia libre de progresión - HIPEC con cualquier platino (3 ECA)",
  pfs_cisplatino    = "Supervivencia libre de progresión - HIPEC con cisplatino (2 ECA)",
  ea_g3_todos       = "Eventos adversos postoperatorios grado >=3 (3 ECA)"
)

forest_figuras_lista <- imap(forest_titulos_vector, function(titulo_texto, analisis_id) {
  ajuste_lista   <- c(hazard_ajustes_lista, binario_ajustes_lista)[[analisis_id]]
  etiqueta_efecto <- if (analisis_id %in% names(hazard_ajustes_lista)) "Hazard ratio" else
    "Riesgo relativo"
  forest_graficar(forest_datos_preparar(ajuste_lista, etiqueta_efecto), titulo_texto)
})

iwalk(forest_figuras_lista, ~ ggsave(
  file.path(hipec_rutas_lista$figuras, paste0("02_forest_", .y, ".png")),
  .x, width = 8, height = 3.8, dpi = 320, bg = "white"
))


# __INDEX___ 8. Escritura de resultados ------------------------------------------------------------
frecuentista_resumen_tbl <- bind_rows(
  map_dfr(hazard_ajustes_lista, "resumen"),
  map_dfr(binario_ajustes_lista, "resumen")
) %>%
  left_join(hipec_mapa_analisis_tbl %>% select(analisis_id, desenlace_id, tipo),
            by = "analisis_id")

saveRDS(
  list(hazard = hazard_ajustes_lista, binario = binario_ajustes_lista),
  file.path(hipec_rutas_lista$modelos, "02_modelos_frecuentistas.rds")
)
saveRDS(frecuentista_resumen_tbl, file.path(hipec_rutas_lista$derivados,
                                             "frecuentista_resumen.rds"))
saveRDS(diferencia_riesgo_mh_tbl,
        file.path(hipec_rutas_lista$derivados, "diferencia_riesgo_mh.rds"))

write_csv(frecuentista_resumen_tbl,
          file.path(hipec_rutas_lista$tablas, "02_frecuentista_resumen.csv"))
write_csv(diferencia_riesgo_mh_tbl,
          file.path(hipec_rutas_lista$tablas, "02_diferencia_riesgo_mh.csv"))
write_csv(influencia_omitir_uno_tbl, file.path(hipec_rutas_lista$tablas, "02_omitir_uno.csv"))


# __INDEX___ 9. Resumen ----------------------------------------------------------------------------
frecuentista_resumen_tbl %>%
  filter(modelo %in% c("REML (principal)", "ensayo único")) %>%
  mutate(resultado = sprintf("%-18s %.2f (%.2f-%.2f)  I2 = %s", analisis_id, estimador, ic_inf,
                             ic_sup, ifelse(is.na(i2), "NA", sprintf("%.0f%%", i2)))) %>%
  pull(resultado) %>%
  walk(~ message("[02] ", .x))

message("[02] Listo | forest plots escritos: ", length(forest_figuras_lista))

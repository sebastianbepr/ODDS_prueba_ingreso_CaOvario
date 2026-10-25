# ==================================================================================================
# 04_umbrales_certeza_grade.R
# ==================================================================================================
#
# AVISO: PRUEBA TÉCNICA. Datos extraídos manualmente de tres artículos publicados; los resultados
# ilustran el método y NO constituyen un análisis formal ni una recomendación clínica.
#
# __INDEX___
#   1. Configuración e insumos
#   2. Efectos absolutos (punto, IC de efectos aleatorios, IC de efecto fijo) por análisis y basal
#   3. Definición de targets y regla de cruce de umbrales
#   4. Calificación de la certeza por análisis y target
#   5. Comparación con MID único (enfoque mínimamente contextualizado)
#   6. Sensibilidad de umbrales y mapa de calor
#   7. Tamaño óptimo de información
#   8. Escritura de resultados (incluye un borrador de tabla de resumen de hallazgos)
#   9. Resumen
# ==================================================================================================


# __INDEX___ 1. Configuración e insumos ------------------------------------------------------------
source(here::here("R", "00_configuracion.R"))

frecuentista_resumen_tbl <- readRDS(file.path(hipec_rutas_lista$derivados,
                                              "frecuentista_resumen.rds"))
diferencia_riesgo_mh_tbl <- readRDS(file.path(hipec_rutas_lista$derivados,
                                              "diferencia_riesgo_mh.rds"))
hazard_curado_tbl        <- readRDS(file.path(hipec_rutas_lista$derivados, "hazard_curado.rds"))
binario_curado_tbl       <- readRDS(file.path(hipec_rutas_lista$derivados, "binario_curado.rds"))

certeza_etiquetas_vector <- c("1" = "Muy baja", "2" = "Baja", "3" = "Moderada", "4" = "Alta")
certeza_simbolos_vector  <- c("1" = "⊕◯◯◯", "2" = "⊕⊕◯◯",
                              "3" = "⊕⊕⊕◯", "4" = "⊕⊕⊕⊕")

message("[04] Insumos cargados.")


# __INDEX___ 2. Efectos absolutos ------------------------------------------------------------------
absoluto_desde_ratio <- function(ratio_vector, basal, tipo) {
  if (tipo == "hazard") 1000 * (basal - basal^ratio_vector) else 1000 * basal * (ratio_vector - 1)
}

ratio_modelos_tbl <- frecuentista_resumen_tbl %>%
  filter(modelo %in% c("REML (principal)", "Efecto fijo", "ensayo único")) %>%
  mutate(modelo_corto = case_when(modelo == "Efecto fijo" ~ "fijo", TRUE ~ "aleatorio")) %>%
  select(analisis_id, desenlace_id, tipo, modelo_corto, estimador, ic_inf, ic_sup) %>%
  pivot_wider(names_from = modelo_corto, values_from = c(estimador, ic_inf, ic_sup)) %>%
  # Los análisis de un solo ensayo tienen un único modelo: el IC fijo es igual al aleatorio
  mutate(
    estimador_fijo = coalesce(estimador_fijo, estimador_aleatorio),
    ic_inf_fijo    = coalesce(ic_inf_fijo, ic_inf_aleatorio),
    ic_sup_fijo    = coalesce(ic_sup_fijo, ic_sup_aleatorio)
  )

basal_escenarios_tbl <- bind_rows(
  hipec_riesgo_basal_tbl %>%
    transmute(desenlace_id, escenario, basal = supervivencia_control, fuente),
  # Basal binario = riesgo control de los ensayos que aportan la razón (renal: solo CHIPOR)
  binario_curado_tbl %>%
    filter(!(desenlace_id == "renal_g3" & estudio != "CHIPOR")) %>%
    group_by(desenlace_id) %>%
    summarise(basal = sum(eventos_control) / sum(n_control), .groups = "drop") %>%
    mutate(escenario = "principal", fuente = "Brazos control combinados")
)

efectos_absolutos_tbl <- ratio_modelos_tbl %>%
  inner_join(basal_escenarios_tbl, by = "desenlace_id", relationship = "many-to-many") %>%
  rowwise() %>%
  mutate(
    absoluto_punto         = absoluto_desde_ratio(estimador_aleatorio, basal, tipo),
    absoluto_ic_inf_aleat  = absoluto_desde_ratio(ic_inf_aleatorio, basal, tipo),
    absoluto_ic_sup_aleat  = absoluto_desde_ratio(ic_sup_aleatorio, basal, tipo),
    absoluto_ic_inf_fijo   = absoluto_desde_ratio(ic_inf_fijo, basal, tipo),
    absoluto_ic_sup_fijo   = absoluto_desde_ratio(ic_sup_fijo, basal, tipo)
  ) %>%
  ungroup() %>%
  # La conversión de supervivencia invierte el orden de los límites; se mantiene inferior < superior
  mutate(
    tmp_inf = pmin(absoluto_ic_inf_aleat, absoluto_ic_sup_aleat),
    tmp_sup = pmax(absoluto_ic_inf_aleat, absoluto_ic_sup_aleat),
    absoluto_ic_inf_aleat = tmp_inf, absoluto_ic_sup_aleat = tmp_sup,
    tmp_inf = pmin(absoluto_ic_inf_fijo, absoluto_ic_sup_fijo),
    tmp_sup = pmax(absoluto_ic_inf_fijo, absoluto_ic_sup_fijo),
    absoluto_ic_inf_fijo = tmp_inf, absoluto_ic_sup_fijo = tmp_sup
  ) %>%
  select(-tmp_inf, -tmp_sup)

# Mortalidad postoperatoria: solo CHIPOR tiene eventos (0 frente a 3), así que la razón sale de un
# ensayo y se aplica al riesgo control combinado. 
efectos_absolutos_tbl <- efectos_absolutos_tbl %>%
  left_join(diferencia_riesgo_mh_tbl %>%
              select(analisis_id, dr_mh_por_mil, dr_mh_ic_inf, dr_mh_ic_sup),
            by = "analisis_id")

message("[04] Efectos absolutos calculados para ", n_distinct(efectos_absolutos_tbl$analisis_id),
        " análisis x escenarios de riesgo basal.")


# __INDEX___ 3. Targets y regla de cruce de umbrales -----------------------------------------------
umbrales_vector_construir <- function(t_pequeno, t_moderado, t_grande) {
  c(-t_grande, -t_moderado, -t_pequeno, t_pequeno, t_moderado, t_grande)
}

# Targets candidatos, expresados como intervalos en la escala absoluta (por 1000)
objetivos_construir <- function(punto, t_pequeno, t_moderado, t_grande) {
  limites_vector <- c(-Inf, -t_grande, -t_moderado, -t_pequeno, t_pequeno, t_moderado,
                      t_grande, Inf)
  zona_indice <- findInterval(punto, limites_vector, rightmost.closed = TRUE)
  tribble(
    ~objetivo_id,                    ~objetivo_inf,           ~objetivo_sup,
    "al_menos_pequeno_beneficio",    -Inf,                    -t_pequeno,
    "al_menos_pequeno_dano",         t_pequeno,               Inf,
    "sin_efecto_importante",         -t_pequeno,              t_pequeno,
    "sin_dano_importante",           -Inf,                    t_pequeno,
    "zona_del_punto",                limites_vector[zona_indice], limites_vector[zona_indice + 1]
  )
}

objetivo_principal_elegir <- function(analisis_id, punto, t_pequeno) {
  case_when(
    analisis_id == "mort_postop_todos" ~ "sin_dano_importante",
    punto <= -t_pequeno                ~ "al_menos_pequeno_beneficio",
    punto >=  t_pequeno                ~ "al_menos_pequeno_dano",
    TRUE                               ~ "sin_efecto_importante"
  )
}

umbrales_cruzados_contar <- function(ic_inf, ic_sup, umbrales_vector, objetivo_inf, objetivo_sup) {
  dentro_ic_logico      <- umbrales_vector > ic_inf & umbrales_vector < ic_sup
  dentro_objetivo_logico <- umbrales_vector > objetivo_inf & umbrales_vector < objetivo_sup
  sum(dentro_ic_logico & !dentro_objetivo_logico)
}


# __INDEX___ 4. Calificación de la certeza por análisis y target -----------------------------------
certeza_calificar <- function(efectos_tbl, umbrales_tbl, escala_umbral = 1) {
  efectos_tbl %>%
    inner_join(umbrales_tbl, by = "desenlace_id") %>%
    inner_join(hipec_juicios_grade_tbl, by = c("analisis_id", "desenlace_id")) %>%
    mutate(across(c(t_pequeno, t_moderado, t_grande), ~ .x * escala_umbral)) %>%
    mutate(objetivos = pmap(list(absoluto_punto, t_pequeno, t_moderado, t_grande),
                            objetivos_construir)) %>%
    unnest(objetivos) %>%
    rowwise() %>%
    mutate(
      ic_imp_inf = if (ic_imprecision == "fijo") absoluto_ic_inf_fijo else absoluto_ic_inf_aleat,
      ic_imp_sup = if (ic_imprecision == "fijo") absoluto_ic_sup_fijo else absoluto_ic_sup_aleat,
      umbrales_cruzados = umbrales_cruzados_contar(
        ic_imp_inf, ic_imp_sup, umbrales_vector_construir(t_pequeno, t_moderado, t_grande),
        objetivo_inf, objetivo_sup
      ),
      imprecision = -min(umbrales_cruzados, 3),
      certeza_num = max(1, 4 + rob + inconsistencia + indirecta + publicacion + imprecision),
      objetivo_es_principal = objetivo_id == objetivo_principal_elegir(analisis_id,
                                                                       absoluto_punto, t_pequeno)
    ) %>%
    ungroup() %>%
    mutate(
      certeza = certeza_etiquetas_vector[as.character(certeza_num)],
      simbolo = certeza_simbolos_vector[as.character(certeza_num)]
    )
}

certeza_completa_tbl <- certeza_calificar(efectos_absolutos_tbl, hipec_umbrales_tbl)

certeza_principal_tbl <- certeza_completa_tbl %>%
  filter(escenario == "principal", objetivo_es_principal)

message("[04] Certeza calificada para ", nrow(certeza_principal_tbl),
        " análisis (target principal).")


# __INDEX___ 5. Comparación con MID único ----------------------------------------------------------
# Enfoque mínimamente contextualizado: un único umbral (el límite de efecto pequeño) a cada lado.
mid_unico_tbl <- efectos_absolutos_tbl %>%
  filter(escenario == "principal") %>%
  inner_join(hipec_umbrales_tbl, by = "desenlace_id") %>%
  inner_join(hipec_juicios_grade_tbl, by = c("analisis_id", "desenlace_id")) %>%
  rowwise() %>%
  mutate(
    ic_imp_inf = if (ic_imprecision == "fijo") absoluto_ic_inf_fijo else absoluto_ic_inf_aleat,
    ic_imp_sup = if (ic_imprecision == "fijo") absoluto_ic_sup_fijo else absoluto_ic_sup_aleat,
    mid_cruzados = sum(c(-t_pequeno, t_pequeno) > ic_imp_inf &
                         c(-t_pequeno, t_pequeno) < ic_imp_sup),
    certeza_mid_num = max(1, 4 + rob + inconsistencia + indirecta + publicacion -
                            min(mid_cruzados, 2))
  ) %>%
  ungroup() %>%
  transmute(analisis_id,
            certeza_mid_unico = certeza_etiquetas_vector[as.character(certeza_mid_num)])


# __INDEX___ 6. Sensibilidad de umbrales -------------------------------------------------
sensibilidad_umbrales_tbl <- map_dfr(hipec_umbrales_escala_vector, function(escala_umbral) {
  certeza_calificar(efectos_absolutos_tbl %>% filter(escenario == "principal"),
                    hipec_umbrales_tbl, escala_umbral) %>%
    filter(objetivo_es_principal) %>%
    transmute(analisis_id, escala_umbral, t_pequeno, certeza_num, certeza)
})

analisis_etiquetas_vector <- c(
  os_todos          = "SG, todos los platinos",
  os_cisplatino     = "SG, cisplatino",
  pfs_todos         = "SLP, todos los platinos",
  pfs_cisplatino    = "SLP, cisplatino",
  ea_g3_todos       = "EA grado >=3",
  renal_g3_chipor   = "Falla renal grado >=3",
  mort_postop_todos = "Mortalidad postoperatoria"
)

sensibilidad_figura <- sensibilidad_umbrales_tbl %>%
  mutate(
    analisis = factor(analisis_etiquetas_vector[analisis_id],
                      levels = rev(analisis_etiquetas_vector)),
    certeza  = factor(certeza, levels = c("Muy baja", "Baja", "Moderada", "Alta"))
  ) %>%
  ggplot(aes(x = factor(escala_umbral), y = analisis, fill = certeza)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = certeza), size = 3) +
  geom_vline(xintercept = which(hipec_umbrales_escala_vector == 1), colour = "black",
             linewidth = 0.6) +
  scale_fill_manual(values = c("Muy baja" = "#F4CCCC", "Baja" = "#FCE5CD",
                               "Moderada" = "#D9EAD3", "Alta" = "#93C47D"), drop = FALSE) +
  labs(
    title = "¿Qué tan robusta es cada calificación de certeza frente a los umbrales elegidos?",
    subtitle = "Todos los umbrales multiplicados por el factor del eje x (1 = por defecto)",
    x = "Factor de escala de los umbrales", y = NULL, fill = "Certeza",
    caption = hipec_aviso_prueba_texto
  )

ggsave(file.path(hipec_rutas_lista$figuras, "04_sensibilidad_umbrales.png"), sensibilidad_figura,
       width = 9.5, height = 4.4, dpi = 320, bg = "white")


# __INDEX___ 7. Tamaño óptimo de información -------------------------------------------------------
ois_eventos_schoenfeld <- function(hr_objetivo, alfa = 0.05, potencia = 0.80) {
  4 * (qnorm(1 - alfa / 2) + qnorm(potencia))^2 / log(hr_objetivo)^2
}
ois_n_binario <- function(p_control, rr_objetivo, alfa = 0.05, potencia = 0.80) {
  p_intervencion <- p_control * rr_objetivo
  2 * (qnorm(1 - alfa / 2) + qnorm(potencia))^2 *
    (p_control * (1 - p_control) + p_intervencion * (1 - p_intervencion)) /
    (p_intervencion - p_control)^2
}

ois_tbl <- bind_rows(
  hazard_curado_tbl %>%
    filter(desenlace_id %in% c("os", "pfs")) %>%
    group_by(desenlace_id) %>%
    summarise(observado = sum(eventos_hipec + eventos_control), .groups = "drop") %>%
    crossing(efecto_objetivo = c(0.70, 0.75, 0.80)) %>%
    mutate(metrica = "eventos", requerido = ceiling(ois_eventos_schoenfeld(efecto_objetivo))),
  binario_curado_tbl %>%
    group_by(desenlace_id) %>%
    summarise(observado = sum(n_hipec + n_control),
              p_control = sum(eventos_control) / sum(n_control), .groups = "drop") %>%
    crossing(efecto_objetivo = c(1.25, 1.50, 2.00)) %>%
    mutate(metrica = "pacientes",
           requerido = ceiling(ois_n_binario(p_control, efecto_objetivo))) %>%
    select(-p_control)
) %>%
  mutate(ois_cumplido = observado >= requerido)


# __INDEX___ 8. Escritura de resultados ------------------------------------------------------------
# Borrador de la tabla de resumen de hallazgos 

sof_borrador_tbl <- certeza_principal_tbl %>%
  left_join(mid_unico_tbl, by = "analisis_id") %>%
  transmute(
    analisis = analisis_etiquetas_vector[analisis_id],
    efecto_relativo = sprintf("%s %.2f (%.2f-%.2f)", if_else(tipo == "hazard", "HR", "RR"),
                              estimador_aleatorio, ic_inf_aleatorio, ic_sup_aleatorio),
    basal_por_mil = round(1000 * if_else(tipo == "hazard", 1 - basal, basal)),
    efecto_absoluto = sprintf("%+.0f por 1000 (%+.0f a %+.0f)", absoluto_punto,
                              absoluto_ic_inf_aleat, absoluto_ic_sup_aleat),
    objetivo = objetivo_id,
    ajustes = sprintf("RoB %d | Inc %d | Ind %d | Pub %d | Imp %d (umbrales cruzados: %d)",
                      rob, inconsistencia, indirecta, publicacion, imprecision, umbrales_cruzados),
    certeza = paste(simbolo, certeza),
    certeza_mid_unico
  )

write_csv(efectos_absolutos_tbl, file.path(hipec_rutas_lista$tablas, "04_efectos_absolutos.csv"))
write_csv(certeza_completa_tbl %>% select(-any_of("objetivos")),
          file.path(hipec_rutas_lista$tablas, "04_certeza_todos_objetivos.csv"))
write_csv(sensibilidad_umbrales_tbl,
          file.path(hipec_rutas_lista$tablas, "04_sensibilidad_umbrales.csv"))
write_csv(ois_tbl,                   file.path(hipec_rutas_lista$tablas, "04_ois.csv"))
write_csv(sof_borrador_tbl,          file.path(hipec_rutas_lista$tablas, "04_sof_borrador.csv"))
saveRDS(certeza_completa_tbl,        file.path(hipec_rutas_lista$derivados, "certeza_completa.rds"))


# __INDEX___ 9. Resumen ----------------------------------------------------------------------------
sof_borrador_tbl %>%
  mutate(resultado = sprintf("%-26s %-26s %s", analisis, efecto_absoluto, certeza)) %>%
  pull(resultado) %>%
  walk(~ message("[04] ", .x))

message(
  "[04] Listo | calificaciones estables en toda la grilla de escalamiento: ",
  sensibilidad_umbrales_tbl %>% group_by(analisis_id) %>%
    summarise(estable = n_distinct(certeza) == 1, .groups = "drop") %>%
    summarise(sum(estable)) %>% pull(),
  " de ", n_distinct(sensibilidad_umbrales_tbl$analisis_id), " análisis"
)

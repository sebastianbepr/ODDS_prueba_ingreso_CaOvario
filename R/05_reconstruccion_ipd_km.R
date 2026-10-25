# ==================================================================================================
# 05_reconstruccion_ipd_km.R
# Reconstrucción de datos individuales aproximados (IPD) a partir de las curvas de Kaplan-Meier
# publicadas (Guyot et al., BMC Med Res Methodol 2012; implementación de IPDfromKM), validación
# contra los resultados publicados y análisis a nivel IPD que los datos agregados no permiten:
#   - modelo de Cox estratificado de un paso y diagnóstico de riesgos proporcionales,
#   - hazard ratios por tramos (en CHIPOR el beneficio aparece hacia los 48 meses),
#   - diferencia del tiempo medio de supervivencia restringido (RMST) a 60 meses, combinada,
#   - modelos paramétricos flexibles (Royston-Parmar) con y sin riesgos no proporcionales.
# ==================================================================================================
#
# AVISO: PRUEBA TÉCNICA. Datos reconstruidos a partir de curvas digitalizadas de artículos
# publicados; los resultados ilustran el método y NO constituyen un análisis formal ni una
# recomendación clínica.
#
# __INDEX___
#   1. Configuración e insumos
#   2. Inventario de curvas digitalizadas (y plantilla si faltan archivos)
#   3. Función auxiliar: reconstruir un brazo con IPDfromKM
#   4. Reconstrucción de todos los pares ensayo x desenlace disponibles
#   5. Validación contra los resultados publicados
#   6. Cox de un paso y diagnóstico de riesgos proporcionales
#   7. Hazard ratios por tramos
#   8. RMST a 60 meses (combinación en dos etapas)
#   9. Modelos paramétricos flexibles
#  10. Figuras
#  11. Escritura de resultados
#  12. Resumen
# ==================================================================================================

# LIMITACIÓN: la IPD reconstruida aproxima los datos originales; no hay covariables, así que no
# se pueden hacer análisis ajustados ni por subgrupos. 

# __INDEX___ 1. Configuración e insumos ------------------------------------------------------------
source(here::here("R", "00_configuracion.R"))
library(survival)

hazard_curado_tbl        <- readRDS(file.path(hipec_rutas_lista$derivados, "hazard_curado.rds"))
hipec_numeros_riesgo_tbl <- readRDS(file.path(hipec_rutas_lista$derivados, "numeros_riesgo.rds"))

rmst_horizonte_meses <- 60
cortes_piecewise_vector <- c(24, 48)

message("[05] Insumos cargados; horizonte del RMST = ", rmst_horizonte_meses, " meses.")


# __INDEX___ 2. Inventario de curvas digitalizadas -------------------------------------------------
curvas_inventario_tbl <- hipec_numeros_riesgo_tbl %>%
  distinct(estudio, desenlace_id, brazo) %>%
  mutate(
    archivo = file.path(hipec_rutas_lista$km_digitizado,
                        sprintf("%s_%s_%s.csv", estudio, desenlace_id, brazo)),
    disponible = file.exists(archivo)
  )

curvas_pares_disponibles_tbl <- curvas_inventario_tbl %>%
  group_by(estudio, desenlace_id) %>%
  summarise(ambos_brazos = all(disponible), .groups = "drop") %>%
  filter(ambos_brazos)

if (nrow(curvas_pares_disponibles_tbl) == 0) {
  plantilla_ruta <- file.path(hipec_rutas_lista$km_digitizado,
                              "PLANTILLA_estudio_desenlace_brazo.csv")
  write_csv(tibble(tiempo_meses = c(0, 1.2, 1.2), supervivencia = c(1, 1, 0.98)), plantilla_ruta)
  message("[05] No se encontraron curvas digitalizadas. Plantilla escrita en: ", plantilla_ruta)
  message("[05] Archivos esperados:\n  ",
          paste(basename(curvas_inventario_tbl$archivo), collapse = "\n  "))
  stop("[05] Se detiene: primero digitalice las curvas de Kaplan-Meier (ver la cabecera).")
}

message("[05] Pares ensayo x desenlace disponibles: ", nrow(curvas_pares_disponibles_tbl), " de ",
        nrow(distinct(curvas_inventario_tbl, estudio, desenlace_id)))


# __INDEX___ 3. Función auxiliar: reconstruir un brazo ---------------------------------------------
brazo_ipd_reconstruir <- function(estudio, desenlace_id, brazo, archivo, arm_id) {
  curva_digitizada_tbl <- read_csv(archivo, show_col_types = FALSE) %>%
    rename_with(~ c("tiempo_meses", "supervivencia")[seq_along(.x)]) %>%
    filter(!is.na(tiempo_meses), !is.na(supervivencia)) %>%
    arrange(tiempo_meses)
  escala_maxima <- if (max(curva_digitizada_tbl$supervivencia) > 1.5) 100 else 1

  riesgo_tbl <- hipec_numeros_riesgo_tbl %>%
    filter(estudio == !!estudio, desenlace_id == !!desenlace_id, brazo == !!brazo) %>%
    arrange(tiempo_meses)

  # Corrección v1: tidyselect no evalúa bien un if() dentro de pull(); se nombra la columna
  eventos_columna <- if (brazo == "HIPEC") "eventos_hipec" else "eventos_control"
  eventos_totales <- hazard_curado_tbl %>%
    filter(estudio == !!estudio, desenlace_id == !!desenlace_id) %>%
    pull(.data[[eventos_columna]])

  # Corrección v1: IPDfromKM exige que todos los tiempos de riesgo estén en el rango digitalizado
  riesgo_tbl <- riesgo_tbl %>%
    filter(tiempo_meses <= max(curva_digitizada_tbl$tiempo_meses))

  preprocesado_obj <- IPDfromKM::preprocess(
    dat      = as.data.frame(curva_digitizada_tbl),
    trisk    = riesgo_tbl$tiempo_meses,
    nrisk    = riesgo_tbl$n_riesgo,
    totalpts = riesgo_tbl$n_riesgo[1],
    maxy     = escala_maxima
  )
  ipd_obj <- IPDfromKM::getIPD(prep = preprocesado_obj, armID = arm_id,
                               tot.events = eventos_totales)
  as_tibble(ipd_obj$IPD) %>%
    transmute(estudio = estudio, desenlace_id = desenlace_id, brazo = brazo,
              tiempo_meses = time, estado = status)
}


# __INDEX___ 4. Reconstrucción de todos los pares ensayo x desenlace -------------------------------
ipd_reconstruida_tbl <- curvas_inventario_tbl %>%
  semi_join(curvas_pares_disponibles_tbl, by = c("estudio", "desenlace_id")) %>%
  mutate(arm_id = if_else(brazo == "HIPEC", 1L, 0L)) %>%
  select(estudio, desenlace_id, brazo, archivo, arm_id) %>%
  pmap_dfr(brazo_ipd_reconstruir) %>%
  mutate(
    hipec   = as.integer(brazo == "HIPEC"),
    estudio = factor(estudio, levels = c("CHIPOR", "MSK", "HORSE")),
    farmaco = if_else(estudio == "MSK", "Carboplatino", "Cisplatino")
  )

message("[05] Registros reconstruidos: ", nrow(ipd_reconstruida_tbl), " en ",
        n_distinct(paste(ipd_reconstruida_tbl$estudio, ipd_reconstruida_tbl$desenlace_id)),
        " pares ensayo x desenlace.")


# __INDEX___ 5. Validación contra los resultados publicados ----------------------------------------
# Reconstrucción aceptable: |diferencia de log HR| < 0,05 y eventos dentro del 2%.
validacion_ipd_tbl <- ipd_reconstruida_tbl %>%
  group_by(estudio, desenlace_id) %>%
  group_modify(function(datos_tbl, clave_tbl) {
    cox_modelo <- coxph(Surv(tiempo_meses, estado) ~ hipec, data = datos_tbl)
    km_ajuste  <- survfit(Surv(tiempo_meses, estado) ~ hipec, data = datos_tbl)
    medianas_vector <- summary(km_ajuste)$table[, "median"]
    tibble(
      hr_reconstruido      = exp(coef(cox_modelo)),
      hr_ic_inf_rec        = exp(confint(cox_modelo))[1],
      hr_ic_sup_rec        = exp(confint(cox_modelo))[2],
      mediana_control_rec  = medianas_vector[1],
      mediana_hipec_rec    = medianas_vector[2],
      eventos_rec          = sum(datos_tbl$estado)
    )
  }) %>%
  ungroup() %>%
  left_join(hazard_curado_tbl %>%
              transmute(estudio, desenlace_id, hr_publicado = hr,
                        mediana_hipec_pub = mediana_hipec, mediana_control_pub = mediana_control,
                        eventos_pub = eventos_hipec + eventos_control),
            by = c("estudio", "desenlace_id")) %>%
  mutate(
    dif_log_hr = abs(log(hr_reconstruido) - log(hr_publicado)),
    reconstruccion_aceptable = dif_log_hr < 0.05 &
      abs(eventos_rec - eventos_pub) / eventos_pub < 0.02
  )

print(validacion_ipd_tbl %>% select(estudio, desenlace_id, hr_publicado, hr_reconstruido,
                                    dif_log_hr, reconstruccion_aceptable))
# Corrección v1: la validación se escribe de inmediato, por si falla una sección posterior
write_csv(validacion_ipd_tbl, file.path(hipec_rutas_lista$tablas, "05_validacion_ipd.csv"))

if (!all(validacion_ipd_tbl$reconstruccion_aceptable)) {
  message("[05] ADVERTENCIA: algunas reconstrucciones superan la tolerancia; redigitalice esas ",
          "curvas antes de interpretar los resultados a nivel IPD.")
}


# __INDEX___ 6. Cox de un paso y diagnóstico de riesgos proporcionales -----------------------------
cox_una_etapa_ajustar <- function(datos_tbl) {
  coxph(Surv(tiempo_meses, estado) ~ hipec + strata(estudio), data = datos_tbl)
}

una_etapa_tbl <- crossing(
  desenlace_id = unique(ipd_reconstruida_tbl$desenlace_id),
  conjunto     = c("todos", "cisplatino")
) %>%
  pmap_dfr(function(desenlace_id, conjunto) {
    datos_tbl <- ipd_reconstruida_tbl %>%
      filter(desenlace_id == !!desenlace_id,
             conjunto == "todos" | farmaco == "Cisplatino")
    if (n_distinct(datos_tbl$estudio) < 2) return(NULL)
    cox_modelo <- cox_una_etapa_ajustar(datos_tbl)
    ph_prueba  <- cox.zph(cox_modelo)
    tibble(
      desenlace_id, conjunto,
      k                 = n_distinct(datos_tbl$estudio),
      hr_una_etapa      = exp(coef(cox_modelo)),
      ic_inf            = exp(confint(cox_modelo))[1],
      ic_sup            = exp(confint(cox_modelo))[2],
      p_schoenfeld      = ph_prueba$table["hipec", "p"]
    )
  })


# __INDEX___ 7. Hazard ratios por tramos -----------------------------------------------------------
supervivencia_global_ipd_tbl <- ipd_reconstruida_tbl %>% filter(desenlace_id == "os")

if (n_distinct(supervivencia_global_ipd_tbl$estudio) < 2) {
  saveRDS(ipd_reconstruida_tbl, file.path(hipec_rutas_lista$derivados, "ipd_reconstruida.rds"))
  stop("[05] Hay menos de dos pares de curvas de SG; se guardó la reconstrucción y se omiten ",
       "los análisis de SG.")
}

piecewise_ajustar <- function(conjunto) {
  datos_tbl <- supervivencia_global_ipd_tbl %>%
    filter(conjunto == "todos" | farmaco == "Cisplatino") %>%
    mutate(id = row_number())
  if (n_distinct(datos_tbl$estudio) < 2) return(NULL)
  dividido_tbl <- survSplit(Surv(tiempo_meses, estado) ~ ., data = datos_tbl,
                            cut = cortes_piecewise_vector, episode = "periodo") %>%
    mutate(periodo = factor(periodo, levels = 1:3, labels = c("0-24 m", "24-48 m", ">48 m")))
  cox_modelo <- coxph(Surv(tstart, tiempo_meses, estado) ~ hipec:periodo + strata(estudio),
                      data = dividido_tbl)
  tibble(conjunto, termino = names(coef(cox_modelo)),
         hr = exp(coef(cox_modelo)),
         ic_inf = exp(confint(cox_modelo))[, 1],
         ic_sup = exp(confint(cox_modelo))[, 2])
}

piecewise_tbl <- map_dfr(c("todos", "cisplatino"), piecewise_ajustar)


# __INDEX___ 8. RMST a 60 meses --------------------------------------------------------------------
rmst_por_estudio_tbl <- ipd_reconstruida_tbl %>%
  group_by(estudio, desenlace_id) %>%
  # Corrección v1: survRM2 exige tau <= el mayor tiempo observado en CADA brazo
  filter(min(tapply(tiempo_meses, hipec, max)) >= rmst_horizonte_meses) %>%
  group_modify(function(datos_tbl, clave_tbl) {
    rmst_resultado <- survRM2::rmst2(time = datos_tbl$tiempo_meses, status = datos_tbl$estado,
                                     arm = datos_tbl$hipec, tau = rmst_horizonte_meses)
    diferencia_fila <- rmst_resultado$unadjusted.result[1, ]
    tibble(rmst_dif = diferencia_fila[1], rmst_ic_inf = diferencia_fila[2],
           rmst_ic_sup = diferencia_fila[3])
  }) %>%
  ungroup() %>%
  mutate(rmst_se = (rmst_ic_sup - rmst_ic_inf) / (2 * qnorm(0.975)),
         farmaco = if_else(estudio == "MSK", "Carboplatino", "Cisplatino"))

rmst_agrupado_tbl <- rmst_por_estudio_tbl %>%
  crossing(conjunto = c("todos", "cisplatino")) %>%
  filter(conjunto == "todos" | farmaco == "Cisplatino") %>%
  group_by(desenlace_id, conjunto) %>%
  filter(n() >= 2) %>%
  group_modify(function(datos_tbl, clave_tbl) {
    modelo_rma <- metafor::rma(yi = datos_tbl$rmst_dif, sei = datos_tbl$rmst_se, method = "REML")
    tibble(k = modelo_rma$k, rmst_dif = as.numeric(modelo_rma$b), ic_inf = modelo_rma$ci.lb,
           ic_sup = modelo_rma$ci.ub, i2 = modelo_rma$I2)
  }) %>%
  ungroup()


# __INDEX___ 9. Modelos paramétricos flexibles -----------------------------------------------------
flexible_tbl <- ipd_reconstruida_tbl %>%
  filter(desenlace_id == "os") %>%
  { list(
      ph     = flexsurv::flexsurvspline(Surv(tiempo_meses, estado) ~ hipec + estudio,
                                        data = ., k = 2, scale = "hazard"),
      no_ph  = flexsurv::flexsurvspline(Surv(tiempo_meses, estado) ~ hipec + estudio,
                                        anc = list(gamma1 = ~ hipec),
                                        data = ., k = 2, scale = "hazard")
    ) } %>%
  imap_dfr(~ tibble(modelo = .y, aic = AIC(.x), loglik = .x$loglik))


# __INDEX___ 10. Figuras ---------------------------------------------------------------------------
km_figura <- ipd_reconstruida_tbl %>%
  filter(desenlace_id == "os") %>%
  group_by(estudio, brazo) %>%
  group_modify(function(datos_tbl, clave_tbl) {
    km_ajuste <- survfit(Surv(tiempo_meses, estado) ~ 1, data = datos_tbl)
    tibble(tiempo = c(0, km_ajuste$time), supervivencia = c(1, km_ajuste$surv),
           ic_inf = c(1, km_ajuste$lower), ic_sup = c(1, km_ajuste$upper))
  }) %>%
  ungroup() %>%
  ggplot(aes(x = tiempo, y = supervivencia, colour = brazo, fill = brazo)) +
  geom_step(linewidth = 0.8) +
  geom_vline(xintercept = rmst_horizonte_meses, linetype = "dotted", colour = "grey40") +
  facet_wrap(~ estudio, nrow = 1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(values = c("Control" = "#7A7A7A", "HIPEC" = "#1F5A8A")) +
  labs(title = "Supervivencia global a partir de la IPD reconstruida",
       subtitle = "Algoritmo de Guyot (IPDfromKM); línea punteada = horizonte del RMST (60 meses)",
       x = "Meses desde la aleatorización", y = "Supervivencia global", colour = NULL,
       fill = NULL, caption = hipec_aviso_prueba_texto) +
  guides(fill = "none")

rmst_figura <- bind_rows(
  rmst_por_estudio_tbl %>% filter(desenlace_id == "os") %>%
    transmute(fila = as.character(estudio), rmst_dif, ic_inf = rmst_ic_inf, ic_sup = rmst_ic_sup,
              tipo = "Ensayo"),
  rmst_agrupado_tbl %>% filter(desenlace_id == "os") %>%
    transmute(fila = paste0("Combinado, ", conjunto), rmst_dif, ic_inf, ic_sup,
              tipo = "Combinado")
) %>%
  mutate(fila = forcats::fct_rev(forcats::fct_inorder(fila))) %>%
  ggplot(aes(x = rmst_dif, y = fila, shape = tipo)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = ic_inf, xmax = ic_sup), size = 0.6) +
  scale_shape_manual(values = c("Ensayo" = 15, "Combinado" = 18)) +
  labs(title = "Tiempo medio de supervivencia restringido ganado con HIPEC a 60 meses",
       x = "Diferencia de RMST (meses); > 0 favorece a HIPEC", y = NULL, shape = NULL,
       caption = hipec_aviso_prueba_texto)

ggsave(file.path(hipec_rutas_lista$figuras, "05_km_reconstruidas_os.png"), km_figura,
       width = 11, height = 4.2, dpi = 320, bg = "white")
ggsave(file.path(hipec_rutas_lista$figuras, "05_rmst_60m.png"), rmst_figura,
       width = 8, height = 3.8, dpi = 320, bg = "white")


# __INDEX___ 11. Escritura de resultados -----------------------------------------------------------
saveRDS(ipd_reconstruida_tbl, file.path(hipec_rutas_lista$derivados, "ipd_reconstruida.rds"))
write_csv(validacion_ipd_tbl,   file.path(hipec_rutas_lista$tablas, "05_validacion_ipd.csv"))
write_csv(una_etapa_tbl,        file.path(hipec_rutas_lista$tablas, "05_cox_una_etapa.csv"))
write_csv(piecewise_tbl,        file.path(hipec_rutas_lista$tablas, "05_hr_piecewise.csv"))
write_csv(rmst_por_estudio_tbl, file.path(hipec_rutas_lista$tablas, "05_rmst_por_estudio.csv"))
write_csv(rmst_agrupado_tbl,    file.path(hipec_rutas_lista$tablas, "05_rmst_agrupado.csv"))
write_csv(flexible_tbl,         file.path(hipec_rutas_lista$tablas, "05_flexible_aic.csv"))


# __INDEX___ 12. Resumen ---------------------------------------------------------------------------
message(
  "[05] Listo | curvas reconstruidas: ", nrow(validacion_ipd_tbl),
  " | aceptables: ", sum(validacion_ipd_tbl$reconstruccion_aceptable),
  " | filas de RMST combinado: ", nrow(rmst_agrupado_tbl),
  " | modelos flexibles comparados: ", nrow(flexible_tbl)
)

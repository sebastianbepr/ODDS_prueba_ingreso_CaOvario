# ==================================================================================================
# 03_metaanalisis_bayesiano.R
# Metaanálisis bayesiano de efectos aleatorios con una prior half-normal débilmente informativa sobre tau. 
# Las muestras de la posterior se convierten en efectos absolutos por 1000 pacientes
# ==================================================================================================
#
# AVISO: PRUEBA TÉCNICA. Datos extraídos manualmente de tres artículos publicados; los resultados
# ilustran el método y NO constituyen un análisis formal ni una recomendación clínica.
#
# __INDEX___
#   1. Configuración e insumos
#   2. Priors (principal y de sensibilidad)
#   3. Ajustar un análisis y muestrear la posterior
#   4. Efectos absolutos y probabilidades por zona de decisión
#   5. Análisis principales
#   6. Sensibilidad a la prior
#   7. Figuras
#   8. Escritura de resultados
#   9. Resumen
# ==================================================================================================

# __INDEX___ 1. Configuración e insumos ------------------------------------------------------------
source(here::here("R", "00_configuracion.R"))

hazard_curado_tbl  <- readRDS(file.path(hipec_rutas_lista$derivados, "hazard_curado.rds"))
binario_curado_tbl <- readRDS(file.path(hipec_rutas_lista$derivados, "binario_curado.rds"))

posterior_n_draws <- 40000

message("[03] Insumos cargados; muestras de la posterior por análisis: ", posterior_n_draws)


# __INDEX___ 2. Priors -----------------------------------------------------------------------------
tau_priores_lista <- list(
  "Half-normal(0,5) [principal]" = function(t) bayesmeta::dhalfnormal(t, scale = 0.5),
  "Half-normal(0,25)"            = function(t) bayesmeta::dhalfnormal(t, scale = 0.25),
  "Half-normal(1,0)"             = function(t) bayesmeta::dhalfnormal(t, scale = 1.0),
  "Half-Cauchy(0,5)"             = function(t) bayesmeta::dhalfcauchy(t, scale = 0.5)
)
tau_prior_principal_etiqueta <- "Half-normal(0,5) [principal]"


# __INDEX___ 3. Ajustar un análisis y muestrear la posterior --------------------------------------
bayesiano_ajustar <- function(log_efecto_vector, se_vector, estudio_vector, tau_prior_funcion) {
  bayesmeta::bayesmeta(
    y             = log_efecto_vector,
    sigma         = se_vector,
    labels        = estudio_vector,
    mu.prior.mean = 0,
    mu.prior.sd   = 4,
    tau.prior     = tau_prior_funcion
  )
}

posterior_muestrear <- function(modelo_bayesmeta, log_efecto_vector, se_vector, n_draws,
                                mu_prior_media = 0, mu_prior_sd = 4, n_grilla = 4000) {
  tau_maximo <- modelo_bayesmeta$qposterior(tau.p = 0.9995)
  tau_grilla_vector <- seq(0, tau_maximo, length.out = n_grilla)
  tau_densidad_vector <- modelo_bayesmeta$dposterior(tau = tau_grilla_vector)
  tau_densidad_vector[!is.finite(tau_densidad_vector)] <- 0
  tau_paso <- tau_grilla_vector[2] - tau_grilla_vector[1]

  tau_draws_vector <- sample(tau_grilla_vector, n_draws, replace = TRUE,
                             prob = tau_densidad_vector) +
    runif(n_draws, -tau_paso / 2, tau_paso / 2)
  tau_draws_vector <- pmax(tau_draws_vector, 0)

  # Posterior condicional de mu dado tau (actualización normal conjugada)
  precision_condicional_vector <- 1 / mu_prior_sd^2 +
    vapply(tau_draws_vector, function(tau_valor) sum(1 / (se_vector^2 + tau_valor^2)),
           numeric(1))
  media_condicional_vector <- (mu_prior_media / mu_prior_sd^2 +
                                 vapply(tau_draws_vector, function(tau_valor) {
                                   sum(log_efecto_vector / (se_vector^2 + tau_valor^2))
                                 }, numeric(1))) / precision_condicional_vector

  mu_draws_vector <- rnorm(n_draws, media_condicional_vector,
                           1 / sqrt(precision_condicional_vector))
  tibble(
    tau             = tau_draws_vector,
    log_mu          = mu_draws_vector,
    log_theta_nuevo = rnorm(n_draws, mu_draws_vector, tau_draws_vector)
  )
}


# __INDEX___ 4. Efectos absolutos y zonas de decisión ------------------------------------------

# Distribución Beta ajustada a una supervivencia a tiempo fijo y su IC 95% (método de momentos)
basal_beta_muestrear <- function(n_draws, media, ic_inf, ic_sup) {
  basal_sd <- (ic_sup - ic_inf) / (2 * qnorm(0.975))
  basal_concentracion <- media * (1 - media) / basal_sd^2 - 1
  rbeta(n_draws, media * basal_concentracion, (1 - media) * basal_concentracion)
}

# Efecto absoluto por 1000: supervivencia con S1 = S0^HR; desenlaces binarios con p0 * (RR - 1)
efecto_absoluto_calcular <- function(log_ratio_vector, basal_vector, tipo) {
  if (tipo == "hazard") {
    1000 * (basal_vector - basal_vector^exp(log_ratio_vector))
  } else {
    1000 * basal_vector * (exp(log_ratio_vector) - 1)
  }
}

# Asigna cada efecto absoluto (negativo = menos eventos = beneficio) a una zona de decisión
zona_decision_asignar <- function(efecto_absoluto_vector, t_pequeno, t_moderado, t_grande) {
  magnitud_vector <- abs(efecto_absoluto_vector)
  sentido_vector  <- ifelse(efecto_absoluto_vector < 0, "Beneficio", "Daño")
  magnitud_etiqueta_vector <- case_when(
    magnitud_vector <  t_pequeno  ~ "Trivial",
    magnitud_vector <  t_moderado ~ "pequeño",
    magnitud_vector <  t_grande   ~ "moderado",
    TRUE                          ~ "grande"
  )
  zona_vector <- ifelse(magnitud_etiqueta_vector == "Trivial", "Trivial",
                        paste(sentido_vector, magnitud_etiqueta_vector))
  factor(zona_vector, levels = names(hipec_paleta_zonas_vector))
}


# __INDEX___ 5. Análisis principales ---------------------------------------------------------------
bayes_mapa_tbl <- tribble(
  ~analisis_id,     ~desenlace_id, ~tipo,     ~estudios,
  "os_todos",       "os",          "hazard",  c("CHIPOR", "MSK", "HORSE"),
  "os_cisplatino",  "os",          "hazard",  c("CHIPOR", "HORSE"),
  "pfs_todos",      "pfs",         "hazard",  c("CHIPOR", "MSK", "HORSE"),
  "pfs_cisplatino", "pfs",         "hazard",  c("CHIPOR", "HORSE"),
  "ea_g3_todos",    "ea_g3",       "binario", c("CHIPOR", "MSK", "HORSE")
)

# Insumos en escala log-razón por análisis
bayes_entrada_obtener <- function(desenlace_id, tipo, estudios) {
  if (tipo == "hazard") {
    hazard_curado_tbl %>%
      filter(desenlace_id == !!desenlace_id, estudio %in% estudios) %>%
      transmute(estudio = as.character(estudio), yi = log_hr, sei = se_log_hr)
  } else {
    metafor::escalc(
      measure = "RR", ai = eventos_hipec, n1i = n_hipec, ci = eventos_control, n2i = n_control,
      data = binario_curado_tbl %>% filter(desenlace_id == !!desenlace_id, estudio %in% estudios)
    ) %>%
      as_tibble() %>%
      transmute(estudio = as.character(estudio), yi = as.numeric(yi), sei = sqrt(as.numeric(vi)))
  }
}

# Muestras del riesgo basal por análisis (escenario principal)
bayes_basal_muestrear <- function(desenlace_id, tipo, n_draws) {
  if (tipo == "hazard") {
    basal_fila_tbl <- hipec_riesgo_basal_tbl %>%
      filter(desenlace_id == !!desenlace_id, escenario == "principal")
    basal_beta_muestrear(n_draws, basal_fila_tbl$supervivencia_control, basal_fila_tbl$ic_inf,
                         basal_fila_tbl$ic_sup)
  } else {
    control_tbl <- binario_curado_tbl %>% filter(desenlace_id == !!desenlace_id)
    rbeta(n_draws, sum(control_tbl$eventos_control) + 0.5,
          sum(control_tbl$n_control) - sum(control_tbl$eventos_control) + 0.5)
  }
}

bayes_analisis_ejecutar <- function(analisis_id, desenlace_id, tipo, estudios, tau_prior_funcion) {
  entrada_tbl <- bayes_entrada_obtener(desenlace_id, tipo, estudios)
  modelo_bayesmeta <- bayesiano_ajustar(entrada_tbl$yi, entrada_tbl$sei, entrada_tbl$estudio,
                                        tau_prior_funcion)
  umbrales_fila_tbl <- hipec_umbrales_tbl %>% filter(desenlace_id == !!desenlace_id)
  draws_tbl <- posterior_muestrear(modelo_bayesmeta, entrada_tbl$yi, entrada_tbl$sei,
                                   posterior_n_draws) %>%
    mutate(
      basal             = bayes_basal_muestrear(desenlace_id, tipo, posterior_n_draws),
      absoluto_mu       = efecto_absoluto_calcular(log_mu, basal, tipo),
      absoluto_prediccion = efecto_absoluto_calcular(log_theta_nuevo, basal, tipo),
      zona_mu           = zona_decision_asignar(absoluto_mu, umbrales_fila_tbl$t_pequeno,
                                                umbrales_fila_tbl$t_moderado,
                                                umbrales_fila_tbl$t_grande)
    )
  list(modelo = modelo_bayesmeta, draws = draws_tbl, analisis_id = analisis_id,
       desenlace_id = desenlace_id, tipo = tipo)
}

bayes_principal_lista <- bayes_mapa_tbl %>%
  { set_names(pmap(., bayes_analisis_ejecutar,
                   tau_prior_funcion = tau_priores_lista[[tau_prior_principal_etiqueta]]),
              .$analisis_id) }

bayes_resumen_tbl <- imap_dfr(bayes_principal_lista, function(ajuste_lista, analisis_id) {
  draws_tbl <- ajuste_lista$draws
  umbrales_fila_tbl <- hipec_umbrales_tbl %>% filter(desenlace_id == ajuste_lista$desenlace_id)
  tibble(
    analisis_id         = analisis_id,
    ratio_mediana       = exp(median(draws_tbl$log_mu)),
    ratio_cri_inf       = exp(quantile(draws_tbl$log_mu, 0.025)),
    ratio_cri_sup       = exp(quantile(draws_tbl$log_mu, 0.975)),
    ratio_prediccion_inf = exp(quantile(draws_tbl$log_theta_nuevo, 0.025)),
    ratio_prediccion_sup = exp(quantile(draws_tbl$log_theta_nuevo, 0.975)),
    tau_mediana         = median(draws_tbl$tau),
    absoluto_mediana    = median(draws_tbl$absoluto_mu),
    absoluto_cri_inf    = quantile(draws_tbl$absoluto_mu, 0.025),
    absoluto_cri_sup    = quantile(draws_tbl$absoluto_mu, 0.975),
    p_beneficio_cualquiera = mean(draws_tbl$absoluto_mu < 0),
    p_beneficio_al_menos_pequeno = mean(draws_tbl$absoluto_mu <= -umbrales_fila_tbl$t_pequeno),
    p_beneficio_al_menos_moderado = mean(draws_tbl$absoluto_mu <= -umbrales_fila_tbl$t_moderado),
    p_dano_al_menos_pequeno = mean(draws_tbl$absoluto_mu >= umbrales_fila_tbl$t_pequeno),
    p_dano_al_menos_moderado = mean(draws_tbl$absoluto_mu >= umbrales_fila_tbl$t_moderado)
  )
})

bayes_zonas_tbl <- imap_dfr(bayes_principal_lista, ~ .x$draws %>%
                              count(zona_mu, .drop = FALSE) %>%
                              mutate(analisis_id = .y, probabilidad = n / sum(n)))

message("[03] Análisis bayesianos principales ajustados: ", length(bayes_principal_lista))


# __INDEX___ 6. Sensibilidad a la prior ------------------------------------------------------------
bayes_sensibilidad_priores_tbl <- imap_dfr(
  tau_priores_lista,
  function(tau_prior_funcion, prior_etiqueta) {
    bayes_mapa_tbl %>%
      pmap_dfr(function(analisis_id, desenlace_id, tipo, estudios) {
        entrada_tbl <- bayes_entrada_obtener(desenlace_id, tipo, estudios)
        modelo_bayesmeta <- bayesiano_ajustar(entrada_tbl$yi, entrada_tbl$sei, entrada_tbl$estudio,
                                              tau_prior_funcion)
        tibble(
          analisis_id   = analisis_id,
          prior_tau     = prior_etiqueta,
          ratio_mediana = exp(modelo_bayesmeta$qposterior(mu.p = 0.5)),
          ratio_cri_inf = exp(modelo_bayesmeta$qposterior(mu.p = 0.025)),
          ratio_cri_sup = exp(modelo_bayesmeta$qposterior(mu.p = 0.975)),
          p_ratio_menor_1 = modelo_bayesmeta$pposterior(mu = 0)
        )
      })
  })

message("[03] Corridas de sensibilidad a la prior: ", nrow(bayes_sensibilidad_priores_tbl))


# __INDEX___ 7. Figuras ----------------------------------------------------------------------------
analisis_etiquetas_vector <- c(
  os_todos       = "SG, todos los platinos (3 ECA)",
  os_cisplatino  = "SG, cisplatino (2 ECA)",
  pfs_todos      = "SLP, todos los platinos (3 ECA)",
  pfs_cisplatino = "SLP, cisplatino (2 ECA)",
  ea_g3_todos    = "EA grado >=3 (3 ECA)"
)

# 7a. Probabilidad posterior apilada de cada zona de decisión: el target GRADE de umbrales múltiples
zonas_figura <- bayes_zonas_tbl %>%
  mutate(analisis = factor(analisis_etiquetas_vector[analisis_id],
                           levels = rev(analisis_etiquetas_vector))) %>%
  ggplot(aes(x = probabilidad, y = analisis, fill = zona_mu)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(probabilidad >= 0.06, scales::percent(probabilidad, 1), "")),
            position = position_stack(vjust = 0.5), size = 3.2) +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = hipec_paleta_zonas_vector, drop = FALSE, name = NULL) +
  labs(
    title = "¿Dónde está el efecto real? Probabilidad posterior por zona de decisión GRADE",
    subtitle = "Efectos aleatorios bayesianos, prior half-normal(0,5) sobre tau; efectos por 1000",
    x = "Probabilidad posterior", y = NULL, caption = hipec_aviso_prueba_texto
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# 7b. Densidad posterior del efecto absoluto sobre la mortalidad, con zonas sombreadas
densidad_os_tbl <- bind_rows(
  bayes_principal_lista$os_todos$draws %>% mutate(analisis = analisis_etiquetas_vector["os_todos"]),
  bayes_principal_lista$os_cisplatino$draws %>%
    mutate(analisis = analisis_etiquetas_vector["os_cisplatino"])
)
umbrales_os_tbl <- hipec_umbrales_tbl %>% filter(desenlace_id == "os")
umbrales_os_vector <- with(umbrales_os_tbl, c(-t_grande, -t_moderado, -t_pequeno,
                                              t_pequeno, t_moderado, t_grande))

densidad_os_figura <- ggplot(densidad_os_tbl, aes(x = absoluto_mu, fill = analisis)) +
  annotate("rect", xmin = -Inf, xmax = -umbrales_os_tbl$t_grande, ymin = -Inf, ymax = Inf,
           fill = hipec_paleta_zonas_vector["Beneficio grande"], alpha = 0.10) +
  annotate("rect", xmin = -umbrales_os_tbl$t_grande, xmax = -umbrales_os_tbl$t_moderado,
           ymin = -Inf, ymax = Inf, fill = hipec_paleta_zonas_vector["Beneficio moderado"],
           alpha = 0.10) +
  annotate("rect", xmin = -umbrales_os_tbl$t_moderado, xmax = -umbrales_os_tbl$t_pequeno,
           ymin = -Inf, ymax = Inf, fill = hipec_paleta_zonas_vector["Beneficio pequeño"],
           alpha = 0.15) +
  geom_density(alpha = 0.45, colour = NA, adjust = 1.2) +
  geom_vline(xintercept = umbrales_os_vector, linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = 0, colour = "grey20") +
  coord_cartesian(xlim = c(-300, 200)) +
  scale_fill_manual(values = set_names(c("#7A7A7A", "#1F5A8A"),
                                       analisis_etiquetas_vector[c("os_todos", "os_cisplatino")]),
                    name = NULL) +
  labs(
    title = "Distribución posterior del efecto sobre la mortalidad a 5 años",
    subtitle = paste("Muertes por 1000 (negativo = menos muertes con HIPEC);",
                     "incertidumbre del riesgo basal propagada"),
    x = "Diferencia de riesgo por 1000 a 5 años", y = "Densidad posterior",
    caption = hipec_aviso_prueba_texto
  )

# 7c. Sensibilidad a la prior
priores_figura <- bayes_sensibilidad_priores_tbl %>%
  mutate(analisis = factor(analisis_etiquetas_vector[analisis_id],
                           levels = rev(analisis_etiquetas_vector))) %>%
  ggplot(aes(x = ratio_mediana, y = analisis, colour = prior_tau)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = ratio_cri_inf, xmax = ratio_cri_sup),
                  position = position_dodge(width = 0.6)) +
  scale_x_log10() +
  labs(title = "Sensibilidad a la prior de heterogeneidad", x = "Razón combinada (escala log)",
       y = NULL, colour = "Prior sobre tau", caption = hipec_aviso_prueba_texto) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file.path(hipec_rutas_lista$figuras, "03_zonas_decision_bayes.png"), zonas_figura,
       width = 10, height = 4.4, dpi = 320, bg = "white")
ggsave(file.path(hipec_rutas_lista$figuras, "03_densidad_mortalidad.png"), densidad_os_figura,
       width = 8.5, height = 4.7, dpi = 320, bg = "white")
ggsave(file.path(hipec_rutas_lista$figuras, "03_sensibilidad_priores.png"), priores_figura,
       width = 8.5, height = 4.7, dpi = 320, bg = "white")


# __INDEX___ 8. Escritura de resultados ------------------------------------------------------------
saveRDS(map(bayes_principal_lista, "modelo"),
        file.path(hipec_rutas_lista$modelos, "03_modelos_bayesmeta.rds"))
saveRDS(bayes_resumen_tbl, file.path(hipec_rutas_lista$derivados, "bayes_resumen.rds"))

write_csv(bayes_resumen_tbl, file.path(hipec_rutas_lista$tablas, "03_bayes_resumen.csv"))
write_csv(bayes_zonas_tbl,                file.path(hipec_rutas_lista$tablas, "03_bayes_zonas.csv"))
write_csv(bayes_sensibilidad_priores_tbl, file.path(hipec_rutas_lista$tablas,
                                                    "03_bayes_sensibilidad_priores.csv"))


# __INDEX___ 9. Resumen ----------------------------------------------------------------------------
bayes_resumen_tbl %>%
  mutate(resultado = sprintf(
    "%-15s razón %.2f (ICr %.2f-%.2f) | abs %+.0f/1000 | P(beneficio >= pequeño) = %.2f",
    analisis_id, ratio_mediana, ratio_cri_inf, ratio_cri_sup, absoluto_mediana,
    p_beneficio_al_menos_pequeno
  )) %>%
  pull(resultado) %>%
  walk(~ message("[03] ", .x))

message("[03] Listo | figuras: 3 | escenarios de prior: ", length(tau_priores_lista))

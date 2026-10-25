# ==================================================================================================
# 06_riesgo_sesgo_rob2.R
# Visualización del riesgo de sesgo (RoB 2) por ensayo y desenlace: gráfico de semáforo y
# resumen por dominio. Los juicios provienen de la tabla de la sección 4 del protocolo.
# El semáforo se dibuja con ggplot
# ==================================================================================================
#
# AVISO: PRUEBA TÉCNICA. Juicios RoB 2 emitidos por un solo revisor sobre artículos publicados;
# ilustran el método y NO constituyen una evaluación formal ni una recomendación clínica.
#
# __INDEX___
#   1. Configuración
#   2. Juicios RoB 2 (dominios D1-D5 y global)
#   3. Gráfico de semáforo
#   4. Resumen por dominio
#   5. Escritura de resultados
#   6. Resumen
# ==================================================================================================
#
# Dominios RoB 2: D1 aleatorización  D2 desviaciones de la intervención D3 datos faltantes 
# D4 medición del desenlace  D5 selección del resultado reportado.


# __INDEX___ 1. Configuración ----------------------------------------------------------------------
source(here::here("R", "00_configuracion.R"))

rob_niveles_vector <- c("Bajo", "Algunas dudas", "Alto")
rob_paleta_vector  <- c("Bajo" = "#4E8F5A", "Algunas dudas" = "#E6B62E", "Alto" = "#C0504D")
rob_simbolos_vector <- c("Bajo" = "+", "Algunas dudas" = "-", "Alto" = "x")
rob_dominios_vector <- c(
  D1 = "D1\nAleatorización",
  D2 = "D2\nDesviaciones",
  D3 = "D3\nDatos faltantes",
  D4 = "D4\nMedición",
  D5 = "D5\nSelección\ndel reporte",
  global = "Global"
)

message("[06] Configuración de RoB 2 lista.")


# __INDEX___ 2. Juicios RoB 2 ----------------------------------------------------------------------
# Códigos: B = bajo, D = algunas dudas, A = alto
rob_codigos_vector <- c(B = "Bajo", D = "Algunas dudas", A = "Alto")

rob_juicios_tbl <- tribble(
  ~estudio, ~desenlace,        ~D1, ~D2, ~D3, ~D4, ~D5, ~global,
  "CHIPOR", "SG",              "B", "B", "B", "B", "B", "B",
  "CHIPOR", "SLP",             "B", "B", "B", "D", "B", "D",
  "CHIPOR", "EA grado >=3",    "B", "B", "B", "D", "B", "D",
  "CHIPOR", "Calidad de vida", "B", "B", "A", "D", "B", "A",
  "MSK",    "SG",              "B", "D", "B", "B", "D", "D",
  "MSK",    "SLP",             "B", "D", "B", "B", "D", "D",
  "MSK",    "EA grado >=3",    "B", "D", "B", "D", "B", "D",
  "HORSE",  "SG",              "D", "B", "B", "B", "D", "D",
  "HORSE",  "SLP",             "D", "B", "B", "D", "D", "D",
  "HORSE",  "EA grado >=3",    "D", "B", "B", "D", "B", "D"
) %>%
  mutate(across(c(D1, D2, D3, D4, D5, global), ~ unname(rob_codigos_vector[.x])))

rob_largo_tbl <- rob_juicios_tbl %>%
  mutate(fila = paste(estudio, "·", desenlace)) %>%
  pivot_longer(c(D1, D2, D3, D4, D5, global), names_to = "dominio", values_to = "juicio") %>%
  mutate(
    juicio  = factor(juicio, levels = rob_niveles_vector),
    dominio = factor(rob_dominios_vector[dominio], levels = rob_dominios_vector),
    fila    = factor(fila, levels = rev(unique(fila)))
  )

stopifnot(!anyNA(rob_largo_tbl$juicio))
message("[06] Juicios cargados: ", nrow(rob_juicios_tbl), " combinaciones ensayo x desenlace.")


# __INDEX___ 3. Gráfico de semáforo ----------------------------------------------------------------
rob_semaforo_figura <- ggplot(rob_largo_tbl, aes(x = dominio, y = fila, colour = juicio)) +
  geom_point(size = 8) +
  geom_text(aes(label = rob_simbolos_vector[as.character(juicio)]), colour = "white",
            fontface = "bold", size = 4.5) +
  geom_vline(xintercept = 5.5, colour = "grey60") +
  scale_colour_manual(values = rob_paleta_vector, drop = FALSE, name = NULL) +
  scale_x_discrete(position = "top") +
  labs(
    title = "Riesgo de sesgo por ensayo y desenlace (RoB 2)",
    x = NULL, y = NULL, caption = hipec_aviso_prueba_texto
  ) +
  theme(panel.grid.major = element_blank(),
        axis.text.x.top  = element_text(vjust = 0, lineheight = 0.9))


# __INDEX___ 4. Resumen por dominio ----------------------------------------------------------------
rob_resumen_tbl <- rob_largo_tbl %>%
  count(dominio, juicio, .drop = FALSE) %>%
  group_by(dominio) %>%
  mutate(proporcion = n / sum(n)) %>%
  ungroup()

rob_resumen_figura <- ggplot(rob_resumen_tbl,
                             aes(x = proporcion, y = forcats::fct_rev(dominio), fill = juicio)) +
  scale_y_discrete(labels = ~ str_replace_all(.x, "\n", " ")) +
  geom_col(width = 0.7, colour = "white") +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = rob_paleta_vector, drop = FALSE, name = NULL) +
  labs(
    title = "Resumen del riesgo de sesgo por dominio (10 combinaciones ensayo x desenlace)",
    x = "Proporción de juicios (sin ponderar)", y = NULL, caption = hipec_aviso_prueba_texto
  )


# __INDEX___ 5. Escritura de resultados ------------------------------------------------------------
ggsave(file.path(hipec_rutas_lista$figuras, "06_rob2_semaforo.png"), rob_semaforo_figura,
       width = 8.5, height = 5.8, dpi = 320, bg = "white")
ggsave(file.path(hipec_rutas_lista$figuras, "06_rob2_resumen.png"), rob_resumen_figura,
       width = 8.5, height = 3.8, dpi = 320, bg = "white")
write_csv(rob_juicios_tbl, file.path(hipec_rutas_lista$tablas, "06_rob2_juicios.csv"))


# __INDEX___ 6. Resumen ----------------------------------------------------------------------------
message(
  "[06] Listo | combinaciones evaluadas: ", nrow(rob_juicios_tbl),
  " | global bajo: ", sum(rob_juicios_tbl$global == "Bajo"),
  " | algunas dudas: ", sum(rob_juicios_tbl$global == "Algunas dudas"),
  " | alto: ", sum(rob_juicios_tbl$global == "Alto")
)

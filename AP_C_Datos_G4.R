# install.packages(c("tidyverse","igraph","tidygraph","ggraph","widyr",
#                    "proxy","cluster","factoextra","ggrepel","scales"))

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(widyr)
library(proxy)
library(cluster)
library(factoextra)
library(ggrepel)
library(scales)

# Ajustá las rutas a donde tengas los archivos
idps      <- read_csv("C:/CienciaDatos/Aplicar_Ciencia_Datos/idps.csv", show_col_types = FALSE)
bloques   <- read_csv("C:/CienciaDatos/Aplicar_Ciencia_Datos/legisladores_bloque_provincia.csv", show_col_types = FALSE)
sesiones  <- read_csv("C:/CienciaDatos/Aplicar_Ciencia_Datos/sesiones_dt.csv", show_col_types = FALSE)

#------------------------------------
#---Módulo 1 – Limpieza y unión------
#------------------------------------

# --- Función para normalizar nombres ---
# Quita comillas (simples y dobles), colapsa espacios, pasa a mayúsculas,
# y normaliza el formato "APELLIDO, NOMBRE"
norm_name <- function(x) {
  x %>%
    str_replace_all('"', "") %>%          # quita comillas dobles
    str_replace_all("'", "") %>%          # quita comillas simples
    str_squish() %>%                       # colapsa espacios múltiples
    str_to_upper()
}

idps <- idps %>%
  mutate(nombre = norm_name(nombre))

bloques <- bloques %>%
  mutate(
    nombre = norm_name(legislador),
    yr     = as.integer(yr)
  ) %>%
  select(nombre, yr, bloque, provincia)

sesiones <- sesiones %>%
  mutate(
    nombre = norm_name(legislador),
    yr     = as.integer(yr)
  )

# --- Unión idps + bloques ---
datos <- idps %>%
  rename(yr = anio) %>%
  left_join(bloques, by = c("nombre", "yr"))

# Chequeo de no-emparejados
faltantes <- datos %>% filter(is.na(bloque)) %>% distinct(nombre, yr)
print(faltantes, n = Inf)

# --- Estandarizar idp_pca POR AÑO (z-score dentro de cada año) ---
datos <- datos %>%
  group_by(yr) %>%
  mutate(idp_pca_z = as.numeric(scale(idp_pca))) %>%
  ungroup()

# --- Tratar a Massa 2022 (presidente de la Cámara) ---
# Lo marcamos para poder excluirlo de análisis de puntos ideales si querés
datos <- datos %>%
  mutate(presidente_camara = str_detect(nombre, "MASSA") & yr == 2022)

#--------------------------------------------------------------
#-----------Módulo 2 – Puntos ideales--------------------------
#--------------------------------------------------------------
#--------2.1 Dispersión por año, coloreado por bloque----------

# Para no saturar, mostramos los 8 bloques más grandes de cada año
top_bloques <- datos %>%
  count(yr, bloque, sort = TRUE) %>%
  group_by(yr) %>%
  slice_max(n, n = 8) %>%
  ungroup()

plot_df <- datos %>%
  semi_join(top_bloques, by = c("yr", "bloque")) %>%
  filter(!presidente_camara)

ggplot(plot_df, aes(x = idp_bayes, y = reorder(bloque, idp_bayes, median),
                    color = bloque)) +
  geom_jitter(width = 0, height = 0.25, alpha = 0.7, size = 1.8) +
  facet_wrap(~ yr, scales = "free_x") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none") +
  labs(x = "Punto ideal (Bayes)", y = NULL,
       title = "Puntos ideales por bloque y año")

#----------------------------------------------------------
#---------2.2 Comparación PCA vs Bayes (por año)-----------
#----------------------------------------------------------

datos %>%
  filter(!is.na(idp_pca), !presidente_camara) %>%
  group_by(yr) %>%
  summarise(cor_pearson = cor(idp_pca_z, idp_bayes, use = "complete.obs"),
            cor_spearman = cor(idp_pca_z, idp_bayes, method = "spearman",
                               use = "complete.obs"))

#----------------------------------------------------------
#---------------2.3 Incertidumbre--------------------------
#----------------------------------------------------------

# Legisladores peor identificados (sd_bayes alta)
datos %>%
  filter(!presidente_camara) %>%
  group_by(yr) %>%
  slice_max(sd_bayes, n = 5) %>%
  select(yr, nombre, bloque, idp_bayes, sd_bayes) %>%
  arrange(yr, desc(sd_bayes)) %>%
  print(n = 30)

#--------------------------------------------------------------
#----------2.4 Trayectorias de legisladores persistentes-------
#--------------------------------------------------------------

persistentes <- datos %>%
  filter(!presidente_camara) %>%
  count(nombre) %>%
  filter(n >= 3) %>% pull(nombre)

datos %>%
  filter(nombre %in% persistentes) %>%
  ggplot(aes(x = yr, y = idp_bayes, group = nombre)) +
  geom_line(alpha = 0.25) +
  geom_point(size = 1, alpha = 0.5) +
  theme_minimal() +
  labs(x = "Año", y = "Punto ideal (Bayes)",
       title = "Trayectorias de legisladores con ≥3 años")


#-------------------------------------------------------------
#--------------Módulo 3 – Distancias--------------------------
#-------------------------------------------------------------
#----------3.1 Matriz de distancias y análisis intra/inter bloque (por año)--------

# Función: distancias euclídeas en idp_bayes y clasificación intra/inter bloque
dist_bloques <- function(d) {
  d <- d %>% filter(!is.na(idp_bayes), !presidente_camara)
  if (nrow(d) < 5) return(NULL)
  m <- dist(d$idp_bayes, method = "euclidean") %>% as.matrix()
  b <- d$bloque
  idx <- expand_grid(i = seq_len(nrow(m)), j = seq_len(nrow(m))) %>%
    filter(i < j) %>%
    mutate(dist = m[cbind(i, j)],
           intra = b[i] == b[j])
  idx
}

distancias_yr <- datos %>%
  group_split(yr) %>%
  map_dfr(~ {
    d <- .x
    y <- unique(d$yr)
    r <- dist_bloques(d)
    if (is.null(r)) return(NULL)
    r %>% mutate(yr = y)
  })

distancias_yr %>%
  group_by(yr, intra) %>%
  summarise(n = n(), dist_media = mean(dist), .groups = "drop") %>%
  arrange(yr, intra)

#-----------------------------------------------------------
#-------------3.2 MDS por año (mapa 2D de afinidades)-------
#-----------------------------------------------------------

mds_plot <- function(d) {
  d <- d %>% filter(!is.na(idp_bayes), !presidente_camara)
  if (nrow(d) < 5) return(NULL)
  mds <- cmdscale(dist(d$idp_bayes), k = 2) %>% as.data.frame()
  bind_cols(d, setNames(mds, c("MDS1", "MDS2")))
}

mds_yr <- datos %>%
  group_split(yr) %>%
  map_dfr(mds_plot)

ggplot(mds_yr, aes(MDS1, MDS2, color = bloque)) +
  geom_point(alpha = 0.7, size = 1.8) +
  facet_wrap(~ yr, scales = "free") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "MDS de puntos ideales por año")

#-------------------------------------------------------------------
#-----------Módulo 4 – Redes de coincidencia------------------------
#-------------------------------------------------------------------
#----------4.1 Construir matriz de coincidencia (Jaccard) por año---

# Jaccard por pares a partir de una matriz legislador x votación (0/1, con NA)
pairwise_jaccard_manual <- function(mat, min_shared = 10) {
  # mat: matriz numérica con filas = legisladores, columnas = votaciones
  #      valores 0/1, NA donde no votó o voto no válido
  # min_shared: mínimo de votaciones en común para considerar el par
  
  n <- nrow(mat)
  if (n < 2) return(tibble())
  
  # Para cada par (i,j):
  #   a = ambos 1
  #   b = i=1, j=0
  #   c = i=0, j=1
  #   d = ambos 0
  #   Jaccard = a / (a + b + c)   (ignorando d, como es habitual en votos)
  #   shared = a + b + c          (votos válidos en ambos)
  
  res <- vector("list", n * (n - 1) / 2)
  k <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      xi <- mat[i, ]
      xj <- mat[j, ]
      ok <- !is.na(xi) & !is.na(xj)
      shared <- sum(ok)
      if (shared < min_shared) next
      a <- sum(xi[ok] == 1 & xj[ok] == 1)
      b <- sum(xi[ok] == 1 & xj[ok] == 0)
      c <- sum(xi[ok] == 0 & xj[ok] == 1)
      jac <- if ((a + b + c) > 0) a / (a + b + c) else NA_real_
      k <- k + 1
      res[[k]] <- tibble(
        item1    = rownames(mat)[i],
        item2    = rownames(mat)[j],
        jaccard  = jac,
        shared   = shared
      )
    }
  }
  bind_rows(res)
}

validos <- c("AFIRMATIVO", "NEGATIVO")

sesiones_v <- sesiones %>%
  filter(voto %in% validos) %>%
  mutate(voto_bin = if_else(voto == "AFIRMATIVO", 1L, 0L))

jaccard_yr <- function(d, min_shared = 10) {
  mat <- d %>%
    select(nombre, id_votacion, voto_bin) %>%
    distinct() %>%
    pivot_wider(names_from = id_votacion, values_from = voto_bin) %>%
    column_to_rownames("nombre") %>%
    as.matrix()
  if (nrow(mat) < 5) return(NULL)
  pairwise_jaccard_manual(mat, min_shared = min_shared)
}

redes <- sesiones_v %>%
  group_split(yr) %>%
  map_dfr(~ {
    y <- unique(.x$yr)
    r <- jaccard_yr(.x, min_shared = 10)
    if (is.null(r) || nrow(r) == 0) return(NULL)
    r %>% mutate(yr = y)
  })

head(redes)

#----------------------------------------------------------
#-----------4.2 Grafo por año y detección de comunidades---
#----------------------------------------------------------

# ---- 4.2 Red de coincidencia 2024 ----

grafo_yr <- function(d, umbral = 0.5) {
  d <- d %>% filter(!is.na(jaccard), jaccard >= umbral)
  cat("Aristas tras filtro:", nrow(d), "\n")
  if (nrow(d) < 5) {
    warning("Menos de 5 aristas: devuelvo NULL")
    return(NULL)
  }
  g <- as_tbl_graph(
    d %>% rename(from = item1, to = item2, weight = jaccard),
    directed = FALSE
  ) %>%
    mutate(comunidad = as.character(group_louvain(weight)))
  g
}

# Diagnóstico
redes %>% filter(yr == 2024) %>%
  summarise(n = n(),
            med = median(jaccard, na.rm = TRUE),
            max = max(jaccard, na.rm = TRUE))

# Loop por umbrales
g2024 <- NULL; u_ok <- NA
for (u in c(0.6, 0.5, 0.4, 0.3, 0.2)) {
  g <- redes %>% filter(yr == 2024) %>% grafo_yr(umbral = u)
  if (!is.null(g) && vcount(g) > 0) {
    cat("Grafo OK con umbral", u, "- nodos:", vcount(g),
        "aristas:", ecount(g), "\n")
    g2024 <- g; u_ok <- u
    break
  }
}

# Dibujo (con print() y u_ok seguro)
if (!is.null(g2024)) {
  u_txt <- if (is.na(u_ok)) "?" else as.character(u_ok)
  print(
    ggraph(g2024, layout = "fr") +
      geom_edge_link(aes(alpha = weight), color = "grey70") +
      geom_node_point(aes(color = comunidad), size = 2) +
      theme_void() +
      labs(title = paste("Red de coincidencia de votos 2024 (Jaccard >=", u_txt, ")"))
  )
} else {
  message("g2024 es NULL: revisá el umbral o min_shared")
}
#----------------------------------------------------------
#-------------4.3 Centralidad y puentes entre bloques------
#----------------------------------------------------------

cent <- g2024 %>%
  mutate(degree = centrality_degree(),
         between = centrality_betweenness()) %>%
  as_tibble() %>%
  arrange(desc(between)) %>%
  slice_head(n = 20)

cent

#----------------------------------------------------------
#---------Módulo 5 – Comparar red vs. puntos ideales-------
#----------------------------------------------------------

# Unir puntos ideales 2024 con centralidad de la red
pi_2024 <- datos %>%
  filter(yr == 2024, !presidente_camara) %>%
  select(nombre, bloque, idp_bayes)

comp <- pi_2024 %>%
  inner_join(cent %>% select(name, degree, between), by = c("nombre" = "name"))

ggplot(comp, aes(idp_bayes, between, color = bloque)) +
  geom_point() +
  geom_text_repel(aes(label = nombre), size = 2.2, max.overlaps = 30) +
  theme_minimal() +
  labs(x = "Punto ideal (Bayes)", y = "Betweenness en red de votos",
       title = "Puentes entre bloques vs. ubicación ideológica (2024)")

#-------------------------------------------------------------
#-------------------------------------------------------------

# ============================================================
# Guardado de objetos y figuras para la presentación
# ============================================================
dir.create("output", showWarnings = FALSE)

# Calcular tabla_com si no existe
if (!exists("tabla_com")) {
  library(tidygraph)
  tabla_com <- g2024 %>%
    as_tibble() %>%
    left_join(
      datos %>% filter(yr == 2024) %>% select(nombre, bloque),
      by = c("name" = "nombre")
    ) %>%
    count(comunidad, bloque) %>%
    pivot_wider(names_from = bloque, values_from = n, values_fill = 0)
}

saveRDS(datos,     "output/datos.rds")
saveRDS(redes,     "output/redes.rds")
saveRDS(g2024,     "output/g2024.rds")
saveRDS(cent,      "output/cent.rds")
saveRDS(tabla_com, "output/tabla_com.rds")

cat("Archivos guardados en output/:\n")
print(list.files("output"))

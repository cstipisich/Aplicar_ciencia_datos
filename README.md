# Análisis de Votaciones - Diputados Argentina (2022-2026)

Este repositorio contiene el análisis de puntos ideales, distancias y redes de coincidencia de los diputados nacionales.

## Archivos principales
* `AP_C_Datos_G4.R`: Script principal con el procesamiento de datos y generación de gráficos.
* `Presentacion_Votaciones_2022_2026.qmd`: Código fuente de la presentación en Quarto (RevealJS).
* `Presentacion_Votaciones_2022_2026.html`: Presentación final compilada.
* Carpeta `output/`: Contiene los objetos guardados (.rds) y las figuras (.png) generadas.

## Metodología
1. Limpieza y unión de datasets (puntos ideales, bloques y sesiones).
2. Análisis de puntos ideales (Bayes y PCA).
3. Cálculo de distancias intra e inter bloque.
4. Redes de coincidencia mediante índice de Jaccard y detección de comunidades (Louvain).

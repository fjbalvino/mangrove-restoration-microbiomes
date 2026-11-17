# Análisis 611 de tendencias funcionales y LOLO

Paquete reproducible para las 51 muestras y 17 perfiles sedimentarios de Carmen, Cozumel y Tuxpan. Implementaciones independientes en Python y R base, ejecutadas y comparadas. Incluye todas las funciones KO y PFAM disponibles, no solo el core o las principales contribuciones al CLR.

## Resultado de esta ejecución

- **Análisis principal CLR:** ninguna función con q BH ≤ 0.05 en ninguna profundidad. No demuestra equivalencia ni ausencia de respuesta.
- **Sensibilidad con proporciones:** K05825 creciente a 40 cm, S = 1, P exacta = 1.286008230452675e-7, q BH = 0.003238040123456790, q BY ≈ 0.034683. Mantiene dirección en los tres LOLO y las tres localidades excluidas; no alcanza selección FDR en ningún entrenamiento. No es un hallazgo principal CLR ni una validación externa.
- **Exploración descriptiva:** 141 registros KO × profundidad (139 KO únicos) y 117 PFAM × profundidad (116 PFAM únicos) con medias monótonas y dirección estable. Ninguno es un descubrimiento principal FDR.
- El UpSet inferencial está vacío. El UpSet descriptivo está separado, sin filtro de P, y no reemplaza el criterio inferencial.

## Ejecutar en landalab

Descomprima este ZIP. La carpeta `611_tendencias_monotonicas_LOLO` contiene el lanzador. Desde ella:

```bash
bash REPLICAR_611_EN_R.sh /home/fjbalvino/Tipping_points/resultados_finales
```

El lanzador usa `Rscript` del entorno activo. No necesita `conda run`, no instala paquetes, no cierra la terminal y no sobrescribe carpetas. Ejecuta secuencialmente el análisis exacto, todas las sensibilidades, todos los LOLO y las figuras. Genera una carpeta con sufijo único, `611_log.txt` y `results/STATUS.txt`.

Si el entorno activo no ofrece Rscript, especifique una ruta existente:

```bash
RSCRIPT_611=/ruta/real/a/Rscript bash REPLICAR_611_EN_R.sh /home/fjbalvino/Tipping_points/resultados_finales
```

También puede ejecutar directamente:

```bash
Rscript code/611_run_all.R --inputs inputs --out /ruta/nueva/611_resultados
```

`--out` debe ser una carpeta que todavía no exista. Solo requiere R base; fue verificado con R 4.6.0 mediante WebR en esta sesión. WebR y Node no son necesarios en landalab. El PDF de R es `results/figures/611_figures_R.pdf`; reproduce datos, selecciones y paneles con un estilo de R base. Los originales de alta resolución de esta ejecución están en `results_python/figures/`.

El código está configurado para este diseño (cinco estadios, tres localidades, perfiles 6/6/5). Sus validaciones detienen diseños incompatibles. No es una herramienta general para cualquier diseño experimental sin revisión metodológica.

## Ejecutar en Python

Requiere Python 3.10 o posterior con `numpy`, `pandas` y `matplotlib`. Versiones usadas en `verification/python_environment.json`.

```bash
python code/611_run_all.py --inputs inputs --out /ruta/nueva/611_resultados_python
```

Verificar resultados suministrados de ambos lenguajes:

```bash
python code/611_verify.py
```

Las pruebas independientes de la distribución exacta pueden ejecutarse con `--unit-only`. El verificador completo usa por defecto las carpetas `results_python/` y `results_R/` de este paquete. Para comparar una réplica, coloque sus resultados en una copia del paquete con esos nombres, preservando los originales.

## Entradas y unidades

| Archivo | Contenido |
|---|---|
| `inputs/metadata.tsv` | 51 muestras; IDs de muestra y perfil, localidad, condición y profundidad canónica |
| `inputs/KO.tsv.gz` | 51 filas × 8,393 KO, más ID; abundancias no negativas fraccionales |
| `inputs/PFAM.tsv.gz` | 51 filas × 7,997 PFAM, más ID; abundancias no negativas fraccionales |
| `inputs/feature_annotations.tsv` | 13 KO ilustrados, etiqueta breve, anotación y fuente oficial |
| `protocol/611_protocolo_aprobado.md` | Protocolo conservado antes de la ejecución |

La profundidad se toma de metadata, **nunca del sufijo de muestra**. Hay nueve muestras conservadas, tres por profundidad. Dentro de cada localidad y profundidad hay un perfil por estadio salvo dos intermedios en Carmen y Cozumel. Los tres estratos pertenecen a los mismos perfiles.

El core de referencia exige **proporción ≥ 0.0001, equivalente a 0.01%, en 9/9 conservados**. La prevalencia es 100%. El umbral no es una proporción de 0.01 (= 1%). El core por profundidad exige 3/3. Estas etiquetas se añaden después y no seleccionan el universo analizado.

Se normaliza cada fila dentro de su ontología. CLR principal: `log2(p + 5e-7) - mean(log2(p + 5e-7))`. El centro incluye toda la ontología. Esta regularización de proporciones difiere del pseudoconteo de una unidad sobre abundancias originales utilizado en algunos análisis previos; aquí se declaran y evalúan 0.1, 0.5 y 1 ppm.

## Contraste e inferencia

Orden: degradado, temprano, intermedio, avanzado, conservado. La condición conservada es un referente ecológico, no una fecha futura medida. Los análisis son separados a 5, 20 y 40 cm; cada uno usa 17 perfiles.

Para cada localidad, S es el promedio de diez concordancias entre pares de estadios. Cada concordancia es la media de `sign(valor_estadio_posterior - valor_estadio_anterior)` de todos sus pares de perfiles. S completa es la media de las tres localidades. Es una adaptación ordinal de tipo Jonckheere–Terpstra con ponderaciones explícitas, no el estadístico estándar de una llamada convencional.

Se enumeran todas las asignaciones de etiquetas posibles dentro de cada localidad. La convolución de sus distribuciones obtiene exactamente la distribución conjunta: 15,552,000 asignaciones en cinco estadios, 86,400 en cuatro. P bilateral incluye empates y cuenta `abs(S_perm) >= abs(S_observada)`. No hay muestreo Monte Carlo, semilla aleatoria ni límite de 9,999 permutaciones.

Para ser contrastable se requiere >0 en al menos tres perfiles y variación. Se conservan todos los registros; los no contrastables tienen P/q NA, pero cuentan como P=1 dentro del ajuste. La familia de BH/BY abarca funciones × tres profundidades dentro de cada ontología: 25,179 KO y 23,991 PFAM. Cada sensibilidad se ajusta por separado y no constituye una búsqueda de la menor q entre métodos.

La monotonía descriptiva exige que las cuatro diferencias de las medias ajustadas tengan el signo de S, permitiendo mesetas con tolerancia 1e-10. Las medias ajustadas dan el mismo peso a cada localidad. El ajuste isotónico por PAVA es descriptivo: no prueba que los datos sean monótonos.

LOLO omite una localidad, reevalúa elegibilidad, calcula dirección/prueba/q en las dos restantes y evalúa después la concordancia de signo en la localidad excluida. Se mantiene el universo completo en cada familia FDR de entrenamiento. Una localidad excluida con S=0 o no contrastable no cuenta como confirmación. Los tres entrenamientos se solapan. Se informa por separado estabilidad del signo y selección FDR transferida; no se exige significación en todos los pliegues para describir estabilidad direccional.

Sensibilidades: CLR 0.1 ppm, CLR 1 ppm, proporciones originales, CLR principal sin conservados. El análisis sin conservados **no** es una prueba de proporciones sin conservados. Todas incluyen LOLO de todo el universo. Las comparaciones entre profundidades son descriptivas; no se ha contrastado una interacción vertical.

## Archivos de resultados

| Archivo en cada carpeta de resultados | Uso |
|---|---|
| `611_main_results.tsv.gz` | 49,170 registros función × profundidad, completos, incluyendo NA |
| `611_sensitivity_results.tsv.gz` | 196,680 registros de cuatro sensibilidades |
| `611_LOLO_results.tsv.gz` | 147,510 resultados de entrenamiento y localidad excluida |
| `611_sensitivity_LOLO.tsv.gz` | 590,040 resultados LOLO de sensibilidades |
| `611_stage_means.tsv.gz` | Cinco medias ajustadas y ajuste isotónico por función y profundidad |
| `611_selected_functions.tsv` | Selección BH principal; en esta ejecución solo tiene encabezados |
| `611_summary.tsv` | Resumen de contraste y selección |
| `report_tables/611_all_results_annotated.tsv.gz` | Tabla principal más descripciones disponibles, sin inventar las faltantes |
| `report_tables/611_descriptive_monotone_candidates.tsv` | Candidatos puramente descriptivos, no descubrimientos FDR |
| `report_tables/611_UpSet_*.tsv` | Conjuntos, intersecciones exclusivas e integrantes; patrones 5/20/40 cm |
| `report_tables/611_direction_changes_across_depths_descriptive.tsv` | S positiva y negativa entre profundidades contrastables; no prueba de interacción |
| `report_tables/611_K05825_*.tsv` | Todas las observaciones, escalas y LOLO de K05825 |

Columnas clave: `eligible`, `n_positive`, `zero_fraction`; `S` y `S_Carmen/Cozumel/Tuxpan`; `p_exact`, `q_BH`, `q_BY`; `mean_CLR_0` a `mean_CLR_4`; `adjacent_difference_1` a `adjacent_difference_4`; `endpoint_difference` (C−D en la escala del análisis); `endpoint_relative_difference` (C−D en proporciones, no porcentajes); `core_reference_9of9`, `core_depth_3of3`; `monotone_observed`, `direction_stable`, `stable_monotone_candidate`. Las sensibilidades de proporciones usan columnas genéricas de efecto y los registros por variante; las medias de la tabla principal son CLR.

`training_selected_n` cuenta entrenamientos con q BH ≤0.05. `training_selected_transferred_n` añade concordancia direccional en la localidad excluida. `lolo_sign_matches_n` y `heldout_matches_training_n` describen signos, no significación. `zero_selection_stable=TRUE` puede significar ausencia de selección en todas las regularizaciones; no debe confundirse con detección robusta.

## Figuras e interpretación

Se entregan mapas de calor KO/PFAM, comparación de escalas de K05825, ejemplos LOLO y UpSet inferencial/descriptivo. Los mapas ilustran dos menores P por signo y profundidad, sin filtro de significación; LOLO ilustra la menor P por ontología y profundidad. Las tablas conservan la selección exacta y las anotaciones core. Los dominios PFAM mantienen sus nombres originales. No se realizó enriquecimiento de rutas.

No se dibujaron intervalos mediante remuestreo ingenuo de celdas con un solo perfil. Las líneas conectan medias de perfiles distintos, no individuos seguidos en el tiempo. S no es porcentaje de recuperación, pendiente temporal ni velocidad.

## Límites para el manuscrito

Las pruebas exactas dependen de intercambiabilidad bajo la hipótesis nula e independencia entre perfiles condicionada al diseño. Los datos no permiten verificar una unidad de intervención independiente distinta del perfil. Si varios perfiles comparten parcela o intervención, habrá que revisar la inferencia a esa escala. Bloquear por localidad no elimina otros factores de confusión.

No detectar funciones tras FDR no equivale a demostrar estabilidad, redundancia o falta de recuperación. Un core presente en degradados no informa por sí mismo de actividad metabólica. Redundancia requiere enlaces taxón–función; recuperación de procesos exige mediciones funcionales; desfases temporales exigen seguimiento. Los párrafos propuestos están en el Word, sin modificar el manuscrito original.

## Verificación y procedencia

`verification/Python_R_validation.json` documenta concordancia de todas las pruebas, clasificadores y medias; máxima diferencia numérica 5.51e-13. `verification/figure_table_parity.json` confirma las tablas gráficas. Pruebas independientes contrastaron fórmulas, empates, casos constantes, convolución y FDR. El redondeo para ordenar es de diez decimales CLR y quince proporciones.

Fuente de entradas: `610_resultados_core100_20260923.zip`, versión 1. Protocolo aprobado, versión 0. La integridad de entradas se contrastó con el manifiesto original. `SHA256SUMS.txt` permite verificar este paquete. Los registros de ejecución están en `verification/logs/`.

Las fuentes metodológicas y anotaciones oficiales están en el protocolo, el informe y `feature_annotations.tsv`. No se atribuye al contraste estándar la adaptación ponderada implementada en este paquete.

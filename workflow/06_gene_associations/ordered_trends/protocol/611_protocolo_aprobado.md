# Protocolo 611: tendencias funcionales ordenadas y estabilidad entre localidades

Fecha: 23 de septiembre de 2026. Estado: diseño previo a ejecución. No contiene nuevos resultados biológicos ni modifica el manuscrito.

## 1. Pregunta y relación con el análisis 610

¿Qué funciones muestran una tendencia relativa creciente o decreciente en el orden Degraded → Early restoration → Intermediate restoration → Advanced restoration → Conserved, a cada profundidad, y cuáles conservan esa dirección al excluir una localidad?

El análisis 610 muestra que gran parte del repertorio funcional del core conservado también se detecta en degradados. Las medianas de cobertura son 99,52% para KO y 99,70% para PFAM; la abundancia representada es 80,07% y 82,57%, respectivamente. Eso motiva estudiar redistribución de abundancias, además de presencia. La distinta proximidad funcional al referente entre localidades es descriptiva y motiva evaluar consistencia geográfica; no equivale todavía a una interacción localidad × estadio demostrada.

Una tendencia positiva significa aumento de representación relativa según el orden ecológico propuesto. Una negativa significa disminución. Ninguna dirección es, por sí sola, beneficiosa o perjudicial. Los datos transversales no identifican velocidad, desfase temporal, causalidad de la intervención, actividad metabólica o recuperación del funcionamiento ecosistémico.

Conservado es una referencia ecológica, no un estadio temporal observado después de avanzada. El análisis principal incluirá los cinco estados solicitados y una sensibilidad D–E–I–A indicará cuánto depende la señal de añadir la referencia conservada.

## 2. Diseño que debe respetarse

Por cada profundidad hay 17 perfiles independientes, distribuidos así:

| Localidad | D | E | I | A | C | Total |
|---|---:|---:|---:|---:|---:|---:|
| Carmen | 1 | 1 | 2 | 1 | 1 | 6 |
| Cozumel | 1 | 1 | 2 | 1 | 1 | 6 |
| Tuxpan | 1 | 1 | 1 | 1 | 1 | 5 |
| Total | 3 | 3 | 5 | 3 | 3 | 17 |

Se ejecutarán análisis separados a 5, 20 y 40 cm. Dentro de cada uno, cada perfil aporta una sola muestra. La localidad será un bloque fijo; no se permutarán muestras entre localidades. Las tres profundidades pertenecen a los mismos perfiles y no constituyen replicaciones independientes de una asociación encontrada a varias profundidades. Se usarán las profundidades de la metadata canónica.

La intercambiabilidad de los perfiles dentro de localidad bajo la hipótesis nula es un supuesto observacional; la restauración no fue asignada aleatoriamente. Antes de ejecutar, comprobar correspondencia entre profile_id y sitio de intervención. Si varios perfiles son submuestras del mismo sitio, el sitio deberá sustituir al perfil como unidad y se recalcularán el diseño y las permutaciones. El diseño propuesto es condicional a la independencia de los 17 perfiles.

## 3. Universo y tratamiento de las abundancias

KO será la capa principal; PFAM será una comprobación con otra anotación de los mismos metagenomas. El universo recibido contiene 8393 KO y 7997 PFAM. No se seleccionarán funciones por contribución al CLR, por P previo, por dirección observada o por pertenencia al core.

En cada profundidad, una función debe tener abundancia original positiva en al menos tres perfiles para entrar al contraste de abundancia. Este filtro es previo al uso de los estadios y no exige presencia en varias localidades: una función geográficamente restringida debe poder aparecer como tal. Las demás permanecerán en la tabla, marcadas como insuficientemente observadas. Los elementos constantes no se contrastarán. Para ajustar multiplicidad se mantendrá la familia completa función × profundidad, asignando P=1 internamente a las no contrastables y conservando NA como resultado visible.

El filtro de tres observaciones no es una definición de core ni usa el umbral de 0,01%. La etiqueta core se añadirá después: abundancia relativa ≥0,01% en 9/9 conservados. También se podrá mostrar core por profundidad, 3/3 conservados, como anotación secundaria. Seleccionar previamente solo funciones que cumplen un criterio en conservados utilizaría el estadio para seleccionar hipótesis. Por ello no se recalculará FDR únicamente dentro del core a posteriori.

Transformación principal, fijada antes de observar asociaciones:

1. Dividir cada fila por su suma dentro de KO o PFAM usando todo el universo disponible de esa capa.
2. Para trabajar con proporciones y abundancias fraccionales, usar `z_ij = log2(p_ij + delta) − mean_k[log2(p_ik + delta)]`, con `delta = 5e−7`, equivalente a 0,5 partes por millón de abundancia relativa.
3. Mantener el mismo universo denominador en todas las profundidades y pliegues. No renormalizar únicamente las funciones seleccionadas, significativas o del core. La transformación es por muestra y no aprende parámetros de la localidad excluida.
4. Repetir como sensibilidad con `delta = 1e−7` y `1e−6`. Estos son valores operativos para evaluar la dependencia del tratamiento de ceros, no límites de detección validados. Auditar cuánto cambia el signo, el tamaño del efecto y la selección.
5. Mostrar abundancia relativa original y proporción de ceros junto al CLR. No convertir las abundancias funcionales fraccionales en conteos enteros ficticios.

Esta especificación no es idéntica al pseudoconteo 1 sobre la escala original del análisis 610. La nueva escala se declara explícitamente para no confundir efectos de unidades con respuestas biológicas. Un cambio en CLR es un cambio en relación con la media geométrica de las funciones, no una estimación de abundancia absoluta por gramo de sedimento.

## 4. Contraste principal de tendencia

Usar un contraste no paramétrico bilateral de concordancia ordinal, construido a partir de comparaciones tipo Jonckheere–Terpstra, con permutaciones exactas dentro de localidad. No aplicar directamente una prueba asintótica que ignore el bloqueo.

Para cada función y profundidad, y cada pareja de estadios a < b en una localidad l:

`C_l,ab = mean[sign(z_j − z_i)]`, con i perteneciente a a y j a b.

Los empates aportan cero. Se promedian todas las parejas de muestras de esos dos estadios. Luego:

`S_l = (1/10) × sum_{a<b} C_l,ab`

`S = (S_Carmen + S_Cozumel + S_Tuxpan)/3`

S está entre −1 y +1. Un valor positivo indica predominio de incrementos entre estadios ordenados; uno negativo, predominio de decrementos. Cada una de las diez parejas de estadios tiene igual peso y cada localidad aporta un tercio, evitando que las réplicas intermedias de Carmen/Cozumel dominen la señal. Esta ponderación es una adaptación explícita al diseño; no debe etiquetarse como el estadístico estándar sin modificar de un paquete.

La hipótesis nula de permutación es ausencia de asociación entre estadio y distribución de la abundancia transformada dentro de localidad. El contraste bilateral usa `P = Pr_perm(|S_perm| >= |S_observado|)`; así se buscan ambas direcciones sin escoger primero el P unilateral más pequeño. La interpretación como desplazamiento de abundancia requiere considerar cambios de dispersión y forma de distribución, difíciles de estimar con la replicación disponible.

### Distribución exacta

Con los cinco estadios y el diseño completo, hay 6!/2! = 360 asignaciones distintas de etiquetas en Carmen, 360 en Cozumel y 5! = 120 en Tuxpan. El espacio conjunto tiene 15 552 000 asignaciones. No es necesario evaluarlas una por una por función: se calcula la distribución discreta local y se combinan sus frecuencias mediante convolución.

En este diseño, `20 × S_l` es entero, incluso con empates, porque las comparaciones entre estadios tienen uno o dos pares. Las distribuciones de ese entero pueden convolucionarse sin aproximar P por una distribución normal ni introducir un límite de Monte Carlo de 9999 permutaciones. La asignación observada está incluida; no se añade la corrección +1 propia del muestreo Monte Carlo a una enumeración exhaustiva. La implementación debe verificar que las frecuencias suman exactamente el número de asignaciones y que el resultado coincide con enumeración directa en ejemplos pequeños, con y sin empates.

## 5. Separar tendencia y forma monótona

Una concordancia ordinal significativa no prueba que cada transición sea creciente o decreciente. Para cada función se calculará también la media CLR de cada estadio con igual peso por localidad. Esta equivale al efecto ajustado por localidad de un modelo aditivo con pesos iguales para las celdas localidad × estadio.

Se conservarán:

- Las cinco medias ajustadas sin restricciones y sus cuatro diferencias adyacentes.
- Un ajuste isotónico creciente o decreciente según la dirección global, permitiendo mesetas. Es una descripción de forma; su curva forzada no constituye una segunda demostración de monotonía.
- La desviación de las medias libres respecto al ajuste isotónico y la mayor inversión adyacente.
- El contraste descriptivo conservado − degradado en CLR y en proporción original.

La etiqueta «perfil ajustado monótono observado» exige que las cuatro diferencias de las medias libres tengan el signo previsto o sean cero, con tolerancia exclusivamente numérica. No exige cuatro pruebas pareadas significativas ni convierte esa condición descriptiva en prueba de monotonía poblacional. Si hay una inversión, se conserva la asociación como «tendencia ordinal con inversión»; si aparece un pico o valle, no se fuerza su interpretación como recuperación monótona.

## 6. LOLO: influencia y transferencia de dirección

| Pliegue | Localidades de entrenamiento | n por profundidad | Localidad excluida | n validación |
|---|---|---:|---|---:|
| Sin Carmen | Cozumel + Tuxpan | 11 | Carmen | 6 |
| Sin Cozumel | Carmen + Tuxpan | 11 | Cozumel | 6 |
| Sin Tuxpan | Carmen + Cozumel | 12 | Tuxpan | 5 |

Se calcularán dos productos distintos:

1. **Influencia LOLO:** recalcular S y su P exacto con las dos localidades restantes. Registrar signo y cambio de magnitud respecto al ajuste completo. Las tres exclusiones se realizan para todas las funciones, no solo para las significativas.
2. **Evaluación de dirección excluida:** fijar la dirección y cualquier selección usando únicamente las dos localidades de entrenamiento; evaluar S_l en la localidad excluida sin reajustar esa dirección. La localidad excluida no escoge pseudoconteo, funciones, umbrales ni orientación.

El criterio mínimo de presencia se aplica a los datos de entrenamiento para una selección dentro del pliegue. Una función no contrastable en entrenamiento se marca así, no se rescata mirando la localidad excluida. Si se presenta etiqueta core en un análisis de entrenamiento, se define con sus seis conservados (o dos conservados a la profundidad analizada), no con los nueve de todas las localidades.

Las tablas mostrarán signo completo, tres signos LOLO, tres signos locales excluidos, P y q de entrenamiento y cambios de tamaño de efecto. No se exigirá q<0,05 en todos los pliegues como condición general, dado el tamaño de 11–12 perfiles, pero se informará la frecuencia con que la función se selecciona en entrenamiento. Un signo local cero no se contará como validación positiva.

Los pliegues se solapan y solo existen tres localidades. Este procedimiento evalúa influencia y transferencia interna de dirección; no es una validación externa ni permite combinar los tres P como si fueran independientes. Para hablar de descubrimiento transferido, deberá distinguirse la selección realizada en entrenamiento de la selección con los datos completos.

## 7. Multiplicidad, etiquetas y tamaños de efecto

Una familia principal incluirá todas las combinaciones KO × tres profundidades (25179 entradas antes de marcar no contrastables). PFAM tendrá su familia propia (23991 entradas), explícitamente secundaria. Se calculará BH, umbral principal q≤0,05, y BY como sensibilidad a dependencia entre funciones. PFAM no constituye replicación biológica independiente.

La anotación core no crea una familia nueva más pequeña. Los P LOLO no sustituyen al P global ni se usa su mínimo. Las sensibilidades no rescatan significación del análisis principal. La clasificación tendrá campos separados para evidencia, forma, estabilidad y selección en entrenamiento; no un único indicador que oculte desacuerdos.

| Clasificación | Requisitos operativos |
|---|---|
| Tendencia global creciente/decreciente | q principal≤0,05 y signo de S |
| Con monotonía descriptiva | La anterior y cuatro transiciones ajustadas sin inversiones |
| Dirección estable entre localidades | Signo concordante en tres LOLO y en las tres localidades excluidas; todas contrastables |
| Transferencia interna de selección | Informar 0–3 pliegues con selección realizada únicamente en entrenamiento y dirección concordante al evaluar la localidad excluida |
| Dependiente del contexto | Signos locales opuestos, o inversión de signo al excluir una localidad; descripción, no interacción formal demostrada |
| Sensible al tratamiento de ceros | Cambios de signo o selección entre pseudoconteos; conservar magnitudes y proporción de ceros |
| Sin evidencia suficiente | q principal>0,05; no significa función estable o equivalente |

La FDR corresponde a las pruebas definidas, no es una garantía nueva para cualquier subconjunto posterior basado en etiquetas. No se impondrá un corte arbitrario de magnitud después de ver los datos. Se informarán S, contrastes CLR y abundancia relativa para valorar relevancia biológica.

## 8. Sensibilidades y límites

1. Repetir D–E–I–A sin conservados. Esta prueba secundaria permite distinguir ordenación apoyada por los restaurados de separación impulsada por la referencia conservada; no sustituye el análisis solicitado de cinco estadios.
2. Contrastar ordenación de proporciones originales como sensibilidad al denominador CLR; interpretar explícitamente las discrepancias entre representación relativa y cociente respecto al centro de la composición.
3. Evaluar pseudoconteos 0,1/0,5/1 parte por millón. No exigir idéntico P, sí documentar cambio de efecto y evidencia.
4. Auditar funciones con muchos ceros y funciones presentes solo en una localidad. No excluirlas retrospectivamente porque contradigan la hipótesis.
5. Comparar resultados entre profundidades descriptivamente. Significativo a 20 cm y no a 5 cm no demuestra que el efecto difiera entre horizontes. Una prueba formal vertical posterior deberá conservar perfil como unidad y probar la diferencia de tendencias explícitamente.
6. No estimar como modelo principal una interacción completa localidad × cinco estadios: con 17 observaciones y 15 celdas quedarían solo dos grados de libertad residuales. Tampoco añadir indiscriminadamente covariables ambientales con esta replicación.
7. No construir intervalos de precisión por remuestreo ingenuo dentro de celdas con un solo perfil. Mostrar los datos, los tres efectos locales y la variación LOLO; cualquier intervalo posterior debe declarar un modelo de incertidumbre compatible con la escasa replicación.

ANCOM-BC2 incluye análisis de patrones ordenados y ajustes por covariables, por lo que sería una alternativa si se recupera una matriz de conteos apropiada y se comprueba su adecuación a este tamaño muestral. Sus supuestos no se trasladan automáticamente a estas abundancias funcionales fraccionales. No redondear datos para hacerlos pasar por conteos.

## 9. Figuras y tablas previstas

- **Figura 611A:** mapas de calor de efectos por estadio, separados por profundidad; bloques crecientes y decrecientes, nombres funcionales legibles, anotación core/no core. Identificar cualquier selección visual y conservar la tabla completa.
- **Figura 611B:** abundancias observadas y medias por estadio de funciones representativas, con color por localidad y panel por profundidad. Líneas solo para resúmenes, sin presentar perfiles distintos como individuos seguidos en el tiempo. Mostrar CLR y proporción original.
- **Figura 611C:** S global, S de los tres entrenamientos LOLO y S de las localidades excluidas. Mostrar discordancias y ceros, no solo ejemplos favorables.
- **Figura 611D:** dos UpSet, uno para tendencias crecientes y otro para decrecientes, comparando conjuntos de 5/20/40 cm. Informar inversiones de dirección entre profundidades por separado. Aplicar los mismos criterios de selección a cada conjunto.

Tabla completa por función × profundidad: ontología, ID, descripción, n positivo, proporción de ceros, core9/9, core3/3, S global y locales, P exacto, qBH/qBY, cinco medias ajustadas, cuatro cambios adyacentes, forma observada, tres resultados LOLO, frecuencia de selección en entrenamiento, sensibilidad de ceros, resultado sin conservados y estado de contraste. Mantener funciones sin evidencia y no contrastables.

## 10. Decisión para la discusión del manuscrito

Después de ejecutar, combinar los hallazgos de 610 y 611, sin anticipar el resultado:

- Core ampliamente presente + tendencias dentro del core: evidencia compatible con redistribución de funciones compartidas.
- Tendencias fuera del core: posibles contribuciones de funciones menos ubicuas, que se perderían al limitar el universo al core.
- Tendencias consistentes en LOLO: candidatos a respuestas compartidas entre estas localidades.
- Signos discordantes: respuestas contextualizadas y limitación de un indicador universal.
- Ninguna tendencia tras FDR: falta de evidencia bajo este diseño; revisar estimabilidad, precisión y tamaño de efecto, sin afirmar equivalencia funcional.

No describir funciones detectadas como activas ni tendencias relativas como tasas de procesos. Redundancia y sustitución de portadores requieren datos taxón–función. Los párrafos del manuscrito se redactarán tras revisar los resultados completos de 611.

## Referencias y procedencia

- Diseño y resultados previos: `610_informe_core100_20260923.docx`, versión 1 consultada el 23 de septiembre de 2026.
- Jonckheere AR. 1954. A distribution-free k-sample test against ordered alternatives. Biometrika 41:133–145. https://doi.org/10.1093/biomet/41.1-2.133
- Documentación del contraste estándar: https://search.r-project.org/CRAN/refmans/PMCMRplus/html/jonckheereTest.html . El contraste estándar usa aproximación normal; aquí se propone una adaptación ponderada y exacta con bloqueo.
- Holt CA, Sullivan SP. 2023. Permutation tests for experimental data. https://doi.org/10.1007/s10683-023-09799-6 . Sustento general para inferencia por permutación y restricciones de intercambio; no describe específicamente la adaptación de este protocolo.
- Lin H, Peddada SD. 2024. Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. Nature Methods 21:83–91. https://doi.org/10.1038/s41592-023-02092-7
- Tutorial oficial ANCOM-BC2 consultado: https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html . Documenta datos de conteos, análisis de patrones y sensibilidad a pseudoconteos.

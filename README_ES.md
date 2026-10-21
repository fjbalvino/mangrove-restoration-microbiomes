# Microbiomas de manglares en restauración

Repositorio mínimo: código, configuración y entradas analíticas disponibles.
Las figuras, tablas estadísticas, modelos ajustados y remuestreos se generan
localmente; no forman parte de los archivos publicados.

La sección 05 usa las matrices corregidas de fragmentos mapeados: conteos y
valores ajustados por longitud para KO y PFAM. Reproduce composición, proximidad
a referencias, núcleo funcional y las figuras desde un solo comando, una vez
aportadas las cuatro matrices externas. Esta actualización publica el código;
las matrices de lecturas se suministran por separado. Sus nombres, tamaños y
SHA-256 están en [FUNCTIONAL_INPUTS_REQUIRED.tsv](docs/FUNCTIONAL_INPUTS_REQUIRED.tsv).
El programa comprueba las entradas antes de iniciar el análisis. Se incluye
una tabla pequeña de valores esperados para comprobar la reproducción. Las imágenes
y PDF se generan localmente y no se suben.

Los análisis de genes, tendencias de KO y redes aún requieren su actualización.
Las matrices 082A y la entrada de redes son históricas y están identificadas como
tales. No alimentan el flujo primario corregido. Los detalles están en
[FUNCTIONAL_REFERENCE.md](docs/FUNCTIONAL_REFERENCE.md).

Se incluyen los metadatos de 51 muestras, las calibraciones fijas de HI/MHI,
las matrices históricas KO/PFAM y la matriz CLR de los 200 KO usados en redes.
Cada entrada está documentada con su SHA-256 en
[INPUTS_MINIMAL.tsv](docs/INPUTS_MINIMAL.tsv).

Faltan el phyloseq canónico, la matriz proteica completa y cuatro archivos
auxiliares de metadatos/QC que el código 601b comprueba antes de transformar
los conteos. Sus rutas están en [MISSING_INPUTS.tsv](docs/MISSING_INPUTS.tsv).
La anotación de genes requiere además su catálogo funcional en el servidor.

Las matrices de los 100.000 genes seleccionados son salidas de 601b y no se
consideran entradas independientes. Los objetos previos 802/803 tampoco son
necesarios: se reconstruyen el orden de aristas y las permutaciones desde los
inputs y semillas originales.

Consulta [README.md](README.md) para los comandos, dependencias y alcance de
cada componente. El flujo corregido de composición utiliza rutas relativas al
repositorio, conserva los perfiles de tres profundidades en las permutaciones y
reporta progreso cada 30 segundos. Sus salidas quedan en `local_runs/functional_reference/`.

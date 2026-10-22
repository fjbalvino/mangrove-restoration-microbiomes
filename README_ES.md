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

Ya se incorporaron los modelos de genes basados en lecturas y el enriquecimiento
funcional por eje ambiental y signo. El código de las figuras 4 y 5 usa el gradiente
azul aprobado y genera PNG, PDF y SVG localmente. Los módulos se mantienen como
candidatos funcionales potenciales, sin afirmar que las rutas estén completas o
activas. Esta actualización contiene solo código, configuración y documentación;
no incorpora bases de datos, matrices ni figuras.

Los comandos, las entradas externas y los alcances de la validación están en
[GENE_ENVIRONMENT.md](docs/GENE_ENVIRONMENT.md). El flujo actual de genes comienza
en `06_05_run_read_based.sh`; los programas `06_01`–`06_04` y el paquete anterior
de la figura 4 se conservan únicamente para reproducción histórica.
La actualización del código de tendencias ordenadas y redes queda separada.
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

El flujo corregido selecciona 100.000 genes después de centrar el CLR sobre
12.540.513 genes retenidos. Para reconstruirlo se necesita el exportado externo
por muestra y su diccionario de anotaciones. El enriquecimiento puede comenzar
directamente con las tablas completas y verificadas de los modelos; no basta
con una lista de genes significativos porque también se necesita el universo
de comparación. Los objetos previos 802/803 tampoco son
necesarios: se reconstruyen el orden de aristas y las permutaciones desde los
inputs y semillas originales.

Consulta [README.md](README.md) para los comandos, dependencias y alcance de
cada componente. El flujo corregido de composición utiliza rutas relativas al
repositorio, conserva los perfiles de tres profundidades en las permutaciones y
reporta progreso cada 30 segundos. Sus salidas quedan en `local_runs/functional_reference/`.

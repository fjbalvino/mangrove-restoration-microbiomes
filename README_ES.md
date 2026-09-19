# Microbiomas de manglares en restauración

Repositorio mínimo: código, configuración y entradas analíticas disponibles.
Las figuras, tablas estadísticas, modelos ajustados y remuestreos se generan
localmente; no forman parte de los archivos publicados.

Se incluyen los metadatos de 51 muestras, las calibraciones fijas de HI/MHI,
las matrices completas KO/PFAM y la matriz CLR de los 200 KO usados en redes.
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
cada componente. El código R conserva parámetros y rutas históricas; deben
indicarse explícitamente las entradas y salidas de cada ejecución.

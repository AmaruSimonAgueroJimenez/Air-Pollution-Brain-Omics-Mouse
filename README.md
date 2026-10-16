# Air Pollution Brain Omics — Mouse

Repositorio local de datos ómicos abiertos de cerebro de ratón (*Mus musculus*) relacionados con PM2.5 y ozono (O₃). Reúne los archivos originales disponibles, datos procesados, metadatos de muestras y un manifiesto reproducible de descarga. La selección de exposición comprende **4 estudios, 64 muestras cerebrales y 55 bibliotecas de RNA-seq**. Se añade un atlas espacial celular de referencia de cuatro ratones sin exposición experimental a contaminantes, descrito más abajo.

## Estudios incluidos

| Estudio | Contaminante | Tejido | Técnica | Muestras seleccionadas | Archivos FASTQ | FASTQ comprimidos |
|---|---|---|---|---:|---:|---:|
| [GSE297423](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297423) | PM2.5 | Corteza cingulada anterior | RNA-seq | 20 | 40 | 121,016 GB |
| [GSE145840](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145840) | PM2.5 | Hipotálamo | RNA-seq | 11 | 22 | 66,797 GB |
| [GSE93967](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93967) | PM2.5 | Hipocampo | Microarray de miRNA | 9 | 0 | No corresponde |
| [E-MTAB-14917](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14917) | O₃ | Hipocampo | RNA-seq | 24 | 48 | 74,098 GB |

Los **110 FASTQ ocupan 261,911 GB**, usando GB decimales (10⁹ bytes). Se añaden matrices procesadas, archivos originales del microarray —aproximadamente 74,864 MB—, metadatos y suplementos. El espacio requerido para descomprimir y analizar las secuencias será adicional.

La selección conserva todos los grupos experimentales cerebrales de estos estudios, incluidos controles y comparadores. GSE145840 contiene otros tejidos y ATAC-seq de hígado; la descarga se limita a sus 11 muestras de hipotálamo. Los cuatro estudios seleccionados son transcriptómica convencional o microarrays y no contienen transcriptómica espacial. Los candidatos de NO₂ están documentados en [docs/no2_access.md](docs/no2_access.md); no se ha confirmado una descarga pública de sus datos.

## Estado de descarga

La descarga de los cuatro estudios de exposición se inició el **14 de octubre de 2026** y terminó el **16 de octubre de 2026**: **136/136 archivos verificados, aproximadamente 261,996 GB en total**, sin fallidos. El estado actualizado se consulta desde la carpeta del repositorio:

```bash
python3 scripts/status.py
```

`logs/download_status.json` registra el estado de cada archivo, y `logs/completion.json` registra el resultado final. La descarga está completa cuando el resultado final indica `complete: true`, sin archivos fallidos y con todos los archivos del manifiesto verificados. Los archivos con sufijo `.part` todavía están en descarga y no deben usarse como datos completos.

## Descargar o reanudar

Requisitos: Python 3 y `curl`, conexión a internet y espacio disponible para los datos. Desde la carpeta del repositorio:

```bash
python3 scripts/prepare_manifest.py
python3 scripts/start_download.py --workers 8
python3 scripts/status.py
```

El primer paso consulta las fuentes originales, selecciona las muestras cerebrales y prepara el manifiesto. El segundo inicia una descarga independiente del chat, escribe el progreso en `logs/download.log` y descarga primero los archivos pequeños y luego los FASTQ. Repetirlo después de que el proceso termine o se interrumpa permite reanudar los archivos parciales. El repositorio bloquea una segunda instancia simultánea. Una suspensión del equipo o un reinicio puede interrumpir la transferencia; el comando permite reanudarla.

Para ejecutarla en primer plano, usar `python3 scripts/download_data.py --workers 8`. La verificación de los archivos procesados se ejecuta con `python3 scripts/verify_repository.py`; al completar toda la descarga, el proceso ejecuta automáticamente `python3 scripts/verify_repository.py --all` y genera `metadata/download_receipts.csv` con las comprobaciones de integridad.

Para descargar únicamente datos procesados y archivos originales de microarray:

```bash
python3 scripts/download_data.py --workers 8 --kind processed --kind raw_array
```

## Organización

```text
metadata/
  datasets.json                  Inventario de estudios
  download_manifest.csv          URL, destino, tamaño y MD5 esperado
  download_manifest.json         El mismo manifiesto para el descargador
  <estudio>/                     Metadatos originales y muestras seleccionadas
data/
  <estudio>/fastq/                Lecturas de RNA-seq comprimidas
  <estudio>/processed/            Matrices o resultados publicados
  <estudio>/raw_array/            Archivos originales de microarray
  <estudio>/supplements/          Suplementos de la publicación
scripts/                         Preparación, descarga y estado
logs/                            Eventos, errores y resultado de descarga
docs/                            Diseño experimental y límites de interpretación
```

Los datos voluminosos y registros de descarga permanecen en el almacenamiento local; Git conserva el código, los metadatos y la documentación. El repositorio local no implica una publicación en GitHub.

## Publicaciones y biblioteca de artículos

La biblioteca principal se limita a **cerebro**, conforme a la selección solicitada. Las cuatro publicaciones originales y dos candidatos cerebrales de NO₂ están en [docs/publications.md](docs/publications.md). La búsqueda distingue el uso explícito de las matrices de una cita al artículo. No se identificó una publicación posterior que reutilice específicamente las muestras de hipotálamo de GSE145840 en las fuentes revisadas al 14 de octubre de 2026; el artículo original ya analizó ese tejido.

Los seis artículos cerebrales de exposición obtenidos se guardan en `papers/pdf/`, incluida una copia de lectura de Li et al., 2026 exportada desde el PDF editorial completo mediante acceso UC. El inventario y las comprobaciones de cada copia están en `metadata/papers/paper_manifest.json`; las referencias para importar a un gestor bibliográfico están en [papers/bibliography.bib](papers/bibliography.bib) y [papers/bibliography.ris](papers/bibliography.ris). El inventario diferencia un PDF obtenido de un depósito de datos abierto confirmado. Los PDF permanecen en la biblioteca local y no se incluyen en Git. Los metadatos bibliográficos se contrastaron con PubMed, los depósitos de datos y las páginas editoriales. Los antecedentes de reutilización de tejidos periféricos están separados de la selección cerebral en [docs/peripheral_reuse_audit.md](docs/peripheral_reuse_audit.md).

Los suplementos cerebrales abiertos de Yan et al., 2020 se guardan en `data/NO2_Yan_2020/supplements/`. Son resultados derivados de NO₂: tablas de enriquecimiento y asociaciones lncRNA–mRNA, junto con genes seleccionados y figuras. No incluyen FASTQ ni una matriz completa de expresión por muestra y no se cuentan entre los cuatro depósitos, las 64 muestras o los 110 FASTQ del manifiesto principal.

## Atlas espacial celular de referencia sin exposición

El **16 de octubre de 2026** se completó la descarga de los cuatro conjuntos **Zhuang-ABCA-1 a Zhuang-ABCA-4** del estudio de [Zhang et al., Nature 2023](https://doi.org/10.1038/s41586-023-06808-9), desde el [Allen Brain Cell Atlas](https://alleninstitute.github.io/abc_atlas_access/notebooks/zhuang_merfish_tutorial.html): **49/49 archivos verificados, sin fallidos**. Es un atlas MERFISH de cerebro adulto de ratón sin intervención experimental de PM2.5, NO₂ u O₃. Las matrices medidas corresponden a un panel de **1.122 genes**. Los valores `raw` son conteos de moléculas detectadas por imagen, no lecturas FASTQ.

El paquete descargable comprende **49 archivos, 10,744 GB**: ocho matrices H5AD (conteos y valores `log2` para cada animal), metadatos celulares y de genes, etiquetas celulares, coordenadas registradas al marco anatómico CCF, taxonomía celular y volúmenes anatómicos NIfTI. Se conserva la estructura de versiones original bajo `data/reference_atlases/Zhuang_2023/`. Los datos medidos no incluyen las imágenes microscópicas originales ni una ampliación del transcriptoma mediante imputación.

El tutorial de Allen distingue aproximadamente 8,4 millones de células en matrices, 5,8 millones con las anotaciones celulares filtradas y 5,4 millones con coordenadas CCF registradas. Estos conjuntos no tienen exactamente las mismas filas: deben cruzarse por `cell_label`, conservando los filtros de confianza. Las coordenadas `x/y` proceden de los cortes alineados; la posición del corte en la tercera dimensión se estima durante el registro. El atlas permite estudiar la expresión y localización de firmas génicas; un efecto espacial de contaminación proyectado sobre él es una inferencia y no una exposición medida en estos ratones.

El manifiesto seleccionado está en [metadata/spatial_atlas/download_manifest.json](metadata/spatial_atlas/download_manifest.json), junto con una copia del manifiesto público de Allen. Se fijó el manifiesto `releases/20261031`, disponible al consultarlo el 16 de octubre de 2026. Las matrices conservan su versión `20230830` y los metadatos celulares su versión `20241115`. Los datos Zhuang se distribuyen bajo **CC BY 4.0**, según las fichas de Allen; aplicar las condiciones de cada recurso anatómico complementario.

El [artículo completo en PDF](papers/pdf/Zhang_2023_whole_mouse_brain_MERFISH.pdf) está descargado desde el depósito institucional Harvard DASH. Su procedencia, tamaño y SHA256 se documentan en [metadata/papers/spatial_atlas_paper.json](metadata/papers/spatial_atlas_paper.json), y la referencia se añadió a las bibliografías BibTeX y RIS.

También se guardaron los artículos de **De Filippo y Schmitz (2024)** y **Xiong et al. (2025)**, que reutilizan expresamente estas matrices. Las referencias y la evidencia están en [docs/publications.md](docs/publications.md), y las comprobaciones de sus PDF están en [metadata/papers/spatial_atlas_reuse.json](metadata/papers/spatial_atlas_reuse.json). La biblioteca cerebral suma nueve PDF.

Este atlas utiliza su propio registro de progreso. Para consultar o reanudar la descarga:

```bash
python3 scripts/download_spatial_atlas.py --status
python3 scripts/download_spatial_atlas.py --background --workers 4
```

`logs/spatial_atlas/completion.json` registra el resultado final y [metadata/spatial_atlas/download_receipts.json](metadata/spatial_atlas/download_receipts.json) guarda las comprobaciones de integridad. Los 49 archivos coincidieron en tamaño y MD5 con el manifiesto publicado. Se revisaron además las dimensiones e identificadores de las ocho matrices y los valores de sus primeras 1.000 células; el alcance de esa inspección está en [metadata/spatial_atlas/h5ad_inspection.json](metadata/spatial_atlas/h5ad_inspection.json). `--verify-only` permite repetir la comprobación completa sin descargar archivos nuevos. Las matrices y los metadatos se mantienen separados de los estudios de exposición y no se suman a sus 64 muestras ni a sus 110 FASTQ.

## Integridad y uso analítico

El descargador comprueba tamaño y MD5 publicados cuando están disponibles, identifica respuestas HTML inesperadas y comprueba la firma de los archivos gzip. También registra el MD5 calculado de cada archivo descargado. Un MD5 local sin una suma publicada permite seguir la integridad del archivo guardado, pero no constituye una comprobación independiente contra la fuente.

Antes de analizar las matrices, revisar [docs/datasets.md](docs/datasets.md): el orden de las columnas de GSE297423 necesita una correspondencia por nombre de biblioteca, GSE93967 contiene valores normalizados que pueden ser negativos y el CSV de E-MTAB-14917 es un resultado de expresión diferencial entre genotipos, **no una matriz de conteos por muestra**. Para rehacer el análisis de O₃ se necesitan los FASTQ y los metadatos experimentales.

## Procedencia y condiciones de uso

Datos obtenidos de NCBI GEO/SRA y EMBL-EBI ENA/BioStudies/ArrayExpress. Las URL originales están en el manifiesto y las descripciones depositadas se conservan en `metadata/`. Citar los accession y las publicaciones correspondientes al reutilizar los datos. Este repositorio no asigna una licencia nueva a datos de terceros ni presupone una licencia común para todos los estudios; aplicar las condiciones de cada fuente.

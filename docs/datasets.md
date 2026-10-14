# Estudios y selección de muestras

Fecha de preparación: **14 de octubre de 2026**. El alcance es cerebro de ratón (*Mus musculus*), con exposición documentada a PM2.5 u O₃. Se incluyen 64 muestras de cuatro estudios; 55 muestras corresponden a bibliotecas de RNA-seq y nueve a microarrays de miRNA. La descarga conserva controles y comparadores necesarios para interpretar las exposiciones.

## GSE297423 — PM2.5 y corteza cingulada anterior

- Fuente: [GEO GSE297423](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297423).
- Publicación vinculada en GEO: [PMID 40906862](https://pubmed.ncbi.nlm.nih.gov/40906862/).
- BioProject: [PRJNA1263837](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1263837).
- Selección: GSM8990943–GSM8990962; 20 muestras y 20 runs, con 40 FASTQ paired-end (121,016 GB comprimidos).
- Modelo: ratones con α-sinucleína humana, de dos meses al inicio; tratamientos durante dos meses.

| Grupo | Muestras |
|---|---:|
| PBS intranasal | 4 |
| PM2.5 intranasal | 4 |
| Inyección de fibrillas preformadas de α-sinucleína (PFF) | 4 |
| Inyección de PM-PFF | 4 |
| Inyección de PM-PFF seguida de PM2.5 intranasal | 4 |

El contraste directo de contaminación es PM2.5 frente a PBS (4 por grupo). Los brazos de PFF permiten investigar la interacción con patología de α-sinucleína y deben identificarse por separado. La administración intranasal debe conservarse en la descripción experimental; no equivale a una concentración ambiental inhalada.

El depósito incluye `GSE297423_mouse_gene_count.txt.gz`: conteos génicos sin normalizar, según el procesamiento declarado en GEO (mm10, GENCODE M25). Se conservan el archivo SOFT, `samples.csv` y las tablas ENA.

**Correspondencia de columnas:** no asignar los grupos a partir del orden de las muestras en GEO. Las columnas de la matriz se identifican mediante el nombre de biblioteca registrado en `!Sample_description`, cuyo orden difiere del orden de los GSM. Para **GSM8990962**, el nombre depositado es `PMPFFPM2778` y la matriz contiene `PMPFFPMA2778`. La correspondencia entre ambos es un alias inferido y debe mantenerse marcada como tal; los archivos originales no deben corregirse silenciosamente. Confirmar esa muestra con los autores o los metadatos de biblioteca antes de una conclusión que dependa de ella.

## GSE145840 — PM2.5 e hipotálamo

- Fuente: [GEO GSE145840](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145840).
- Publicación vinculada en GEO: [PMID 32780721](https://pubmed.ncbi.nlm.nih.gov/32780721/).
- BioProject: [PRJNA608484](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA608484).
- Selección: **GSM4337093–GSM4337103**, las 11 muestras de hipotálamo; 11 runs y 22 FASTQ paired-end (66,797 GB comprimidos).
- Modelo: ratones C57BL/6J machos; exposición a PM2.5 ambiental concentrado durante 14 semanas, descrita en GEO como aproximadamente 60–120 µg/m³.

| Grupo | Muestras | Accession |
|---|---:|---|
| Dieta alta en grasa (HFD) | 4 | GSM4337093–GSM4337096 |
| Aire filtrado, dieta estándar (FA) | 4 | GSM4337097–GSM4337100 |
| PM2.5, dieta estándar (PM) | 3 | GSM4337101–GSM4337103 |

El contraste de PM2.5 es PM frente a FA (3 frente a 4 muestras). HFD es un comparador metabólico, no un control adicional de aire filtrado. Cada muestra cuenta con un archivo procesado `GSM..._Hypothalamus_....txt.gz`; el depósito declara alineación a mm10 y extracción a nivel de gen con featureCounts.

**Filtro de tejido:** el proyecto completo también contiene hígado, otros tejidos, grupos de cese de exposición y ensayos ATAC-seq/OMNI. Los ATAC-seq depositados son de **hígado**, no de cerebro. Las tablas `ena_all_project_runs.tsv` conservan el inventario fuente completo para trazabilidad; `ena_selected_brain_runs.tsv` contiene solo los runs cerebrales incluidos en la descarga. No mezclar ambos inventarios al construir un análisis cerebral.

**Lista de genes:** los archivos FA y PM contienen 50.417 genes y los HFD 50.508. La intersección es de 50.417 y la unión de 50.508. Se preservaron los originales; los genes ausentes de un archivo no deben convertirse automáticamente en conteos cero. Esta diferencia está registrada en `metadata/validation_report.json` y debe resolverse antes de combinar HFD con los otros grupos.

## GSE93967 — PM2.5 y miRNA de hipocampo

- Fuente: [GEO GSE93967](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93967).
- Publicación vinculada en GEO: [PMID 28851397](https://pubmed.ncbi.nlm.nih.gov/28851397/).
- Selección: GSM2465786–GSM2465794; nueve muestras de ratones machos.
- Técnica: microarray de miRNA; **no genera FASTQ**.
- Exposición: aspiración orofaríngea de solución salina o PM2.5 cada dos días durante cuatro semanas, según el diseño depositado.

| Grupo | Muestras |
|---|---:|
| Solución salina | 3 |
| PM2.5, 1 mg/kg | 3 |
| PM2.5, 5 mg/kg | 3 |

La descarga incluye `GSE93967_RAW.tar` con archivos originales del microarray y `GSE93967_series_matrix.txt.gz` con expresión procesada. El archivo original del microarray ocupa aproximadamente **74,864 MB**. Los metadatos SOFT incluyen la anotación de plataforma y los valores de cada muestra.

**Escala de expresión:** GEO describe los valores `VALUE` como normalizados por cuantiles con GeneSpring. La matriz depositada puede contener valores negativos; esto no implica por sí mismo corrupción del archivo. Estos valores no son conteos de RNA-seq: no eliminarlos por ser negativos ni aplicar directamente modelos que requieren conteos enteros no negativos. Conservar la escala publicada o reprocesar los archivos originales con un método adecuado para microarrays.

## E-MTAB-14917 — O₃, APOE e hipocampo

- Fuente: [BioStudies/ArrayExpress E-MTAB-14917](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14917).
- Publicación: [Nakamya et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11942628/), DOI [10.3390/ijms26062407](https://doi.org/10.3390/ijms26062407).
- Proyecto ENA: [PRJEB86782](https://www.ebi.ac.uk/ena/browser/view/PRJEB86782).
- Selección: 24 muestras de hipocampo; 24 runs y 48 FASTQ paired-end (74,098 GB comprimidos).
- Modelo: ratones machos con reemplazo del gen por APOE3 o APOE4 humano; edades de 3 y 17 meses al inicio, según el protocolo depositado.
- Exposición: ocho ciclos; cada ciclo incluye cinco días de aire filtrado/O₃ a 0,8 ppm durante siete horas por día, seguidos de nueve días de aire filtrado.

El diseño tiene ocho combinaciones de edad × genotipo × exposición, con **tres bibliotecas de RNA-seq por combinación**. Los protocolos describen más animales por grupo experimental; el número de bibliotecas depositadas para transcriptómica es 24 y no debe sustituirse por el número total de animales.

Se conservan los archivos IDF/SDRF, el registro BioStudies, las tablas de runs ENA, el XML de la publicación y sus suplementos. El SDRF tiene dos filas por biblioteca paired-end: **48 filas de archivos no representan 48 muestras biológicas**. Usar `Source Name` y `Comment[ENA_RUN]` para identificar cada biblioteca y sus dos mates.

**Límite del archivo procesado:** `E4_vs_E3_005_total.csv` es una tabla de resultados de expresión diferencial **APOE4 frente a APOE3**. No es una matriz genes × muestras ni proporciona conteos por biblioteca para rehacer los contrastes de O₃. Para estimar efectos de O₃, edad y genotipo desde el principio es necesario procesar los FASTQ y generar una matriz de conteos, usando el SDRF como diseño experimental. La descarga de secuencias no implica que esa cuantificación o ese análisis se hayan ejecutado.

## Alcance y trazabilidad

Estos conjuntos no incluyen datos humanos, de rata, ni transcriptómica espacial. Las dosis, rutas, regiones cerebrales, edades y técnicas difieren entre estudios; conservar esas diferencias al compararlos. El repositorio recopila datos y metadatos, sin atribuir de antemano cambios de expresión o efectos causales que requieran análisis.

El manifiesto contiene URL, ruta local, tamaño esperado y MD5 cuando existe una suma disponible. La tabla `datasets.json` resume la selección; los SOFT/IDF/SDRF preservan la descripción fuente. El estado de integridad y descarga se registra en `logs/`, separado de los metadatos experimentales.

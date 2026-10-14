# NO₂: estudios cerebrales pendientes de acceso a datos

Fecha de revisión: **14 de octubre de 2026**. Se identificaron publicaciones de transcriptómica cerebral de ratón con exposición a NO₂, pero no se confirmó un accession GEO/SRA/ENA u otra descarga pública de sus secuencias o matrices completas por muestra. Sí se obtuvieron resultados derivados públicos en los suplementos de Yan et al., 2020. Estos se guardan por separado y no se añaden al manifiesto de FASTQ ni al total de 64 muestras de los cuatro depósitos principales.

| Publicación | Modelo y tejido descritos | Tipo de datos de interés | Estado de acceso |
|---|---|---|---|
| [Yan et al., 2020 — PMID 32203807](https://pubmed.ncbi.nlm.nih.gov/32203807/) | Ratón; corteza de crías tras exposición prenatal a NO₂ | Perfiles de mRNA y lncRNA | Depósito público descargable no confirmado |
| [Li et al., 2026 — PMID 42526541](https://pubmed.ncbi.nlm.nih.gov/42526541/) | Ratón; hipocampo de crías tras exposición gestacional materna a NO₂ | Transcriptómica | Depósito público descargable no confirmado |

## Artículos obtenidos y revisión del texto completo

Se obtuvo el [PDF editorial abierto de Yan et al., 2020](../papers/pdf/NO2_Yan_2020_Environment_International.pdf) y una [copia de lectura del PDF editorial de Li et al., 2026 mediante acceso UC](../papers/pdf/NO2_Li_2026_Neurotoxicology_reading_copy.pdf). Las copias tienen 11 y 12 páginas respectivamente y sus DOI fueron verificados. El guardado directo del PDF original de Li permaneció deshabilitado; la copia se exportó desde el visor editorial y se revisaron sus páginas inicial, intermedia y final.

Yan et al. describe mRNA-seq de corteza en días posnatales 1, 7 y 21 y lncRNA-seq en día 21, con tres crías por edad, sexo y grupo (página 3, métodos 2.2 y 2.4). El artículo y sus hipervínculos no identifican un accession ni una sección formal de disponibilidad. Se obtuvieron los suplementos oficiales [DOC](../data/NO2_Yan_2020/supplements/Yan2020_mmc2.doc) y [XLS](../data/NO2_Yan_2020/supplements/Yan2020_mmc3.xls). El XLS contiene cuatro hojas: enriquecimiento GO y asociaciones lncRNA–mRNA, separadas por sexo. El DOC añade genes seleccionados, primers y figuras. La inspección completa no identificó expresión por muestra, conteos/FPKM completos, FASTQ ni identificadores de bibliotecas/runs. Son **resultados cerebrales derivados abiertos**, y no datos primarios completos. La estructura, procedencia y las comprobaciones se conservan en [no2_supplement_inspection.json](../metadata/papers/no2_supplement_inspection.json).

Li et al. declara en la página 10, sección «Data availability», **“Data will be made available on request.”** Su RNA-seq es de hipocampo de crías **macho**, con tres crías por grupo según la Figura 2; los ensayos conductuales y otras mediciones incluyen ambos sexos. No se identifican bibliotecas/runs de un depósito público. Li et al. cita a Yan et al. pero describe datos propios; no se confirmó reutilización de los perfiles de corteza de Yan.

**No confirmar el acceso no demuestra que los datos no existan.** Los autores pueden disponer de los datos o de un depósito que aún no se haya vinculado de forma verificable a la publicación. Tampoco se ha confirmado un conjunto abierto de transcriptómica espacial cerebral con exposición específica a NO₂.

Para incorporar un estudio, el acceso debe permitir vincular los archivos al trabajo y a las muestras cerebrales de ratón: accession estable o enlace de descarga verificable, tabla de muestras, exposición/control y diseño experimental. En datos de secuenciación se requieren las bibliotecas/runs y sus archivos originales; una publicación o una tabla de genes destacados no sustituye una matriz por muestra.

Si el depósito no se identifica, el siguiente paso es solicitar a los autores las matrices completas, archivos originales y metadatos de muestra, junto con sus condiciones de reutilización. Este repositorio no registra ninguna solicitud enviada ni presume autorización para redistribuir datos recibidos de forma privada. Una respuesta verificable podrá documentarse aquí y añadirse al manifiesto con su procedencia.

# Publicaciones de datos ómicos de cerebro de ratón

Revisión al **14 de octubre de 2026**, limitada a **cerebro**. La biblioteca contiene seis artículos: cuatro publicaciones originales vinculadas a los conjuntos descargables y dos estudios cerebrales de NO₂ cuyos datos completos por muestra aún no están disponibles de forma pública confirmada. Se verificaron las seis copias PDF, sus DOI y sus páginas. La copia de Li et al., 2026 se exportó para lectura desde el PDF editorial completo: el guardado directo del original permanecía deshabilitado. La descarga de matrices y secuencias tiene su propio estado y no se considera completada por obtener un artículo.

| Conjunto | Contaminante | Región seleccionada | Publicación original | Artículo obtenido |
|---|---|---|---|---|
| [GSE297423](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297423) | PM2.5 | Corteza cingulada anterior | Zhang et al., 2025 | [PDF](../papers/pdf/GSE297423_Zhang_2025_Science.pdf) |
| [GSE145840](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145840) | PM2.5 | Hipotálamo | Rajagopalan et al., 2020 | [PDF](../papers/pdf/GSE145840_Rajagopalan_2020_JCI.pdf) |
| [GSE93967](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93967) | PM2.5 | Hipocampo | Ku et al., 2017 | [PDF](../papers/pdf/GSE93967_Ku_2017_Particle_Fibre_Toxicology.pdf) |
| [E-MTAB-14917](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14917) | O₃ | Hipocampo | Nakamya et al., 2025 | [PDF](../papers/pdf/E-MTAB-14917_Nakamya_2025_IJMS.pdf) |

## Originales: quién generó y publicó los datos

### GSE297423 — Zhang et al., 2025

**Lewy body dementia promotion by air pollutants.** *Science*. 2025;389(6764):eadu4132. DOI: [10.1126/science.adu4132](https://doi.org/10.1126/science.adu4132). [PMID 40906862](https://pubmed.ncbi.nlm.nih.gov/40906862/).

Xiaodi Zhang, Haiqing Liu, Xiao Wu y Longgang Jia comparten primera autoría; último autor **Xiaobo Mao**, Johns Hopkins. Correspondencia: Mao, Shizhong Han, Ted M. Dawson y Wu. **Evidencia:** GEO vincula ese PMID; la disponibilidad del [texto original](https://pmc.ncbi.nlm.nih.gov/articles/PMC12459341/) declara RNA-seq original/procesado en **GSE297423 / PRJNA1263837**. Son 20 muestras cerebrales, con PM2.5 intranasal y brazos de α-sinucleína. El PDF local corresponde a la versión editorial obtenida mediante acceso UC.

### GSE145840 — Rajagopalan et al., 2020

**Metabolic effects of air pollution exposure and reversibility.** *Journal of Clinical Investigation*. 2020;130(11):6034–6040. DOI: [10.1172/JCI137315](https://doi.org/10.1172/JCI137315). [PMID 32780721](https://pubmed.ncbi.nlm.nih.gov/32780721/).

**Sanjay Rajagopalan y Bongsoo Park** comparten primera autoría; último autor **Shyam Biswal**. Correspondencia: Rajagopalan y Biswal, Case Western Reserve/Johns Hopkins. **Evidencia:** GEO y «Data availability» del [artículo](https://www.jci.org/articles/view/137315) identifican **GSE145840**. Aquí se seleccionan 11 muestras de hipotálamo. El original analiza este tejido en su RNA-seq y remite sus vías funcionales a **Figura suplementaria 9**; los otros tejidos del depósito quedan fuera de la selección.

### GSE93967 — Ku et al., 2017

**NF-κB-regulated microRNA-574-5p underlies synaptic and cognitive impairment in response to atmospheric PM2.5 aspiration.** *Particle and Fibre Toxicology*. 2017;14:34. DOI: [10.1186/s12989-017-0215-3](https://doi.org/10.1186/s12989-017-0215-3). [PMID 28851397](https://pubmed.ncbi.nlm.nih.gov/28851397/).

Primera autora **Tingting Ku**; último autor y correspondencia **Nan Sang**, Shanxi University. **Evidencia:** GEO vincula ese PMID; «MicroRNA microarray analysis» y disponibilidad del [original](https://pmc.ncbi.nlm.nih.gov/articles/PMC5575838/) declaran **GSE93967**. Son nueve microarrays de miRNA de hipocampo tras aspiración orofaríngea de PM2.5; no genera FASTQ.

### E-MTAB-14917 — Nakamya et al., 2025

**Age- and ApoE Genotype-Dependent Transcriptomic Responses to O₃ in the Hippocampus of Mice.** *International Journal of Molecular Sciences*. 2025;26(6):2407. DOI: [10.3390/ijms26062407](https://doi.org/10.3390/ijms26062407). [PMID 40141051](https://pubmed.ncbi.nlm.nih.gov/40141051/).

**Mary F. Nakamya y Kaili Hu** comparten primera autoría; último autor **Rui-Ming Liu**. Correspondencia: Liu y **Zechen Chong**, University of Alabama at Birmingham. **Evidencia:** métodos RNA-seq del [original](https://pmc.ncbi.nlm.nih.gov/articles/PMC11942628/) identifican **E-MTAB-14917**, enlazado a **PRJEB86782**. Son 24 bibliotecas de hipocampo con edad, APOE3/APOE4 y O₃/aire filtrado.

## Quién ha publicado después usando datos cerebrales

El original de Rajagopalan **ya analizó el hipotálamo**. No se verificó una publicación posterior que reanalice las muestras cerebrales de ninguno de estos cuatro conjuntos en la búsqueda por accession, alias y texto completo realizada hasta la fecha indicada. El uso posterior de tejidos periféricos de GSE145840 no equivale a un análisis cerebral; su evidencia queda en [auditoría periférica](peripheral_reuse_audit.md), fuera de la selección. Esta búsqueda no demuestra ausencia universal de reutilizaciones.

## NO₂: artículos cerebrales y acceso a datos primarios

**Yan et al., 2020. Prenatal NO₂ exposure and neurodevelopmental disorders in offspring mice: Transcriptomics reveals sex-dependent changes in cerebral gene expression.** *Environment International*. 138:105659. DOI: [10.1016/j.envint.2020.105659](https://doi.org/10.1016/j.envint.2020.105659). [PMID 32203807](https://pubmed.ncbi.nlm.nih.gov/32203807/). Primera autoría **Wei Yan**; último autor/correspondencia **Nan Sang**, Shanxi University. Corteza de crías tras exposición prenatal; perfiles mRNA/lncRNA. El resumen primario confirma tejido y diseño; **accession descargable no confirmado**.

[PDF editorial abierto, 11 páginas](../papers/pdf/NO2_Yan_2020_Environment_International.pdf). La revisión del texto completo y sus enlaces no identificó un accession GEO/SRA/ENA ni una sección formal de disponibilidad. Los métodos 2.2 y 2.4 describen RNA-seq de corteza en diferentes edades posnatales y ambos sexos; el artículo original ya contiene esos análisis. Se descargaron el [documento](../data/NO2_Yan_2020/supplements/Yan2020_mmc2.doc) y la [tabla suplementaria](../data/NO2_Yan_2020/supplements/Yan2020_mmc3.xls) de la editorial. Contienen resultados derivados cerebrales —enriquecimiento GO, asociaciones lncRNA–mRNA y genes seleccionados—; la inspección completa no identificó una matriz de expresión por muestra ni FASTQ. Se registran por separado del manifiesto de datos primarios.

**Li et al., 2026. Gestational NO₂ exposure induces anxiety- and depression-like behaviors in mice via hippocampal PTGDS-PGD₂ activation and tryptophan metabolism reprogramming.** *Neurotoxicology*. 116:103527. DOI: [10.1016/j.neuro.2026.103527](https://doi.org/10.1016/j.neuro.2026.103527). [PMID 42526541](https://pubmed.ncbi.nlm.nih.gov/42526541/). Primer autor **Wei Li**; última autora **Fang Xiao**; correspondencia registrada **Di Ji**, Central South University. Hipocampo de crías en día posnatal 22; alteraciones transcriptómicas tras exposición gestacional. El resumen primario confirma diseño; **accession descargable no confirmado**.

[Copia de lectura del PDF editorial obtenida mediante UC, 12 páginas](../papers/pdf/NO2_Li_2026_Neurotoxicology_reading_copy.pdf). El RNA-seq se realizó en **crías macho**, con tres crías por grupo según la Figura 2; los experimentos conductuales y otras mediciones incluyen ambos sexos. La página 10, sección «Data availability», declara: **“Data will be made available on request.”** No indica un accession público. La cita a Yan et al. no demuestra reutilización de sus datos: Li et al. describe su propio RNA-seq de hipocampo.

Obtener los artículos de NO₂ no incorpora automáticamente sus matrices o secuencias. Consultar [no2_access.md](no2_access.md), [metadatos originales](../metadata/papers/originals_pubmed.json), [registros PubMed de NO₂](../metadata/papers/no2_pubmed.xml) y [datasets.md](datasets.md) para trazabilidad y selección de muestras.

Las referencias se pueden importar desde [BibTeX](../papers/bibliography.bib) o [RIS](../papers/bibliography.ris). El [inventario de PDF](../metadata/papers/paper_manifest.json) registra procedencia, tipo de copia, páginas y SHA-256; las comprobaciones locales permiten identificar cambios del archivo y no sustituyen una suma publicada por la editorial.

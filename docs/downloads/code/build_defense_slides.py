#!/usr/bin/env python3
"""Author paired Quarto/Reveal.js defense decks from the frozen project outputs.

This script reorganizes completed analyses. It does not refit a scientific model.
Visible slides, detailed speaker notes, provenance and EN/ES IDs are generated
together so that language switching preserves the same scientific slide.
"""
from pathlib import Path
import argparse
import hashlib
import html
import json
import re
import subprocess
import zipfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/defense_slides"
OUT.mkdir(parents=True, exist_ok=True)
BASE = "https://amarusimonaguerojimenez.github.io/Air-Pollution-Brain-Omics-Mouse/"
GH = "https://github.com/AmaruSimonAgueroJimenez/Air-Pollution-Brain-Omics-Mouse/blob/main/"
DATE = "2026-10-29"
SLIDES = []
SOURCES = set()

PAPERS = {
    "HY": ("Rajagopalan et al. 2020", "https://doi.org/10.1172/JCI137315"),
    "ACA": ("Zhang et al. 2025", "https://doi.org/10.1126/science.adu4132"),
    "miRNA": ("Ku et al. 2017", "https://doi.org/10.1186/s12989-017-0215-3"),
    "O3": ("Nakamya et al. 2025", "https://doi.org/10.3390/ijms26062407"),
    "atlas": ("Zhang et al. 2023", "https://doi.org/10.1038/s41586-023-06808-9"),
    "DESeq2": ("Love et al. 2014", "https://doi.org/10.1186/s13059-014-0550-8"),
    "limma": ("Ritchie et al. 2015", "https://doi.org/10.1093/nar/gkv007"),
    "Salmon": ("Patro et al. 2017", "https://doi.org/10.1038/nmeth.4197"),
    "tximport": ("Soneson et al. 2015", "https://doi.org/10.12688/f1000research.7563.2"),
    "PyDESeq2": ("Muzellec et al. 2023", "https://doi.org/10.1093/bioinformatics/btad547"),
    "apeglm": ("Zhu et al. 2019", "https://doi.org/10.1093/bioinformatics/bty895"),
    "GSEA": ("Subramanian et al., PNAS 2005", "https://doi.org/10.1073/pnas.0506580102"),
    "GO": ("Gene Ontology Consortium", "https://geneontology.org/docs/go-enrichment-analysis/"),
    "GPML": ("Rasmussen & Williams, 2006, ch. 2", "https://gaussianprocess.org/gpml/chapters/RW2.pdf"),
    "Murphy": ("Murphy, 2007, Gaussian conjugacy", "https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf"),
    "NCBI": ("NCBI Gene: Tll1", "https://www.ncbi.nlm.nih.gov/gene/21892"),
    "UniProt": ("UniProtKB: Tll1", "https://www.uniprot.org/uniprotkb/Q62381/entry"),
    "Beta": ("SciPy, beta-binomial distribution", "https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.betabinom.html"),
}
PAPERS.update({
    "HY": ("Rajagopalan et al., JCI 2020", PAPERS["HY"][1]),
    "ACA": ("Zhang et al., Science 2025", PAPERS["ACA"][1]),
    "miRNA": ("Ku et al., Part. Fibre Toxicol. 2017", PAPERS["miRNA"][1]),
    "O3": ("Nakamya et al., IJMS 2025", PAPERS["O3"][1]),
    "atlas": ("Zhang et al., Nature 2023", PAPERS["atlas"][1]),
    "DESeq2": ("Love et al., Genome Biol. 2014", PAPERS["DESeq2"][1]),
    "limma": ("Ritchie et al., NAR 2015", PAPERS["limma"][1]),
})

def J(path):
    SOURCES.add(path)
    return json.loads((ROOT/path).read_text())

def T(path):
    SOURCES.add(path)
    return pd.read_csv(ROOT/path, sep="\t" if ".tsv" in path else ",")

def table(headers, rows):
    return "\n".join(["| " + " | ".join(map(str, headers)) + " |",
                      "|" + "|".join("---" for _ in headers) + "|"] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])

def block(text, cls):
    return f"\n::: {{.{cls}}}\n{text}\n:::\n"

def add(sid, en_title, es_title, en, es, en_notes, es_notes, sources=(), cls=""):
    for src in sources:
        if src not in PAPERS and not src.startswith("http"):
            SOURCES.add(src)
    SLIDES.append(dict(id=sid, en_title=en_title, es_title=es_title, en=en.strip(),
                       es=es.strip(), en_notes=en_notes.strip(), es_notes=es_notes.strip(),
                       sources=list(sources), cls=cls))

def fig(sid, en_title, es_title, stem, en_caption, es_caption,
        en_notes, es_notes, sources=()):
    add(sid, en_title, es_title, "", "", en_notes, es_notes, sources, "figure-slide")
    SLIDES[-1].update(figure=stem, en_caption=en_caption, es_caption=es_caption)

def section(sid, en_title, es_title, en, es, en_notes, es_notes):
    add(sid,en_title,es_title,en,es,en_notes,es_notes,cls="section-slide")

def content():
    hy=J("results/hypothalamus/summary.json")
    aca=J("results/cortex/summary.json")
    mir=J("results/mirna/summary.json")
    o3=J("results/ozone/summary.json")
    rep=J("results/replication_pm/cortex_refit_summary.json")
    adv=J("results/advanced_omics/summary.json")
    gauss=J("results/spatial_bayes/summary.json")
    det=J("results/spatial_detection_bayes/summary.json")
    sections=J("results/spatial_sections/summary.json")
    genes=T("results/gene_priorities/prioritized_genes.tsv")
    enrich=T("results/enrichment/summary_counts.tsv")

    add("scope","Scientific scope and levels of evidence","Alcance científico y niveles de evidencia",
        """- **Question:** which regional expression patterns accompany the original mouse exposure experiments?
- **Reproduction:** recover deposited author results using the same source animals.
- **Project extensions:** new primary contrasts, enrichment, interactions and spatial-reference models.
- **Spatial interpretation:** measured unexposed cells, fitted reference estimates and an illustrative PM₂.₅ scenario remain separately identified.""",
        """- **Pregunta:** ¿qué patrones regionales de expresión acompañan los experimentos originales de exposición en ratón?
- **Reproducción:** recuperar resultados de los autores usando los mismos animales fuente.
- **Extensiones del proyecto:** contrastes principales nuevos, enriquecimiento, interacciones y modelos de referencia espacial.
- **Interpretación espacial:** se distinguen células medidas sin exposición, estimaciones ajustadas y un escenario ilustrativo con PM₂.₅.""",
        """This is a defense-style presentation of completed computational work, not a claim of a new animal experiment or independent biological replication. “New” means added by this project, not proven biological novelty. The analyses and reference outputs were completed on 18–20 October 2026; this presentation was prepared on 28 October 2026. The data do not contain matched pollutant-exposed spatial or single-cell brains. No NO₂ expression matrix enters the four analyzed studies. Use the Notes button for the complete methodological discussion, arrow keys to navigate, Esc for overview and S for the optional speaker window. The two languages share explicit slide identifiers.""",
        """Esta presentación organiza trabajo computacional terminado, sin afirmar un experimento animal nuevo ni replicación biológica independiente. «Nuevo» significa incorporado por el proyecto, no novedad biológica demostrada. Los análisis y referencias se completaron el 18–20 de octubre de 2026; la presentación se preparó el 28 de octubre de 2026. Los datos no contienen cerebros espaciales o single-cell expuestos y emparejados. Ninguna matriz de expresión de NO₂ entra en los cuatro estudios analizados. El botón Notas amplía la discusión metodológica, las flechas navegan, Esc abre la vista general y S permite la ventana opcional del presentador. Los idiomas comparten identificadores explícitos de diapositiva.""",
        ["results/replication_pm/summary.json","results/spatial_bayes/summary.json"])

    add("intro","Brain regions and exposure designs define the question","La región cerebral y el diseño de exposición definen la pregunta",
        """Regional RNA measurements can identify expression associations under a particular route, dose, genotype and age.

The project addresses three linked questions:

1. Can the deposited author analyses be reconstructed?
2. What do explicitly specified reanalyses and response-interaction tests show?
3. Where are prioritized transcripts detected in an independent measured brain reference?""",
        """Las mediciones regionales de ARN permiten estudiar asociaciones de expresión bajo una ruta, dosis, genotipo y edad concretos.

El proyecto aborda tres preguntas relacionadas:

1. ¿Se pueden reconstruir los análisis depositados por los autores?
2. ¿Qué muestran los reanálisis explícitos y las pruebas de interacción de respuesta?
3. ¿Dónde se detectan los transcritos priorizados en una referencia cerebral medida e independiente?""",
        """The background is deliberately brief. The four primary papers motivate hypothalamic, cortical and hippocampal questions involving particulate matter and ozone, but their routes and units are not interchangeable. An exposure association in bulk tissue may reflect regulation, tissue composition, source-quality variation or combinations of these. Spatial reference anatomy helps locate baseline detection, but cannot identify which individual cells responded in the original experiment. This presentation prioritizes the computational evidence and the limits needed to defend each inference.""",
        """La introducción es deliberadamente breve. Los cuatro artículos originales motivan preguntas hipotalámicas, corticales e hipocampales sobre material particulado y ozono, pero sus rutas y unidades no son intercambiables. Una asociación en tejido bulk puede reflejar regulación, composición, calidad de la fuente o una combinación. La anatomía de referencia permite localizar detección basal, pero no identificar qué células individuales respondieron en el experimento original. La presentación prioriza la evidencia computacional y los límites necesarios para defender cada inferencia.""",
        ["HY","ACA","miRNA","O3","atlas"])

    add("studies-rna","Two PM₂.₅ RNA-seq studies","Dos estudios RNA-seq de PM₂.₅",
        table(["Feature","Hypothalamus","Anterior cingulate cortex"],[
            ["Deposit","GSE145840","GSE297423"],
            ["Primary animals","4 FA + 3 PM","4 PBS + 4 PM"],
            ["Model / context","Male C57BL/6J, chow-fed","Humanized hWT α-synuclein, sex unreported"],
            ["Exposure route","Concentrated ambient aerosol inhalation","Repeated intranasal PM suspension"],
            ["Starting material","Seven deposited gene-count files","Eight selected libraries from a 20-library matrix"],
            ["Primary inference","PyDESeq2, PM versus FA","PyDESeq2, PM versus PBS"]]),
        table(["Característica","Hipotálamo","Corteza cingulada anterior"],[
            ["Depósito","GSE145840","GSE297423"],
            ["Animales principales","4 FA + 3 PM","4 PBS + 4 PM"],
            ["Modelo / contexto","C57BL/6J macho, dieta estándar","α-sinucleína hWT humanizada, sexo no informado"],
            ["Ruta de exposición","Inhalación de aerosol ambiental concentrado","Suspensión de PM intranasal repetida"],
            ["Material de partida","Siete archivos de conteos génicos","Ocho bibliotecas de una matriz de 20"],
            ["Inferencia principal","PyDESeq2, PM frente a FA","PyDESeq2, PM frente a PBS"]]),
        """FA means filtered air, whereas PBS is the cortex vehicle control. The four HFD hypothalamus libraries are not additional FA controls. The cortex RNA-seq arm begins at two months in humanized hWT α-synuclein mice; it must not be conflated with the paper's distinct older conventional-mouse arm. Both primary reanalyses start from deposited counts, so read alignment and read-level QC were not independently reconstructed. Detailed exposure protocols and platform information appear in the appendix.""",
        """FA significa aire filtrado y PBS es el vehículo de control de corteza. Las cuatro bibliotecas hipotalámicas HFD no son controles FA adicionales. El brazo RNA-seq cortical comienza a los dos meses en ratones con α-sinucleína hWT humanizada; no debe confundirse con el brazo distinto de ratones convencionales mayores del artículo. Ambos reanálisis parten de conteos depositados, sin reconstruir alineamiento ni QC de lecturas de forma independiente. El apéndice detalla protocolos y plataformas.""",
        ["HY","ACA","results/hypothalamus/methods.json","results/cortex/methods.json"])

    add("studies-hippocampus","Hippocampal miRNA and ozone studies","Estudios hipocampales de miRNA y ozono",
        table(["Feature","PM₂.₅ miRNA","Ozone RNA-seq"],[
            ["Deposit","GSE93967","E-MTAB-14917"],
            ["Sequenced / assayed animals","9 arrays, 3 per group","24 animals, 3 in each of 8 strata"],
            ["Design","Saline, 1 or 5 mg/kg PM","2 ages × 2 APOE genotypes × 2 exposures"],
            ["Assay","Agilent miRNA V19 microarray","Paired-end left-hippocampus RNA-seq"],
            ["Input","Deposited quantile-normalized signal","Original FASTQ pairs"],
            ["Scale","Original signal units","Estimated gene fragments + length offsets"]]),
        table(["Característica","miRNA con PM₂.₅","RNA-seq de ozono"],[
            ["Depósito","GSE93967","E-MTAB-14917"],
            ["Animales medidos","9 arrays, 3 por grupo","24 animales, 3 en cada uno de 8 estratos"],
            ["Diseño","Salina, PM de 1 o 5 mg/kg","2 edades × 2 genotipos APOE × 2 exposiciones"],
            ["Ensayo","Microarray Agilent miRNA V19","RNA-seq paired-end de hipocampo izquierdo"],
            ["Entrada","Señal con normalización cuantílica depositada","Pares FASTQ originales"],
            ["Escala","Unidades originales de señal","Fragmentos génicos estimados + offsets de longitud"]]),
        """The miRNA experiment uses oropharyngeal aspiration, not an ambient airborne dose. The ozone ages 3 and 17 months describe exposure onset. Its paper-wide animal numbers are not the transcriptomic sample size: only three sequenced mice per age/genotype/exposure stratum enter this analysis. APOE3/APOE4 are human targeted-replacement genotypes. Mouse-reference quantification does not establish allele-specific human APOE expression. The miRNA scale has signed values and an uncertain log2 preprocessing history, so no second log transform or count likelihood was imposed.""",
        """El experimento de miRNA utiliza aspiración orofaríngea, no dosis ambiental aérea. Las edades de ozono de 3 y 17 meses corresponden al inicio de exposición. El número total de animales del artículo no es el tamaño transcriptómico: entran tres ratones secuenciados por estrato edad/genotipo/exposición. APOE3/APOE4 son genotipos humanos de reemplazo dirigido. Cuantificar con referencia de ratón no establece expresión alélica humana de APOE. La señal miRNA tiene valores negativos y un historial log₂ no confirmado; no se impuso otra transformación logarítmica ni una verosimilitud de conteos.""",
        ["miRNA","O3","results/mirna/methods.json","results/ozone/methods.json"])

    add("units","Experimental units and analytical populations","Unidades experimentales y poblaciones analíticas",
        """- **Bulk and arrays:** the sequenced or assayed animal is the biological replicate.
- **Author reproduction:** the same deposited animals supply computational agreement.
- **Spatial reference:** 52,393 ACA/HY cells come from one unexposed donor.
- **Spatial validation:** reserving sections measures within-donor prediction.

No pooled pollutant dose, cross-study effect or independent-cell exposure test is estimated.""",
        """- **Bulk y arrays:** el animal secuenciado o medido es la réplica biológica.
- **Reproducción de autores:** los mismos animales depositados permiten evaluar concordancia computacional.
- **Referencia espacial:** 52,393 células ACA/HY proceden de un donante sin exposición.
- **Validación espacial:** reservar cortes evalúa predicción dentro del donante.

No se estima una dosis conjunta de contaminantes, un efecto combinado entre estudios ni una prueba de exposición con células independientes.""",
        """This slide prevents pseudoreplication. Thousands of cells do not replace independent exposed and control animals. The spatial folds are useful for checking prediction at unsampled reference sections, but they cannot validate across donors or exposures. A repeated gene–cell evaluation in a ten-gene pooled score is not ten new cells. Similarly, fitting another software implementation to the same seven animals does not create a new experiment.""",
        """Esta diapositiva previene pseudorreplicación. Miles de células no sustituyen animales expuestos y controles independientes. Los bloques espaciales permiten comprobar predicción en cortes de referencia reservados, pero no validar entre donantes o exposiciones. Evaluar una célula para diez genes no crea diez células nuevas. Ajustar otra implementación de software a los mismos siete animales tampoco crea un experimento nuevo.""",
        ["results/hypothalamus/summary.json","results/spatial_sections/summary.json"])

    fig("workflow","Project workflow and provenance","Flujo del proyecto y procedencia","context_workflow",
        "Four experimental studies and a separate spatial reference feed explicitly labeled analytical branches.",
        "Cuatro estudios experimentales y una referencia espacial independiente alimentan ramas analíticas explícitas.",
        """Inputs, sample mappings, models and results have retained identifiers, hashes and downloadable tables. The ozone branch begins with original reads; PM bulk starts with supplied counts and miRNA with deposited normalized arrays. Author reproduction and primary reanalysis remain distinct. Gene prioritization connects regional results to measured panel coverage, and the atlas supports observed sections, illustrative scaling and baseline models. The arrows describe processing order, not a biological causal diagram.""",
        """Las entradas, correspondencias muestrales, modelos y resultados conservan identificadores, hashes y tablas descargables. La rama de ozono parte de lecturas originales; bulk PM, de conteos suministrados; miRNA, de arrays normalizados depositados. La reproducción y el reanálisis principal permanecen diferenciados. La priorización vincula resultados regionales con cobertura del panel medido; el atlas sustenta cortes observados, escalamiento ilustrativo y modelos basales. Las flechas representan procesamiento, no un diagrama causal biológico.""",
        ["results/ozone/methods.json","results/replication_pm/methods.json","metadata/spatial_sections/methods.json"])

    section("methods","Methods","Metodología",
        "Original-study reconstruction, new statistical analyses and measured spatial reference",
        "Reconstrucción de los estudios originales, análisis estadísticos nuevos y referencia espacial medida",
        "The next block defines what each model estimates before showing results. Detailed assumptions and source links remain available in notes and the appendix.",
        "El bloque siguiente define qué estima cada modelo antes de presentar resultados. Los supuestos detallados y las fuentes permanecen en las notas y el apéndice.")

    add("reproduction-criteria","Criteria for computational reproduction","Criterios de reproducción computacional",
        table(["Layer","Comparison","Interpretation"],[
            ["Inputs","Identifiers, counts, samples and source hashes","Verify the same analytical starting point"],
            ["Model","Design, normalization, filter and multiplicity family","Reconstruct the stated procedure"],
            ["Outputs","Feature identities, effects, p values and selected sets","Counts alone are insufficient"],
            ["Missing specification","Unrecovered parameters and altered quantification","Report sensitivity or concordance explicitly"]]),
        table(["Nivel","Comparación","Interpretación"],[
            ["Entradas","Identificadores, conteos, muestras y hashes","Verificar el mismo punto de partida"],
            ["Modelo","Diseño, normalización, filtro y familia de multiplicidad","Reconstruir el procedimiento declarado"],
            ["Salidas","Identidades, efectos, valores p y conjuntos seleccionados","El número de genes por sí solo es insuficiente"],
            ["Especificación faltante","Parámetros no recuperados o cuantificación distinta","Declarar sensibilidad o concordancia"]]),
        """The analysis distinguishes exact output reconstruction from a transparent substitute workflow. Exact recovery applies only to the output layer actually checked. It does not imply recovery of raw-read alignment, histology, protein experiments or every figure. Missing parameters are not selected after inspecting which setting maximizes agreement. This project preserves disagreements and source-quality questions alongside successful numerical reconstruction.""",
        """Se distingue reconstrucción exacta de salidas de un flujo sustituto transparente. La recuperación exacta se limita al nivel efectivamente cotejado; no implica reconstruir alineamiento, histología, proteínas ni todas las figuras. Los parámetros ausentes no se eligen buscando maximizar concordancia. El proyecto conserva discrepancias y preguntas de calidad junto con la reproducción numérica exitosa.""",
        ["results/replication_pm/paper_claim_audit.csv","results/replication_mirna_ozone/paper_replication_audit.csv"])

    add("cortex-author-method","Cortex author-aligned count analysis","Análisis cortical de conteos alineado con los autores",
        """1. Match the original Zenodo count matrix to the GEO deposit.
2. Retain all **20 libraries in five treatment groups**.
3. Keep genes with ≥10 counts in ≥4 libraries: **21,388 genes**.
4. Reconstruct R DESeq2 with design <code>~group</code>.
5. Preserve default independent-filtering α = 0.1, then select **BH q < 0.05**.
6. Reconstruct the code-defined PCA using TMM log₂CPM, top 500 variable genes and variance scaling.""",
        """1. Cotejar la matriz original de Zenodo con el depósito GEO.
2. Conservar las **20 bibliotecas de cinco grupos de tratamiento**.
3. Retener genes con ≥10 conteos en ≥4 bibliotecas: **21,388 genes**.
4. Reconstruir R DESeq2 con diseño <code>~group</code>.
5. Conservar α = 0.1 del filtrado independiente por defecto y seleccionar **BH q < 0.05**.
6. Reconstruir el PCA definido por código: TMM log₂CPM, 500 genes más variables y escalamiento por varianza.""",
        """The five arms are PBS, nasal PM, PFF inoculation, PM-PFF inoculation and PM-PFF plus nasal PM. They do not constitute a standard complete PM×PFF factorial because the inoculated materials differ. The repository reimplements the relevant code rather than running an arbitrary downloaded script wholesale. The source prose mentions filterByExpr but the deposited executable workflow uses an explicit count filter; this discrepancy is documented. R DESeq2 alpha controls independent-filter threshold selection; the final gene-set gate remains q<.05. Do not replace this 20-animal reconstruction with the eight-animal project primary result.""",
        """Los cinco brazos son PBS, PM nasal, inoculación PFF, inoculación PM-PFF y PM-PFF más PM nasal. No forman un factorial PM×PFF completo estándar, porque el material inoculado difiere. El repositorio reimplementa el código pertinente, sin ejecutar íntegramente un script descargado arbitrario. El texto menciona filterByExpr, pero el código depositado usa un filtro explícito de conteos; se documenta la discrepancia. El α de R DESeq2 controla la elección del filtrado independiente; el umbral final sigue siendo q<.05. Esta reconstrucción de 20 animales no se sustituye por el resultado principal del proyecto con ocho animales.""",
        ["ACA","https://doi.org/10.5281/zenodo.15374632","results/replication_pm/methods.json"],"small")

    add("hypothalamus-author-method","Hypothalamic RUVg specification remains incomplete","La especificación RUVg hipotalámica sigue incompleta",
        """The paper describes lane normalization, empirical controls and RUVg before downstream testing.

- The exact **5,000 negative-control genes** and unwanted-factor count **k** are unavailable.
- A declared sensitivity uses all 11 FA/PM/HFD libraries with TMM/voom and limma, **without invented RUVg factors**.
- Shared source universe: 50,417 features. Exclude 91 HFD-only rows rather than invent zeros.
- GENCODE M13 coding filter + CPM filter yields 13,398 tested genes.""",
        """El artículo describe normalización por carril, controles empíricos y RUVg antes de las pruebas.

- No están disponibles los **5,000 genes de control negativo** exactos ni el número **k** de factores.
- Una sensibilidad declarada utiliza las 11 bibliotecas FA/PM/HFD con TMM/voom y limma, **sin inventar factores RUVg**.
- Universo común: 50,417 características. Se excluyen 91 filas exclusivas HFD sin inventar ceros.
- El filtro codificante GENCODE M13 y el filtro CPM dejan 13,398 genes evaluados.""",
        """The sensitivity keeps CPM>2 in at least three libraries, then requires BH<.05, absolute log2FC>.8 and AveExpr>1. The paper's main and supplementary abundance wording differ, and the original negative-control selection is ambiguous. No substitute can uniquely recover the unavailable RUVg model. Therefore the sensitivity is a methodological comparison, not a failed exact reconstruction and not proof that the biology is absent. The separate seven-animal primary PyDESeq2 analysis does not include HFD.""",
        """La sensibilidad retiene CPM>2 en al menos tres bibliotecas y exige BH<.05, |log₂FC|>.8 y AveExpr>1. La redacción de abundancia principal y suplementaria difiere, y la selección original de controles negativos es ambigua. Ningún sustituto recupera de forma única el RUVg ausente. Por eso se trata de una comparación metodológica, no una reproducción exacta fallida ni evidencia de ausencia biológica. El análisis PyDESeq2 principal de siete animales excluye HFD.""",
        ["HY","results/replication_pm/methods.json"])

    add("mirna-author-method","miRNA target reconstruction and testing gates","Reconstrucción del objetivo miRNA y criterios de prueba",
        """- Recover the prespecified Figure 6 target **mmu-miR-574-5p**.
- Preserve all nine arrays and the original normalized signal.
- Reconstruct ANOVA and pooled-MSE LSD pairwise comparisons.
- Distinguish an **unprotected pairwise test** from a **protected LSD procedure requiring omnibus significance**.
- Compare the target with the new moderated limma analysis and its array-wide BH family.""",
        """- Recuperar el objetivo preespecificado de la Figura 6: **mmu-miR-574-5p**.
- Conservar los nueve arrays y su señal normalizada original.
- Reconstruir ANOVA y comparaciones LSD con MSE agrupado.
- Distinguir la **prueba pareada no protegida** del **LSD protegido que exige significación ómnibus**.
- Comparar el objetivo con el análisis limma moderado nuevo y su familia BH de todo el array.""",
        """The paper gives general ANOVA/LSD language but not the complete GeneSpring array-specific gates or original heatmap feature selection. Consequently alternative readings of the testing sequence are shown. “Pairwise” refers to comparing independent groups, not paired animals. The original heatmap cannot be exactly reconstructed from an unknown row-selection rule. The arrays do not independently reproduce qPCR, ChIP, rescue, protein, BACE1 or behavior experiments.""",
        """El artículo describe ANOVA/LSD de forma general, sin aportar toda la configuración de GeneSpring ni la selección original del heatmap. Se muestran las lecturas alternativas de la secuencia de pruebas. «Entre pares de grupos» se refiere a grupos independientes, no animales pareados. No se puede reconstruir exactamente el heatmap si se desconoce la selección de filas. Los arrays no reproducen independientemente qPCR, ChIP, rescate, proteínas, BACE1 ni comportamiento.""",
        ["miRNA","results/replication_mirna_ozone/methods.json"])

    add("ozone-author-method","Ozone table reconstruction and cross-pipeline comparison","Reconstrucción de tablas de ozono y comparación entre flujos",
        """**Exact selection reconstruction**

The deposited E4/E3 table contains genotype effects. Apply its nominal p < 0.05 and |log₂FC| > 0.3 rule.

**Conditional effect concordance**

Compare exact-symbol finite effects from the new Salmon/tximport counts with published results. Reestimate candidate-model dispersions under their declared designs.

The original HISAT2/featureCounts pipeline, historical annotation and complete testing specification are not fully reconstructed.""",
        """**Reconstrucción exacta de selección**

La tabla depositada E4/E3 contiene efectos de genotipo. Aplicar su regla de p nominal < 0.05 y |log₂FC| > 0.3.

**Concordancia condicional de efectos**

Comparar efectos finitos con símbolo exacto obtenidos de los conteos Salmon/tximport nuevos con resultados publicados. Reestimar dispersiones en los modelos candidatos declarados.

No se reconstruyen íntegramente HISAT2/featureCounts, la anotación histórica ni toda la especificación de pruebas.""",
        """The original genotype CSV is not an exposure count matrix. Other original comparisons use an absolute effect threshold .1, so thresholds are matched to each actual contrast. Two candidate genotype models are reported rather than selecting the one with the closest paper count. The 17-month ozone E4/E3 S3 comparison is conditional on the published selected genes and the same source animals. Neither a high correlation on that list nor an UpSet overlap supplies independent replication or an interaction test.""",
        """El CSV original de genotipo no es una matriz de conteos de exposición. Otras comparaciones usan un umbral absoluto .1; se respetan los umbrales de cada contraste. Se informan dos modelos candidatos de genotipo, sin seleccionar por cercanía al número publicado. La comparación S3 E4/E3 bajo ozono a 17 meses está condicionada a los genes seleccionados y a los mismos animales fuente. Ni una correlación alta en esa lista ni un UpSet aportan replicación independiente o una prueba de interacción.""",
        ["O3","results/replication_mirna_ozone/methods.json"])

    add("count-model","Negative-binomial model for RNA counts","Modelo binomial negativo para conteos de ARN",
        r"""$$K_{gi}\sim\mathrm{NB}(\mu_{gi},\alpha_g),\qquad
\operatorname{Var}(K_{gi})=\mu_{gi}+\alpha_g\mu_{gi}^{2}$$

$$\log\mu_{gi}=\log o_{gi}+\mathbf x_i^\top\boldsymbol\beta_g$$

- $K_{gi}$: assigned counts or rounded estimated fragments.
- $\alpha_g$: gene-specific dispersion, estimated with shrinkage.
- $o_{gi}$: normalization factor. Ozone additionally carries effective-length information.
- The design matrix defines the comparison before the test.""",
        r"""$$K_{gi}\sim\mathrm{NB}(\mu_{gi},\alpha_g),\qquad
\operatorname{Var}(K_{gi})=\mu_{gi}+\alpha_g\mu_{gi}^{2}$$

$$\log\mu_{gi}=\log o_{gi}+\mathbf x_i^\top\boldsymbol\beta_g$$

- $K_{gi}$: conteos asignados o fragmentos estimados redondeados.
- $\alpha_g$: dispersión por gen, estimada con regularización.
- $o_{gi}$: factor de normalización. En ozono incorpora información de longitud efectiva.
- La matriz de diseño define la comparación antes de la prueba.""",
        """This slide uses the usual natural-log GLM notation. Reported DESeq2 coefficients and effects are converted to log2 units. For the two simple PM contrasts, o_gi reduces to a sample size factor estimated from the retained counts. For ozone, DESeqDataSetFromTximport constructs gene/sample normalization from effective lengths plus depth, so treating every offset as only library size would misdescribe the implemented model. The likelihood accommodates overdispersion but does not automatically control unmeasured confounding or infer cell composition.""",
        """La notación del GLM utiliza logaritmo natural. Los coeficientes y efectos reportados por DESeq2 se expresan en log₂. En los contrastes PM simples, o_gi se reduce al factor muestral estimado con los conteos retenidos. En ozono, DESeqDataSetFromTximport construye normalización gen/muestra a partir de longitudes efectivas y profundidad; describirla como mero tamaño de biblioteca sería incorrecto. La verosimilitud admite sobredispersión, pero no controla automáticamente confusión no medida ni infiere composición celular.""",
        ["DESeq2","PyDESeq2","results/ozone/methods.json"],"formula-slide")

    add("filters","Count filtering and normalization populations","Filtrado de conteos y poblaciones de normalización",
        table(["Analysis","Count filter","Retained genes","Fitted samples"],[
            ["HY primary","≥10 in ≥3","19,611","7"],
            ["ACA primary","≥10 in ≥4","19,235","8"],
            ["ACA author reconstruction","≥10 in ≥4","21,388","20"],
            ["O₃ primary","≥10 in ≥3","18,537","24"]])+
        "\nPM primary analyses use median-of-ratios size factors. Different fitted sample populations change the filter and dispersion estimates.",
        table(["Análisis","Filtro de conteos","Genes retenidos","Muestras ajustadas"],[
            ["HY principal","≥10 en ≥3","19,611","7"],
            ["ACA principal","≥10 en ≥4","19,235","8"],
            ["ACA reconstrucción autores","≥10 en ≥4","21,388","20"],
            ["O₃ principal","≥10 en ≥3","18,537","24"]])+
        "\nLos análisis PM principales usan factores por mediana de cocientes. Cambiar las muestras ajustadas modifica el filtro y las dispersiones.",
        """These filters are applied to their specified analytical sample sets. The same numerical count threshold can retain different gene universes. Hypothalamus input has 50,417 symbols, cortex 56,748 Ensembl IDs, and ozone 54,347 gene summaries. Exact source identifiers remain intact. Retention does not guarantee a valid p value: Cook filtering, independent filtering or numerical validity may still withhold inference. Never equate the retained-gene denominator with the finite BH-tested family without checking the specific table.""",
        """Los filtros se aplican a la población analítica especificada. El mismo umbral puede dejar universos distintos. La entrada hipotalámica tiene 50,417 símbolos; corteza, 56,748 IDs Ensembl; ozono, 54,347 resúmenes génicos. Se mantienen identificadores exactos. Retener un gen no garantiza p válido: Cook, filtrado independiente o validez numérica pueden retirar inferencia. No se equipara el denominador retenido con la familia finita BH sin revisar la tabla correspondiente.""",
        ["results/hypothalamus/methods.json","results/cortex/methods.json","results/replication_pm/methods.json","results/ozone/methods.json"])

    add("testing","Wald evidence, shrinkage and multiple testing","Evidencia Wald, regularización y pruebas múltiples",
        r"""$$z_g=\frac{\widehat\beta_g}{\operatorname{SE}(\widehat\beta_g)},
\qquad \mathrm{CI}_{95\%}=\widehat\beta_g\pm1.959964\,\operatorname{SE}(\widehat\beta_g)$$

$$q_{(i)}=\min\!\left(1,\min_{j\ge i}\frac{m\,p_{(j)}}{j}\right)$$

- Original MLEs, SEs and p values define the Wald/BH evidence.
- apeglm-style shrunken effects support visualization and ranking stability.
- Pointwise intervals are not simultaneous or selection-adjusted.
- Each declared test family has its own multiplicity correction.""",
        r"""$$z_g=\frac{\widehat\beta_g}{\operatorname{SE}(\widehat\beta_g)},
\qquad \mathrm{IC}_{95\%}=\widehat\beta_g\pm1.959964\,\operatorname{SE}(\widehat\beta_g)$$

$$q_{(i)}=\min\!\left(1,\min_{j\ge i}\frac{m\,p_{(j)}}{j}\right)$$

- MLE, SE y valores p originales definen la evidencia Wald/BH.
- Los efectos regularizados tipo apeglm apoyan visualización y estabilidad.
- Los intervalos puntuales no son simultáneos ni ajustados por selección.
- Cada familia de pruebas declarada tiene su corrección de multiplicidad.""",
        """The BH expression describes the monotone adjusted p value among the finite valid tests in that family. A shrunken mode is not substituted into the original Wald statistic. Primary PM PyDESeq2 uses its documented apeGLM-style adaptive prior, whereas ozone and author-aligned R workflows use apeglm in R. These implementations are described separately rather than claimed identical. A nominal interval excluding zero can coexist with BH q>.05. Increasing or changing the family can change the correction in nontrivial ways; the across-family outputs remain separate analyses.""",
        """La expresión BH describe el valor p ajustado monótono entre las pruebas finitas válidas de esa familia. No se sustituye la moda regularizada en el estadístico Wald original. PyDESeq2 usa su prior adaptativo tipo apeGLM documentado; ozono y los flujos R alineados con autores usan apeglm en R. No se declaran idénticas estas implementaciones. Un intervalo nominal que excluye cero puede coexistir con BH q>.05. Cambiar la familia modifica la corrección de forma no trivial; las salidas entre familias son análisis separados.""",
        ["DESeq2","apeglm","results/hypothalamus/methods.json"],"formula-slide")

    add("qc","Quality control and animal sensitivity","Control de calidad y sensibilidad a animales",
        """- For RNA-count inputs, verify nonnegative integers, exact IDs and sample assignments.
- Inspect library size, normalization, PCA and sample correlations.
- Preserve Cook-filtered inference and convergence flags.
- Retain all originally selected animals. PCA separation alone does not justify exclusion.
- Refit the HY primary model leaving one animal out at a time.
- Compare the fixed HY hypothesis across Python and R implementations.""",
        """- Para entradas de conteos de ARN, verificar enteros no negativos, IDs exactos y asignación muestral.
- Inspeccionar tamaño de biblioteca, normalización, PCA y correlaciones.
- Conservar indicadores de inferencia filtrada por Cook y convergencia.
- Retener los animales seleccionados. El PCA por sí solo no justifica exclusión.
- Reajustar HY omitiendo un animal por vez.
- Cotejar la hipótesis HY fija entre implementaciones Python y R.""",
        """QC plots diagnose rather than certify absence of batch or identity problems. HY uses design-aware VST whereas ACA uses blind VST; both PCAs use 500 variable genes with mean centering and no variance scaling. These differ from the original cortex code's TMM logCPM PCA with scaling. The HY Cook status withholds 159 gene tests, not 159 animals. Count replacement requires at least seven replicates in an eligible group and is unavailable in the small primary PM groups. Leave-one-out variation is reported without choosing an omission that improves significance.""",
        """Los gráficos QC diagnostican, sin certificar ausencia de lote o problemas de identidad. HY usa VST con diseño y ACA VST blind; ambos PCA usan 500 genes variables centrados, sin escalar por varianza. Difieren del PCA cortical original TMM logCPM escalado. Cook retira 159 pruebas génicas de HY, no 159 animales. Reemplazar conteos exige al menos siete réplicas en un grupo elegible, condición ausente en los grupos PM principales. Se informa la variación leave-one-out sin elegir una omisión favorable.""",
        ["results/hypothalamus/methods.json","results/hypothalamus/independent_deseq2_check/summary.json","results/cortex/methods.json"])

    add("mirna-model","Linear modeling of deposited miRNA signals","Modelado lineal de señales miRNA depositadas",
        r"""$$y_{fi}=\mathbf x_i^\top\boldsymbol\beta_f+\epsilon_{fi},
\qquad \Delta_f=\beta_{f,\mathrm{PM}}-\beta_{f,\mathrm{Control}}$$

- 1,268 input rows → 21 controls/unannotated + 2 unusable rows excluded → **1,245 tested features**.
- limma design <code>~0+condition</code>, two pairwise contrasts and a separate omnibus moderated F test.
- Empirical Bayes: trend = TRUE, robust = FALSE.
- Keep 823 low-detection flags and negative signal values.
- Effects remain **mean differences in deposited signal units**.""",
        r"""$$y_{fi}=\mathbf x_i^\top\boldsymbol\beta_f+\epsilon_{fi},
\qquad \Delta_f=\beta_{f,\mathrm{PM}}-\beta_{f,\mathrm{Control}}$$

- 1,268 filas → se excluyen 21 controles/no anotadas + 2 inutilizables → **1,245 características evaluadas**.
- Diseño limma <code>~0+condition</code>, dos contrastes y prueba F moderada ómnibus separada.
- Bayes empírico: trend = TRUE, robust = FALSE.
- Se conservan 823 indicadores de baja detección y las señales negativas.
- Los efectos son **diferencias medias en unidades de señal depositada**.""",
        """The deposited GeneSpring signal is quantile-normalized but its log2 preprocessing is not explicitly established. Negative values are retained because they are part of that supplied scale. A second transformation or a negative-binomial count model would invent a measurement model. Feature-level moderated SEs and t intervals use limma posterior residual variance and total degrees of freedom. BH runs separately across all retained features for each contrast and separately for the omnibus test. Detection flags provide context and do not remove features from the test family.""",
        """La señal GeneSpring depositada tiene normalización cuantílica, pero su transformación log₂ no está establecida explícitamente. Los negativos se conservan porque pertenecen a esa escala. Otra transformación o un modelo binomial negativo inventaría el modelo de medición. Las SE e intervalos t moderados usan varianza residual posterior y grados de libertad totales de limma. BH se aplica a todas las características retenidas por contraste y de forma separada al ómnibus. Los indicadores de detección contextualizan, sin excluir características de la familia.""",
        ["miRNA","limma","results/mirna/methods.json"],"formula-slide small")

    add("ozone-raw-qc","Ozone processing begins with all original reads","Ozono se procesa desde todas las lecturas originales",
        """1. Match 24 biological SourceNames to paired FASTQ runs and eight strata.
2. Inspect the first 200,000 read pairs per library with fastp.
3. Apply consistent adapter trimming to full libraries.
4. Retain full before/after QC. Disable quality, length, poly-G and duplicate filtering.
5. Process **532,781,354 paired fragments** across all 24 libraries.

Sampled Q30: 92.37–96.49%. All sequenced animals remain in the analysis.""",
        """1. Vincular 24 SourceNames biológicos con pares FASTQ y ocho estratos.
2. Inspeccionar los primeros 200,000 pares por biblioteca con fastp.
3. Recortar adaptadores de forma consistente en las bibliotecas completas.
4. Conservar QC antes/después. Desactivar filtros de calidad, longitud, poli-G y duplicados.
5. Procesar **532,781,354 fragmentos paired-end** de las 24 bibliotecas.

Q30 muestreado: 92.37–96.49%. Se conservan todos los animales secuenciados.""",
        """The SDRF has two rows per animal because of paired files; biological counts deduplicate SourceName, not FASTQ rows. The mapping is independently checked against deposited metadata. A strain-string exception occurs for Sample1/ERR14708304 and is retained as a source-quality flag, not represented as a complete genotype confound. Initial QC is a prefix inspection rather than a random sample. Its adapter findings motivate full-library trimming. The fastp adapter-dimer behavior is documented even though other filtering is disabled.""",
        """SDRF tiene dos filas por animal debido a los archivos pareados; el número biológico deduplica SourceName, no filas FASTQ. La correspondencia se coteja independientemente con metadatos depositados. Sample1/ERR14708304 tiene una cadena de cepa distinta; se conserva como indicador de calidad, sin declararla confusión completa por genotipo. El QC inicial inspecciona un prefijo, no una muestra aleatoria. Sus adaptadores motivan recorte de bibliotecas completas. Se documenta la conducta de fastp ante dímeros aunque los demás filtros estén desactivados.""",
        ["results/ozone/methods.json","results/ozone/summary.json","results/replication_mirna_ozone/sample_mapping_validation.json"])

    add("ozone-quant","Selective alignment and effective-length offsets","Alineamiento selectivo y offsets de longitud efectiva",
        """**Salmon 2.5.1**

GENCODE M25 transcriptome + whole GRCm38 primary-genome decoys. Automatic library type, sequence and GC bias models.

**tximport**

Estimated fragments and sample-specific abundance-weighted effective gene lengths. <code>countsFromAbundance = "no"</code>.

**DESeqDataSetFromTximport**

Gene/sample length normalization offsets enter inference. Mapping ranges from **72.31% to 88.33%**.""",
        """**Salmon 2.5.1**

Transcriptoma GENCODE M25 + genoma primario GRCm38 completo como decoys. Tipo de biblioteca automático y modelos de sesgo de secuencia y GC.

**tximport**

Fragmentos estimados y longitudes génicas efectivas ponderadas por abundancia por muestra. <code>countsFromAbundance = "no"</code>.

**DESeqDataSetFromTximport**

Los offsets de longitud gen/muestra entran en la inferencia. El mapeo varía de **72.31% a 88.33%**.""",
        """This workflow differs deliberately from the original HISAT2/featureCounts pipeline. Selective alignment with genomic decoys reduces spurious transcript assignments but is not full reproduction of historical gene counts. TPM values are not used as a substitute for negative-binomial input counts. Effective gene lengths vary with transcript usage and enter the offsets. Mouse-reference counts cannot resolve human APOE allele-specific expression. Versions, reference hashes and all quantification checkpoints are preserved in the repository.""",
        """Este flujo difiere deliberadamente del HISAT2/featureCounts original. El alineamiento selectivo con decoys genómicos limita asignaciones transcriptómicas espurias, sin reproducir íntegramente los conteos históricos. TPM no sustituye los conteos de entrada binomial negativa. La longitud efectiva varía con el uso de transcritos y entra en los offsets. La referencia de ratón no resuelve expresión alélica humana de APOE. Se conservan versiones, hashes de referencia y puntos de control de cuantificación.""",
        ["Salmon","tximport","results/ozone/methods.json"])

    add("ozone-factorial","Eight-parameter ozone factorial model","Modelo factorial de ozono con ocho parámetros",
        r"""$$\log\mu=\log o+\beta_0+\beta_A A+\beta_G G+\beta_E E
+\beta_{AG}AG+\beta_{AE}AE+\beta_{GE}GE+\beta_{AGE}AGE$$

Define $A=1$ for 17-month onset, $G=1$ for APOE4 and $E=1$ for O₃.

- Full design: age × genotype × exposure.
- 24 animals, 3 in each of 8 strata. Design rank = 8.
- 18,537 retained genes.
- Fixed prefilter, Cook checks and validity masks. Adaptive independent filtering is disabled.""",
        r"""$$\log\mu=\log o+\beta_0+\beta_A A+\beta_G G+\beta_E E
+\beta_{AG}AG+\beta_{AE}AE+\beta_{GE}GE+\beta_{AGE}AGE$$

$A=1$: inicio a 17 meses; $G=1$: APOE4; $E=1$: O₃.

- Diseño completo: edad × genotipo × exposición.
- 24 animales, 3 en cada uno de 8 estratos. Rango del diseño = 8.
- 18,537 genes retenidos.
- Prefiltro fijo, verificaciones Cook y máscaras de validez. Sin filtrado independiente adaptativo.""",
        """The model has residual information from 24 animals rather than eight group means. The mathematical GLM uses natural-log coefficients; reported DESeq2 response contrasts are in log2 units. Genes all-zero within a required stratum have the corresponding effect and p value withheld before BH, even if a numerical coefficient exists. No batch or RIN covariates were invented because deposited annotations do not supply them. The factorial analysis is a project reanalysis, not an exact reconstruction of the original authors' model.""",
        """El modelo utiliza información de 24 animales, no solo ocho medias. El GLM matemático emplea coeficientes en logaritmo natural; los contrastes DESeq2 se reportan en log₂. Si todos los conteos son cero en un estrato necesario, se retiran su efecto y p antes de BH, aunque exista un coeficiente numérico. No se inventan covariables de lote o RIN ausentes. El factorial es un reanálisis del proyecto, no una reconstrucción exacta del modelo original.""",
        ["O3","results/ozone/methods.json"],"formula-slide")

    add("ozone-contrasts","Stratum effects and covariance-aware differences","Efectos por estrato y diferencias con covarianza",
        r"""| O₃ versus FA effect | Linear combination |
|---|---|
| 3M / E3 | $\beta_E$ |
| 3M / E4 | $\beta_E+\beta_{GE}$ |
| 17M / E3 | $\beta_E+\beta_{AE}$ |
| 17M / E4 | $\beta_E+\beta_{AE}+\beta_{GE}+\beta_{AGE}$ |

$$\widehat\Delta=\mathbf c^\top\widehat{\boldsymbol\beta},
\qquad \operatorname{SE}(\widehat\Delta)
=\sqrt{\mathbf c^\top\widehat V\mathbf c}$$

A difference between response effects requires its own contrast.""",
        r"""| Efecto O₃ frente a FA | Combinación lineal |
|---|---|
| 3M / E3 | $\beta_E$ |
| 3M / E4 | $\beta_E+\beta_{GE}$ |
| 17M / E3 | $\beta_E+\beta_{AE}$ |
| 17M / E4 | $\beta_E+\beta_{AE}+\beta_{GE}+\beta_{AGE}$ |

$$\widehat\Delta=\mathbf c^\top\widehat{\boldsymbol\beta},
\qquad \operatorname{SE}(\widehat\Delta)
=\sqrt{\mathbf c^\top\widehat V\mathbf c}$$

Comparar respuestas exige un contraste específico.""",
        """The covariance term matters: subtracting two estimated effects and treating their SEs as independent would be incorrect because both derive from the same fitted model. An equivalent eight-group parameterization supplies the single within-stratum coefficients used for apeglm display. New secondary response-difference tests preserve the original fitted dispersions and effective-length offsets. Their pointwise Wald intervals are not simultaneous and are not intervals for atlas cells.""",
        """La covarianza importa: restar dos efectos suponiendo SE independientes sería incorrecto, pues ambos proceden del mismo modelo. Una parametrización equivalente de ocho grupos aporta coeficientes individuales para el despliegue apeglm. Los análisis secundarios de diferencias preservan dispersiones y offsets de longitud originales. Sus intervalos Wald puntuales no son simultáneos ni corresponden a células del atlas.""",
        ["results/ozone/methods.json","metadata/advanced_omics/methods.json"],"formula-slide")

    add("ozone-lrt","Global exposure and response-heterogeneity tests","Pruebas globales de exposición y heterogeneidad",
        r"""$$D=2(\ell_{\mathrm{full}}-\ell_{\mathrm{reduced}})
\ \overset{H_0}{\approx}\ \chi^2_{\Delta df}$$

| Test | Removed terms under the null | df |
|---|---|---|
| Any exposure effect | E, AE, GE, AGE | 4 |
| Genotype-dependent exposure | GE, AGE | 2 |
| Age-dependent exposure | AE, AGE | 2 |
| Three-way interaction | AGE | 1 |

Separate BH families are reported, plus supplementary pooled-family corrections.""",
        r"""$$D=2(\ell_{\mathrm{completo}}-\ell_{\mathrm{reducido}})
\ \overset{H_0}{\approx}\ \chi^2_{\Delta df}$$

| Prueba | Términos retirados bajo la nula | gl |
|---|---|---|
| Cualquier efecto de exposición | E, AE, GE, AGE | 4 |
| Respuesta dependiente de genotipo | GE, AGE | 2 |
| Respuesta dependiente de edad | AE, AGE | 2 |
| Interacción triple | AGE | 1 |

Se informan familias BH separadas y correcciones suplementarias conjuntas.""",
        """The global ozone LRT asks whether any exposure-associated coefficient is needed. The three secondary LRTs ask narrower response-heterogeneity questions. Their hypotheses, degrees of freedom and multiplicity families differ. A significant response in one stratum and a nonsignificant response in another do not demonstrate interaction. The signed single-contrast Wald tests answer different alternatives from multi-parameter LRTs and therefore need not yield nested discovery sets.""",
        """El LRT global pregunta si se necesita algún coeficiente asociado con exposición. Los tres LRT secundarios plantean preguntas más específicas sobre heterogeneidad. Sus hipótesis, grados de libertad y familias de multiplicidad difieren. Significación en un estrato y ausencia en otro no demuestra interacción. Los Wald con signo responden alternativas distintas de los LRT multiparámetro; no se exigen conjuntos de hallazgos anidados.""",
        ["results/ozone/methods.json","metadata/advanced_omics/methods.json"],"formula-slide")

    add("ora","GO over-representation uses a testable universe","La sobrerrepresentación GO utiliza un universo evaluable",
        r"""$$P(X\ge k)=\sum_{x=k}^{\min(K,n)}
\frac{\binom Kx\binom{N-K}{n-x}}{\binom Nn}$$

- $N$: unique mapped genes with finite adjusted DE p values.
- $n$: BH-significant foreground, up/down separately by MLE sign.
- $K$: GO-BP term members in that universe; $k$: overlap.
- Test 15–500-member terms. Keep zero-overlap terms at p = 1.
- BH across all eligible terms within each contrast/direction family.""",
        r"""$$P(X\ge k)=\sum_{x=k}^{\min(K,n)}
\frac{\binom Kx\binom{N-K}{n-x}}{\binom Nn}$$

- $N$: genes únicos mapeados con p ajustado DE finito.
- $n$: conjunto BH significativo, separado por signo MLE.
- $K$: miembros GO-BP del universo; $k$: intersección.
- Evaluar términos de 15–500 miembros. Conservar los de intersección cero con p = 1.
- BH entre todos los términos elegibles por contraste/dirección.""",
        """Official frozen mouse GO Biological Process annotations use taxon10090, exclude NOT and propagate only is_a and part_of. Exact symbols or unique aliases map to canonical genes. Duplicate features choose the highest baseMean representative before looking at significance. The universe includes testable mapped genes without specific BP annotation, rather than restricting to well-annotated genes. These choices do not remove all detectability or gene-length bias. No miRNA target enrichment is claimed without an independently supported target map.""",
        """Las anotaciones oficiales GO Biological Process de ratón usan taxón10090, excluyen NOT y propagan solo is_a y part_of. Se mapean símbolos exactos o alias únicos. Los duplicados eligen el representante de mayor baseMean antes de examinar significación. El universo incluye genes mapeados evaluables sin BP específico, sin limitarse a genes bien anotados. Esto no elimina todo sesgo de detectabilidad o longitud. No se infiere enriquecimiento de blancos miRNA sin un mapa respaldado independientemente.""",
        ["results/enrichment/methods.json"],"formula-slide")

    add("gsea","Ranked enrichment preserves weak coordinated signals","El enriquecimiento por rango conserva señales coordinadas débiles",
        r"""$$RS(t)=
\frac{\sum_{j\le t,\ j\in S}|r_j|}{\sum_{j\in S}|r_j|}
-\frac{\#\{j\le t:j\notin S\}}{N-|S|}$$

- Rank $r_j$: unshrunken Wald statistic with finite original p value.
- No gene-level BH or effect-size preselection.
- GSEApy prerank, weight 1, 1,000 gene-set permutations, fixed seed.
- NES and **algorithmic GSEA FDR** differ from ORA/BH.
- Shared genes and GO ancestry make enriched terms correlated.""",
        r"""$$RS(t)=
\frac{\sum_{j\le t,\ j\in S}|r_j|}{\sum_{j\in S}|r_j|}
-\frac{\#\{j\le t:j\notin S\}}{N-|S|}$$

- Rango $r_j$: estadístico Wald no regularizado con p original finito.
- Sin preselección por BH génico o tamaño de efecto.
- GSEApy prerank, peso 1, 1,000 permutaciones de conjuntos y semilla fija.
- NES y **FDR algorítmica GSEA** difieren de ORA/BH.
- Compartir genes y ancestros GO correlaciona los términos.""",
        """The displayed formula gives the weight1 running enrichment score before sign-dependent extremum selection and normalization. Original p values must be finite even when a Wald statistic remains finite, to exclude Cook-filtered unreliable tests. Genes with a finite original p but missing adjusted p from independent filtering can remain in the rank. These are gene-set/label permutations, not phenotype permutations, and do not propagate full sample-level uncertainty or preserve gene–gene correlation. Permutation q=0 is a finite-simulation estimate, not a zero probability of error. Enrichment supports exploratory ranked associations, not demonstrated pathway activation or suppression.""",
        """La fórmula muestra el puntaje acumulado con peso1 antes de seleccionar su extremo con signo y normalizar. Se exige p original finito aunque persista estadístico Wald, para excluir pruebas no fiables filtradas por Cook. Un gen con p original finito y p ajustado ausente por filtrado independiente puede seguir en el rango. Son permutaciones de conjuntos/etiquetas génicas, no de fenotipos; no propagan toda la incertidumbre muestral ni preservan correlación entre genes. q=0 es una estimación de simulación finita, no probabilidad de error nula. El enriquecimiento apoya asociaciones exploratorias de rango, sin demostrar activación o supresión de vías.""",
        ["results/enrichment/methods.json"],"formula-slide")

    add("atlas","A measured spatial reference from one donor","Una referencia espacial medida de un donante",
        table(["Reference feature","ACA","HY"],[
            ["Measured cells","14,582","37,811"],
            ["Physical sections","18","20"],
            ["Cell classes","16","26"],
            ["Anatomical subdivisions","10","49"],
            ["Occupied 100 µm bins","1,656","3,772"]])+
        "\nZhuang-ABCA-2: one unexposed adult male, one registered hemisphere, **1,122 measured MERFISH genes**.",
        table(["Característica","ACA","HY"],[
            ["Células medidas","14,582","37,811"],
            ["Cortes físicos","18","20"],
            ["Clases celulares","16","26"],
            ["Subdivisiones anatómicas","10","49"],
            ["Casillas ocupadas de 100 µm","1,656","3,772"]])+
        "\nZhuang-ABCA-2: un macho adulto sin exposición, un hemisferio registrado y **1,122 genes MERFISH medidos**.",
        """The reference is selected from a larger atlas, but the project uses this specific donor to preserve cell identity, imaging counts and coordinates. All 52,393 ACA/HY cells are retained without subsampling. ACA includes ACAd and ACAv, each with layers1,2/3,5,6a,6b. HY has 49 nuclei/subregion labels, including HY-unassigned and ZI-unassigned, so it is incorrect to call all49 named nuclei. The atlas annotations originate from the source workflow; reference expression and metadata are not matched controls for the source pollution animals.""",
        """La referencia proviene de un atlas mayor, pero el proyecto usa este donante específico para conservar identidad celular, conteos de imagen y coordenadas. Se mantienen las 52,393 células ACA/HY sin submuestreo. ACA incluye ACAd y ACAv con capas1,2/3,5,6a,6b. HY tiene49 etiquetas de núcleos/subregiones, incluidas HY-unassigned y ZI-unassigned; no son49 núcleos nombrados. Las anotaciones provienen del flujo fuente y la referencia no constituye controles emparejados de los animales expuestos.""",
        ["atlas","results/spatial_sections/summary.json","results/spatial_cortex/summary.json"])

    add("coordinates","Deposited section coordinates and registered anatomy","Coordenadas depositadas del corte y anatomía registrada",
        """**Within-section maps**

Deposited experimental x/y after source tissue rotation/alignment. Experimental z estimates acquisition AP.

**3D anatomy and Gaussian basis**

Registered CCF AP/DV/ML in millimeters. The 3D display uses ML, AP and −DV.

Nominal coronal spacing is about 200 µm. Surface mesh steps of 50–150 µm describe rendering, not expression acquisition.""",
        """**Mapas dentro del corte**

x/y experimentales depositadas tras rotación/alineación del tejido. z experimental estima AP de adquisición.

**Anatomía 3D y base gaussiana**

CCF AP/DV/ML registradas en milímetros. La visualización usa ML, AP y −DV.

La separación coronal nominal es aproximadamente 200 µm. Los pasos de malla de 50–150 µm describen representación, no adquisición de expresión.""",
        """Experimental x/y are not interchangeable with CCF ML/DV. The 3D mesh derives from the Allen CCF25µm annotation via marching cubes, with rendering steps wholebrain.15mm,HY.075mm,ACA.05mm. Increasing mesh density does not add measured RNA resolution. Cell coordinates do not use the centered RAS NIfTI affine when computing source CCF distances. Coordinate-parcellation agreement is94.05% in10,000 HY diagnostic cells and93.54% acrossACA; this is a consistency diagnostic, not ground-truth registration accuracy. No missing hemisphere is mirrored.""",
        """Las x/y experimentales no equivalen a CCF ML/DV. La malla3D deriva de la anotación AllenCCF25µm con marching cubes y pasos de representación .15mm cerebro,.075mmHY,.05mmACA. Una malla más densa no agrega resolución de ARN medida. No se usa la transformación RAS centrada del NIfTI para distancias CCF fuente. La concordancia coordenada-parcelación es94.05% en10,000 célulasHY diagnósticas y93.54% deACA; es consistencia, no exactitud frente a una verdad de referencia. No se refleja el hemisferio ausente.""",
        ["metadata/spatial_sections/methods.json","results/spatial_viewer/summary.json","results/spatial_cortex/summary.json"])

    add("gene-selection","Gene priorities and measured-panel coverage","Priorización génica y cobertura del panel medido",
        """1. Preserve the original primary regional BH q and effect estimates.
2. Rank candidates by q, then source symbol or exact ID.
3. Map HY by exact source symbol and ACA by exact Ensembl ID.
4. Annotate products and functions with frozen official NCBI/UniProt records.

**Panel coverage:** 1 of 11 HY primary hits and 76 of 1,167 ACA primary hits are physically measured.

The HY displayed top ten include nine exploratory candidates. All ten ACA candidates satisfy primary BH q < 0.05.""",
        """1. Conservar BH q y efectos regionales principales originales.
2. Ordenar por q y después símbolo fuente o ID exacto.
3. Mapear HY por símbolo fuente exacto y ACA por ID Ensembl exacto.
4. Anotar productos y funciones con registros oficiales NCBI/UniProt congelados.

**Cobertura:** se mide físicamente 1 de 11 hallazgos HY y 76 de 1,167 hallazgos ACA.

El top diez HY incluye nueve candidatos exploratorios. Los diez ACA cumplen BH q < 0.05 principal.""",
        """This is transparent evidence ranking, not a causal-importance score. No current alias is silently substituted for the original mapping. Current annotation is shown alongside source identity. HY lacks original Ensembl IDs, which limits cross-resource identity resolution. Six HY measured candidates show strong shrinkage. The gene-product table includes noncoding features and pseudogenes explicitly; absent atlas genes are not imputed. Forty-six distinct genes produce51 selection rows because genes can enter more than one regional/selection group.""",
        """Es un ordenamiento transparente de evidencia, no un puntaje de importancia causal. Los alias actuales no reemplazan silenciosamente la identidad original. La anotación actual se muestra junto a la fuente. HY carece de IDs Ensembl originales, lo que limita vinculación entre recursos. Seis candidatosHY medidos presentan regularización fuerte. Se identifican ARN no codificante y pseudogenes; los genes ausentes del atlas no se imputan. Cuarenta y seis genes distintos producen51 filas por la pertenencia a grupos regionales/de selección.""",
        ["metadata/gene_priorities/methods.json","results/gene_priorities/summary.json"])

    add("detection-method","Detection fractions retain measured zeros","Las fracciones de detección conservan ceros medidos",
        r"""$$d_{ig}=\mathbf1(c_{ig}>0),\qquad
\widehat f_{gA}=\frac{\sum_{i\in A}d_{ig}}{n_A}$$

$$b_x=\lfloor x/0.1\rfloor,\qquad b_y=\lfloor y/0.1\rfloor$$

- $A$: selected physical section, anatomy, cell class or occupied bin.
- The denominator includes every selected measured cell, including zeros.
- Bins are 0.1 × 0.1 mm and remain separate by physical section.
- Empty bins are unobserved. No spatial interpolation or density claim is added.""",
        r"""$$d_{ig}=\mathbf1(c_{ig}>0),\qquad
\widehat f_{gA}=\frac{\sum_{i\in A}d_{ig}}{n_A}$$

$$b_x=\lfloor x/0.1\rfloor,\qquad b_y=\lfloor y/0.1\rfloor$$

- $A$: corte, anatomía, clase celular o casilla ocupada seleccionados.
- El denominador incluye todas las células medidas, incluidos los ceros.
- Las casillas de 0.1 × 0.1 mm se separan por corte físico.
- Las casillas vacías no están observadas. No se agrega interpolación ni una afirmación de densidad.""",
        """A raw imaging zero means no detected count in this assay, not proven biological absence. Bins use half-open intervals at global experimental-coordinate multiples. A one-cell edge section remains visible with its denominator. An empty filter combination has unavailable detection, not zero detection. Static section figures and3D use log2(1+count), whereas the interactive count map and Gaussian working response use natural ln(1+count). Every transform is explicitly labeled. Descriptive detection summaries do not receive an independent-cell exposure p value.""",
        """Cero de imagen significa conteo no detectado en este ensayo, no ausencia biológica demostrada. Las casillas usan intervalos semiabiertos en múltiplos globales experimentales. El corte de borde con una célula conserva su denominador. Una combinación vacía de filtros tiene detección no disponible, no0%. Las figuras estáticas y3D usan log₂(1+conteo); el mapa interactivo y la respuesta gaussiana usan ln(1+conteo). Cada transformación se rotula. La detección descriptiva no recibe p de exposición con células independientes.""",
        ["metadata/spatial_sections/methods.json"],"formula-slide")

    add("scenario-method","Uniform regional PM₂.₅ scenario","Escenario regional uniforme con PM₂.₅",
        r"""$$E_{\mathrm{scenario}}(c,g)
=E_{\mathrm{reference}}(c,g)\,
2^{\widehat\beta^{\,\mathrm{shrunk}}_g}$$

| Example | Regional shrunken log₂FC | Multiplier |
|---|---:|---:|
| HY / Tll1 | −1.22405 | 0.42808 |
| ACA / Fosb | +0.99711 | 1.99600 |

Cell positions and numbers remain fixed. Shared per-gene color limits support comparison.

This is an illustrative propagation of a regional effect into an unexposed reference.""",
        r"""$$E_{\mathrm{escenario}}(c,g)
=E_{\mathrm{referencia}}(c,g)\,
2^{\widehat\beta^{\,\mathrm{regularizada}}_g}$$

| Ejemplo | log₂FC regional regularizado | Multiplicador |
|---|---:|---:|
| HY / Tll1 | −1.22405 | 0.42808 |
| ACA / Fosb | +0.99711 | 1.99600 |

Se mantienen posiciones y número de células. Una escala común por gen permite comparar.

Es una propagación ilustrativa del efecto regional a una referencia sin exposición.""",
        """Uniform scaling is underidentified biologically: the same bulk effect could arise from cell-specific regulation or composition. This display does not estimate which cells respond, cell loss, migration or a causal mechanism. Wald MLE±1.96SE provides alternative group-level scenarios, separate from the shrunken default, not per-cell confidence intervals. A shared optional99th-percentile color cap affects the display only; numeric values are retained. An arbitrary smoother or mirrored hemisphere would add assumptions, not observations. The atlas donor is not a matched experimental control.""",
        """El escalamiento uniforme no está identificado biológicamente: el mismo efecto bulk puede surgir de regulación específica o composición. No estima células respondedoras, pérdida, migración ni mecanismo causal. MLE Wald±1.96SE aporta escenarios grupales alternativos separados del valor regularizado por defecto, no intervalos celulares. El límite opcional común del percentil99 afecta solo el color; conserva los valores. Suavizar o reflejar un hemisferio agregaría supuestos, no observaciones. El donante del atlas no es control experimental emparejado.""",
        ["results/spatial_viewer/summary.json","results/spatial_cortex/summary.json"],"formula-slide")

    add("gaussian-features","Gaussian baseline with fixed spatial basis functions","Referencia gaussiana con bases espaciales fijas",
        r"""$$y_{ig}=\ln(1+c_{ig})
=\mathbf x_i^\top\boldsymbol\beta_g+\epsilon_{ig}$$

$$\phi_j(\mathbf r)=
\exp\!\left[-\frac12\sum_k\frac{(r_k-m_{jk})^2}{\ell_k^2}\right],
\qquad\boldsymbol\ell=(0.60,0.45,0.45)\ \mathrm{mm}$$

Features: intercept, anatomical indicators, three standardized CCF coordinates, 36 anisotropic RBFs and section nuisance intercepts.

No expression-derived cell-class predictors enter the fit.""",
        r"""$$y_{ig}=\ln(1+c_{ig})
=\mathbf x_i^\top\boldsymbol\beta_g+\epsilon_{ig}$$

$$\phi_j(\mathbf r)=
\exp\!\left[-\frac12\sum_k\frac{(r_k-m_{jk})^2}{\ell_k^2}\right],
\qquad\boldsymbol\ell=(0.60,0.45,0.45)\ \mathrm{mm}$$

Variables: intercepto, indicadores anatómicos, tres coordenadas CCF estandarizadas, 36 RBF anisótropas e interceptos de corte.

No se usan predictores de clase celular derivados de expresión.""",
        """RBF means radial basis function. The36centers form4×3×3 training-coordinate quantiles: AP.05/.35/.65/.95 andDV/ML.1/.5/.9. Lengths correspond to AP,DV,ML and remain fixed before held-out evaluation. Standardization, RBF centering and anatomy vocabularies use training data only. Section effects absorb nuisance variation but cannot distinguish every technical source from biology. Excluding expression-derived cell labels reduces one route for target-transcript leakage. This model is finite-basis Bayesian regression, not a fitted full Gaussian process or an implementation of Tangram, BayesSpace or cell2location.""",
        """RBF significa función de base radial. Los36centros forman cuantiles de entrenamiento4×3×3: AP.05/.35/.65/.95 yDV/ML.1/.5/.9. Las longitudes AP,DV,ML se fijan antes de evaluar. Estandarización, centrado RBF y vocabulario anatómico usan solo entrenamiento. Los efectos de corte absorben variación, sin separar toda fuente técnica de biología. Excluir etiquetas celulares derivadas de expresión reduce filtración del transcrito objetivo. Es regresión bayesiana con bases finitas, no un proceso gaussiano completo ni Tangram, BayesSpace o cell2location.""",
        ["metadata/spatial_bayes/methods.json"],"formula-slide")

    add("gaussian-priors","Proper priors and analytical posterior","Priors propios y posterior analítica",
        r"""$$\boldsymbol\beta\mid\sigma^2\sim
\mathcal N(0,\sigma^2V_0),\qquad \sigma^2\sim IG(3,1)$$

$$V_n=(X^\top X+V_0^{-1})^{-1},\quad
m_n=V_nX^\top y,\quad a_n=3+n/2$$

$$b_n=1+\frac12\left(y^\top y-m_n^\top V_n^{-1}m_n\right)$$

Prior variance multipliers: intercept 25, anatomy 1, linear coordinates 0.25, each RBF 1/36, each section 0.25.

Cholesky computation checks the normal equations and posterior scale.""",
        r"""$$\boldsymbol\beta\mid\sigma^2\sim
\mathcal N(0,\sigma^2V_0),\qquad \sigma^2\sim IG(3,1)$$

$$V_n=(X^\top X+V_0^{-1})^{-1},\quad
m_n=V_nX^\top y,\quad a_n=3+n/2$$

$$b_n=1+\frac12\left(y^\top y-m_n^\top V_n^{-1}m_n\right)$$

Multiplicadores de varianza previa: intercepto25, anatomía1, coordenadas0.25, cada RBF1/36 y cada corte0.25.

El cálculo Cholesky verifica ecuaciones normales y escala posterior.""",
        """The inverse-Gamma convention uses shape a and scale b in a density proportional to (sigma²)^(-a-1)exp(-b/sigma²). Section random effects are represented as proper Gaussian-prior coefficients with variance .25sigma². All priors are modeling assumptions, not empirically validated neural priors. The implementation checks the residual-plus-prior quadratic representation of b_n to reduce numerical cancellation concerns. It also checks covariance symmetry/inversion. Posterior integration is analytic conditional on this model, not proof that the Gaussian likelihood fits sparse nonnegative counts.""",
        """La inversa-Gamma usa forma a y escala b con densidad proporcional a (sigma²)^(-a-1)exp(-b/sigma²). Los efectos aleatorios de corte son coeficientes con prior gaussiano propio de varianza.25sigma². Los priors son supuestos de modelado, no priors neuronales validados empíricamente. La implementación verifica b_n mediante la suma cuadrática residual más prior para revisar cancelación numérica, además de simetría e inversión de covarianza. La integración es analítica condicional al modelo, sin demostrar adecuación gaussiana a conteos no negativos dispersos.""",
        ["metadata/spatial_bayes/methods.json","results/spatial_bayes/validation.json"],"formula-slide")

    add("gaussian-intervals","Credible and new-section predictive intervals","Intervalos creíbles y predictivos de un corte nuevo",
        r"""Section-adjusted latent mean:

$$f_*\mid y\sim t_{2a_n}\left(x_*^\top m_n,\
\sqrt{\frac{b_n}{a_n}x_*^\top V_nx_*}\right)$$

New cell in a new physical section:

$$y_*\mid y\sim t_{2a_n}\left(x_*^\top m_n,\
\sqrt{\frac{b_n}{a_n}\left(1.25+x_*^\top V_nx_*\right)}\right)$$

The second arguments are **Student-t scales**, not variances.

Maps show the full-reference fit at measured locations. Held-out prediction is evaluated separately.""",
        r"""Media latente ajustada por corte:

$$f_*\mid y\sim t_{2a_n}\left(x_*^\top m_n,\
\sqrt{\frac{b_n}{a_n}x_*^\top V_nx_*}\right)$$

Célula nueva en un corte físico nuevo:

$$y_*\mid y\sim t_{2a_n}\left(x_*^\top m_n,\
\sqrt{\frac{b_n}{a_n}\left(1.25+x_*^\top V_nx_*\right)}\right)$$

Los segundos argumentos son **escalas t de Student**, no varianzas.

Los mapas muestran ajuste con toda la referencia en posiciones medidas. La predicción reservada se evalúa aparte.""",
        """Latent maps average the nuisance section offset to zero rather than reporting each fitted section shift. The predictive1.25 combines residual variance1 and independent new-section variance.25. An unseen test anatomy adds prior variance1 and mean0, explicitly rather than silently substituting a known anatomy. Negative means, bounds and predictive support are retained as diagnostics of the Gaussian working approximation. Estimates are not exponentiated into expected molecule counts. Uncertainty is conditional on the basis, priors, anatomy and registration and omits unmodeled donor and registration uncertainty.""",
        """Los mapas latentes promedian a cero el desplazamiento de corte, sin informar el corrimiento ajustado de cada uno. El1.25 predictivo suma varianza residual1 y de corte nuevo.25. Una anatomía no vista añade varianza previa1 y media0 explícitas, sin sustituirse por una conocida. Se conservan medias, límites y soporte negativos como diagnóstico de la aproximación. No se aplica exponencial para afirmar conteos esperados. La incertidumbre depende de bases, priors, anatomía y registro, y omite variación no modelada entre donantes y registro.""",
        ["metadata/spatial_bayes/methods.json"],"formula-slide")

    fig("guarded-cv","Whole-section validation with neighboring guards","Validación por cortes completos con separación vecina","spatial_guarded_folds",
        "Five contiguous test blocks. Adjacent represented sections are guards, and every measured cell is tested once.",
        "Cinco bloques contiguos de prueba. Cortes representados vecinos actúan como separación y cada célula se evalúa una vez.",
        """Sections are ordered by experimental acquisition AP, while the Gaussian predictors use registered CCF coordinates. The nearest represented section on each side of a test block is removed from training. If the region is absent from intervening acquisition cuts, this guard can span more than one physical acquisition step. Anatomy vocabulary, scales and basis centers are learned only from remaining training cells. Fixed priors and basis lengths are not selected from test performance. Both Bayesian models use these identical folds. This is within-donor reference prediction, not exposure or cross-donor validation.""",
        """Los cortes se ordenan por AP experimental de adquisición; los predictores gaussianos usan CCF registrado. Se retira del entrenamiento el corte representado más próximo de cada lado. Si la región está ausente en cortes intermedios, la separación puede abarcar más de un paso de adquisición. Vocabulario anatómico, escalas y centros se aprenden con el entrenamiento restante. Priors y longitudes no se eligen por rendimiento de prueba. Ambos modelos usan los mismos bloques. Es predicción de referencia dentro del donante, sin validar exposición ni generalización entre donantes.""",
        ["results/spatial_bayes/ACA/section_folds.tsv","results/spatial_bayes/HY/section_folds.tsv"])

    add("beta-model","Hierarchical model of bounded detection probability","Modelo jerárquico de probabilidad de detección acotada",
        r"""$$Y_{sag}\mid p_{sag}\sim
\operatorname{Binomial}(n_{sa},p_{sag})$$

$$p_{sag}\mid\mu_{ag},\kappa_{ag}\sim
\operatorname{Beta}(\mu_{ag}\kappa_{ag},(1-\mu_{ag})\kappa_{ag})$$

$$\mu_{ag}\sim\operatorname{Beta}(0.5,0.5),\qquad
\ln\kappa_{ag}\sim\mathcal N(\ln20,1.5^2)$$

Genes/anatomies have separate fits. Section propensities are exchangeable within each fit. Probabilities stay between 0 and 1.""",
        r"""$$Y_{sag}\mid p_{sag}\sim
\operatorname{Binomial}(n_{sa},p_{sag})$$

$$p_{sag}\mid\mu_{ag},\kappa_{ag}\sim
\operatorname{Beta}(\mu_{ag}\kappa_{ag},(1-\mu_{ag})\kappa_{ag})$$

$$\mu_{ag}\sim\operatorname{Beta}(0.5,0.5),\qquad
\ln\kappa_{ag}\sim\mathcal N(\ln20,1.5^2)$$

Se ajustan genes/anatomías por separado. Las propensiones de corte son intercambiables dentro de cada ajuste. Las probabilidades permanecen entre 0 y 1.""",
        """Y is positive measured cells and n includes all measured zeros in a section×anatomy group. Mu is the mean of exchangeable section propensities, not a pooled-cell weighted detection fraction. Kappa controls section heterogeneity; the beta-binomial intraclass correlation is1/(kappa+1) conditional on the hyperparameters. Priors are proper and fixed before held-out outcomes. The model does not include an AP trend, local spatial process, class-specific response, or exposure effect. Exchangeability does not separate biological AP variation from technical section variation. Conditional binomial assumptions can miss residual cellular dependence.""",
        """Y cuenta células positivas y n incluye ceros en un grupo corte×anatomía. Mu es la media de propensiones de corte intercambiables, no la fracción agrupada ponderada por células. Kappa controla heterogeneidad; la correlación intraclase beta-binomial condicional es1/(kappa+1). Los priors son propios y se fijan antes de evaluar. No incluye tendencia AP, proceso espacial local, respuesta por clase ni exposición. La intercambiabilidad no separa variación AP biológica de variación técnica. La binomial condicional puede omitir dependencia residual entre células.""",
        ["metadata/spatial_detection_bayes/methods.json"],"formula-slide")

    add("beta-integration","Finite-grid integration and uncertainty computation","Integración por cuadrícula y cálculo de incertidumbre",
        """- Main grid: **64 μ × 48 κ prior quantiles**, equal prior-quantile weights.
- Normalize log beta-binomial likelihood weights with log-sum-exp.
- Hypermean μ intervals: weighted discrete-grid quantiles.
- Known-section p and held-out positive-count intervals: **8,192 seeded mixture draws**.
- Full-reference new-section p intervals: numerical inversion of the weighted Beta-mixture CDF.
- Check 32 × 24 and 96 × 64 grids without selecting a model from test results.""",
        """- Cuadrícula principal: **64 μ × 48 κ cuantiles previos**, con pesos iguales.
- Normalizar pesos de log-verosimilitud beta-binomial con log-sum-exp.
- Intervalos de μ: cuantiles ponderados de la cuadrícula discreta.
- Intervalos de p conocido y conteo reservado: **8,192 sorteos de mezcla con semilla**.
- Intervalos de p de un corte nuevo con referencia completa: inversión numérica de la CDF de mezcla Beta.
- Revisar 32 × 24 y 96 × 64 sin seleccionar por desempeño reservado.""",
        """In prior-quantile coordinates, prior masses are already represented by equally weighted nodes; adding another prior density/Jacobian would double-count the prior. The discrete posterior is an approximation to a continuous posterior, not an MCMC convergence result. Known-section mixtures update each beta component with that section's y/n using full-reference hyperposterior weights. Predictive count draws additionally sample binomial counts at the actual held-out denominator. The fixed seed is20261020 with reproducible keyed streams. Sensitivity evaluates numerical precision separately from scientific prediction performance.""",
        """En coordenadas de cuantiles previos, los nodos equiponderados ya representan el prior; multiplicar de nuevo por densidad/Jacobiano lo contaría dos veces. La posterior discreta aproxima una continua, sin afirmar convergencia MCMC. Las mezclas conocidas actualizan cada Beta con y/n del corte y pesos de hiperpósterior de referencia completa. La predicción de conteos sortea además binomiales con el denominador real reservado. La semilla20261020 usa flujos reproducibles por clave. La sensibilidad numérica se distingue del rendimiento científico.""",
        ["metadata/spatial_detection_bayes/methods.json","results/spatial_detection_bayes/scientific_diagnostics.json"])

    add("beta-targets","Three distinct detection uncertainty targets","Tres objetivos distintos de incertidumbre de detección",
        r"""| Target | Conditioning | Interval interpretation |
|---|---|---|
| Known-section $p_{sa}$ | Its observed positives and denominator | Credible interval of that fitted section |
| Anatomy mean $\mu_a$ | Exchangeable section information | Credible interval of the hypermean |
| New-section $p_*$ or $Y_*$ | No counts from that new section | Predictive distribution of its propensity or count |

An unseen held-out anatomy stays **prior-only** and remains in evaluation.

A fitted known-section interval is not a held-out prediction.""",
        r"""| Objetivo | Condicionamiento | Interpretación |
|---|---|---|
| $p_{sa}$ de corte conocido | Positivas y denominador observados | Intervalo creíble del corte ajustado |
| Media anatómica $\mu_a$ | Información de cortes intercambiables | Intervalo creíble de la hipermedia |
| $p_*$ o $Y_*$ de corte nuevo | Sin conteos de ese corte nuevo | Distribución predictiva de propensión o conteo |

Una anatomía reservada ausente del entrenamiento queda **solo con prior** y se conserva en evaluación.

Un intervalo ajustado conocido no es una predicción reservada.""",
        """The three intervals estimate different random quantities. Do not assert that a new-section interval is always numerically wider than a known-section interval, since location and parameter mixtures can differ. Held-out positive-count coverage refers to Y at its observed n, not to coverage of an unobserved true cell probability or a latent Gaussian intensity. In the interface, selecting a cell class suppresses this all-class fit because it has not estimated class-specific probabilities. Prior-only mu has mean.5 under the fixed symmetric prior, which can perform poorly for rare genes; that failure is scored rather than hidden.""",
        """Los intervalos estiman cantidades aleatorias diferentes. No se afirma que el intervalo de corte nuevo siempre sea más ancho, pues posición y mezcla paramétrica pueden diferir. La cobertura reservada de conteos corresponde a Y dado n observado, no a una probabilidad verdadera no observada ni intensidad gaussiana. En la interfaz, filtrar una clase oculta este ajuste de todas las clases porque no estima probabilidades específicas. Mu solo con prior tiene media.5 bajo el prior simétrico fijo, lo que puede predecir mal genes raros; se conserva ese fallo en la evaluación.""",
        ["metadata/spatial_detection_bayes/methods.json"],"formula-slide")

    add("metrics","Prediction metrics and baseline comparisons","Métricas predictivas y comparación con referencias",
        r"""$$\operatorname{RMSE}=
\sqrt{\frac1N\sum_i(y_i-\widehat y_i)^2}$$

$$\operatorname{Brier}=
\frac{\sum_j[Y_j(1-\widehat p_j)^2+
(n_j-Y_j)\widehat p_j^2]}{\sum_j n_j}$$

Gaussian baseline: compare anatomy-only regression and 30-neighbor prediction.

Detection baseline: collapse anatomy within each permitted training section.

Report held-out interval coverage **together with width**. Lower prediction loss is better.""",
        r"""$$\operatorname{RMSE}=
\sqrt{\frac1N\sum_i(y_i-\widehat y_i)^2}$$

$$\operatorname{Brier}=
\frac{\sum_j[Y_j(1-\widehat p_j)^2+
(n_j-Y_j)\widehat p_j^2]}{\sum_j n_j}$$

Referencia gaussiana: comparar regresión de solo anatomía y predicción con30 vecinos.

Referencia de detección: colapsar anatomía en cada corte permitido de entrenamiento.

Informar cobertura reservada **junto con amplitud**. Menor pérdida predictiva es mejor.""",
        """Brier is descriptive cell-weighted loss reconstructed from section×anatomy counts, not a significance test with independent cells. When pooling ten genes, group and cell evaluation counts multiply by ten; actual cells and sections do not. kNN uses30 nearest training cells within anatomy when available, with explicit broader fallback otherwise. Both baselines share the same held-out and guard exclusions. Gaussian sensitivity changes fixed basis lengths by.5/2 and spatial prior variances by.25/4. Overall95% prediction coverage is marginal: conditioning on the realized outcome being positive creates a different diagnostic without a nominal95% target.""",
        """Brier es pérdida descriptiva ponderada por células, reconstruida con grupos corte×anatomía, no una prueba con células independientes. Agrupar diez genes multiplica evaluaciones de grupos/células, no células o cortes reales. kNN usa30 vecinos de entrenamiento de la misma anatomía cuando existen y una alternativa más amplia explícita. Las referencias respetan prueba y separación idénticas. La sensibilidad gaussiana cambia longitudes por.5/2 y varianzas previas espaciales por.25/4. La cobertura95% es marginal: condicionar por resultado positivo crea un diagnóstico distinto, sin objetivo nominal95%.""",
        ["metadata/spatial_bayes/methods.json","metadata/spatial_detection_bayes/methods.json"],"formula-slide")

    section("replicated-results","Results reconstructed from the original studies","Resultados reconstruidos de los estudios originales",
        "Exact output recovery, partial methodological sensitivity and conditional concordance",
        "Recuperación exacta de salidas, sensibilidad metodológica parcial y concordancia condicional",
        "This results block uses the authors' original source animals. Each slide names the level of reproduction rather than grouping all forms of agreement under biological replication.",
        "Este bloque utiliza los animales fuente originales. Cada diapositiva nombra el nivel de reproducción, sin agrupar toda concordancia como replicación biológica.")

    fig("cortex-reproduction","Exact cortex significant-gene sets are recovered","Se recuperan exactamente los conjuntos corticales significativos","rep_cortex_reproduction",
        "Four author-aligned contrasts recover 109 / 1,205 / 1,222 / 639 genes. Every FDR-set Jaccard equals 1.",
        "Cuatro contrastes alineados recuperan 109 / 1,205 / 1,222 / 639 genes. Todos los Jaccard de conjuntos FDR son 1.",
        "The contrast order is PFF, PM, PM-PFF and PM-PFF plus nasal PM, each versus PBS. There are zero discordant significant genes in every comparison. The maximum MLE effect discrepancy across contrasts is about 8.12×10⁻⁵, and the maximum nominal-p discrepancy is 1.94×10⁻⁶. All coefficient fits converge. Code-defined PCA recovers 29.2014% and 9.7606%, matching the rounded values 29.2% and 9.8%. Also recovered are 582 PM-exclusive genes and 323 genes exclusive to the PM/PM-PFF pair. Original GO table row counts were checked, without independently recomputing the original GSEA. This reconstructs count-based inference and selected outputs, rather than raw-read alignment or every experiment in the paper.",
        "El orden de los contrastes es PFF, PM, PM-PFF y PM-PFF más PM nasal, cada uno frente a PBS. No hay genes significativos discordantes. La discrepancia MLE máxima es aproximadamente 8.12×10⁻⁵, y la de p nominal es 1.94×10⁻⁶. Todos los ajustes convergen. El PCA definido por el código recupera 29.2014% y 9.7606%, que coinciden con 29.2% y 9.8% al redondear. Se recuperan 582 genes exclusivos de PM y 323 exclusivos del par PM/PM-PFF. Se cotejaron las cantidades de filas GO originales sin repetir su GSEA. Se reconstruyen la inferencia de conteos y salidas seleccionadas, sin afirmar reproducción del alineamiento original ni de todos los experimentos.",
        ["ACA","results/replication_pm/cortex_gene_level_metrics.csv","results/replication_pm/cortex_refit_summary.json"])

    fig("cortex-source-quality","Reproducible output retains a source-quality question","La salida reproducible conserva una pregunta de calidad","rep_cortex_y_markers",
        "All four PBS controls have zero raw Ddx3y, Uty and Kdm5d counts. Other treatment groups are heterogeneous.",
        "Los cuatro controles PBS tienen cero conteos Ddx3y, Uty y Kdm5d. Los otros grupos son heterogéneos.",
        "The plot preserves all 20 animals and their source treatment labels. These Y-linked markers raise a possible sex-composition or sample-identity question, while the source RNA-seq sex metadata remain unreported. They do not establish sex assignments or complete confounding by themselves. No animals were relabeled or removed, and no invented sex covariate entered the reproduction. Exact numerical agreement can coexist with an unresolved source-design concern. Independent metadata verification is the appropriate next step, rather than post hoc omission of inconvenient samples.",
        "Se preservan los 20 animales y sus tratamientos originales. Los marcadores ligados a Y plantean una pregunta sobre composición por sexo o identidad de las muestras, mientras el sexo de los animales RNA-seq sigue sin informarse. Por sí solos no establecen el sexo individual ni una confusión completa. No se reasignaron ni excluyeron animales, y no se inventó una covariable de sexo. La reproducción numérica exacta puede coexistir con una preocupación de diseño. El siguiente paso es verificar los metadatos de manera independiente, sin omitir muestras inconvenientes a posteriori.",
        ["results/replication_pm/cortex_marker_expression_by_animal.csv","results/replication_pm/paper_claim_audit.csv"])

    fig("hy-reproduction","Hypothalamus comparison remains method-sensitive","La comparación hipotalámica conserva sensibilidad metodológica","rep_hypothalamus_sensitivity",
        "Paper HY lists imply 17 PM and 761 HFD genes. The declared non-RUVg sensitivity yields 0 and 3,600 under its stated gates.",
        "Las listas HY implican 17 genes PM y 761 HFD. La sensibilidad declarada sin RUVg produce 0 y 3,600 bajo sus criterios.",
        "The different numbers cannot be attributed uniquely to the unavailable RUVg adjustment. The empirical controls, number of factors, abundance-gate details and complete executable analysis are underdetermined. The project's independent primary analysis uses only 7 chow-fed FA/PM animals and returns 11 BH discoveries, without the author sensitivity analysis's coding/effect gates. These are different estimands and testing families. The defensible claim is transparent methodological sensitivity with unresolved reconstruction, rather than proof of biological nonreplication.",
        "Las diferencias no pueden atribuirse de manera única al ajuste RUVg ausente. Se desconocen los controles empíricos, la cantidad de factores, los detalles del criterio de abundancia y el código completo. El análisis principal independiente del proyecto usa siete animales FA/PM con dieta estándar y obtiene 11 hallazgos BH, sin los criterios de genes codificantes y efecto de la sensibilidad metodológica. Son estimandos y familias de pruebas distintos. La afirmación defendible es sensibilidad metodológica transparente con reconstrucción incompleta, sin demostrar una ausencia de replicación biológica.",
        ["HY","results/replication_pm/hypothalamus_refit_summary.json"])

    fig("mirna-reproduction","miR-574-5p direction agrees but BH confirmation does not","La dirección de miR-574-5p coincide sin confirmación BH","rep_mirna574",
        "High-dose difference −0.26945 signal units. limma p = 0.0253, BH q = 0.312. ANOVA p = 0.0788 prevents protected LSD.",
        "Diferencia a dosis alta −0.26945 unidades. limma p = 0.0253, BH q = 0.312. ANOVA p = 0.0788 impide LSD protegido.",
        "Unprotected high-dose LSD gives p = .0368, and nominal intervals exclude zero, while array-wide BH does not confirm this target. Protected LSD requires a significant omnibus test, which fails at .05. These statements are compatible because they use different gates and multiplicity rules. The target's lower signal at both doses reproduces direction, without reproducing independent qPCR or a mechanistic rescue result. The new descriptive heatmap has a disclosed selection rule and does not claim the unknown original row order or selection.",
        "LSD no protegido a dosis alta da p = .0368, e intervalos nominales que excluyen cero, pero BH de todo el array no confirma el objetivo. LSD protegido exige un ómnibus significativo, que no supera el criterio de .05. Estas afirmaciones son compatibles porque usan criterios y reglas de multiplicidad distintos. La menor señal a ambas dosis reproduce dirección, sin reproducir qPCR independiente ni un rescate mecanístico. El mapa de calor descriptivo nuevo tiene una selección declarada y no afirma recuperar el orden o la selección original desconocidos.",
        ["miRNA","results/replication_mirna_ozone/mirna_author_targets.csv"])

    fig("ozone-reproduction","Ozone genotype effects show conditional concordance","Los efectos de genotipo en ozono muestran concordancia condicional","rep_ozone_concordance",
        "Published genotype selection: 2,501 genes. Whole-effect comparison ρ = 0.865; selected S3 comparison ρ = 0.962 with 98.8% direction agreement.",
        "Selección publicada de genotipo: 2,501 genes. Comparación global de efectos ρ = 0.865; lista S3 seleccionada ρ = 0.962 y 98.8% de dirección concordante.",
        "The source table rule exactly recovers 1490 higher and 1011 lower genes. The new pooled genotype candidate has 15,559 common finite exact-symbol effects, 2362 nominal selections and a selected-set Jaccard of .612 in the common testable universe. The 17M/O₃ E4-versus-E3 S3 list contains 959 genes: 858 have finite matched effects, 848 retain direction, 543 retain nominal selection and 71 retain BH < .05. Selection-conditioned correlation cannot estimate unbiased transcriptome-wide agreement. The common-universe UpSet uses 15,661 genes, a different denominator from the 15,559-effect comparison. Hemoglobin concordance cannot establish blood–brain barrier function from bulk RNA alone.",
        "La regla de la tabla fuente recupera exactamente 1490 genes mayores y 1011 menores. El candidato conjunto de genotipo tiene 15,559 efectos finitos con símbolo exacto, 2362 selecciones nominales y Jaccard de .612 en el universo común evaluable. La lista S3 de E4 frente a E3 bajo O₃ a 17M contiene 959 genes: 858 tienen efectos finitos cotejados, 848 conservan dirección, 543 conservan selección nominal y 71 conservan BH < .05. La correlación condicionada por selección no estima concordancia transcriptómica global no sesgada. El UpSet común usa 15,661 genes, un denominador distinto. La concordancia de hemoglobinas no demuestra función de barrera hematoencefálica con ARN bulk.",
        ["O3","results/replication_mirna_ozone/summary.json"])

    add("reproduction-summary","Reproduction conclusions differ across studies","Las conclusiones de reproducción difieren entre estudios",
        table(["Study","Defensible result","Remaining limit"],[
            ["Cortex","Exact four significant-gene sets and code PCA","Read-level pipeline and source sex/identity unresolved"],
            ["Hypothalamus","Declared TMM/voom sensitivity","RUVg controls and k unavailable"],
            ["miRNA","Target direction and unprotected nominal agreement","Omnibus gate and array-wide BH do not confirm target"],
            ["Ozone","Exact deposited threshold selection and matched effects","Different quantifier and selected-list conditioning"]]),
        table(["Estudio","Resultado defendible","Límite pendiente"],[
            ["Corteza","Cuatro conjuntos y PCA del código exactos","Sin flujo de lecturas completo ni sexo/identidad resueltos"],
            ["Hipotálamo","Sensibilidad TMM/voom declarada","Sin controles RUVg ni k"],
            ["miRNA","Dirección y evidencia nominal no protegida","Ómnibus y BH del array no confirman el objetivo"],
            ["Ozono","Selección exacta de tabla y efectos cotejados","Otro cuantificador y condicionamiento por selección"]]),
        """This is the conclusion of the reproduction block, not a ranking of biological plausibility. All comparisons reuse source animals. Exact output agreement is strongest for cortex and still requires source-quality follow-up. Partial reproduction should not be silently upgraded into full mechanism validation. The next block presents new project analyses, including models whose validation reveals limitations rather than successful reconstruction.""",
        """Se cierra el bloque de reproducción, sin ordenar plausibilidad biológica. Todos reutilizan animales fuente. Corteza ofrece la concordancia de salidas más fuerte y aún requiere aclarar calidad. Una reproducción parcial no se eleva a validación mecanística. El bloque siguiente contiene análisis nuevos del proyecto, incluidos modelos cuya evaluación identifica limitaciones en vez de reconstrucción exitosa.""",
        ["results/replication_pm/paper_claim_audit.csv","results/replication_mirna_ozone/paper_replication_audit.csv"])

    section("new-results","New analyses within this project","Análisis nuevos dentro del proyecto",
        "Primary reanalysis, response heterogeneity, gene sets and spatial-reference modeling",
        "Reanálisis principal, heterogeneidad de respuesta, conjuntos génicos y modelado de referencia espacial",
        "These results add analyses to the repository, not new animals. “New” denotes project contribution and does not assert publication novelty or a confirmed biological mechanism.",
        "Estos resultados incorporan análisis al repositorio, no animales nuevos. «Nuevo» denota aporte del proyecto, sin afirmar novedad publicable o mecanismo biológico confirmado.")

    fig("hy-primary","Primary HY analysis identifies 11 BH discoveries","El análisis HY principal identifica 11 hallazgos BH","rep_hypothalamus_primary",
        "4 FA versus 3 PM animals. Two higher and nine lower genes at BH q < 0.05. Intervals are pointwise MLE Wald intervals.",
        "4 animales FA frente a 3 PM. Dos genes mayores y nueve menores con BH q < 0.05. Intervalos Wald MLE puntuales.",
        "The primary analysis retains 19,611 genes and has 19,452 finite Wald p values after 159 Cook-filtered tests. Tll1 has MLE −1.343708, q = .003497 and shrunken effect −1.224053. All coding primary hits decrease under the frozen current annotation, while the two increases are noncoding/pseudogene features. Leave-one-out produces 0–825 BH discoveries and preserves the direction of all 11 primary hits, but not significance in every fit. Tll1 retains q < .05 in 5/7 omissions. This small-n regional association does not identify a cell-specific or protein result.",
        "Se retienen 19,611 genes y se obtienen 19,452 valores p Wald finitos tras 159 pruebas filtradas por Cook. Tll1 tiene MLE = −1.343708, q = .003497 y efecto regularizado −1.224053. Todos los hallazgos codificantes disminuyen bajo la anotación actual congelada, mientras los dos aumentos corresponden a una característica no codificante y un pseudogén. Leave-one-out produce 0–825 hallazgos BH y conserva la dirección de los 11, sin conservar significación en todos los ajustes. Tll1 mantiene q < .05 en 5/7 omisiones. La asociación regional con n pequeño no identifica un resultado celular específico ni proteico.",
        ["results/hypothalamus/summary.json","results/hypothalamus/leave_one_out_primary_hits.tsv"])

    fig("aca-primary","Primary cortex analysis uses eight PM/PBS animals","El análisis cortical principal usa ocho animales PM/PBS","rep_cortex_primary",
        "1,167 BH discoveries: 658 higher and 509 lower. This population and model differ from the 20-animal, 1,205-gene author reconstruction.",
        "1,167 hallazgos BH: 658 mayores y 509 menores. La población y el modelo difieren de la reconstrucción de autores con 20 animales y 1,205 genes.",
        "The primary matrix retains 19,235 genes and has 17,360 finite adjusted p values, with 10 Cook-filtered tests and additional adaptive independent filtering. Forty-five genes satisfy the exploratory joint gate BH < .05 and absolute shrunken LFC ≥ 1. The volcano displays the primary evidence, without choosing another threshold or claiming formal reproduction of the author's volcano. Fosb has MLE = 1.01816 and q = 2.265×10⁻¹⁹. Hspa5 is the top primary q-ranked example but is absent from the measured atlas panel. The unresolved study-level sex/identity signal limits attribution without implying that every selected gene is sex-linked.",
        "La matriz principal retiene 19,235 genes y tiene 17,360 valores p ajustados finitos, con 10 pruebas filtradas por Cook y filtrado independiente adicional. Cuarenta y cinco genes cumplen el criterio exploratorio conjunto BH < .05 y |LFC regularizado| ≥ 1. El gráfico volcano representa la evidencia principal, sin elegir otro umbral ni afirmar una reproducción formal del volcano original. Fosb tiene MLE = 1.01816 y q = 2.265×10⁻¹⁹. Hspa5 es el ejemplo mejor ordenado por q principal, pero no se mide en el panel. La señal de sexo/identidad del estudio limita la atribución, sin implicar que todo gen seleccionado esté ligado al sexo.",
        ["results/cortex/summary.json","results/cortex/DE_PM_vs_PBS_MLE.tsv"])

    add("mirna-primary","New miRNA analysis retains the full test family","El análisis miRNA nuevo conserva toda la familia de pruebas",
        table(["Contrast","Tested features","BH discoveries","Effect scale"],[
            ["PM 1 mg/kg versus control","1,245","1 higher","Original signal difference"],
            ["PM 5 mg/kg versus control","1,245","0","Original signal difference"],
            ["Separate omnibus F","1,245","1","Joint condition difference"]])+
        "\nThe low-dose discovery is mmu-miR-5119. Lower detection flags remain visible, and no array is removed.",
        table(["Contraste","Características","Hallazgos BH","Escala"],[
            ["PM 1 mg/kg frente a control","1,245","1 mayor","Diferencia original de señal"],
            ["PM 5 mg/kg frente a control","1,245","0","Diferencia original de señal"],
            ["F ómnibus separada","1,245","1","Diferencia conjunta de condición"]])+
        "\nEl hallazgo de dosis baja es mmu-miR-5119. Se muestran indicadores de baja detección y no se excluye ningún array.",
        "There are 823 low-detection features and 7722 negative signal values retained. The lower-bound pile-up near −3.3 and small n of 3/group caution against treating every continuous value as a precise abundance measurement. The absence of clear PCA separation does not justify reselecting animals. No high-dose BH discoveries does not establish equivalence or a monotonic dose-response reversal. Pairwise and omnibus families are corrected separately and have different null hypotheses. No miRNA-target network or mRNA GO analysis is inferred from names alone.",
        "Se conservan 823 características con baja detección y 7722 valores negativos. La acumulación cerca de −3.3 y n de 3/grupo limitan su interpretación como abundancia precisa. La ausencia de separación clara en PCA no justifica seleccionar animales nuevamente. Cero hallazgos BH a dosis alta no establece equivalencia ni una inversión monotónica de dosis-respuesta. Las familias por contraste y ómnibus se corrigen por separado y plantean hipótesis nulas distintas. No se infiere una red de blancos ni GO de mRNA a partir de nombres de miRNA.",
        ["results/mirna/summary.json","results/mirna/methods.json"])

    fig("ozone-primary","Ozone discovery counts depend on stratum and test","Los hallazgos de ozono dependen del estrato y la prueba","rep_ozone_primary",
        "Within-stratum BH discoveries: 0 / 399 / 2 / 3. The separate 4-df global exposure LRT identifies 24 genes.",
        "Hallazgos BH por estrato: 0 / 399 / 2 / 3. El LRT global separado de exposición, con 4 gl, identifica 24 genes.",
        "The order is 3M/E3, 3M/E4, 17M/E3 and 17M/E4, with 3 FA and 3 O₃ animals in each stratum. Young E4 has 86 higher and 313 lower genes, old E3 has two higher genes, and old E4 has three lower genes. The global LRT has 18,511 valid tests. Supplementary BH across all four Wald families gives 1/43/4/6 and defines a different correction family, so its selected lists need not nest inside the per-contrast BH lists. We do not choose the largest discovery count as the preferred correction. Unequal stratum discoveries motivate direct heterogeneity tests, rather than proving age/genotype effect modification by themselves.",
        "El orden es 3M/E3, 3M/E4, 17M/E3 y 17M/E4, con 3 FA y 3 O₃ por estrato. E4 joven tiene 86 genes mayores y 313 menores, E3 mayor tiene dos mayores y E4 mayor tiene tres menores. El LRT global tiene 18,511 pruebas válidas. La BH suplementaria de las cuatro familias Wald da 1/43/4/6 y define otra familia de corrección, por lo que sus listas no tienen que anidarse en las de BH por contraste. No se escoge el mayor número como corrección preferida. Los recuentos distintos motivan pruebas directas de heterogeneidad, sin demostrar por sí solos modificación por edad/genotipo.",
        ["results/ozone/summary.json"])

    add("interaction-results","Age-dependent exposure has the clearest joint signal","La respuesta dependiente de edad tiene la señal conjunta más clara",
        table(["Secondary joint LRT","df","Within-family BH hits","Across-three-family BH hits"],[
            ["Genotype-dependent exposure","2","0","0"],
            ["Age-dependent exposure","2","8","1"],
            ["Three-way interaction","1","0","0"]])+
        "\nAge-response genes: Ryr3, Rps24, Pkig, Rpl26, Inpp5d, Stk11, Rpl34 and Gm21586.",
        table(["LRT conjunto secundario","gl","Hallazgos BH por familia","BH conjunta tres familias"],[
            ["Respuesta dependiente de genotipo","2","0","0"],
            ["Respuesta dependiente de edad","2","8","1"],
            ["Interacción triple","1","0","0"]])+
        "\nGenes de respuesta por edad: Ryr3, Rps24, Pkig, Rpl26, Inpp5d, Stk11, Rpl34 y Gm21586.",
        "All secondary LRTs preserve the original 18537-gene universe, fitted dispersions and length offsets. The age LRT tests AE and AGE jointly. Its 8 within-family hits reduce to 1 under the supplementary BH family across three LRTs. Signed genotype-response differences give 6 hits at 3M and 4 at 17M. Age-response differences give 2 in E3 and 4 in E4, while the triple signed contrast gives 0. These tests use different contrasts, validity masks and multiplicity families, so their discoveries need not agree with the joint LRTs. Zero discoveries in another family do not prove identical effects.",
        "Los LRT secundarios conservan el universo original de 18537 genes, las dispersiones ajustadas y los offsets. El LRT de edad prueba AE y AGE conjuntamente. Sus 8 hallazgos por familia pasan a 1 con la BH suplementaria entre tres LRT. Las diferencias de respuesta por genotipo dan 6 a 3M y 4 a 17M. Las diferencias por edad dan 2 en E3 y 4 en E4, mientras el contraste triple da 0. Son contrastes, máscaras de validez y familias de multiplicidad diferentes, por lo que sus hallazgos no tienen que coincidir con los LRT conjuntos. Cero hallazgos en otra familia no demuestra igualdad de efectos.",
        ["results/advanced_omics/summary.json","results/advanced_omics/age_response_FDR_genes.tsv"])

    fig("interaction-animals","Animal-level profiles contextualize response differences","Los perfiles animales contextualizan diferencias de respuesta","context_ozone_animals",
        "All 24 animals remain visible for each gene. Points encode animals and short horizontal segments show descriptive means.",
        "Se muestran los 24 animales por gen. Los puntos representan animales y los segmentos cortos muestran medias descriptivas.",
        "Ryr3 and Rps24 come from the fixed age-response LRT selection, rather than a new search for visually impressive genes. Values are log2(1 + existing normalized count), with fixed tiny horizontal offsets and no pairing lines or inferential error bars. Model-based interaction p values come from the negative-binomial fit with covariance, rather than the visual separation of these points. Transformed means in the plot differ conceptually from count-scale GLM means. Three animals per stratum limit generalization.",
        "Ryr3 y Rps24 proceden de la selección fija del LRT de edad, sin buscar nuevos gráficos llamativos. Los valores son log₂(1 + conteo normalizado existente), con pequeños desplazamientos horizontales fijos y sin líneas de pareamiento ni barras inferenciales. Los valores p de interacción derivan del ajuste binomial negativo con covarianza, sin inferirse de la separación visual. Las medias transformadas difieren de las medias GLM en escala de conteos. Tres animales por estrato limitan la generalización.",
        ["results/advanced_omics/ozone_interaction_animal_profiles.tsv"])

    add("enrichment-results","ORA and ranked GSEA describe different evidence","ORA y GSEA por rango describen evidencia diferente",
        table(["Contrast","ORA up / down terms","GSEA algorithmic-FDR terms"],[
            [str(row.dataset).replace("hypothalamus","HY").replace("cortex","ACA").replace("ozone_","O₃ "),
             f"{row.ORA_up_significant_terms} / {row.ORA_down_significant_terms}",
             row.GSEA_terms_algorithmic_FDR_0_05] for row in enrich.itertuples()])+
        "\nThese correlated terms are not counts of independent biological mechanisms.",
        table(["Contraste","Términos ORA mayores / menores","Términos FDR algorítmica GSEA"],[
            [str(row.dataset).replace("hypothalamus","HY").replace("cortex","ACA").replace("ozone_","O₃ "),
             f"{row.ORA_up_significant_terms} / {row.ORA_down_significant_terms}",
             row.GSEA_terms_algorithmic_FDR_0_05] for row in enrich.itertuples()])+
        "\nEstos términos correlacionados no cuentan mecanismos biológicos independientes.",
        "HY has no ORA terms but 916 ranked GSEA terms because a thresholded foreground and coordinated ranks ask different questions. The large count does not imply 916 independent pathways or resolve the seven-animal uncertainty. ACA has 86 up and 2 down ORA terms. All ozone ORA families are empty at their thresholds, while GSEA finds 4, 368, 85 and 23 terms. Mapping, finite-test universes, propagated ancestry and algorithmic GSEA FDR rather than BH all matter. Enrichment is exploratory and inherits source-design limitations.",
        "HY no tiene términos ORA, pero sí 916 términos GSEA, porque el conjunto umbralizado y los rangos coordinados responden preguntas distintas. No implica 916 vías independientes ni resuelve la incertidumbre de siete animales. ACA presenta 86 términos ORA positivos y 2 negativos. Las familias ORA de ozono están vacías con sus umbrales, mientras GSEA obtiene 4, 368, 85 y 23 términos. Importan el mapeo, los universos finitos, los ancestros propagados y la FDR algorítmica de GSEA, distinta de BH. El enriquecimiento es exploratorio y hereda las limitaciones del diseño.",
        ["results/enrichment/summary_counts.tsv","results/enrichment/methods.json"])

    fig("go-patterns","Selected GO patterns across the six RNA contrasts","Patrones GO seleccionados en seis contrastes de ARN","context_go_patterns",
        "Fixed display selection, shared NES scale. Dot: algorithmic GSEA FDR < 0.05. Gray: ineligible or unavailable.",
        "Selección gráfica fija y escala NES común. Punto: FDR algorítmica de GSEA < 0.05. Gris: no elegible o no disponible.",
        "The presentation redraws the first 8 terms from the existing 22-term project display selection for projected readability, without rerunning enrichment. The full 22-term figure and 132 source rows remain in the repository. The selected terms are not an exhaustive or unbiased summary of all GO. A positive NES means enrichment toward the positive contrast rank, rather than demonstrated pathway activation. The separate 20-node GO overlap network has 15 edges at Jaccard ≥ .25 and depicts shared annotations, without establishing regulation, coexpression or spatial connectivity.",
        "Se redibujan los primeros 8 términos de la selección gráfica existente de 22 para facilitar la lectura proyectada, sin repetir el enriquecimiento. La figura completa y las 132 filas fuente permanecen en el repositorio. Los términos seleccionados no constituyen un resumen exhaustivo ni no sesgado de todo GO. NES positivo significa enriquecimiento hacia el rango positivo, sin demostrar activación de vías. La red GO separada de 20 nodos tiene 15 aristas con Jaccard ≥ .25 y representa anotaciones compartidas, sin establecer regulación, coexpresión ni conectividad espacial.",
        ["results/advanced_omics/GO_multi_contrast_NES_source.tsv","metadata/advanced_omics/methods.json"])

    selected=["Tll1","Fosb","S1pr1","Hspa5","Ksr2","C730002L08Rik"]
    ann=genes.drop_duplicates("gene").set_index("gene")
    product_sources=T("results/gene_priorities/annotations.tsv").drop_duplicates("gene").set_index("gene")
    def gene_link(g):
        return f'[{g}]({product_sources.loc[g,"source_url"]})'
    add("gene-products","Gene products add biological context without measuring protein","Los productos génicos contextualizan sin medir proteínas",
        table(["Gene","Annotated product / role","Measured panel"],[
            [gene_link(g),{"Tll1":"Tolloid-like extracellular processing protease",
                "Fosb":"FosB, AP-1 transcription-factor subunit",
                "S1pr1":"Sphingosine-1-phosphate receptor",
                "Hspa5":"BiP/GRP78 endoplasmic-reticulum chaperone",
                "Ksr2":"RAF–MEK/MAPK scaffold",
                "C730002L08Rik":"Noncoding RNA, no annotated protein"}[g],
             "Yes" if bool(ann.loc[g,"atlas_measured"]) else "No"] for g in selected]),
        table(["Gen","Producto / función anotada","Panel medido"],[
            [gene_link(g),{"Tll1":"Proteasa extracelular de procesamiento tolloid-like",
                "Fosb":"FosB, subunidad del factor transcripcional AP-1",
                "S1pr1":"Receptor de esfingosina-1-fosfato",
                "Hspa5":"Chaperona BiP/GRP78 del retículo endoplásmico",
                "Ksr2":"Andamio RAF–MEK/MAPK",
                "C730002L08Rik":"ARN no codificante, sin proteína anotada"}[g],
             "Sí" if bool(ann.loc[g,"atlas_measured"]) else "No"] for g in selected]),
        "Product definitions summarize the frozen NCBI/UniProt evidence, with qualifiers preserved in the full annotation table. RNA expression does not establish enzymatic activity, protein abundance or an exposure mechanism. Fosb gene-level RNA cannot resolve FosB versus ΔFosB isoforms. Ksr2 is interpreted as a scaffold: reported low in-vitro kinase activity does not establish active catalysis here. Noncoding C730002L08Rik has no established protein product. Hspa5 and Ksr2 remain useful bulk candidates even though no measured spatial map is available. No absent gene is replaced with a similarly named atlas gene.",
        "Las definiciones resumen registros NCBI/UniProt congelados y conservan sus calificadores en la tabla completa. La expresión de ARN no establece actividad enzimática, abundancia proteica ni un mecanismo de exposición. El ARN de Fosb no distingue las isoformas FosB/ΔFosB. Ksr2 se interpreta como andamio: la baja actividad cinasa in vitro informada no demuestra catálisis activa aquí. C730002L08Rik no tiene producto proteico establecido. Hspa5 y Ksr2 siguen siendo candidatos bulk útiles aunque carezcan de un mapa medido. No se reemplaza un gen ausente por otro de nombre parecido.",
        ["results/gene_priorities/annotations.tsv","metadata/gene_priorities/methods.json"])

    fig("aca-sections","Measured cortical layers organize the reference cells","Las capas corticales medidas organizan las células de referencia","spatial_aca_measured",
        "Largest ACA section selected by cell count: 1,940 cells. Anatomy colors and measured Fosb detection use the deposited experimental plane.",
        "Corte ACA de mayor n: 1,940 células. Los colores anatómicos y la detección de Fosb usan el plano experimental depositado.",
        "The section is Zhuang-ABCA-2.023, selected by maximal regional n independently of expression. All measured cells, including zeros, contribute to bin denominators. Layer structure follows the original spatial positions and CCF annotations, without random placement or a simulated point cloud. Empty bins remain unobserved. The complete 18-section ACA reference includes 14582 cells and 10 layers/subdivisions, with the one-cell edge section preserved. These observations describe one unexposed male donor, without identifying layer-specific PM effects.",
        "El corte es Zhuang-ABCA-2.023, seleccionado por mayor n regional sin usar la expresión. Todas las células medidas, incluidos los ceros, aportan a los denominadores. La organización laminar procede de las posiciones originales y las anotaciones CCF, sin colocación aleatoria ni una nube simulada. Las casillas vacías quedan sin observación. La referencia ACA completa contiene 14582 células, 18 cortes y 10 subdivisiones, incluido el corte de una célula. Describe un macho sin exposición, sin identificar efectos PM por capa.",
        ["results/spatial_sections/ACA/sections.tsv","results/spatial_sections/ACA/detection_by_bins.tsv.gz"])

    fig("hy-sections","Measured hypothalamic nuclei and Tll1 detection","Núcleos hipotalámicos medidos y detección Tll1","spatial_hy_measured",
        "Largest HY section: 3,565 cells. Tll1 is positive in 667 of 37,811 cells across the complete reference HY region.",
        "Corte HY de mayor n: 3,565 células. Tll1 es positivo en 667 de las 37,811 células de la región HY completa.",
        "The example physical section is .034, selected without using Tll1 expression. The complete HY reference has 20 sections and 49 anatomical labels, including unassigned labels. Overall Tll1 detection is 1.764%; a section-specific fraction can differ and retains its own denominator. A low raw count does not distinguish biological absence from limited assay detection. No new intermediate sections or bilateral copies are generated. Descriptive detection by nucleus is not an exposure association test.",
        "El corte físico de ejemplo es .034, seleccionado sin usar Tll1. La referencia HY tiene 20 cortes y 49 etiquetas anatómicas, incluidas etiquetas no asignadas. La detección global de Tll1 es 1.764%; una fracción por corte puede diferir y conserva su propio denominador. Un conteo bajo no distingue ausencia biológica de detectabilidad limitada. No se generan cortes intermedios ni copias bilaterales. La detección descriptiva por núcleo no constituye una prueba de asociación con la exposición.",
        ["results/spatial_sections/HY/sections.tsv","results/spatial_reference/summary.json"])

    add("scenario-results","Illustrative regional scaling preserves the measured geometry","El escalamiento ilustrativo preserva la geometría medida",
        """**ACA / Fosb:** the shrunken bulk effect gives an approximately twofold scenario.

**HY / Tll1:** the shrunken bulk effect gives a 0.428-fold scenario.

- Every positive reference value receives the same regional multiplier for its gene.
- Measured zeros stay zero, and positions and cell counts stay fixed.
- The interface links measured candidates, gene products, source animal results and uncertainty targets.

[Open the paired 3D cortical viewer](https://amarusimonaguerojimenez.github.io/Air-Pollution-Brain-Omics-Mouse/cortex_viewer.html).""",
        "**ACA / Fosb:** el efecto bulk regularizado produce un escenario cercano al doble.\n\n**HY / Tll1:** el efecto bulk regularizado produce un escenario de 0.428 veces.\n\n- Cada valor positivo recibe el mismo multiplicador regional de su gen.\n- Los ceros, las posiciones y la cantidad de células se mantienen.\n- La interfaz vincula candidatos medidos, productos génicos, resultados de animales fuente e incertidumbre.\n\n[Abrir el visor cortical 3D emparejado](https://amarusimonaguerojimenez.github.io/Air-Pollution-Brain-Omics-Mouse/es/cortex_viewer.html).",
        "The interactive side-by-side display is a project deliverable, without independently validating a spatial effect. The static full-resolution reference/scenario figure remains linked in the results. A shared scale helps interpret uniform scaling but cannot create cell-specific response information. Relative spatial pattern and detection support remain unchanged under a positive multiplier. Bayesian baseline uncertainty is separate from the bulk Wald alternatives and must not be presented as an exposure interval.",
        "El visor lado a lado es un aporte del proyecto, sin validar de manera independiente un efecto espacial. La figura estática completa de referencia/escenario permanece en los resultados. La escala común ayuda a interpretar el escalamiento uniforme, sin crear información de respuesta celular específica. Un multiplicador positivo conserva el patrón relativo y el soporte de detección. La incertidumbre basal bayesiana es distinta de las alternativas Wald del bulk y no representa un intervalo de exposición.",
        ["results/spatial_cortex/summary.json","results/spatial_viewer/summary.json"])
    SLIDES[-1]["demo_figure"]="spatial_cortex/reference_vs_simulation"

    fig("gaussian-rmse","Spatial Gaussian regression rarely improves the baseline","La regresión gaussiana espacial mejora pocas veces la referencia","spatial_gaussian_rmse",
        "RMSE improves over anatomy-only prediction for 4/10 ACA genes and 1/10 HY genes under the same guarded section folds.",
        "El RMSE mejora respecto a solo anatomía en 4/10 genes ACA y 1/10 HY, con los mismos bloques de cortes.",
        "Tll1 RMSE is .109713 for the spatial model, .107582 for anatomy only and .107175 for kNN. Fosb RMSE is .297687, .293223 and .301802, respectively. All losses use natural ln(1 + raw imaging count). Ratios below 1 favor the spatial model, but these are descriptive comparisons within one donor, without a significance test assuming independent cells. Four prespecified sensitivity configurations inspect length and prior choices without test-selected tuning. Smooth maps alone do not establish predictive gain.",
        "El RMSE de Tll1 es .109713 para el modelo espacial, .107582 para solo anatomía y .107175 para kNN. Para Fosb es .297687, .293223 y .301802, respectivamente. Las pérdidas usan ln(1 + conteo de imagen). Cocientes menores que 1 favorecen el modelo espacial, pero son comparaciones descriptivas de un donante, sin una prueba que suponga células independientes. Cuatro configuraciones predefinidas revisan longitudes y priors, sin selección por los datos de prueba. Un mapa suave no demuestra ganancia predictiva.",
        ["results/spatial_bayes/summary.json"])

    fig("gaussian-diagnostics","Overall coverage hides important Gaussian limitations","La cobertura global oculta limitaciones gaussianas importantes","spatial_gaussian_diagnostics",
        "Overall coverage, positive-only inclusion and impossible negative predictive support are different diagnostics.",
        "La cobertura global, la inclusión solo en positivos y el soporte predictivo negativo imposible son diagnósticos distintos.",
        "Median marginal coverage is 97.27% in HY and 94.00% in ACA. Positive-tail inclusion medians are 0% and 66.46%; the 95% nominal marginal target does not apply after conditioning on a positive realized outcome. Nevertheless, median predictive mass below zero is 41.06% in HY and 38.19% in ACA, which is physically impossible for ln(1 + count). Tll1 overall coverage of 98.24% coexists with 0/667 positive inclusion and 41.24% negative predictive mass. These findings and limited baseline gains prevent a validated cellular RNA reconstruction claim. Negative means and bounds remain visible rather than being clipped.",
        "Las medianas de cobertura marginal son 97.27% en HY y 94.00% en ACA. Las medianas de inclusión en la cola positiva son 0% y 66.46%; el objetivo marginal nominal de 95% no se aplica al condicionar por el resultado positivo. Aun así, la masa predictiva bajo cero es 41.06% en HY y 38.19% en ACA, imposible para ln(1 + conteo). Tll1 tiene 98.24% de inclusión global, 0/667 en positivos y 41.24% de masa negativa. Junto con las ganancias basales limitadas, esto impide afirmar una reconstrucción celular del ARN validada. Se conservan las medias y los límites negativos sin recortarlos.",
        ["results/spatial_bayes/summary.json"])

    fig("beta-performance","Bounded detection improves modestly in ACA but fails in HY","La detección acotada mejora modestamente en ACA y falla en HY","spatial_beta_performance",
        "ACA Brier 0.13725 versus 0.14051 regional. HY Brier 0.09246 versus 0.03621. High count-interval coverage accompanies wide intervals.",
        "Brier ACA de 0.13725 frente a 0.14051 regional. HY de 0.09246 frente a 0.03621. La cobertura alta de conteos acompaña intervalos amplios.",
        "Pooled ten-gene count PI coverage is 1253/1280 = 97.89% for ACA and 1848/1860 = 99.35% for HY. Mean count PI width divided by n is 41.9 and 61.5 percentage points. These are gene×section×anatomy evaluations, without independent animals. In HY, 79/186 unique groups per gene and 10416/37811 actual cells are prior-only under the held-out folds. The fixed symmetric prior can predict detection poorly for rare genes. Trained-anatomy-only HY Brier of .032557 versus .032877 is a small restricted improvement and does not replace full-support failure. The model targets detection, rather than count magnitude.",
        "La cobertura agrupada de diez genes es 1253/1280 = 97.89% para ACA y 1848/1860 = 99.35% para HY. La amplitud media del intervalo dividida por n abarca 41.9 y 61.5 puntos porcentuales. Son evaluaciones gen×corte×anatomía, sin animales independientes. En HY, 79/186 grupos únicos por gen y 10416/37811 células reales quedan solo con prior en los bloques reservados. El prior simétrico fijo puede predecir mal la detección de genes raros. El Brier HY restringido a anatomías entrenadas, .032557 frente a .032877, mejora ligeramente, sin reemplazar el fallo con soporte completo. El objetivo es detección, no magnitud del conteo.",
        ["results/spatial_detection_bayes/summary.json"])

    fig("beta-aca","ACA posterior detection separates uncertainty targets","La detección posterior ACA separa objetivos de incertidumbre","spatial_beta_aca_posterior",
        "Fosb: observed fractions, fitted known-section probabilities and anatomy/new-section uncertainty remain separately labeled.",
        "Fosb: fracciones observadas, probabilidades del corte conocido e incertidumbre anatómica/de un corte nuevo, rotuladas por separado.",
        "At the selected known section, p conditioned on its own observed counts is a fitted summary, without being an out-of-fold prediction. Anatomy μ is a mean over exchangeable section propensities, whereas the measured pooled fraction weights actual cells. The new-section p distribution models section heterogeneity and differs from a count prediction with denominator n. Small layers can have wider or differently shaped uncertainty, without resolving exposure effects. ACA/Fosb held-out Brier of .08745 versus .08803 shows only a small gain, with count PI coverage of 99.2% and width fraction .345.",
        "La p del corte conocido se condiciona a sus propios conteos: es un resumen ajustado, sin ser una predicción reservada. μ promedia propensiones de cortes intercambiables, mientras la fracción agrupada pondera células reales. La distribución de p de un corte nuevo modela heterogeneidad y difiere de una predicción de conteos con denominador n. Las capas pequeñas pueden tener incertidumbre amplia o de forma distinta, sin resolver efectos de exposición. El Brier de Fosb, .08745 frente a .08803, indica una ganancia pequeña, con cobertura de 99.2% y ancho/n de .345.",
        ["results/spatial_detection_bayes/ACA/known_section_posterior.tsv.gz","results/spatial_detection_bayes/ACA/anatomy_posterior.tsv.gz"])

    fig("beta-hy","HY fitted anatomy does not guarantee held-out prediction","La anatomía HY ajustada no garantiza predicción reservada","spatial_beta_hy_posterior",
        "Tll1: fitted probabilities stay bounded, while unseen held-out anatomies expose dependence on the stated prior.",
        "Tll1: probabilidades ajustadas acotadas. Las anatomías reservadas no vistas revelan dependencia del prior.",
        "HY/Tll1 all-group held-out Brier of .08329 is worse than the regional value .01772. Count coverage is 98.9%, but mean width/n is .600 and 79 of 186 groups are prior-only. The displayed full-reference posterior uses all reference sections and can therefore look reasonable even when guarded prediction fails. Showing both prevents conflating fit and generalization. A class filter in the interface hides the all-class estimates instead of pretending to refit a class-specific model. The posterior forest adds neither spatial resolution nor PM-responsive cells.",
        "El Brier reservado de Tll1, .08329, es peor que el regional, .01772. La cobertura es 98.9%, pero la amplitud media/n es .600 y 79 de 186 grupos quedan solo con prior. La posterior de referencia completa utiliza todos los cortes y puede parecer razonable aunque falle la predicción con separación. Mostrar ambas evita confundir ajuste y generalización. Filtrar por clase oculta las estimaciones de todas las clases, sin simular un reajuste específico de clase. El gráfico de intervalos no agrega resolución ni células respondedoras a PM.",
        ["results/spatial_detection_bayes/HY/known_section_posterior.tsv.gz","results/spatial_detection_bayes/HY/gene_metrics.tsv"])

    add("numerics","Numerical stability and predictive usefulness are separate","Estabilidad numérica y utilidad predictiva son distintas",
        table(["64×48 versus 96×64 grid","ACA","HY"],[
            ["Maximum mean μ change","0.000192","0.000537"],
            ["Maximum μ interval endpoint change","0.0203","0.0181"],
            ["Maximum new-section p endpoint change","0.00254","0.00358"],
            ["Held-out count coverage","97.89% → 98.05%","99.35% unchanged"]])+
        "\nFinite-grid and seeded simulation checks support computational review. They do not repair weak prediction or validate the priors biologically.",
        table(["Cuadrícula 64×48 frente a 96×64","ACA","HY"],[
            ["Cambio máximo de media μ","0.000192","0.000537"],
            ["Cambio máximo de extremo IC de μ","0.0203","0.0181"],
            ["Cambio máximo de extremo p nuevo","0.00254","0.00358"],
            ["Cobertura reservada de conteos","97.89% → 98.05%","99.35% sin cambio"]])+
        "\nLas verificaciones numéricas y de simulación reproducible apoyan la revisión computacional. No reparan una predicción débil ni validan biológicamente los priors.",
        "OOF Brier differences are at most 10⁻⁶ in the fine-grid comparison. Coverage status changes for four ACA groups and no HY groups. Discrete μ interval quantiles can move more than means, so reporting mean stability alone would overstate uncertainty convergence. Deterministic new-p CDF inversion removes Monte Carlo noise from that particular diagnostic; known-section p and held-out counts still use 8192 seeded draws. An independent Gaussian posterior simulation check verifies analytic uncertainty transport without proving likelihood adequacy. Implementation checks do not establish blanket scientific success.",
        "Las diferencias OOF de Brier son como máximo 10⁻⁶. Cambia el estado de cobertura de cuatro grupos ACA y ninguno HY. Los cuantiles discretos de μ pueden moverse más que las medias, por lo que informar solo estabilidad media exageraría la convergencia de la incertidumbre. La inversión determinista de la CDF elimina el ruido Monte Carlo de ese diagnóstico particular; p del corte conocido y conteos reservados mantienen 8192 sorteos. Un chequeo independiente de la posterior gaussiana verifica el transporte analítico de incertidumbre, sin validar la verosimilitud. Las verificaciones de implementación no equivalen a éxito científico general.",
        ["results/spatial_detection_bayes/scientific_diagnostics.json","results/spatial_bayes/validation.json"])

    add("limits","Limits that change the scientific interpretation","Límites que cambian la interpretación científica",
        table(["Evidence layer","Supported inference","Boundary"],[
            ["Author reproduction","Numerical recovery using source data","No new biological replicates"],
            ["Bulk / miRNA contrasts","Study-specific regional associations","Composition, small n and source metadata"],
            ["GO and response interactions","Exploratory coordinated patterns / defined contrasts","Multiplicity, correlated terms, model assumptions"],
            ["Measured atlas","Detection at original cell positions","One donor, targeted panel, section gaps"],
            ["Spatial Bayesian models","Conditional reference estimates and held-out diagnostics","Limited baseline gain, prior support, model misspecification"],
            ["PM₂.₅ scenario","Uniform group-effect illustration","No identified cellular exposure response"]]),
        table(["Nivel","Inferencia respaldada","Límite"],[
            ["Reproducción","Recuperación numérica con datos fuente","Sin réplicas biológicas nuevas"],
            ["Bulk / miRNA","Asociaciones regionales de cada estudio","Composición, n pequeño y metadatos"],
            ["GO e interacciones","Patrones exploratorios / contrastes definidos","Multiplicidad, correlación y supuestos"],
            ["Atlas medido","Detección en posiciones originales","Un donante, panel dirigido y espacios entre cortes"],
            ["Bayesianos espaciales","Estimaciones condicionales y diagnósticos reservados","Ganancia limitada, prior y modelo inadecuado"],
            ["Escenario PM₂.₅","Ilustración uniforme de efecto grupal","Sin respuesta celular identificada"]]),
        "These operational limits determine the claims the figures support. One donor cannot establish population-level spatial prevalence. A finer plotting grid cannot increase acquisition resolution. Direct pollutant-exposed spatial RNA or matched single-cell data would be needed to estimate cellular responses rather than assume uniform scaling. Additional independent donors are needed to assess the transfer of spatial priors and calibrate uncertainty beyond one reference specimen.",
        "Estos límites operativos determinan qué afirmaciones pueden respaldar las figuras. Un donante no establece prevalencia espacial poblacional. Una cuadrícula de visualización más fina no eleva la resolución de adquisición. Para estimar la respuesta celular se requieren mediciones espaciales de ARN con exposición o datos single-cell emparejados, sin asumir escalamiento uniforme. Donantes independientes permitirían evaluar la transferencia de priors y calibrar la incertidumbre más allá de una referencia.",
        ["results/replication_pm/summary.json","metadata/spatial_bayes/methods.json","metadata/spatial_detection_bayes/methods.json"],"small")

    add("conclusions","Conclusions supported by the completed work","Conclusiones respaldadas por el trabajo completado",
        """1. Cortex author-aligned reconstruction recovers all four significant-gene sets, alongside an unresolved source-quality signal.
2. New reanalyses expose model, testing-family and small-sample sensitivity. Age-dependent ozone response receives the clearest joint interaction signal.
3. Measured sections, gene products and bounded detection models improve interpretation of the reference atlas.
4. Spatial validation reports limited prediction gains and model failures openly. The work does not establish an exposed cellular brain reconstruction.""",
        "1. La reconstrucción cortical recupera los cuatro conjuntos significativos y conserva una pregunta sobre la calidad de la fuente.\n2. Los reanálisis muestran sensibilidad al modelo, a la familia de pruebas y al tamaño muestral. La respuesta de ozono dependiente de edad ofrece la señal conjunta más clara.\n3. Los cortes medidos, los productos génicos y los modelos de detección acotada mejoran la interpretación del atlas.\n4. La validación espacial informa ganancias limitadas y fallos. El trabajo no establece una reconstrucción celular de cerebro expuesto.",
        "A defense should preserve the strengths and failures together. The strongest exact computational result is cortex count-level reproduction. The most direct newly tested exposure-heterogeneity evidence is the age-dependent joint LRT, with multiplicity-family sensitivity explicit. The spatial contribution is a reproducible measured-reference explorer and evaluated conditional models, without resolution enhancement. “New” denotes a project-specific contribution. Independent experimental validation and complete source metadata remain future requirements.",
        "La defensa conserva las fortalezas y los fallos. El resultado exacto más fuerte es la reproducción cortical a nivel de conteos. La evidencia nueva más directa de diferencias de respuesta es el LRT conjunto de edad, con sensibilidad explícita a la familia de multiplicidad. El aporte espacial es un explorador medido reproducible y modelos condicionales evaluados, sin aumento de resolución. «Nuevo» se refiere al proyecto. La verificación experimental independiente y los metadatos completos siguen pendientes.",
        ["results/replication_pm/summary.json","results/advanced_omics/summary.json","results/spatial_detection_bayes/summary.json"])

    add("future","Future work requires independent biological information","El trabajo futuro requiere información biológica independiente",
        """- Recover source sex/identity metadata and missing author RUVg specifications.
- Validate bulk candidates in independent animals and matched cell populations.
- Evaluate sparse-count or hurdle models and cross-donor spatial prediction.
- Acquire pollutant-exposed spatial/single-cell measurements before inferring cell-specific responses.
- Explore UCSC Genome Browser, Ensembl and ENCODE for compatible brain epigenomic evidence.""",
        "- Recuperar los metadatos fuente de sexo/identidad y las especificaciones RUVg faltantes.\n- Validar candidatos en animales y poblaciones celulares independientes.\n- Evaluar modelos de conteos dispersos o hurdle y predicción entre donantes.\n- Obtener mediciones espaciales/single-cell con exposición antes de inferir respuestas específicas.\n- Explorar UCSC Genome Browser, Ensembl y ENCODE para evidencia epigenómica cerebral compatible.",
        "The epigenomics item records the user's proposed next direction, without completed retrieval or analysis. This deck claims no epigenetic track, enhancer assignment, chromatin accessibility result or NO₂ exposure dataset. Compatibility would require matching species, genome build, brain region, cell type, age and assay. A compatible reference chromatin track would provide baseline regulatory context, rather than a measured pollution response unless the experimental design supports that contrast. Alternative statistical models should be judged using prespecified held-out units and uncertainty calibration, rather than selecting models for more visually detailed maps.",
        "La epigenómica recoge la línea propuesta por el usuario, sin una búsqueda o análisis completado. No se afirman capas epigenéticas, asignaciones de enhancers, resultados de accesibilidad ni un dataset de NO₂. La integración exigiría compatibilidad de especie, ensamblaje, región, tipo celular, edad y ensayo. Una capa de referencia aportaría contexto regulatorio basal, sin demostrar una respuesta medida a la contaminación sin un contraste experimental. Los modelos alternativos deben juzgarse con unidades reservadas y calibración preespecificadas, sin seleccionarlos por crear mapas visualmente detallados.",
        ["results/report_validation/experimental_design_review.json"])

    add("discussion","Reproducibility and discussion","Reproducibilidad y discusión",
        """**Complete evidence:** source mappings, analytical populations, model code, frozen references, tables and bilingual figures.

**Presentation controls:** arrow keys navigate, Esc opens overview, Notes opens the detailed explanation, and the language switch keeps the same slide.

[Scientific report](https://amarusimonaguerojimenez.github.io/Air-Pollution-Brain-Omics-Mouse/) · [Repository](https://github.com/AmaruSimonAgueroJimenez/Air-Pollution-Brain-Omics-Mouse)

The following appendix contains protocols, formulas, sensitivity checks, candidate lists and references.""",
        """**Evidencia completa:** correspondencias con las fuentes, poblaciones analíticas, código, referencias congeladas, tablas y figuras bilingües.

**Controles:** flechas para navegar, Esc para la vista general, Notas para la explicación detallada y cambio de idioma en la misma diapositiva.

[Informe científico](https://amarusimonaguerojimenez.github.io/Air-Pollution-Brain-Omics-Mouse/es/) · [Repositorio](https://github.com/AmaruSimonAgueroJimenez/Air-Pollution-Brain-Omics-Mouse)

El apéndice contiene protocolos, fórmulas, análisis de sensibilidad, listas de candidatos y referencias.""",
        """This is a convenient end of the main presentation. The appendix supports methodological questioning without requiring every technical table during the oral defense. The HTML embeds its figures, stylesheet, runtime and native math for offline display. Online source/report links still require access to the repository or published site. Notes are available through the in-page drawer even when the standard external speaker window cannot load offline. No new model fitting occurs when opening the slides.""",
        """Aquí puede cerrar la presentación principal. El apéndice apoya las preguntas metodológicas sin necesidad de proyectar cada tabla técnica. El HTML incluye figuras, estilos, el motor de presentación y fórmulas nativas para visualizarlo sin conexión. Los enlaces a fuentes y al informe requieren acceso al repositorio o a la página. Las notas están disponibles en el panel de la propia diapositiva, aunque la ventana externa del presentador no pueda cargar sin conexión. Abrir las diapositivas no ejecuta ajustes nuevos.""",
        ["metadata/analysis_environment.json"])

    section("appendix","Technical appendix","Apéndice técnico",
        "Protocols, testing families, prior definitions and reproducibility",
        "Protocolos, familias de pruebas, definición de priors y reproducibilidad",
        "The appendix supports detailed examination. It is not necessary to project all appendix slides during the main presentation.",
        "El apéndice apoya el examen detallado. No es necesario proyectarlo completo durante la presentación principal.")

    add("a-pm-protocols","PM₂.₅ exposure protocols remain distinct","Los protocolos de PM₂.₅ permanecen diferenciados",
        table(["Arm","Starting context","Protocol","Timing caution"],[
            ["HY / GSE145840","Male C57BL/6J, 3 weeks","About 60–120 µg/m³, mean 80; 6 h/day, 5 days/week, 14 weeks","Terminal age/harvest lag not independently reported"],
            ["ACA / GSE297423","hWT αSyn, 2 months, sex unreported","2 µL of 5 mg/mL PM intranasal, every 2 days for 2 months","Separate from older conventional-mouse arm"],
            ["miRNA / GSE93967","Male C57BL/6, about 8 weeks","1 or 5 mg/kg oropharyngeal aspiration, every other day for 4 weeks","Harvest 18 h after final dose"]]),
        table(["Brazo","Contexto inicial","Protocolo","Precaución temporal"],[
            ["HY / GSE145840","C57BL/6J macho, 3 semanas","60–120 µg/m³, media 80; 6 h/día, 5 días/semana, 14 semanas","Edad terminal e intervalo hasta la toma de tejido no informados"],
            ["ACA / GSE297423","hWT αSyn, 2 meses, sexo no informado","2 µL de 5 mg/mL intranasal, cada 2 días por 2 meses","Distinto del brazo convencional de mayor edad"],
            ["miRNA / GSE93967","C57BL/6 macho, aprox. 8 semanas","1 o 5 mg/kg por aspiración orofaríngea, días alternos por 4 semanas","Toma de tejido 18 h tras la dosis final"]]),
        """These protocols come from the saved source-paper/deposited metadata review, with per-arm evidence locators in experimental_design_review.json. Atmosphere concentration is not converted into an instilled body-weight dose. Arithmetic starting-age plus duration is not substituted for a reported terminal age. The cortex and miRNA liquid administration routes should not be described as chamber inhalation. Exact platform and sample IDs remain linked in the accession metadata.""",
        """Los protocolos proceden de la revisión guardada de artículos y metadatos, con localizadores por brazo en experimental_design_review.json. No se transforma concentración atmosférica en dosis instilada por peso. La suma aritmética de edad inicial y duración no reemplaza una edad terminal informada. Las administraciones líquidas cortical y miRNA no se describen como inhalación en cámara. Las plataformas e IDs permanecen en los metadatos de acceso.""",
        ["HY","ACA","miRNA","results/report_validation/experimental_design_review.json"],"appendix small")

    add("a-ozone-protocol","Ozone protocol and sequenced sample size","Protocolo de ozono y tamaño muestral secuenciado",
        """- Male human APOE3/APOE4 targeted-replacement mice.
- Start exposure at 3 or 17 months.
- **0.8 ppm O₃, 7 h/day for five days**, followed by nine days of filtered air.
- Repeat eight cycles. Left hippocampus enters RNA-seq.
- Paper-wide 14–16 animals per group differs from **three sequenced animals per stratum**.
- One SDRF strain string differs for Sample1 / ERR14708304.""",
        """- Ratones machos de reemplazo dirigido APOE3/APOE4 humano.
- Inicio de exposición a los 3 o 17 meses.
- **0.8 ppm O₃, 7 h/día por 5 días**, seguidos de 9 días con aire filtrado.
- Ocho ciclos. RNA-seq del hipocampo izquierdo.
- Los 14–16 animales/grupo del artículo difieren de **tres secuenciados por estrato**.
- Una cadena de cepa SDRF difiere en Sample1 / ERR14708304.""",
        """Eight 14-day cycles imply 112 days, a calendar calculation rather than an individually recorded collection date. Chamber geometry and the exact lag from the last ozone exposure to tissue collection are not independently established in the reviewed source, so no whole-body/nose-only label or precise terminal age is added. The one strain-string exception remains an annotation caveat, without an assertion that it confounds the entire genotype comparison. The sequenced subset defines the inferential n.""",
        """Ocho ciclos de 14 días implican 112 días como cálculo de calendario, no como fecha individual de toma de tejido. La fuente revisada no establece independientemente la geometría de la cámara ni el intervalo exacto desde la última exposición a ozono hasta la toma de tejido. Por ello, no se agrega una etiqueta de exposición de cuerpo completo o solo por nariz, ni una edad terminal precisa. La excepción de cepa se conserva como cautela, sin declararla confusión completa por genotipo. El subconjunto secuenciado define el n inferencial.""",
        ["O3","results/report_validation/experimental_design_review.json"],"appendix")

    add("a-normalization","Median-of-ratios normalization and length offsets","Mediana de cocientes y offsets de longitud",
        r"""For count genes with a usable geometric mean:

$$G_g=\left(\prod_i K_{gi}\right)^{1/m},\qquad
s_i=\operatorname{median}_{g}\left(\frac{K_{gi}}{G_g}\right)$$

Normalized PM counts are $K_{gi}/s_i$. Size factors use the retained analysis universe.

Ozone uses sample/gene effective-length normalization through tximport. TPM is not substituted as negative-binomial input.

VST/PCA and row-standardized heatmaps are descriptive transforms, separate from the count likelihood.""",
        r"""Para genes con media geométrica utilizable:

$$G_g=\left(\prod_i K_{gi}\right)^{1/m},\qquad
s_i=\operatorname{median}_{g}\left(\frac{K_{gi}}{G_g}\right)$$

Los conteos PM normalizados son $K_{gi}/s_i$. Los factores se estiman con el universo retenido.

Ozono usa normalización gen/muestra por longitud efectiva con tximport. TPM no sustituye la entrada binomial negativa.

VST/PCA y heatmaps estandarizados por fila son transformaciones descriptivas, separadas de la verosimilitud.""",
        """This is the standard ratio definition on usable positive geometric means, not a promise that every retained gene contributes to every factor. The project saves estimated size factors and unchanged source matrices. Within-gene heatmap z scores remove magnitude and are not cross-gene expression levels. A PCA dimension is not a biological cell type or an adjusted exposure effect. Those visual transformations do not alter the original likelihood input or multiplicity family.""",
        """Es la definición de cocientes con medias geométricas positivas utilizables, sin afirmar que todos los genes retenidos aportan a cada factor. Se guardan factores y matrices fuente intactas. Los z de heatmap por gen eliminan magnitud y no equivalen a expresión comparable entre genes. Un eje PCA no es un tipo celular ni efecto ajustado. Las transformaciones gráficas no alteran entrada o familia de pruebas.""",
        ["DESeq2","tximport","results/hypothalamus/methods.json"],"appendix formula-slide")

    add("a-limma","Empirical-Bayes moderation in limma","Moderación Bayes empírica en limma",
        r"""$$s_{f,\mathrm{post}}^2=
\frac{d_0s_{0f}^2+d_fs_f^2}{d_0+d_f}$$

$$t_{f,\mathrm{mod}}=
\frac{\widehat\Delta_f}{s_{f,\mathrm{post}}\sqrt{v_f}}$$

The prior variance trend depends on signal level. Moderation borrows information across features.

Pointwise moderated t intervals and BH across features answer different questions.

No verified pairing or batch variable is invented for the nine arrays.""",
        r"""$$s_{f,\mathrm{post}}^2=
\frac{d_0s_{0f}^2+d_fs_f^2}{d_0+d_f}$$

$$t_{f,\mathrm{mod}}=
\frac{\widehat\Delta_f}{s_{f,\mathrm{post}}\sqrt{v_f}}$$

La tendencia de varianza previa depende del nivel de señal. La moderación comparte información entre características.

Los intervalos t puntuales y BH entre características responden preguntas diferentes.

No se inventa pareamiento o lote para los nueve arrays.""",
        """The formula expresses the residual variance moderation used by limma, with v_f the unscaled contrast variance. trend=TRUE permits a signal-dependent prior scale. The model uses robust=FALSE. Total moderated degrees of freedom enter t tests and intervals, subject to limma's implementation conventions. This is an empirical-Bayes array analysis, distinct from the fixed-prior hierarchical spatial model. It does not resolve uncertainty in the supplied measurement scale.""",
        """La fórmula representa la moderación residual limma, con v_f como varianza no escalada del contraste. trend=TRUE permite escala previa dependiente de señal; robust=FALSE. Los grados de libertad moderados entran en t e intervalos bajo convenciones de implementación limma. Es Bayes empírico de arrays, distinto del modelo espacial con priors fijos. No resuelve la incertidumbre de escala depositada.""",
        ["limma","results/mirna/methods.json"],"appendix formula-slide")

    add("a-multiplicity","Multiplicity families and non-nested discovery lists","Familias de multiplicidad y listas no anidadas",
        table(["Example","Family","BH discoveries"],[
            ["O₃ stratum tests","Each stratum separately","0 / 399 / 2 / 3"],
            ["O₃ stratum tests","All four valid p-value families","1 / 43 / 4 / 6"],
            ["Age-dependent exposure LRT","Its own gene family","8"],
            ["Three secondary LRTs","All valid tests pooled","1 age-response hit"]])+
        "\nPooled BH adjusts a different p-value collection. Neither list is selected to maximize discoveries.",
        table(["Ejemplo","Familia","Hallazgos BH"],[
            ["Pruebas O₃ por estrato","Cada estrato separado","0 / 399 / 2 / 3"],
            ["Pruebas O₃ por estrato","Cuatro familias válidas juntas","1 / 43 / 4 / 6"],
            ["LRT exposición por edad","Familia propia de genes","8"],
            ["Tres LRT secundarios","Pruebas válidas agrupadas","1 hallazgo por edad"]])+
        "\nBH conjunta ajusta otra colección de valores p. No se elige la lista que maximiza los hallazgos.",
        """A pooled BH family can alter both m and the ranks of a stratum's p values, so selected sets need not nest within a smaller separately corrected family. All denominator and validity rules are preserved. The primary hypothesis and supplementary sensitivity remain labeled. GSEA algorithmic FDR is not inserted into these BH comparisons. A q threshold is a decision rule under a model, not a discontinuous biological boundary.""",
        """BH conjunta cambia tanto m como el rango de los p de un estrato; los conjuntos no tienen por qué anidarse en la corrección separada. Se conservan denominadores y criterios de validez. Se identifican la hipótesis principal y la sensibilidad suplementaria. La FDR algorítmica GSEA no entra en estas comparaciones BH. El umbral q es regla de decisión bajo un modelo, no frontera biológica discontinua.""",
        ["results/ozone/summary.json","results/advanced_omics/summary.json"],"appendix")

    software=J("results/hypothalamus/independent_deseq2_check/summary.json")
    loo=T("results/hypothalamus/leave_one_out_summary.tsv")
    add("a-software-loo","Software and leave-one-animal-out sensitivity","Sensibilidad a software y omisión de animales",
        table(["HY check","Saved result"],[
            ["PyDESeq2 versus R BH hits","11 versus 12"],
            ["Effect Pearson correlation",f"{software['numerical_agreement_metrics']['common_universe_log2FC_Pearson']:.10f}"],
            ["Cook status agreement","All 19,611 retained genes"],
            ["Discordant BH threshold gene","Tnrc6c: Python 0.05211, R 0.03992"],
            ["Leave-one-out BH hit range",f"{loo.genes_fdr_below_0_05.min()}–{loo.genes_fdr_below_0_05.max()}"],
            ["Direction of 11 primary hits","Preserved in all seven omissions"]]),
        table(["Verificación HY","Resultado guardado"],[
            ["Hallazgos BH PyDESeq2 frente a R","11 frente a 12"],
            ["Correlación Pearson de efectos",f"{software['numerical_agreement_metrics']['common_universe_log2FC_Pearson']:.10f}"],
            ["Acuerdo de estado Cook","Los 19,611 genes retenidos"],
            ["Gen discordante en umbral BH","Tnrc6c: Python 0.05211, R 0.03992"],
            ["Rango de hallazgos leave-one-out",f"{loo.genes_fdr_below_0_05.min()}–{loo.genes_fdr_below_0_05.max()}"],
            ["Dirección de 11 hallazgos principales","Conservada en las siete omisiones"]]),
        """Tnrc6c effects differ by only about −3.3×10⁻⁶ log2 units, while small SE differences cross the q decision boundary. This is implementation sensitivity, not evidence that one software is biologically correct. Neither the software choice nor an animal omission was selected to increase hits. No primary hit retains q < .05 in every omission, even when every direction agrees. These checks illustrate why effect continuity, uncertainty and threshold decisions should be discussed together.""",
        """Los efectos de Tnrc6c difieren solo aproximadamente −3.3×10⁻⁶ unidades log₂, pero pequeñas diferencias del error estándar cruzan el umbral q. Es sensibilidad de implementación, sin definir qué software es biológicamente correcto. No se eligieron el software o las omisiones para aumentar hallazgos. Ningún hallazgo principal conserva q < .05 en todas las omisiones, aunque la dirección siempre coincida. Estas comprobaciones permiten discutir conjuntamente continuidad del efecto, incertidumbre y decisiones basadas en umbrales.""",
        ["results/hypothalamus/independent_deseq2_check/summary.json","results/hypothalamus/leave_one_out_summary.tsv"],"appendix")

    for region,group,name_en,name_es in [
        ("HY","hypothalamus_measured_top10","Hypothalamic measured candidates","Candidatos hipotalámicos medidos"),
        ("ACA","cortex_measured_top10","Cortical measured candidates","Candidatos corticales medidos")]:
        sub=genes[genes.selection_group==group].sort_values("selection_rank")
        rows=[[r.gene,f"{r.padj:.3g}",f"{r.log2FoldChange:+.3f}",f"{r.log2FoldChange_shrunk:+.3f}",
               "q<.05" if r.primary_significant else "q≥.05"] for r in sub.itertuples()]
        note_en="Only Tll1 is primary-BH significant in HY; nine displayed candidates remain exploratory." if region=="HY" else "All ten displayed ACA candidates satisfy the primary regional BH threshold."
        note_es="Solo Tll1 cumple BH principal en HY; nueve candidatos permanecen exploratorios." if region=="HY" else "Los diez candidatos ACA cumplen el umbral BH regional principal."
        add("a-genes-"+region.lower(),name_en,name_es,
            table(["Gene","Primary q","MLE log₂FC","Shrunken log₂FC","Evidence"],rows)+"\n"+note_en,
            table(["Gen","q principal","log₂FC MLE","log₂FC regularizado","Evidencia"],rows)+"\n"+note_es,
            """The rows preserve the frozen measured-panel priority rule and original regional inference. The source table also provides pointwise 95% intervals, gene products, biotypes, source URLs and QC flags. These estimates are regional bulk effects, not individual-cell changes. The Bayesian atlas models use baseline measurements of these candidate transcripts and do not refit the regional contrast.""",
            """Se conservan el orden congelado del panel medido y la inferencia regional original. La tabla fuente aporta además intervalos puntuales del 95%, productos, biotipos, URLs e indicadores de control de calidad. Son efectos regionales bulk, no cambios de células individuales. Los modelos del atlas usan mediciones basales de estos transcritos sin reajustar el contraste regional.""",
            ["results/gene_priorities/prioritized_genes.tsv"],"appendix tiny")

    add("a-gaussian-assumptions","Gaussian assumptions and prespecified sensitivity","Supuestos gaussianos y sensibilidad preespecificada",
        table(["Component","Primary specification","Sensitivity / limitation"],[
            ["Residual","Normal, variance σ²","Sparse nonnegative response may violate support"],
            ["Section effect","Normal, variance 0.25σ²","New section contributes 0.25σ²"],
            ["RBF lengths AP/DV/ML","0.60 / 0.45 / 0.45 mm","Multiply all by 0.5 or 2"],
            ["Spatial prior variances","Fixed diagonal entries","Multiply by 0.25 or 4"],
            ["Unknown anatomy","Prior mean 0, extra variance 1×σ²","Report unavailable training anatomy"],
            ["Validation","5 guarded whole-section blocks","Within one donor only"]]),
        table(["Componente","Especificación principal","Sensibilidad / límite"],[
            ["Residual","Normal, varianza σ²","La respuesta dispersa puede violar el soporte"],
            ["Efecto de corte","Normal, varianza 0.25σ²","Un corte nuevo aporta 0.25σ²"],
            ["Longitudes RBF AP/DV/ML","0.60 / 0.45 / 0.45 mm","Multiplicar todas por 0.5 o 2"],
            ["Varianzas espaciales previas","Diagonal fija","Multiplicar por 0.25 o 4"],
            ["Anatomía desconocida","Media previa 0, varianza extra 1×σ²","Informar anatomía ausente del entrenamiento"],
            ["Validación","5 bloques con separación","Dentro de un donante"]]),
        """The implementation performs 60 fold-level Gaussian fits and two full-reference fits across the two regions. Each solve processes ten response columns while retaining gene-specific coefficients and residual scales. The fold-level fits cover the primary spatial model, the anatomy-only baseline and four prespecified spatial sensitivity configurations. kNN is evaluated separately. A successful linear algebra check or posterior simulation test verifies computation, without validating the Gaussian likelihood. Hyperparameters remain fixed, with no held-out tuning.""",
        """La implementación realiza 60 ajustes gaussianos por bloque y dos ajustes con la referencia completa entre las dos regiones. Cada resolución procesa diez columnas de respuesta, conservando coeficientes y escalas residuales específicos por gen. Los ajustes por bloque cubren el modelo espacial principal, la referencia de solo anatomía y cuatro configuraciones espaciales de sensibilidad preespecificadas. kNN se evalúa aparte. Verificar álgebra lineal o transporte de incertidumbre por simulación comprueba el cálculo, sin validar la verosimilitud gaussiana. Los hiperparámetros son fijos, sin ajustar con resultados reservados.""",
        ["metadata/spatial_bayes/methods.json","results/spatial_bayes/summary.json"],"appendix")

    add("a-beta-marginal","Beta-binomial likelihood and new-section prediction","Verosimilitud beta-binomial y predicción de corte nuevo",
        r"""$$P(Y=y\mid n,\mu,\kappa)=\binom ny
\frac{B(y+\mu\kappa,n-y+(1-\mu)\kappa)}
{B(\mu\kappa,(1-\mu)\kappa)}$$

$$\rho=\frac1{\kappa+1}$$

Known-section components update to

$$p_s\mid y,n,\mu,\kappa\sim
\operatorname{Beta}(\mu\kappa+y,(1-\mu)\kappa+n-y).$$

Hyperposterior-weighted mixtures preserve parameter uncertainty.""",
        r"""$$P(Y=y\mid n,\mu,\kappa)=\binom ny
\frac{B(y+\mu\kappa,n-y+(1-\mu)\kappa)}
{B(\mu\kappa,(1-\mu)\kappa)}$$

$$\rho=\frac1{\kappa+1}$$

Los componentes del corte conocido se actualizan:

$$p_s\mid y,n,\mu,\kappa\sim
\operatorname{Beta}(\mu\kappa+y,(1-\mu)\kappa+n-y).$$

Las mezclas ponderadas por hiperpósterior conservan incertidumbre paramétrica.""",
        """For a new section, draw or integrate p from the unupdated Beta(muκ, (1−mu)κ) mixture, then sample a binomial count at the desired n if the target is Y. The known-section posterior conditions on the observed section itself, which is legitimate for fitted estimation but not for an out-of-fold prediction. Conditional intraclass correlation rho summarizes the beta-binomial overdispersion model, not an experimentally measured cell correlation. Fits are independent by anatomy/gene and do not borrow a spatially adjacent anatomy's data as if it were the same category.""",
        """Para un corte nuevo se integra o sortea p de la mezcla Beta no actualizada y después un conteo binomial con n si el objetivo es Y. La posterior conocida condiciona al propio corte, algo válido para estimación ajustada pero no para predicción reservada. Rho resume la sobredispersión bajo el modelo, no una correlación celular medida experimentalmente. El ajuste es independiente por anatomía/gen, sin tratar anatomías vecinas como una misma categoría.""",
        ["metadata/spatial_detection_bayes/methods.json"],"appendix formula-slide")

    add("a-coverage","Marginal coverage and outcome-conditioned diagnostics","Cobertura marginal y diagnósticos condicionados por resultado",
        r"""A nominal prediction interval targets a marginal statement:

$$P\{Y_*\in I_{0.95}(X_*,D)\}\approx0.95.$$

It does not impose

$$P\{Y_*\in I_{0.95}(X_*,D)\mid Y_*>0\}=0.95.$$

Positive-only inclusion is a descriptive tail-recovery diagnostic.

Impossible negative predictive support is a separate likelihood limitation.""",
        r"""Un intervalo nominal busca una propiedad marginal:

$$P\{Y_*\in I_{0.95}(X_*,D)\}\approx0.95.$$

No exige

$$P\{Y_*\in I_{0.95}(X_*,D)\mid Y_*>0\}=0.95.$$

La inclusión solo positiva describe recuperación de la cola.

El soporte predictivo negativo imposible es otra limitación de verosimilitud.""",
        """This distinction prevents the misleading conclusion that positive-only inclusion below 95% by itself disproves marginal calibration. The Gaussian model can have high overall coverage because zeros dominate. Poor positive-tail recovery remains scientifically informative for the goal of reconstructing positive RNA observations. Combined with negative support and weak baseline performance, it limits usefulness. Beta count-interval coverage uses a different discrete target and must be interpreted alongside its substantial widths and prior-only strata.""",
        """Esta distinción evita afirmar que la inclusión positiva bajo el 95% por sí sola refuta la calibración marginal. La abundancia de ceros puede producir una cobertura global alta. La mala recuperación de positivos sigue siendo informativa para reconstruir observaciones de ARN. Junto al soporte negativo y al rendimiento basal limitado, restringe la utilidad del modelo. La cobertura beta de conteos tiene otro objetivo discreto y se interpreta junto con sus amplitudes y estratos sostenidos solo por el prior.""",
        ["metadata/spatial_bayes/methods.json","results/spatial_bayes/summary.json"],"appendix formula-slide")

    add("a-outstanding","Data and analyses outside the completed scope","Datos y análisis fuera del alcance completado",
        """- NO₂ brain papers remain in the literature inventory, without a complete expression matrix entering these four models.
- No exposed spatial or single-cell brain dataset is synthesized as measured evidence.
- Genes outside the targeted MERFISH panel remain unmeasured.
- UCSC, Ensembl and ENCODE epigenomic integration remains proposed future work.
- Independent cross-donor and exposure validation remains outstanding.""",
        """- Los artículos cerebrales NO₂ permanecen en el inventario, sin matriz completa incorporada a estos cuatro modelos.
- No se sintetiza un dataset cerebral espacial/single-cell expuesto como evidencia medida.
- Los genes fuera del panel MERFISH permanecen no medidos.
- La integración epigenómica con UCSC, Ensembl y ENCODE es trabajo futuro.
- Siguen pendientes validación entre donantes y exposición independientes.""",
        """A defense deck should not turn an earlier search, downloaded paper or proposed future direction into a completed analysis. The original scope includes PM2.5 and O3 brain experiments and a separate unexposed spatial reference. An unavailable matrix cannot be replaced by supplementary hit lists for sample-level modeling. Cross-study differences in routes, ages, genotypes and measurement scales preclude a pooled pollutant coefficient in this project.""",
        """No se convierte una búsqueda, un artículo descargado o una propuesta en un análisis terminado. El alcance incluye experimentos cerebrales con PM₂.₅/O₃ y una referencia espacial sin exposición. Las listas suplementarias no sustituyen una matriz para modelar muestras. Las diferencias de ruta, edad, genotipo y escala impiden un coeficiente conjunto de contaminantes en este proyecto.""",
        ["results/report_validation/experimental_design_review.json","metadata/papers/reuse_search_summary.json"],"appendix")

    primary_refs=[
        "**Rajagopalan et al. (2020).** Metabolic effects of air pollution exposure and reversibility. J Clin Invest. [doi:10.1172/JCI137315](https://doi.org/10.1172/JCI137315)",
        "**Zhang et al. (2025).** Lewy body dementia promotion by air pollutants. Science. [doi:10.1126/science.adu4132](https://doi.org/10.1126/science.adu4132)",
        "**Ku et al. (2017).** NF-κB-regulated microRNA-574-5p underlies synaptic and cognitive impairment in response to atmospheric PM₂.₅ aspiration. Particle and Fibre Toxicology. [doi:10.1186/s12989-017-0215-3](https://doi.org/10.1186/s12989-017-0215-3)",
        "**Nakamya et al. (2025).** Age- and ApoE Genotype-Dependent Transcriptomic Responses to O₃ in the Hippocampus of Mice. Int J Mol Sci. [doi:10.3390/ijms26062407](https://doi.org/10.3390/ijms26062407)",
        "**Zhang et al. (2023).** Molecularly defined and spatially resolved cell atlas of the whole mouse brain. Nature. [doi:10.1038/s41586-023-06808-9](https://doi.org/10.1038/s41586-023-06808-9)"]
    add("a-primary-refs","Primary studies and spatial reference","Estudios originales y referencia espacial",
        "\n\n".join(primary_refs),"\n\n".join(primary_refs),
        "These are the primary source studies, not five matched exposure cohorts. Exact full bibliographic records, article-reading copies and source locators are retained locally. Public material links to original publications rather than redistributing restricted full text.",
        "Estas son fuentes originales, no cinco cohortes de exposición emparejadas. Se conservan registros bibliográficos completos, copias de lectura y localizadores. El material público enlaza los artículos originales en vez de redistribuir textos restringidos.",
        ["papers/bibliography.bib"],"appendix small")

    method_refs=[
        "**DESeq2:** Love et al. 2014. [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)",
        "**PyDESeq2:** Muzellec et al. 2023. [10.1093/bioinformatics/btad547](https://doi.org/10.1093/bioinformatics/btad547)",
        "**limma:** Ritchie et al. 2015. [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)",
        "**Salmon:** Patro et al. 2017. [10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197)",
        "**tximport:** Soneson et al. 2015. [10.12688/f1000research.7563.2](https://doi.org/10.12688/f1000research.7563.2)",
        "**apeglm:** Zhu et al. 2019. [10.1093/bioinformatics/bty895](https://doi.org/10.1093/bioinformatics/bty895)"]
    add("a-method-refs","Statistical software references","Referencias de software estadístico",
        "\n\n".join(method_refs),"\n\n".join(method_refs),
        "The presentation describes the installed frozen versions in the repository, rather than claiming that every contemporary release has identical behavior. Versioned R session and Python environment records accompany analyses. Software citation does not establish biological validity.",
        "La presentación describe versiones congeladas del repositorio, sin afirmar que toda versión actual se comporte igual. Los análisis incluyen registros de sesiones R y del entorno Python. Citar software no establece validez biológica.",
        ["metadata/analysis_environment.json"],"appendix small")

    add("a-bayes-refs","Bayesian derivation and implementation sources","Fuentes de derivación e implementación bayesiana",
        """- Fixed-basis Bayesian regression and prediction: [Rasmussen and Williams, chapter 2](https://gaussianprocess.org/gpml/chapters/RW2.pdf).
- Conjugate Gaussian analysis: [Murphy, technical note](https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf).
- Beta-binomial distribution: [official SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.betabinom.html).
- Project formulas, priors and validation: frozen method files and source implementations.

The reference models are project implementations with explicitly evaluated assumptions.""",
        """- Regresión bayesiana con bases y predicción: [Rasmussen y Williams, capítulo 2](https://gaussianprocess.org/gpml/chapters/RW2.pdf).
- Análisis gaussiano conjugado: [nota técnica de Murphy](https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf).
- Distribución beta-binomial: [documentación oficial de SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.betabinom.html).
- Fórmulas, priors y validación del proyecto: métodos congelados e implementaciones.

Los modelos de referencia son implementaciones del proyecto con supuestos evaluados.""",
        "These references support the mathematical foundations. They are not validations of the particular fixed length scales, priors, anatomy categories or Gaussian count approximation used here. The complete source code defines the implementation.",
        "Estas fuentes fundamentan la matemática, sin validar las longitudes, los priors, las categorías anatómicas o la aproximación gaussiana particulares. El código completo define la implementación.",
        ["metadata/spatial_bayes/methods.json","metadata/spatial_detection_bayes/methods.json"],"appendix")

    add("a-rebuild","Rebuilding the presentation without refitting the science","Reconstruir la presentación sin reajustar los análisis",
        """The builder reads saved result tables and produces presentation-sized figures, paired QMD sources and HTML.

    .venv-analysis/bin/python scripts/build_defense_slides.py --figures
    .venv-analysis/bin/python scripts/render_docs.py --page defense.qmd
    .venv-analysis/bin/python scripts/build_defense_slides.py --package
    .venv-analysis/bin/python scripts/validate_defense_slides.py

Source hashes, equal EN/ES slide IDs, figure provenance and browser review accompany the deck.""",
        """El constructor lee tablas guardadas y genera figuras para proyectar, fuentes QMD pareadas y HTML.

    .venv-analysis/bin/python scripts/build_defense_slides.py --figures
    .venv-analysis/bin/python scripts/render_docs.py --page defense.qmd
    .venv-analysis/bin/python scripts/build_defense_slides.py --package
    .venv-analysis/bin/python scripts/validate_defense_slides.py

Las presentaciones conservan hashes, IDs EN/ES iguales, procedencia de las figuras y revisión en el navegador.""",
        "Rebuilding this presentation does not rerun the original RNA-seq/array models, raw-read pipeline or Bayesian fits. Their separate reproducibility commands remain in the website and source scripts. The presentation audit verifies scope, source consistency, slide structure, assets and layout within checked limits, not universal scientific validity.",
        "Reconstruir esta presentación no repite modelos RNA-seq/arrays, el procesamiento de lecturas o los ajustes bayesianos. Los comandos científicos separados permanecen en la web y los scripts. La auditoría revisa el alcance, la consistencia, la estructura, los recursos y el diseño dentro de sus límites, sin certificar validez científica universal.",
        ["metadata/analysis_environment.json"],"appendix small")


def source_url(source, lang):
    if source in PAPERS:
        return PAPERS[source]
    if source.startswith("http"):
        return ("Source" if lang=="en" else "Fuente",source)
    p=Path(source)
    if not (ROOT/p).exists():
        raise FileNotFoundError(source)
    if p.parts[0]=="results":
        target=BASE+"downloads/"+str(Path(*p.parts[1:]))
    elif p.parts[0]=="metadata":
        target=BASE+"downloads/provenance/"+str(Path(*p.parts[1:]))
    elif p.parts[0]=="papers":
        target=BASE+"downloads/papers/"+p.name
    else:
        target=GH+source
    label="/".join(p.parts[-2:])
    return label,target

def scholarly_references(s):
    """Scientific citations are publications/resources, never internal JSON logs."""
    sid=s["id"]
    if sid in {"discussion","a-rebuild","future","a-outstanding"}: return []
    if "gaussian" in sid:
        return ["atlas","Murphy","GPML"][:2 if "diagnostics" in sid or "rmse" in sid else 3]
    if sid in {"a-coverage","metrics","numerics"}:
        return ["atlas","Murphy","Beta"]
    if "beta" in sid or sid in {"detection-method","guarded-cv"}:
        return ["atlas","Beta"] if "beta" in sid else ["atlas"]
    if sid in {"ora","gsea","enrichment-results","go-patterns"}:
        return ["GSEA","GO"]
    if sid=="gene-selection": return ["atlas","HY","ACA"]
    if sid=="gene-products": return ["NCBI","UniProt","atlas"]
    if sid in {"scope","intro","units","workflow","reproduction-criteria","reproduction-summary","limits","conclusions"}:
        return ["HY","ACA","O3"]
    keys=[p for p in s["sources"] if p in PAPERS]
    for source in s["sources"]:
        if "/hypothalamus/" in source: keys.append("HY")
        if "/cortex/" in source: keys.append("ACA")
        if "/mirna/" in source: keys.append("miRNA")
        if "/ozone/" in source or "advanced_omics" in source: keys.append("O3")
        if "spatial_" in source: keys.append("atlas")
        if "replication_pm" in source: keys.append("HY" if "hy-" in sid or "hypothalamus" in sid else "ACA")
        if "replication_mirna_ozone" in source: keys.append("miRNA" if "mirna" in sid else "O3")
    if sid=="count-model": keys=["DESeq2","PyDESeq2"]
    if sid in {"testing","filters","qc","a-normalization","a-multiplicity","ozone-factorial","ozone-contrasts","ozone-lrt","interaction-results"}:
        keys=list(keys)+["DESeq2"]
    if "mirna" in sid or sid=="a-limma": keys=list(keys)+["limma"]
    if sid=="a-software-loo": keys=["DESeq2","PyDESeq2","HY"]
    if sid=="a-bayes-refs": keys=["GPML","Murphy","Beta"]
    return list(dict.fromkeys(keys))[:3]

def phase(s):
    index=SLIDES.index(s)
    mark={v["id"]:j for j,v in enumerate(SLIDES)}
    if index<mark["methods"]: return ("01  /  Context","01  /  Contexto")
    if index<mark["replicated-results"]: return ("02  /  Methods","02  /  Metodología")
    if index<mark["new-results"]: return ("03  /  Author reproduction","03  /  Reproducción de autores")
    if index<mark["limits"]: return ("04  /  Project results","04  /  Resultados del proyecto")
    if index<mark["appendix"]: return ("05  /  Interpretation","05  /  Interpretación")
    return ("06  /  Technical appendix","06  /  Apéndice técnico")

def composition(s,lang):
    """Compose the slide instead of treating every page as the same text layout."""
    body=s[lang]
    classes=[]
    if s["cls"]=="section-slide": return body,classes
    if ".formula" in s["cls"] or "formula-slide" in s["cls"]:
        formulas=re.findall(r"\$\$[\s\S]*?\$\$",body)
        if not formulas: return body,classes
        rest=re.sub(r"\$\$[\s\S]*?\$\$","",body).strip()
        wide=s["id"] in {"testing","mirna-model","ozone-factorial","ozone-contrasts","ozone-lrt","gaussian-priors","gaussian-intervals","scenario-method","a-coverage","gaussian-features"}
        core="\n\n".join(formulas)
        if s["id"]=="gaussian-intervals":
            labels=("Section-adjusted latent mean","New cell in a new physical section") if lang=="en" else ("Media latente ajustada por corte","Célula nueva en un corte físico nuevo")
            core="\n\n".join(f"**{label}**\n\n{formula}" for label,formula in zip(labels,formulas))
            rest=re.sub(r"Section-adjusted latent mean:|New cell in a new physical section:|Media latente ajustada por corte:|Célula nueva en un corte físico nuevo:","",rest).strip()
        inside=block(core,"equation-core")+block(rest,"equation-detail")
        classes=["layout-equation"]
        if wide:
            classes.append("wide-equation")
            return inside,classes
        return block(inside,"equation-grid"),classes
    hero={
        "scope":("04","independent exposure studies","estudios de exposición independientes"),
        "units":("01","unexposed spatial donor","donante espacial sin exposición"),
        "mirna-primary":("1 / 0","BH discoveries at low / high dose","hallazgos BH a dosis baja / alta"),
        "interaction-results":("08","age-dependent response genes; 1 under pooled BH","genes de respuesta por edad; 1 con BH conjunta"),
        "conclusions":("04","conclusions supported by the evidence","conclusiones respaldadas por evidencia"),
        "future":("NEXT","independent biological validation","validación biológica independiente"),
    }
    if s["id"] in hero:
        number,en,es=hero[s["id"]]
        label=en if lang=="en" else es
        left=f'<div class="statement-number">{number}</div><div class="statement-label">{label}</div>'
        return block(block(left,"statement-metric")+block(body,"statement-copy"),"statement-grid"),["layout-statement"]
    if s["id"] in {"studies-rna","studies-hippocampus","a-pm-protocols"}:
        classes.append("layout-study")
    if "results" in s["id"] or s["id"] in {"reproduction-summary","numerics","limits"}:
        classes.append("layout-result")
    return body,classes

def write_deck(lang):
    prefix="../" if lang=="es" else ""
    title=("Air pollution and mouse brain omics" if lang=="en" else "Contaminación atmosférica y ómica cerebral de ratón")
    subtitle=("Detailed methods, original-study reproduction and project analyses" if lang=="en" else "Metodología detallada, reproducción de estudios y análisis del proyecto")
    yaml=f"""---
title: "{title}"
subtitle: "{subtitle}"
date: '{DATE}'
date-format: D MMMM YYYY
lang: {lang}
page-layout: full
toc: false
format:
  revealjs:
    theme: simple
    css: {prefix}assets/defense-slides/defense.css
    width: 1600
    height: 900
    margin: 0.055
    center: false
    transition: none
    background-transition: none
    slide-number: c/t
    controls: true
    controls-layout: bottom-right
    progress: true
    hash: true
    navigation-mode: linear
    keyboard: true
    overview: true
    touch: true
    html-math-method: mathml
    embed-resources: true
    pdf-separate-fragments: false
    link-external-newwindow: true
    show-notes: false
    title-slide-attributes:
      data-background-color: "#103d49"
    include-after-body:
      text: |
        <script src="{prefix}assets/defense-slides/defense.js"></script>
execute:
  echo: false
  eval: false
---

"""
    chunks=[yaml]
    for s in SLIDES:
        body,more_classes=composition(s,lang)
        classes=(s["cls"]+" "+" ".join(more_classes)).strip()
        attrs=f"#{s['id']} "+("."+classes.replace(" "," .") if classes else "")
        if s["cls"]=="section-slide":
            attrs+=' data-background-color="#103d49"'
        chunks.append(f"## {s[lang+'_title']} {{{attrs.strip()}}}\n\n")
        chunks.append(block(phase(s)[0 if lang=="en" else 1],"slide-eyebrow"))
        if s.get("figure"):
            stem=s["figure"]+("_es" if lang=="es" else "")
            image_path=f"{prefix}assets/figures/defense_slides/{stem}.png"
            chunks.append(f'![]({image_path}){{fig-alt="{html.escape(s[lang+"_caption"],quote=True)}" width=100%}}\n')
            chunks.append(block(s[lang+"_caption"],"figure-caption"))
        else:
            chunks.append(body+"\n\n")
        links=[]
        for src in scholarly_references(s):
            label,url=PAPERS[src]
            links.append(f"[{label}]({url})")
        if links:
            tag=("Analysis: this project" if lang=="en" else "Análisis: este proyecto") if "figure" in s or phase(s)[0].startswith("04") else ("Scientific references" if lang=="en" else "Referencias científicas")
            chunks.append(block(f'<span class="analysis-tag">{tag}</span>'+" · ".join(links),"source"))
        notes=s[lang+"_notes"]
        notes=re.sub(r"[\w./-]+\.json","reproducibility record" if lang=="en" else "registro de reproducibilidad",notes)
        if s.get("demo_figure"):
            path=s["demo_figure"]+("_es" if lang=="es" else "")+".png"
            label="Full reference/scenario figure" if lang=="en" else "Figura completa referencia/escenario"
            notes+=f"\n\n[{label}]({BASE}assets/figures/{path})."
        chunks.append(block(notes,"notes"))
    path=ROOT/("docs/defense.qmd" if lang=="en" else "docs/es/defense.qmd")
    path.write_text("\n".join(chunks))
    return path

def package():
    files=[ROOT/"docs/defense.html",ROOT/"docs/es/defense.html"]
    if not all(p.exists() for p in files):
        raise RuntimeError("Render both HTML presentations before packaging")
    target=ROOT/"docs/downloads/defense-slides-bilingual.zip"
    with zipfile.ZipFile(target,"w",zipfile.ZIP_DEFLATED) as z:
        for p in files:
            z.write(p,p.relative_to(ROOT/"docs"))
        z.writestr("READ-ME.txt",
            "Open defense.html for English or es/defense.html for Spanish.\n"
            "Each HTML embeds its figures, styling, formulas and presentation runtime.\n"
            "Keep the folder structure for the language switch. Online evidence links require internet.\n"
            "Arrow keys: navigate. Esc: overview. Notes/Notas: detailed explanation.\n")
    print("Portable bilingual package:",target.relative_to(ROOT),target.stat().st_size)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figures",action="store_true")
    parser.add_argument("--render",action="store_true")
    parser.add_argument("--package",action="store_true")
    args=parser.parse_args()
    if args.figures:
        for script in ("build_defense_context_figures.py","build_defense_replication_figures.py","build_defense_spatial_figures.py"):
            subprocess.run([str(ROOT/".venv-analysis/bin/python"),str(ROOT/"scripts"/script)],check=True)
    content()
    ids=[s["id"] for s in SLIDES]
    if len(ids)!=len(set(ids)):
        raise ValueError("Duplicate slide IDs")
    for lang in ("en","es"):
        print("Source:",write_deck(lang).relative_to(ROOT))
    source_files=[]
    for path in sorted(SOURCES):
        p=ROOT/path
        if not p.exists(): raise FileNotFoundError(path)
        source_files.append({"path":path,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()})
    receipt={"presentation_date":DATE,"analysis_snapshot":"Completed project analyses, 18–20 October 2026",
             "scientific_models_refitted":False,"same_source_animals_not_independent_replication":True,
             "slides_per_language":len(SLIDES)+1,"explicit_ids":ids,"title_slide_id":"title-slide",
             "main_content_slides":next(i+1 for i,s in enumerate(SLIDES) if s["id"]=="appendix"),
             "source_files":source_files,"slides":[{k:v for k,v in s.items() if k not in {"en","es","en_notes","es_notes"}} for s in SLIDES],
             "presentation_code_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/"scripts"/name for name in ("build_defense_slides.py","build_defense_context_figures.py","build_defense_replication_figures.py","build_defense_spatial_figures.py")]},
             "source_commit_before_presentation":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()}
    (OUT/"presentation_manifest.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
    print(f"{len(SLIDES)+1} slides per language, including cover and technical appendix")
    if args.render:
        subprocess.run([str(ROOT/".venv-analysis/bin/python"),str(ROOT/"scripts/render_docs.py"),"--page","defense.qmd"],check=True,cwd=ROOT)
    if args.package: package()

if __name__=="__main__": main()

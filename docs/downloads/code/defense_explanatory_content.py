"""Visible reading guides for the bilingual defense, using frozen project facts.

This module edits presentation content only. It neither fits a model nor changes
the scientific figure data. ``augment`` mutates and returns a supplied SLIDES
list and is idempotent, so repeated HTML builds use the same content.
"""


CAPTIONS = {
    "workflow": (
        "Follow the upper branch from exposure data to author reproduction and new study-specific analyses. The lower branch starts with an independent, unexposed MERFISH donor. Its measured anatomy supports reference modeling. The dashed regional bulk-effect transfer enters only the illustrative scenario, so the diagram does not imply measured cellular exposure responses.",
        "La rama superior sigue los datos de exposición hasta la reproducción y los análisis nuevos por estudio. La inferior parte de un donante MERFISH independiente sin exposición. Su anatomía medida sustenta modelos de referencia. La transferencia bulk discontinua entra solo al escenario ilustrativo, sin implicar respuestas celulares medidas a la exposición.",
    ),
    "guarded-cv": (
        "Each row marks one validation fold along the ordered physical sections. Teal sections supply the held-out test, rust sections are excluded neighbors, and pale sections train the model. Every measured cell is tested once. Reserving whole sections reduces local leakage, while all folds still describe the same donor.",
        "Cada fila muestra un bloque de validación sobre cortes físicos ordenados. Los cortes verdeazulados forman la evaluación, los rojizos son vecinos excluidos y los claros entrenan el modelo. Cada célula medida se evalúa una vez. Reservar cortes completos reduce filtración local, aunque todos los bloques siguen perteneciendo al mismo donante.",
    ),
    "cortex-reproduction": (
        "Matching bars recover 109 / 1,205 / 1,222 / 639 significant genes in the four author-aligned contrasts. The effect comparison lies on the identity line. Every significant-set Jaccard equals 1, confirming numerical recovery at the deposited-count level. The same source animals are reused, so this agreement does not add independent biological replication.",
        "Las barras recuperan 109 / 1,205 / 1,222 / 639 genes significativos en los cuatro contrastes originales. Los efectos coinciden con la línea de identidad. Todos los Jaccard son 1: se recuperan las salidas numéricas desde conteos depositados. Se reutilizan los animales fuente, por lo que la concordancia no añade replicación biológica independiente.",
    ),
    "cortex-source-quality": (
        "Points retain all 20 source animals. All four PBS controls have zero raw Ddx3y, Uty and Kdm5d counts, whereas other treatment groups vary. The Xist comparison provides another descriptive marker view. These patterns motivate checking source sex and sample identity, without assigning sex or proving confounding from expression alone.",
        "Los puntos conservan los 20 animales fuente. Los cuatro controles PBS tienen cero conteos Ddx3y, Uty y Kdm5d, mientras los otros grupos varían. La comparación con Xist añade otra vista descriptiva. El patrón motiva verificar sexo e identidad originales, sin asignar sexo ni demostrar confusión solo mediante expresión.",
    ),
    "hy-reproduction": (
        "The original lists contain 17 PM and 761 HFD genes. The declared TMM/voom sensitivity returns 0 and 3,600 under different specified gates. The primary seven-animal analysis and its R implementation check are shown separately. Missing RUVg specifications prevent exact reconstruction, so these counts describe methodological sensitivity rather than biological nonreplication.",
        "Las listas originales contienen 17 genes PM y 761 HFD. La sensibilidad TMM/voom devuelve 0 y 3,600 bajo otros criterios declarados. Se muestran aparte el análisis principal de siete animales y su comprobación en R. Faltan especificaciones RUVg para reconstruir exactamente: los conteos describen sensibilidad metodológica, sin demostrar ausencia de replicación biológica.",
    ),
    "mirna-reproduction": (
        "Animal points show lower miR-574-5p signal at high dose. Its estimated difference is −0.26945 deposited signal units. The nominal limma interval excludes zero, but p = 0.0253 becomes BH q = 0.312 across the array. ANOVA p = 0.0788 also prevents protected LSD. Directional agreement therefore exceeds the supported multiplicity-adjusted claim.",
        "Los puntos muestran menor señal miR-574-5p a dosis alta, con diferencia de −0.26945 unidades depositadas. El intervalo limma nominal excluye cero, pero p = 0.0253 pasa a BH q = 0.312 en el array. ANOVA p = 0.0788 impide LSD protegido. La coincidencia direccional no alcanza confirmación con ajuste de multiplicidad.",
    ),
    "ozone-reproduction": (
        "The deposited genotype threshold recovers 2,501 genes exactly. Across common finite effects, correlation is ρ = 0.865. The author-selected S3 subset gives ρ = 0.962 and 98.8% directional agreement among matched effects. Because the subset was selected and the source animals are reused, its stronger agreement cannot establish unbiased transcriptome-wide replication.",
        "El umbral depositado de genotipo recupera exactamente 2,501 genes. La correlación de efectos finitos comunes es ρ = 0.865. La lista S3 seleccionada alcanza ρ = 0.962 y 98.8% de concordancia direccional entre efectos cotejados. Su selección previa y los mismos animales impiden interpretar esa mayor concordancia como replicación transcriptómica global no sesgada.",
    ),
    "hy-primary": (
        "For 4 filtered-air and 3 PM animals, the primary model identifies 11 BH discoveries: two higher and nine lower. Points show MLE log₂ fold changes, and bars are pointwise Wald 95% intervals. Direction indicates the regional contrast. These intervals are neither simultaneous confidence bands nor evidence that each spatial cell changes similarly.",
        "Con 4 animales de aire filtrado y 3 PM, el modelo identifica 11 hallazgos BH: dos mayores y nueve menores. Los puntos son log₂ fold changes MLE y las barras, intervalos Wald puntuales del 95%. La dirección corresponde al contraste regional, sin intervalos simultáneos ni evidencia de un cambio equivalente en cada célula espacial.",
    ),
    "aca-primary": (
        "The volcano separates effect magnitude from nominal Wald evidence, while colors identify BH discoveries. The seven smallest-q genes have pointwise MLE intervals. The eight-animal primary comparison yields 1,167 discoveries, 658 higher and 509 lower. Its population and model differ from the 20-animal reconstruction that recovered 1,205 PM genes.",
        "El volcán separa magnitud del efecto y evidencia Wald nominal. Los colores identifican hallazgos BH. Los siete genes con menor q tienen intervalos MLE puntuales. El contraste principal de ocho animales produce 1,167 hallazgos: 658 mayores y 509 menores. Su población y modelo difieren de la reconstrucción de 20 animales con 1,205 genes PM.",
    ),
    "ozone-primary": (
        "Separate within-stratum BH procedures yield 0 / 399 / 2 / 3 genes in starting-age/APOE order 3M/E3, 3M/E4, 17M/E3 and 17M/E4. The joint 4-df exposure LRT finds 24 genes in another testing family. Unequal discovery counts reflect evidence and precision together, and cannot by themselves establish an interaction.",
        "BH por estrato produce 0 / 399 / 2 / 3 genes en el orden edad inicial/APOE 3M/E3, 3M/E4, 17M/E3 y 17M/E4. El LRT conjunto de exposición con 4 grados de libertad encuentra 24 genes en otra familia de pruebas. Los conteos reflejan evidencia y precisión, sin establecer por sí solos una interacción.",
    ),
    "interaction-animals": (
        "Each gene retains all 24 sequenced animals, grouped by starting age, APOE and exposure. Points show log₂(1 + normalized count), and short segments give descriptive means. Ryr3 and Rps24 belong to the eight age-dependent LRT discoveries. Fixed horizontal offsets only separate symbols, while the factorial model supplies the interaction test.",
        "Cada gen conserva los 24 animales secuenciados, agrupados por edad inicial, APOE y exposición. Los puntos representan log₂(1 + conteo normalizado) y los segmentos, medias descriptivas. Ryr3 y Rps24 pertenecen a los ocho hallazgos LRT dependientes de edad. El desplazamiento horizontal separa símbolos, mientras el modelo factorial aporta la prueba de interacción.",
    ),
    "go-patterns": (
        "Color gives the signed normalized enrichment score on a common scale, and dots mark algorithmic GSEA FDR < 0.05. Gray cells are ineligible or unavailable. These are the first eight terms of a fixed 22-term display selection. Correlated annotations and selected display terms prevent reading the matrix as independent mechanisms or exhaustive pathway activity.",
        "El color muestra NES con signo en una escala común y los puntos indican FDR algorítmica GSEA < 0.05. Gris significa no elegible o no disponible. Son los primeros ocho términos de una selección gráfica de 22. La correlación y selección de anotaciones impiden leer la matriz como mecanismos independientes o actividad exhaustiva de vías.",
    ),
    "aca-sections": (
        "The largest ACA section contains 1,940 measured cells. Left, deposited positions and anatomical labels reveal cortical layers. Right, each occupied 100-µm bin shows measured Fosb-positive cells divided by all measured cells, including zeros. Empty bins remain unobserved. The pattern describes one unexposed donor and does not estimate a layer-specific PM effect.",
        "El mayor corte ACA contiene 1,940 células medidas. A la izquierda, posiciones depositadas y etiquetas anatómicas muestran capas corticales. A la derecha, cada casilla ocupada de 100 µm divide células Fosb positivas por todas las medidas, incluidos ceros. Las casillas vacías siguen sin observarse. Describe un donante sin exposición, sin estimar efectos PM por capa.",
    ),
    "hy-sections": (
        "The largest HY section contains 3,565 measured cells. Anatomical labels organize the deposited cell positions, and occupied bins display positive/all-cell detection. Across the complete reference HY, Tll1 is detected in 667 of 37,811 cells, or 1.764%. Sparse detection can reflect biology and assay sensitivity, without resolving either explanation from these data alone.",
        "El mayor corte HY contiene 3,565 células medidas. Las etiquetas anatómicas organizan posiciones depositadas y las casillas ocupadas muestran positivas/todas. En HY completo, Tll1 se detecta en 667 de 37,811 células: 1.764%. La detección escasa puede reflejar biología y sensibilidad del ensayo, sin distinguir ambas explicaciones con estos datos solos.",
    ),
    "gaussian-rmse": (
        "Each point compares spatial-model RMSE with anatomy-only RMSE on identical held-out sections. A ratio below 1 favors the spatial model, which improves 4/10 ACA genes and 1/10 HY genes. Added coordinates and basis functions therefore provide limited benefit here. A smooth full-reference map alone would not reveal this prediction failure.",
        "Cada punto compara RMSE del modelo espacial con RMSE anatómico en los mismos cortes reservados. Una razón menor que 1 favorece al modelo espacial: mejora 4/10 genes ACA y 1/10 HY. Coordenadas y bases aportan beneficio limitado. Un mapa suavizado ajustado con toda la referencia no mostraría por sí solo esta limitación predictiva.",
    ),
    "gaussian-diagnostics": (
        "Left, overall interval inclusion can be high while inclusion among positive observations is poor. Positive-only inclusion is a conditional diagnostic, without a nominal 95% target. Right, the Gaussian predictive distribution assigns mass below zero on log(1 + count), an impossible range. Together these checks expose limitations hidden by marginal coverage alone.",
        "A la izquierda, la inclusión global puede ser alta con mala inclusión entre observaciones positivas. La inclusión solo en positivas es un diagnóstico condicional, sin objetivo nominal del 95%. A la derecha, la predicción gaussiana asigna masa bajo cero en log(1 + conteo), un rango imposible. Ambos diagnósticos revelan límites ocultos por la cobertura marginal.",
    ),
    "beta-performance": (
        "Lower Brier loss favors the model. ACA improves modestly, 0.13725 versus 0.14051, while HY worsens, 0.09246 versus 0.03621. The coverage–width plot shows that high predictive count coverage accompanies wide intervals. These held-out comparisons retain prior-only anatomical groups, so bounded probabilities alone do not establish accurate or transferable prediction.",
        "Menor Brier favorece al modelo. ACA mejora modestamente: 0.13725 frente a 0.14051. HY empeora: 0.09246 frente a 0.03621. El gráfico cobertura–amplitud muestra cobertura predictiva alta con intervalos amplios. Se conservan grupos anatómicos sostenidos solo por el prior, de modo que acotar probabilidades no establece precisión ni transferencia predictiva.",
    ),
    "beta-aca": (
        "For Fosb, the left panel compares observed fractions with known-section posterior probabilities. The right panel separates anatomy-mean credible intervals from new-section propensity predictions. The pale interval expresses section-to-section uncertainty, and its target differs from the dark hypermean interval. These full-reference estimates describe the unexposed donor, while held-out performance is evaluated separately.",
        "Para Fosb, el panel izquierdo compara fracciones observadas y probabilidades posteriores del corte conocido. El derecho separa intervalos creíbles de la media anatómica y predicciones de propensión de un corte nuevo. El intervalo claro expresa variación entre cortes y estima otro objetivo que el oscuro. Describen un donante sin exposición. La validación reserva cortes.",
    ),
    "beta-hy": (
        "For Tll1, fitted section probabilities remain within 0–1 despite sparse detections. The anatomy mean and new-section propensity carry different uncertainty. A fitted map can still look plausible when held-out anatomical labels lack training support. The full HY evaluation retains those groups, making the resulting prior dependence visible in the reported prediction error.",
        "Para Tll1, las probabilidades ajustadas permanecen entre 0 y 1 pese a detecciones escasas. La media anatómica y la propensión de corte nuevo tienen incertidumbres distintas. Un mapa puede parecer plausible aunque falte anatomía en entrenamiento. La evaluación HY conserva esos grupos y hace visible la dependencia del prior mediante el error predictivo.",
    ),
}


def lesson(sid, en_title, es_title, en, es, en_notes, es_notes, sources):
    return dict(
        id=sid,
        en_title=en_title,
        es_title=es_title,
        en=en.strip(),
        es=es.strip(),
        en_notes=en_notes.strip(),
        es_notes=es_notes.strip(),
        sources=list(sources),
        cls="explanation-slide",
    )


INSERTIONS = [
    ("units", lesson(
        "explain-bulk-spatial",
        "Bulk contrasts and spatial cells answer different questions",
        "El contraste bulk y las células espaciales responden preguntas distintas",
        """**Exposure comparison.** Bulk RNA-seq combines RNA from many cells in each sampled region. The animal supplies the independent replicate. A regional difference can reflect expression within cells, changes in cellular composition, or both.

**Spatial reference.** MERFISH counts a targeted gene panel at measured cell positions in one unexposed donor. It shows where a transcript is detected in that reference.

**Connection.** We map regional candidates onto the measured panel. The atlas supplies anatomical context, while the separate exposure experiment supplies the regional contrast.""",
        """**Comparación de exposición.** RNA-seq bulk reúne ARN de muchas células de cada región. El animal aporta la réplica independiente. Una diferencia regional puede reflejar expresión dentro de células, cambios de composición o ambos.

**Referencia espacial.** MERFISH cuenta un panel dirigido en posiciones celulares medidas de un donante sin exposición. Muestra dónde se detecta un transcrito en esa referencia.

**Conexión.** Vinculamos los candidatos regionales con el panel medido. El atlas aporta contexto anatómico y el experimento de exposición independiente aporta el contraste regional.""",
        "The atlas and exposure studies have different animals, assays and estimands. Spatial cells cannot be assigned the pollution label of an animal from another experiment. Increasing the number of reference cells improves description within the donor but does not increase the number of exposed biological replicates. The measured panel also limits which regional candidates can appear on a spatial map. Missing panel genes remain missing rather than receiving invented coordinates or expression.",
        "Atlas y experimentos de exposición tienen animales, ensayos y estimandos distintos. Las células espaciales no reciben la etiqueta de contaminación de un animal de otro experimento. Más células de referencia mejoran la descripción dentro del donante, pero no aumentan las réplicas biológicas expuestas. El panel medido limita qué candidatos regionales aparecen en mapas. Los genes ausentes siguen sin medición, sin inventarles coordenadas o expresión.",
        ["HY", "ACA", "atlas"],
    )),
    ("workflow", lesson(
        "explain-gene-trace",
        "Fosb across the experimental and spatial analyses",
        "Fosb a través de los análisis experimental y espacial",
        """**Regional evidence.** The eight-animal ACA comparison estimates higher Fosb RNA after PM exposure. Its shrunken log₂ fold change is +0.99711, corresponding to a 1.996-fold multiplier.

**Measured reference.** The separate atlas detects Fosb in 1,418 of 14,582 ACA cells. Original coordinates place these detections within the registered anatomy.

**Illustrative display.** Multiplying every reference value by 1.996 changes the color and magnitude while preserving positions and measured zeros. This assumption illustrates the regional estimate without identifying the cells that actually respond to PM.""",
        """**Evidencia regional.** El contraste ACA de ocho animales estima mayor ARN Fosb con PM. Su log₂ fold change regularizado es +0.99711, equivalente a un multiplicador de 1.996.

**Referencia medida.** El atlas independiente detecta Fosb en 1,418 de 14,582 células ACA. Las coordenadas originales ubican las detecciones en la anatomía registrada.

**Visualización ilustrativa.** Multiplicar cada valor de referencia por 1.996 cambia color y magnitud, conservando posiciones y ceros medidos. Este supuesto ilustra la estimación regional sin identificar qué células responden realmente a PM.""",
        "Fosb is a concrete example linking the project's primary count analysis, measured panel and viewer. The primary model uses four PBS and four PM animals. The atlas denominator is all measured ACA cells, including measured zeros, and comes from one unexposed donor. A positive bulk coefficient does not imply activation in every Fosb-positive cell. Protein abundance, AP-1 activity, isoform-specific FosB versus ΔFosB effects and causal mechanisms are outside what the gene-level RNA measurement establishes.",
        "Fosb conecta el análisis principal de conteos, el panel medido y el visor. El modelo utiliza cuatro animales PBS y cuatro PM. El denominador del atlas incluye todas las células ACA medidas, incluidos ceros, de un donante sin exposición. Un coeficiente bulk positivo no implica activación en cada célula Fosb positiva. El ARN a nivel génico no establece abundancia proteica, actividad AP-1, efectos específicos FosB/ΔFosB ni mecanismos causales.",
        ["ACA", "atlas", "DESeq2"],
    )),
    ("reproduction-criteria", lesson(
        "explain-reproduction",
        "Numerical recovery and independent biological validation",
        "Recuperación numérica y validación biológica independiente",
        """**Numerical recovery** asks whether the same deposited inputs and analysis rules recover the authors' effects, probabilities and selected genes. Cortex recovers all four tested significant-gene sets exactly.

**Methodological sensitivity** asks what changes when a model, filter or analytical population changes. Hypothalamic reconstruction remains incomplete because key RUVg specifications are unavailable.

**Independent biological validation** would require additional animals or experiments. Reanalyzing the source animals, even with another software package, checks computation and assumptions while leaving the biological sample size unchanged.""",
        """**Recuperación numérica** pregunta si las mismas entradas y reglas recuperan efectos, probabilidades y genes seleccionados por los autores. Corteza recupera exactamente los cuatro conjuntos significativos cotejados.

**Sensibilidad metodológica** pregunta qué cambia al modificar modelo, filtro o población analítica. La reconstrucción hipotalámica queda incompleta porque faltan especificaciones RUVg.

**Validación biológica independiente** exigiría animales o experimentos adicionales. Reanalizar los animales fuente, incluso con otro software, comprueba cálculos y supuestos sin modificar el tamaño muestral biológico.""",
        "Identical gene counts alone are a weak reproduction criterion because different sets can have the same size. The cortex comparison also checks identities, effect estimates, p values and selected-set Jaccard. This is strong evidence for recovery of the checked count-level output, not validation of every experiment in the source article. Source metadata concerns remain relevant even when calculations match. Transparent nonagreement and explicit missing specifications are scientific results of a reproduction audit rather than reasons to tune an analysis until it matches.",
        "El mismo número de genes constituye un criterio débil porque conjuntos distintos pueden tener igual tamaño. En corteza también se cotejan identidades, efectos, valores p y Jaccard. Es evidencia fuerte de recuperación de salidas de conteos examinadas, sin validar cada experimento del artículo. Las preguntas de metadatos siguen siendo pertinentes aunque coincidan los cálculos. Las discrepancias y especificaciones faltantes se informan, sin ajustar parámetros hasta forzar coincidencia.",
        ["ACA", "HY"],
    )),
    ("count-model", lesson(
        "explain-count-model",
        "Normalization and dispersion in the count model",
        "Normalización y dispersión en el modelo de conteos",
        """**Normalization** accounts for technical differences in the amount of RNA sequenced. Otherwise, a deeply sequenced library could appear higher for many genes even without a biological change. Ozone also uses effective-length information from transcript quantification.

**Dispersion** describes variation among animals beyond ordinary count noise. A variable gene needs a larger, more consistent group difference to supply strong evidence.

**The coefficient** estimates the chosen regional comparison after these adjustments. Neither normalization nor dispersion estimation creates biological replicates or removes an unmeasured source of confounding.""",
        """**Normalización** ajusta diferencias técnicas en la cantidad de ARN secuenciado. Sin ella, una biblioteca profunda podría parecer mayor en muchos genes sin cambio biológico. Ozono incorpora además longitudes efectivas de la cuantificación transcriptómica.

**Dispersión** describe variación entre animales adicional al ruido ordinario de conteo. Un gen variable necesita una diferencia grupal mayor y más consistente para aportar evidencia fuerte.

**El coeficiente** estima el contraste regional elegido tras esos ajustes. Normalizar y estimar dispersión no crea réplicas biológicas ni elimina una fuente de confusión no medida.""",
        "DESeq2 and PyDESeq2 use negative-binomial count models with gene-specific dispersion estimation and information sharing. For the simple PM comparisons the normalization factor reduces to sample-specific size factors. The tximport ozone workflow adds gene/sample normalization from effective lengths, so a single library-size description would be incomplete. These count procedures do not apply to the signed, deposited miRNA continuous signals, which use limma. Low counts generally limit precision, and a retained zero remains an observation rather than a missing value filled by normalization.",
        "DESeq2 y PyDESeq2 emplean modelos binomiales negativos con dispersión por gen e información compartida. En los contrastes PM simples, los factores de normalización se reducen a factores de tamaño de biblioteca específicos de cada muestra. El flujo tximport de ozono incorpora longitudes efectivas por gen y muestra. Estas reglas de conteos no se aplican a las señales continuas miRNA depositadas, analizadas con limma. Conteos bajos suelen limitar precisión y un cero retenido sigue siendo una observación, sin rellenarse por normalización.",
        ["DESeq2", "PyDESeq2", "tximport"],
    )),
    ("testing", lesson(
        "explain-testing",
        "Effect size, uncertainty and the discovery threshold",
        "Magnitud del efecto, incertidumbre y umbral de hallazgo",
        """**Effect size** describes direction and magnitude. A log₂ fold change of +1 corresponds to twice the regional RNA level, while −1 corresponds to half, under the fitted contrast.

**The standard error and Wald test** compare the estimated effect with its sampling uncertainty. A large fold change from inconsistent animals may have weak evidence.

**BH adjustment** evaluates many genes as a declared testing family. A nominal p < 0.05 can become q ≥ 0.05. Shrinkage stabilizes displayed effect estimates, while the stated inferential model determines the reported tests.""",
        """**La magnitud del efecto** expresa dirección y tamaño. Un log₂ fold change de +1 corresponde al doble de ARN regional y −1, a la mitad, bajo el contraste ajustado.

**El error estándar y Wald** comparan el efecto estimado con su incertidumbre de muestreo. Un cambio grande entre animales inconsistentes puede aportar evidencia débil.

**El ajuste BH** evalúa numerosos genes como una familia declarada. Un p nominal < 0.05 puede pasar a q ≥ 0.05. La regularización estabiliza efectos mostrados y el modelo inferencial declarado determina las pruebas.""",
        "BH q values should not be described as the posterior probability that an individual gene is false. Under its assumptions, the procedure controls the expected false-discovery proportion among rejections at the chosen level. A pointwise 95% Wald interval is an approximate interval for one model coefficient, not a simultaneous interval for all genes or a per-cell interval. Different filtering populations and multiplicity families can change the discovery list even when underlying effects are similar. The MLE and shrunken effects are retained separately throughout the project.",
        "Un q BH no es la probabilidad posterior de que un gen individual sea falso. Bajo sus supuestos, el procedimiento controla la proporción esperada de falsos descubrimientos entre rechazos al nivel elegido. Un intervalo Wald puntual del 95% aproxima un coeficiente del modelo, sin simultaneidad para todos los genes ni interpretación por célula. Filtros y familias de multiplicidad diferentes pueden cambiar las listas pese a efectos similares. El proyecto conserva por separado los efectos MLE y regularizados.",
        ["DESeq2", "apeglm"],
    )),
    ("ozone-lrt", lesson(
        "explain-interaction",
        "Interaction tests compare exposure effects directly",
        "Las pruebas de interacción comparan directamente efectos de exposición",
        """**Within a stratum**, the model compares ozone with filtered air at one starting age and APOE genotype. Each stratum has only three animals per exposure group.

**An interaction** asks whether the ozone effect itself differs between ages or genotypes, using coefficient differences and their covariance. One significant stratum and one nonsignificant stratum do not answer that question.

**Here**, the joint age-dependent exposure LRT identifies eight genes under its own BH family, and one under BH pooled across the three joint heterogeneity tests. The conclusion depends on the stated test and family.""",
        """**Dentro de un estrato**, el modelo compara ozono y aire filtrado a una edad inicial y genotipo APOE. Cada estrato tiene solo tres animales por exposición.

**Una interacción** pregunta si el efecto del ozono cambia entre edades o genotipos, usando diferencias de coeficientes y su covarianza. Un estrato significativo y otro no significativo no responden esa pregunta.

**Aquí**, el LRT conjunto de respuesta dependiente de edad identifica ocho genes con su propia familia BH y uno al agrupar las tres pruebas conjuntas de heterogeneidad. La conclusión depende de la prueba y familia declaradas.""",
        "The 24 RNA-seq samples form the eight age×genotype×exposure combinations with three animals each. The full factorial design has eight parameters and estimates all groups jointly. The four-degree-of-freedom exposure LRT and the two-degree-of-freedom age-dependent exposure LRT test different null hypotheses. Signed individual interaction contrasts also differ from the joint LRTs. Starting ages of three and seventeen months are source design labels at exposure onset. These labels should not be relabeled as terminal sampling ages.",
        "Las 24 muestras RNA-seq forman ocho combinaciones edad×genotipo×exposición, con tres animales cada una. El diseño factorial completo tiene ocho parámetros y estima conjuntamente todos los grupos. El LRT global de exposición de cuatro grados de libertad y el LRT de respuesta por edad de dos contrastan hipótesis nulas distintas. Los contrastes individuales con signo también difieren de los LRT conjuntos. Tres y diecisiete meses son edades al inicio de exposición, sin sustituirlas por edades terminales.",
        ["O3", "DESeq2"],
    )),
    ("gsea", lesson(
        "explain-enrichment",
        "ORA and GSEA use different evidence from the same contrast",
        "ORA y GSEA utilizan evidencia distinta del mismo contraste",
        """**ORA** starts with genes passing a declared discovery gate and asks whether an annotation contains more selected genes than expected from the tested background. Its answer depends on the foreground and universe.

**GSEA** uses the signed gene ranking and can detect coordinated shifts even when few individual genes pass BH. The analysis uses the frozen ranked-test procedure and algorithmic FDR.

**Interpretation** concerns enrichment of annotated genes. A positive score follows the positive rank direction, while shared genes make related GO terms correlated. RNA enrichment alone does not measure pathway activity.""",
        """**ORA** parte de genes que superan un criterio declarado y pregunta si una anotación contiene más seleccionados de lo esperado en el fondo evaluado. Depende de la selección y del universo.

**GSEA** utiliza el rango génico con signo y detecta desplazamientos coordinados aunque pocos genes superen BH. Se aplica el procedimiento congelado por rangos y su FDR algorítmica.

**Interpretación** se refiere a enriquecimiento de genes anotados. Un valor positivo sigue la dirección positiva del rango. Los genes compartidos correlacionan términos GO y el enriquecimiento de ARN por sí solo no mide actividad de vías.""",
        "The project's ORA universe consists of mapped genes with finite adjusted p values and retains genes lacking a propagated biological-process annotation in the background. GSEA retains finite original p values and signed Wald statistics, including finite-p genes excluded by independent filtering from adjusted-p reporting. These universes and null procedures differ. Algorithmic preranked GSEA FDR is neither BH q nor a phenotype-permutation result. No enrichment result here establishes causal regulation, protein activity, spatial connectivity or an independent mechanism for each correlated GO term.",
        "El universo ORA incluye genes mapeados con p ajustado finito y conserva en el fondo genes sin anotación propagada de proceso biológico. GSEA retiene p originales y Wald con signo finitos, incluidos genes que no reciben p ajustado por filtrado independiente. Los universos y procedimientos nulos difieren. La FDR algorítmica GSEA por rangos no equivale a q BH ni a permutación de fenotipos. El enriquecimiento no establece regulación causal, actividad proteica, conectividad espacial ni un mecanismo independiente por término GO correlacionado.",
        ["GO", "GSEA"],
    )),
    ("scenario-method", lesson(
        "explain-scenario",
        "The assumption behind the reference–PM comparison",
        "El supuesto detrás de la comparación referencia–PM",
        """**Measured input.** Each cell has its original coordinate and reference expression value. These measurements come from an unexposed animal.

**Scenario rule.** Apply one gene-specific regional multiplier to all reference values. Fosb receives approximately ×2 in ACA. Tll1 receives ×0.428 in HY. Measured zeros stay zero and no cells move, appear or disappear.

**Scientific meaning.** This makes a regional estimate visible in anatomical space. Many cellular response patterns could produce the same bulk effect, so the rule cannot recover cell-specific pollution responses or add spatial resolution.""",
        """**Entrada medida.** Cada célula conserva su coordenada y expresión de referencia originales. Las mediciones proceden de un animal sin exposición.

**Regla del escenario.** Aplicar un multiplicador regional por gen a todos los valores. Fosb recibe aproximadamente ×2 en ACA y Tll1, ×0.428 en HY. Los ceros siguen en cero y ninguna célula cambia de posición, aparece o desaparece.

**Significado científico.** Permite visualizar una estimación regional en anatomía. Múltiples respuestas celulares pueden producir el mismo efecto bulk. La regla no recupera respuestas celulares específicas ni añade resolución espacial.""",
        "The scenario uses the shrunken bulk coefficient as its default multiplier. Alternative group-level scenarios can use the pointwise MLE Wald coefficient limits, but those are not confidence intervals for each atlas cell. Uniform scaling preserves the reference's detection pattern exactly, so it cannot predict newly expressing or silenced cells. The atlas donor is not an experimental control matched to the pollution animals. Shared per-gene color scales and camera settings are necessary for comparison, and the captions identify the second view as a modeled illustration.",
        "El escenario utiliza por defecto el coeficiente bulk regularizado. Los límites Wald MLE puntuales permiten escenarios grupales alternativos, pero no son intervalos de confianza por célula del atlas. El escalamiento uniforme conserva exactamente el patrón de detección de referencia y no predice células que empiezan o dejan de expresar. El donante no es un control emparejado de los animales expuestos. La comparación exige la misma escala por gen y cámara, e identifica la segunda vista como ilustración modelada.",
        ["ACA", "HY", "atlas"],
    )),
    ("gaussian-priors", lesson(
        "explain-bayes",
        "How the prior and reference measurements form a posterior",
        "Cómo el prior y las mediciones forman una posterior",
        """**Before the observations**, a proper prior constrains coefficient magnitudes and unexplained variation. It regularizes sparse anatomical groups and spatial basis functions without supplying new measurements.

**The likelihood** describes how measured reference expression would vary for a given coefficient set. Combining it with the prior gives a posterior distribution of plausible reference parameters.

**Prediction** then propagates that distribution to a defined target. A smooth posterior mean is only a conditional estimate. Its scientific usefulness still requires held-out comparisons, interval checks and sensitivity to the stated assumptions.""",
        """**Antes de las observaciones**, un prior propio restringe magnitudes de coeficientes y variación no explicada. Regulariza grupos anatómicos escasos y bases espaciales sin aportar mediciones nuevas.

**La verosimilitud** describe cómo varía la expresión medida de referencia para unos coeficientes dados. Al combinarla con el prior obtenemos una distribución posterior de parámetros plausibles.

**La predicción** propaga esa distribución a un objetivo definido. Una media posterior suave es una estimación condicional. Su utilidad científica exige comparaciones reservadas, revisión de intervalos y sensibilidad a los supuestos declarados.""",
        "The Gaussian implementation is finite-basis regression with fixed length scales and a conjugate normal–inverse-gamma prior. It is not a fitted full Gaussian process. Prior settings were fixed before held-out evaluation, and numerical sensitivity analyses do not select a setting to maximize test performance. Where anatomy is absent from training, the model must fall back to its stated prior handling. The separate beta-binomial model changes both the observation target and distribution: it models positive/all-cell detection probabilities rather than RNA magnitude.",
        "La implementación gaussiana es regresión con bases finitas, longitudes fijas y prior conjugado normal–gamma inversa. No es un proceso gaussiano completo ajustado. Los priors se fijaron antes de evaluar y la sensibilidad numérica no selecciona el ajuste con mejor rendimiento reservado. Una anatomía ausente de entrenamiento depende del tratamiento a priori declarado. El modelo beta-binomial cambia objetivo y distribución: representa probabilidades positivas/todas las células en lugar de magnitud de ARN.",
        ["Murphy", "GPML", "Beta"],
    )),
    ("guarded-cv", lesson(
        "explain-validation",
        "Why validation reserves sections and their neighbors",
        "Por qué la validación reserva cortes y sus vecinos",
        """**The prediction question** is how the model performs at a physical section it has not seen. Randomly splitting nearby cells could leave very similar local information in training and testing.

**The guard** removes the represented sections adjacent to the test block from training. Feature standardization, anatomical vocabularies and spatial basis preparation use the remaining training data only.

**The comparison** uses identical folds for spatial and simpler baseline models. Better fit to the complete atlas cannot replace better held-out prediction. These folds assess transfer within one donor, while cross-donor and exposure transfer remain untested.""",
        """**La pregunta predictiva** es cómo funciona el modelo en un corte físico no observado durante el ajuste. Dividir células vecinas al azar puede dejar información local muy similar en entrenamiento y evaluación.

**La separación** excluye del entrenamiento los cortes representados adyacentes al bloque evaluado. Estandarización, vocabulario anatómico y bases espaciales usan únicamente el entrenamiento restante.

**La comparación** utiliza los mismos bloques para modelos espaciales y referencias simples. Ajustar mejor el atlas completo no sustituye predecir mejor datos reservados. Se evalúa transferencia dentro de un donante y quedan pendientes otros donantes y exposición.""",
        "The five test blocks contain contiguous whole physical sections ordered by acquisition position. Adjacent represented sections are guards, which does not imply filling an unobserved intermediate section. Every measured cell appears in a test fold once. Unknown test anatomy remains visible in scoring, with the specified prior-only handling rather than post hoc exclusion. Multiple genes evaluated in the same cell are correlated evaluations. Cross-validation protects against selected routes of leakage but cannot create independent donors or validate the uniform pollutant scenario.",
        "Los cinco bloques reservan cortes físicos completos contiguos, ordenados por posición de adquisición. Se excluyen cortes representados vecinos sin rellenar cortes intermedios no observados. Cada célula aparece una vez en evaluación. Las anatomías desconocidas se conservan al puntuar con el tratamiento solo-prior declarado, sin excluirlas a posteriori. Evaluar varios genes en una célula genera observaciones correlacionadas. La validación reduce ciertas filtraciones, sin crear donantes independientes ni validar el escenario uniforme de contaminación.",
        ["atlas", "Murphy", "Beta"],
    )),
    ("beta-targets", lesson(
        "explain-uncertainty",
        "Different intervals answer different detection questions",
        "Distintos intervalos responden distintas preguntas de detección",
        """**A known section:** after observing its positives and all-cell denominator, estimate that section's underlying detection probability. Its credible interval conditions on those observations.

**An anatomical mean:** combine the exchangeable section information to estimate the typical propensity for that anatomy. This hypermean has its own credible interval.

**A new section:** predict its propensity before seeing its counts, or predict a positive-cell count at a specified denominator. Section heterogeneity and, for counts, sampling variation enter this prediction. Always read the target before comparing interval widths or coverage.""",
        """**Un corte conocido:** tras observar sus positivas y el denominador total, estimamos la probabilidad subyacente de ese corte. Su intervalo creíble está condicionado por esas observaciones.

**Una media anatómica:** combinamos información de cortes intercambiables para estimar la propensión típica de esa anatomía. Esta hipermedia tiene su propio intervalo creíble.

**Un corte nuevo:** predecimos su propensión antes de ver sus conteos, o el número de positivas para un denominador definido. La heterogeneidad entre cortes y, en conteos, la variación de muestreo entran en la predicción. Hay que identificar el objetivo antes de comparar amplitud o cobertura.""",
        "The anatomy hypermean is the mean of the exchangeable section propensities rather than necessarily the pooled-cell detection fraction. The probability of a known section and the propensity of a new section both lie in [0,1], but the new positive count has integer support from zero to the chosen denominator. Credible intervals for fitted quantities are not held-out count predictions. The observed group fraction is also noisy, especially when its denominator is small. There is no universal ordering that requires every new-section propensity interval to be wider than every hypermean interval.",
        "La hipermedia anatómica resume propensiones de cortes intercambiables y no necesariamente coincide con la fracción celular agrupada. Probabilidades de corte conocido y nuevo pertenecen a [0,1], mientras el conteo nuevo tiene soporte entero entre cero y el denominador elegido. Un intervalo creíble ajustado no es una predicción de conteo reservado. La fracción observada también tiene ruido, sobre todo con pocos datos. No existe un orden universal que exija que cada intervalo de propensión nueva sea más ancho que cada intervalo de hipermedia.",
        ["atlas", "Beta"],
    )),
    ("beta-performance", lesson(
        "explain-model-failure",
        "Prediction quality depends on error, coverage and support",
        "La calidad predictiva depende de error, cobertura y soporte",
        """**Prediction error.** The bounded model has slightly lower Brier loss than the regional baseline in ACA. In HY, its full evaluation is worse than the baseline.

**Training support.** In each HY gene evaluation, 79 of 186 reserved section–anatomy groups contain anatomical labels unavailable in training. Their predictions depend only on the stated prior and remain in the score.

**Uncertainty.** Wide intervals can cover many outcomes while providing little precision. Bounded probabilities, stable numerical integration and visually plausible maps do not by themselves demonstrate useful prediction or higher spatial resolution.""",
        """**Error predictivo.** El modelo acotado tiene Brier ligeramente menor que la referencia regional en ACA. En HY, la evaluación completa es peor que la referencia.

**Soporte del entrenamiento.** Por cada gen HY, 79 de 186 grupos corte–anatomía reservados tienen etiquetas anatómicas ausentes del entrenamiento. Sus predicciones dependen solo del prior declarado y permanecen en la puntuación.

**Incertidumbre.** Intervalos amplios pueden cubrir numerosos resultados con poca precisión. Acotar probabilidades, integrar de forma estable y mostrar mapas plausibles no demuestra por sí solo predicción útil ni mayor resolución espacial.""",
        "Brier loss is the mean squared error of predicted detection probabilities against cell detection outcomes, here using the specified cell weighting. Lower values are better only for the same evaluation target and cells. The 79/186 figures are unique HY section×anatomy groups per gene, and 10,416 of 37,811 actual reference cells fall in prior-only groups. Repeating evaluation across ten genes does not create new cells. The restricted trained-anatomy HY comparison improves slightly, but excluding unsupported groups would answer a narrower question and cannot replace the full reported failure.",
        "Brier representa error cuadrático medio de probabilidades de detección frente a detección celular, con la ponderación declarada. Menor valor favorece la predicción del mismo objetivo y células. Los 79/186 son grupos HY únicos corte×anatomía por gen y 10,416 de 37,811 células reales quedan en grupos solo-prior. Repetir diez genes no crea células. La comparación HY restringida a anatomía entrenada mejora ligeramente, pero excluir grupos sin soporte responde otra pregunta y no sustituye el fallo completo informado.",
        ["atlas", "Murphy", "Beta"],
    )),
]


def augment(slides):
    """Add paired explanation slides and replace the original19 reading guides.

    Existing figure notes and all source files are preserved. Unknown new figures
    from another builder are left intact. Missing expected anchors raise an error
    rather than quietly placing an explanation in an unrelated section.
    """
    ids = {s["id"] for s in slides}
    missing = {anchor for anchor, _ in INSERTIONS} - ids
    missing |= set(CAPTIONS) - ids
    if missing:
        raise ValueError("Missing explanation anchor(s): " + ", ".join(sorted(missing)))
    for s in slides:
        if s["id"] in CAPTIONS:
            s["en_caption"], s["es_caption"] = CAPTIONS[s["id"]]
            classes = s.get("cls", "").split()
            if "figure-interpretation" not in classes:
                classes.append("figure-interpretation")
            s["cls"] = " ".join(classes)
    after = {anchor: item for anchor, item in INSERTIONS}
    result = []
    for s in slides:
        result.append(s)
        if s["id"] in after and after[s["id"]]["id"] not in ids:
            result.append(dict(after[s["id"]]))
    slides[:] = result
    return slides

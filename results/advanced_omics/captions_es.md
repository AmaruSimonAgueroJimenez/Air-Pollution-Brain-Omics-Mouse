# Pies de figuras

## ozone_effect_uncertainty_es

A: Cambios log2 de ozono frente a aire sin regularizar e intervalos puntuales de Wald del95%; 3 animales por exposición en cada estrato, con dispersión estimada en24 animales. Símbolos rellenos:FDR de BH del estrato<0.05; los intervalos no son simultáneos. B:Efectos regularizados por apeglm de los mismos genes. Gris o ausente indica estimaciones no disponibles, no cero. Seis genes se seleccionaron por la LRT global entre genes codificantes con conteo medio>=20; Fos, Gfap, Aif1, Trem2 son objetivos secundarios de visualización. La selección es posterior al análisis y no establece replicación génica. La regularización se estima por separado en cada contraste y puede llevar el modo casi a cero pese a una prueba de Wald significativa sin regularizar; por ejemplo, Fos en M17/APOE4 tiene MLE=-2.271, q=0.0212, pero modo≈0. La matriz diagnostica la regularización, no demuestra ausencia del efecto original de exposición.

Alt: Cuatro gráficos de efectos e intervalos y un mapa de efectos regularizados de diez genes seleccionados del estudio de ozono.

Source: `results/advanced_omics/figures/ozone_effect_uncertainty_es.png`

## ozone_interaction_evidence_es

A:LRT de modelos anidados con universo génico fijo; cada gen es un punto. Cada prueba usa una familia de BH separada. Se conservan24 animales y los ajustes por longitud efectiva; los rangos completo/reducido son8/6, 8/6 y8/7. B:Ocho genes muestran respuesta dependiente de la edad; ninguno supera FDR<0.05 en la prueba ómnibus de genotipo ni en la interacción triple. Un resultado nulo no implica equivalencia. C:Diferencias sin regularizar de los cambios log2 por ozono en los ocho hallazgos de edad; los puntos negros indican FDR<0.05 en la familia condicional de Wald correspondiente. Las pruebas condicionales de1gl y las ómnibus responden hipótesis distintas; no requieren los mismos hallazgos. Son análisis secundarios exploratorios.

Alt: Tres gráficos de evidencia de interacción, conteos de hallazgos y mapa de diferencias de respuesta con signo.

Source: `results/advanced_omics/figures/ozone_interaction_evidence_es.png`

## ozone_interaction_animals_es

Se muestran los24 animales independientes por gen, 3 por combinación exposición/edad/genotipo. Los conteos usan los factores originales de normalización por longitud génica de DESeq2 y después log2(1+conteo) para visualizar. Los puntos usan rangos verticales específicos por gen que no comienzan en cero, iguales entre grupos dentro de cada panel. Las barras conectadas son medias aritméticas de los valores transformados; no conectan animales pareados ni son intervalos inferenciales. Se muestran los seis menores valores q de la LRT de respuesta dependiente de la edad, como selección secundaria. Las pruebas del modelo de conteos no se ajustaron sobre estos valores logarítmicos. Los nombres génicos y los grupos no establecen expresión específica de células.

Alt: Seis paneles con los valores de los24 animales individuales en genes seleccionados de respuesta por edad.

Source: `results/advanced_omics/figures/ozone_interaction_animals_es.png`

## GO_multi_contrast_heatmap_es

Se muestran NES con signo en una selección determinista con reducción de redundancia:hasta2 términos con FDR<0.05 por dirección y contraste, ordenados por q estimado y |NES|; un término nuevo requiere Jaccard<0.75 respecto de los ya elegidos, usando la unión de miembros elegibles. Los círculos delineados indican FDR algorítmica de GSEA<0.05 en esa celda; no son valores p de BH. Gris× significa no elegible/no analizado, nunca enriquecimiento cero. Las permutaciones de conjuntos no conservan la correlación génica ni la incertidumbre muestral. Los contrastes tienen tejidos/modelos/ensayos distintos y no se combinan; diferencias de NES o q no son pruebas formales entre estudios. Se conservan los nombres originales de GO.

Alt: Matriz de NES con signo en seis contrastes, con elegibilidad ausente explícita y marcas de FDR algorítmica.

Source: `results/advanced_omics/figures/GO_multi_contrast_heatmap_es.png`

## GO_overlap_network_es

Diez términos GO-BP de NES positivo y diez de NES negativo con menor q para ozono M3/APOE4, desempatados por |NES|. Los nodos representan tamaño del conjunto en la lista elegible y NES, no abundancia celular. Las aristas unen conjuntos con Jaccard>=0.25; no representan coexpresión, interacción proteica, regulación causal ni vías independientes. Los componentes se distribuyen en una cuadrícula determinista con posiciones circulares; las distancias son gráficas, no anatómicas ni moleculares. La selección es exploratoria y corresponde a un único contraste; se conservan las estimaciones originales de FDR algorítmica. Los q iguales a cero reflejan permutaciones finitas. Se proporcionan nombres completos, miembros y tablas de nodos/aristas.

Alt: Red de solapamiento de miembros de veinte términos GO seleccionados con tabla numerada.

Source: `results/advanced_omics/figures/GO_overlap_network_es.png`

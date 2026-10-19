# Figure captions

## ozone_effect_uncertainty

A: Unshrunken ozone-versus-air log2 fold changes and pointwise 95% Wald intervals; 3 animals per exposure in each stratum, with dispersion estimation shared across 24 animals. Filled symbols mark within-stratum BH FDR<0.05; intervals are not simultaneous. B: apeglm-regularized effects on the same genes. Gray/missing marks mean unavailable estimates, not zero. Six genes were selected by the global ozone LRT among protein-coding genes with mean normalized count>=20; Fos, Gfap, Aif1, Trem2 are named secondary display targets. Selection is post hoc and does not establish gene-specific replication. Regularization is estimated separately for each contrast and can contract a mode almost to zero despite a significant unshrunken Wald result; for example, Fos in M17/APOE4 has MLE=-2.271, q=0.0212, but mode≈0. The matrix is a regularization diagnostic, not evidence that the original exposure effect is absent.

Alt: Four forest plots and a regularized-effect heatmap for ten selected ozone-study genes.

Source: `results/advanced_omics/figures/ozone_effect_uncertainty.png`

## ozone_interaction_evidence

A: Fixed-universe nested-model likelihood-ratio tests; each gene is one point. A separate BH family is used for each test. All 24 animals and effective-length offsets are retained; the full/reduced ranks are 8/6, 8/6, and8/7. B: Eight age-dependent-response genes, zero genotype-omnibus and zero three-way genes meet BH FDR<0.05. A null test is not equivalence. C: Unshrunken differences of ozone log2 fold changes for all eight age-LRT discoveries; black dots indicate the corresponding conditional Wald-family FDR<0.05. These one-df conditional tests and omnibus tests answer different hypotheses; their discovery sets need not agree. All analyses are secondary and exploratory.

Alt: Three interaction evidence plots, discovery counts and a signed response-difference heatmap.

Source: `results/advanced_omics/figures/ozone_interaction_evidence.png`

## ozone_interaction_animals

All 24 independent animals per gene are shown, 3 per exposure/age/genotype cell. Counts use the original DESeq2 gene-length normalization factors, then log2(1+count) for display. Point-position plots use gene-specific nonzero y-ranges, identical across groups within a panel. Connected bars represent arithmetic means of displayed transformed values; they do not connect paired animals and are not inferential intervals. Six lowest-q age-dependent-exposure LRT genes are displayed, preserving secondary selection. Count-model tests were not fitted on these log-transformed values. Gene names and grouping do not establish cell-specific expression.

Alt: Six panels showing all24 individual animal values for selected age-response genes.

Source: `results/advanced_omics/figures/ozone_interaction_animals.png`

## GO_multi_contrast_heatmap

Signed NES values are displayed for a deterministic, redundancy-reduced selection:up to2 FDR<0.05 terms per direction per contrast, ordered by estimated q then absolute NES; new terms require Jaccard<0.75 with already selected terms using the union of rank-eligible GO members. Outlined circles mark that cell's GSEA algorithmic FDR<0.05; these are not BH p-values. Gray× means the term was not eligible/tested, never zero enrichment. Gene-set permutations do not preserve gene-gene correlation or sample uncertainty. Contrasts have distinct tissues/models/assays and are not pooled; differences in NES or q do not constitute a formal between-study test. Source GO names are retained.

Alt: A signed NES matrix across six contrasts with explicit missing eligibility and algorithmic FDR markers.

Source: `results/advanced_omics/figures/GO_multi_contrast_heatmap.png`

## GO_overlap_network

Ten lowest-q positive and ten lowest-q negative GO-BP terms for ozone M3/APOE4, with |NES| tie-breaks. Nodes encode gene-set size in the eligible rank and NES, not cell abundance. Edges join sets with Jaccard membership overlap>=0.25; they do not represent coexpression, protein interaction, causal regulation or independent pathway discoveries. Connected components use a deterministic grid and circular positions; visual distances are graphical, not anatomical or molecular. Selection is exploratory and restricted to one contrast; the original algorithmic GSEA FDR estimates are retained. Zero q values reflect finite permutations. Full names, memberships, node and edge tables are available.

Alt: A membership-overlap network of twenty selected GO terms with a numbered term table.

Source: `results/advanced_omics/figures/GO_overlap_network.png`

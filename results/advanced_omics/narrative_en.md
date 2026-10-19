# Secondary analyses of ozone response heterogeneity and comparative pathway patterns

## Introduction

Bulk brain RNA-seq can reveal exposure-associated transcriptional changes, but different numbers of differentially expressed genes across age or genotype groups do not, by themselves, demonstrate different exposure responses. The response difference must be tested directly. We therefore extended the existing ozone analysis with explicit age-by-exposure and genotype-by-exposure tests, alongside figures that show effect size, uncertainty, individual animals, and pathway redundancy. These are new secondary exploratory analyses; they are distinct from attempts to reproduce the original authors’ statistical workflow.

The biological unit is the animal. The ozone experiment contains 24 animals in eight age × APOE genotype × exposure groups, with three animals per group: 3 or 17 months (M3/M17), APOE3 or APOE4 (E3/E4), and filtered air or ozone (FA/O3). These measurements are bulk tissue expression, not exposed single-cell or spatial transcriptomics. Cross-study pathway figures additionally show the existing hypothalamus and anterior cingulate cortex PM2.5 analyses as separate contrasts; their animals and controls are never pooled.

## Methods

### Formal interaction models

We reused the stored DESeq2 full-model object (`results/ozone/dds_full_wald.rds`), its fixed 18,537-gene retained universe, sample-specific tximport effective-length normalization factors, and fitted dispersions. The full design was `~ age * genotype * exposure`. Dispersion estimates were preserved exactly in all added likelihood-ratio tests (LRTs). Each full-model design matrix had rank eight, and each reduced model was checked for rank and nesting.

Three model comparisons were declared before these secondary fits: (1) genotype-dependent exposure response, full versus `~ age * genotype + exposure + age:exposure` (two degrees of freedom); (2) age-dependent exposure response, full versus `~ age * genotype + exposure + genotype:exposure` (two degrees of freedom); and (3) the three-way interaction, full versus `~ (age + genotype + exposure)^2` (one degree of freedom). The first two tests jointly evaluate the relevant interaction terms across both levels of the other factor. An LRT is an omnibus test and does not have a single signed fold-change interpretation. This construction follows the parameter-count and interaction principles documented by [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions).

Cook’s filtering remained enabled and adaptive independent filtering was disabled. Nonconverged full or reduced fits would have been withheld before multiplicity adjustment; none occurred. Each LRT retained 18,511 finite p-values, with 26 genes withheld by the existing Cook’s filtering. Benjamini–Hochberg (BH) adjustment was performed separately within each declared LRT family. A supplementary column adjusts all finite p-values across the three LRT families together; this is a more extensive multiplicity family and is reported explicitly rather than substituted silently.

Five signed Wald contrasts estimated the E4–E3 difference in ozone response at each age, the M17–M3 difference in ozone response within each genotype, and the three-way difference of response differences. Numeric linear combinations of the original full-model coefficients provided covariance-aware standard errors. We report unshrunken maximum-likelihood estimates (MLEs) and pointwise 95% Wald confidence intervals, MLE ± 1.959964 × SE. They are not simultaneous intervals. Estimates requiring an entirely zero-count six-animal age/genotype stratum were additionally withheld, and each of the five Wald contrasts had a separate BH family. Joint adjustment across all five is provided as a supplementary sensitivity analysis.

### Figure and term selection

The ozone forest plot displays the six lowest global-exposure-LRT q-values among protein-coding genes with mean normalized count at least 20, plus the named display targets Fos, Gfap, Aif1 and Trem2. This selection is post hoc. The interaction heatmap retains all eight age-LRT discoveries, including the weakly expressed Gm21586. Individual-animal panels show the six strongest age-LRT genes using log2(1 + normalized count). Lines join exposure-group means of the displayed values, not paired animals or confidence limits. Each gene uses its own y-axis range.

The comparative GO Biological Process heatmap uses existing preranked GSEA results for six separate contrasts. We selected up to two significant terms per direction and contrast, ordered by algorithmic GSEA FDR, absolute normalized enrichment score (NES), and GO identifier. A greedy Jaccard filter removed highly redundant candidate sets (membership overlap at least 0.75), yielding 22 terms. The final display then shows each selected term wherever it was eligible in every contrast. Missing eligibility is gray with a cross and is never encoded as zero. Significance circles use the original algorithmic GSEA FDR, not BH-adjusted gene p-values.

The GO overlap network uses the ten most significant positive-NES and ten negative-NES terms in young APOE4 ozone animals. Edges represent Jaccard overlap of at least 0.25 among GO members present in that contrast’s gene ranking. Node area is proportional to eligible set size and node color represents NES. The layout is a deterministic grid of connected components with circular positions within each component; spatial distances have no anatomical or molecular interpretation. This is a term-membership redundancy network, not a regulatory or coexpression network.

## Results

### Age-related response heterogeneity is supported for a small gene set

At within-test BH FDR < 0.05, eight genes showed age-dependent ozone responses: **Ryr3, Rps24, Pkig, Rpl26, Inpp5d, Stk11, Rpl34 and Gm21586**. Ryr3 had the strongest evidence (LRT statistic 28.357, two degrees of freedom; p = 6.96 × 10−7; within-family q = 0.0129). The remaining seven q-values ranged from 0.0396 to 0.0412. When BH adjustment was expanded across all three LRT families, only Ryr3 remained below 0.05 (joint q = 0.0386). This sensitivity matters: the evidence is a focused secondary signal rather than widespread, multiplicity-insensitive age modulation.

No genes crossed within-family FDR < 0.05 for the two-degree-of-freedom genotype-dependent-response omnibus test or the one-degree-of-freedom three-way test. These null findings do not establish equivalent responses. With only three animals per exposure group, small or variable interaction effects may remain imprecise.

Conditional Wald tests yielded 6 and 4 genotype-response differences at M3 and M17, and 2 and 4 age-response differences within E3 and E4, respectively; the three-way Wald test yielded none. Joint adjustment across all five Wald families yielded 3, 3, 2, 3 and 0, respectively. These counts refer to different hypotheses and multiplicity families from the omnibus tests. Their discovery sets therefore need not coincide. For example, the Ryr3 age-response difference was +1.106 log2 units in E3 (pointwise 95% CI 0.436 to 1.776; conditional q = 0.702) and +1.456 in E4 (0.786 to 2.126; q = 0.0591). A pointwise interval excluding zero does not imply genome-wide BH significance.

The individual-animal panels reveal the direction and dispersion underlying these selected age-response signals. They should be read with the formal models, rather than used to infer untested pairwise effects from visual separation. Bulk tissue changes may reflect cell-state changes, cell-composition changes, or both; these data do not separate those explanations.

### Regularization strongly contracts some existing ozone estimates

The forest plot uses MLEs and their intervals for effect interpretation. The adjacent apeglm matrix is a diagnostic of the existing regularized estimates. Three ozone strata show strong contraction of most modes toward zero. For example, Fos in M17/APOE4 has MLE −2.271 log2 units, SE 0.474, and within-stratum q = 0.0212, but its saved apeglm mode is approximately −0.000009 (posterior SD 0.00144). The original prior hyperparameters were not retained, so this report does not assert a verified prior scale or diagnose a computational failure. It documents a material difference between estimands and regularization behavior. A near-zero regularized mode is not evidence that the corresponding original Wald exposure association is absent. Original source results have not been changed.

### Pathway patterns differ across distinct experiments and ozone strata

The six-contrast heatmap highlights shared annotations with differing ranking directions. For oxidative phosphorylation, NES was +1.847 in hypothalamus PM2.5 (algorithmic FDR 0.0326), +3.205 in M3/E4 ozone (estimated algorithmic FDR 0), and −3.054 in M17/E3 ozone (estimated FDR 0). Translation at synapse showed positive NES in both PM2.5 contrasts and M3/E4 ozone, but negative NES in both older ozone strata. These are descriptive comparisons of separately normalized ranking-based enrichment scores, not a statistical test of pathway differences between studies or an assay of pathway activation. Values estimated as zero reflect finite permutations, not probability zero.

The M3/E4 overlap network contains 20 selected terms and 15 edges. Positive-NES mitochondrial energy and translation annotations form overlapping components; several negative-NES developmental or synaptic annotations are less connected at the declared threshold. Multiple related GO labels can therefore summarize some of the same ranked genes. Apparent annotation breadth should not be counted as independent mechanistic evidence. Names such as myoblast fusion, lens development or cardiac hypertrophy are ontology annotations of shared genes and do not establish those tissue phenotypes in brain.

## Interpretation and reproducibility limits

The strongest new finding is limited evidence of age-dependent ozone transcriptional responses, with Ryr3 surviving the broader three-family adjustment. Genotype-wide response heterogeneity and three-way modification are not established by the declared omnibus tests. These analyses extend the current reanalysis and are not a replication of an original author claim or an independent cohort validation.

All gene selections are disclosed, every plotted table is exported, and English/Spanish figures are available as PNG, SVG and PDF. The pathway permutation procedure uses the existing 1,000 gene-set permutations; it does not preserve gene–gene correlation or propagate sample-level uncertainty. Different regions, exposure methods, age structures and eligible gene sets limit cross-study biological comparisons. No exposed single-cell transcriptome, anatomical localization of an exposure effect, or causal network is reconstructed here.

`validation.json` independently checks LRT chi-square tails, Wald normal tails, BH calculations, interval arithmetic, convergence masks and model invariants. This validates computation and provenance, not biological truth. Figure captions, alternative text and export hashes are in `figure_manifest.json`; inputs and methods are in `metadata/advanced_omics/methods.json`; visual inspection is recorded in `visual_QA.json`.

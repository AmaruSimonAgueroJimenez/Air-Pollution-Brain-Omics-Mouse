# Regional gene priorities and encoded products

## Introduction

Regional transcriptomics can prioritize reproducible molecular candidates while a spatial reference atlas shows where their baseline RNA is detected. These two measurements answer different questions: differential expression estimates a regional exposure association, and an unexposed atlas supplies spatial context. Functional annotation identifies an encoded product or noncoding locus; it does not measure that product's abundance or activity.

## Methods

The existing primary hypothalamus and cingulate-cortex differential-expression tables were ranked without refitting models or removing animals. Four selection groups retain all 11 hypothalamic BH q<0.05 hits, the 20 lowest-q cortical protein-coding hits, the ten measured hypothalamic genes with lowest finite q, and the ten lowest-q measured cortical hits. Finite effect estimates and standard errors are required for measured-panel lists. Cortex panel availability uses exact original Ensembl IDs; the symbol-only hypothalamic source uses exact original symbols. The original Pcnx symbol is preserved and its unique current official alias Pcnx1 is disclosed.

Official NCBI mouse gene types and summaries and uniquely mapped reviewed mouse UniProt records were retrieved and frozen. Current primary names take precedence over historical aliases; this prevents the historical Fzd3 alias of Fzd9 from creating a false mapping. English protein names follow the reviewed record where available. Spanish names and both functional summaries are curated paraphrases with prediction/similarity qualifications preserved. Every gene links to its official record. Source URLs, hashes, selected records, software versions and ranking rules accompany the analysis. A default rerun uses this frozen annotation snapshot.

The forest figures show the unchanged MLE log2 fold change and pointwise normal 95% interval (MLE±1.96SE), with the shrunken estimate as a separate diamond. These are distinct estimates; the intervals are not confidence intervals for the shrunken coefficient and do not adjust for multiplicity or candidate selection. Each panel labels its own horizontal scale.

## Results and interpretation

There are 46 unique annotated loci across 51 selection rows: 43 protein-coding loci, two noncoding RNAs and one pseudogene. All nine coding hypothalamic primary hits decrease, including Hipk3, Ksr2, Tll1, Fzd3 and Fstl4. The two increases correspond to Gm13423 (an Npm1 pseudogene) and B230312C02Rik (noncoding RNA); functional NPM1 protein cannot be inferred from the former signal. All eleven directions persist in the seven leave-one-animal-out fits, but no gene passes BH q<0.05 in all seven fits.

The leading cortical genes include Hspa5/BiP, Pdia4, Xbp1, Pdia6, Creld2, Mlec and Manf, whose official annotations converge on protein folding and endoplasmic-reticulum homeostasis. This is descriptive convergence, not a new enrichment test or proof of an activated pathway. Fosb, Egr3 and Fosl2 encode transcription regulators; Homer1 encodes a postsynaptic scaffold; Nr1d1 and Ciart are circadian transcriptional repressors. Gene-level counts do not isolate ΔFosB, Homer1a or spliced XBP1s. Ak4 is described conservatively as a mitochondrial nucleotide/energy-homeostasis protein: its mouse annotation is by similarity to human AK4, whose reviewed record records both inactive and active catalytic studies. Ksr2 is treated as a MAPK scaffold, without assuming active kinase catalysis in these samples. Detailed molecular roles remain unresolved in the retrieved mouse records for Pcnx1 and Lyrm9.

Only Tll1 among the ten measured hypothalamic candidates passes primary BH q<0.05. The other nine are explicitly exploratory; six large MLE effects shrink close to zero. In cortex, all ten measured candidates pass primary BH q<0.05, including noncoding C730002L08Rik. Atlas measurement alone does not validate an exposure effect. Primary hits absent from the panel remain in the tables with atlas availability marked false.

## Limits

Small animal numbers, hypothalamic sensitivity to individual animals, and unresolved cortical sex/identity metadata limit interpretation. Missing cortical leave-one-out values are unavailable, not zero stability. The reference atlas contains unexposed cells; a regional RNA effect projected onto it remains a simulation and cannot identify exposed cell types. No protein abundance, enzyme activity, isoform-specific response or causal mechanism is established here. Statistical rank is an explicit prioritization rule rather than proof of biological importance.


## Official annotation sources

- Ak4 (Ak4): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/11639); [UniProt](https://www.uniprot.org/uniprotkb/Q9WUR9/entry)
- B230312C02Rik (B230312C02Rik): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/320485)
- Bcl11b (Bcl11b): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/58208); [UniProt](https://www.uniprot.org/uniprotkb/Q99PV8/entry)
- C730002L08Rik (C730002L08Rik): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/320155)
- Car12 (Ca12): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/76459); [UniProt](https://www.uniprot.org/uniprotkb/Q8CI85/entry)
- Cbfa2t3 (Cbfa2t3): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/12398); [UniProt](https://www.uniprot.org/uniprotkb/O54972/entry)
- Ciart (Ciart): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/229599); [UniProt](https://www.uniprot.org/uniprotkb/Q3TQ03/entry)
- Col19a1 (Col19a1): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/12823); [UniProt](https://www.uniprot.org/uniprotkb/Q0VF58/entry)
- Cpne9 (Cpne9): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/211232); [UniProt](https://www.uniprot.org/uniprotkb/Q1RLL3/entry)
- Creld2 (Creld2): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/76737); [UniProt](https://www.uniprot.org/uniprotkb/Q9CYA0/entry)
- Crym (Crym): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/12971); [UniProt](https://www.uniprot.org/uniprotkb/O54983/entry)
- Egr2 (Egr2): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/13654); [UniProt](https://www.uniprot.org/uniprotkb/P08152/entry)
- Egr3 (Egr3): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/13655); [UniProt](https://www.uniprot.org/uniprotkb/P43300/entry)
- Fmnl1 (Fmnl1): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/57778); [UniProt](https://www.uniprot.org/uniprotkb/Q9JL26/entry)
- Fosb (Fosb): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/14282); [UniProt](https://www.uniprot.org/uniprotkb/P13346/entry)
- Fosl2 (Fosl2): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/14284); [UniProt](https://www.uniprot.org/uniprotkb/P47930/entry)
- Fstl4 (Fstl4): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/320027); [UniProt](https://www.uniprot.org/uniprotkb/Q5STE3/entry)
- Fzd3 (Fzd3): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/14365); [UniProt](https://www.uniprot.org/uniprotkb/Q61086/entry)
- Gm13423 (Gm13423): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/100416062)
- Hgf (Hgf): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/15234); [UniProt](https://www.uniprot.org/uniprotkb/Q08048/entry)
- Hipk3 (Hipk3): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/15259); [UniProt](https://www.uniprot.org/uniprotkb/Q9ERH7/entry)
- Homer1 (Homer1): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/26556); [UniProt](https://www.uniprot.org/uniprotkb/Q9Z2Y3/entry)
- Hspa5 (Hspa5): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/14828); [UniProt](https://www.uniprot.org/uniprotkb/P20029/entry)
- Iqgap2 (Iqgap2): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/544963); [UniProt](https://www.uniprot.org/uniprotkb/Q3UQ44/entry)
- Ksr2 (Ksr2): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/333050); [UniProt](https://www.uniprot.org/uniprotkb/Q3UVC0/entry)
- Lmbrd2 (Lmbrd2): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/320506); [UniProt](https://www.uniprot.org/uniprotkb/Q8C561/entry)
- Lyrm9 (Lyrm9): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/66274); [UniProt](https://www.uniprot.org/uniprotkb/Q3UN90/entry)
- Manf (Manf): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/74840); [UniProt](https://www.uniprot.org/uniprotkb/Q9CXI5/entry)
- Mlec (Mlec): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/109154); [UniProt](https://www.uniprot.org/uniprotkb/Q6ZQI3/entry)
- Mmp16 (Mmp16): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/17389); [UniProt](https://www.uniprot.org/uniprotkb/Q9WTR0/entry)
- Nr1d1 (Nr1d1): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/217166); [UniProt](https://www.uniprot.org/uniprotkb/Q3UV55/entry)
- Nr4a3 (Nr4a3): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/18124); [UniProt](https://www.uniprot.org/uniprotkb/Q9QZB6/entry)
- P4ha1 (P4ha1): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/18451); [UniProt](https://www.uniprot.org/uniprotkb/Q60715/entry)
- Pcdh8 (Pcdh8): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/18530); [UniProt](https://www.uniprot.org/uniprotkb/Q7TSK3/entry)
- Pcnx (Pcnx1): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/54604); [UniProt](https://www.uniprot.org/uniprotkb/Q9QYC1/entry)
- Pdia4 (Pdia4): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/12304); [UniProt](https://www.uniprot.org/uniprotkb/P08003/entry)
- Pdia6 (Pdia6): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/71853); [UniProt](https://www.uniprot.org/uniprotkb/Q922R8/entry)
- Ptar1 (Ptar1): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/72351)
- Rorb (Rorb): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/225998); [UniProt](https://www.uniprot.org/uniprotkb/Q8R1B8/entry)
- S1pr1 (S1pr1): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/13609); [UniProt](https://www.uniprot.org/uniprotkb/O08530/entry)
- Sgsm1 (Sgsm1): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/52850); [UniProt](https://www.uniprot.org/uniprotkb/Q8BPQ7/entry)
- Tll1 (Tll1): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/21892); [UniProt](https://www.uniprot.org/uniprotkb/Q62381/entry)
- Usp43 (Usp43): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/216835); [UniProt](https://www.uniprot.org/uniprotkb/Q8BUM9/entry)
- Xbp1 (Xbp1): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/22433); [UniProt](https://www.uniprot.org/uniprotkb/O35426/entry)
- Zbtb40 (Zbtb40): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/230848)
- Zkscan2 (Zkscan2): [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/210162)

Ak4 catalytic qualification: [human reviewed AK4](https://www.uniprot.org/uniprotkb/P27144/entry), [Liu et al.](https://pubmed.ncbi.nlm.nih.gov/19130895/), [Panayiotou et al.](https://pubmed.ncbi.nlm.nih.gov/19766732/).

# Related work

Every reference below was checked against PubMed (PMIDs listed) on 2026-09-24.
The descriptions are short summaries for orientation. Read the papers before
citing them for anything specific.

## Studies closest to this project's question

| Study | Disease | Expression data | Disease phenotype actually measured | Relation to this project |
|---|---|---|---|---|
| Freeze et al. 2018, *NeuroImage: Clinical* 18:456-466 (PMID 29868450) | PD | AHBA | Regional atrophy in PD patients | Links a candidate gene set (alpha-synuclein transfer) to measured atrophy. This project uses GWAS-derived genes and has no measured phenotype. |
| Keo et al. 2020, *Communications Biology* 3:101 (PMID 32139796) | PD | AHBA | Braak-stage regions (anatomically defined) | Closest PD comparator. It studies transcriptomic signatures of regions affected at different disease stages. |
| Zheng et al. 2019, *PLoS Biology* 17:e3000495 (PMID 31751329) | PD | AHBA | Atrophy maps from PD patients' MRI, plus connectivity | Models spread through connectivity, weighted by regional gene expression. |
| Rittman et al. 2016, *Neurobiology of Aging* 48:153-160 (PMID 27697694) | PD, PSP | AHBA | Loss of network hubs, cognition | Single-gene (MAPT) regional expression against network measures. |
| Grothe et al. 2018, *Brain* 141:2755-2771 (PMID 30016411) | AD | AHBA | Amyloid and tau PET patterns | Regional molecular properties against measured AD pathology. |
| Anand et al. 2025, *Brain* 148:3679-3693 (PMID 40631882) | AD | AHBA | Tau PET, plus the connectome | Recent AD comparator. AD risk-gene expression and connectivity explain regional tau vulnerability. |
| McColgan et al. 2018, *Biological Psychiatry* 83:456-465 (PMID 29174593) | HD | AHBA | White matter loss in HD | Enrichment of gene categories in regions that degenerate. |
| Liscovitch & French 2014, *PLoS ONE* 9:e115029 (PMID 25493648) | PD | AHBA | None (expression and co-expression only) | Source of the SNCA checkpoint: SNc 2nd of 105 regions, SNr 88th. |
| Kamath et al. 2022, *Nature Neuroscience* 25:588-595 (PMID 35513515) | PD | Single-nucleus RNA-seq from human SN | Loss of a specific dopamine-neuron subtype in PD | Shows vulnerability is cell-type specific, below the resolution of any AHBA region. |

## Methods references

| Reference | Used for |
|---|---|
| Hawrylycz et al. 2012, *Nature* 489:391-399 (PMID 22996553) | The Allen Human Brain Atlas dataset |
| Arnatkeviciute et al. 2019, *NeuroImage* 189:353-367 (PMID 30648605) | AHBA processing steps (probe reannotation, intensity filtering, probe selection, normalization) |
| Markello et al. 2021, *eLife* 10:e72129 (PMID 34783653) | abagen, and evidence that processing choices change downstream results |
| Burt et al. 2020, *NeuroImage* 220:117038 (PMID 32585343) | BrainSMASH spatial null model |
| Fulcher et al. 2021, *Nature Communications* 12:2669 (PMID 33976144) | False positives in gene-set enrichment on spatial brain data |
| Wu & Smyth 2012, *Nucleic Acids Research* 40:e133 (PMID 22638577) | CAMERA competitive gene-set test with an inter-gene correlation correction |
| Pauli et al. 2018, *Scientific Data* 5:180063 (PMID 29664465) | CIT168 subcortical atlas (SNc and SNr) |

## What this project can add

Most of the studies above relate AHBA expression to a phenotype measured in
patients (atrophy, PET, Braak regions). This project has no patient
phenotype. Its targets are the regions a disease is known to damage early.

A defensible contribution is therefore narrow. It tests how stable GWAS
gene-set enrichment in healthy tissue is across:

- how the substantia nigra is defined (anatomical sample labels versus atlas
  coordinates)
- normalization choices
- donors
- gene-mapping rules
- two different null models

It also includes a direct replication check (SNCA) and one external dataset
(GTEx). Whether that is new has to be checked against the literature,
especially Keo et al. 2020 for PD and Anand et al. 2025 for AD, before any
novelty claim goes into a write-up.

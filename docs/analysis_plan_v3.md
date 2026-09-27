# Analysis plan, v3: cell types inside the damaged region

Frozen 2026-09-27, 11:25 CDT. At that point no expression value of any gene
had been computed from the datasets below. Only their metadata (donors,
disease labels, cell-type labels, cell counts) and the storage format (raw
integer counts) had been checked.

## What was already known

- v2 (`docs/v2_results.md`): PD and AD GWAS gene sets show no large
  enrichment in the region each disease damages, in bulk AHBA tissue.
- Amendment 1: cell-type marker genes are detected easily, but the design
  could only have detected large regional effects (PD would have needed
  z of about 4.75; AD about 3.75).
- Published single-cell work locates PD genetic risk in dopamine (DA)
  neurons: Kamath et al. 2022 (Nat Neurosci 25:588, PMID 35513515), and
  Agarwal et al. 2020 (Nat Commun 11:4183, PMID 32826893), which also
  implicates oligodendrocytes. Bryois et al. 2020 (Nat Genet 52:482,
  PMID 32341526) reached similar conclusions from mouse data. AD risk
  variants sit in microglial regulatory elements (Nott et al. 2019,
  Science 366:1134, PMID 31727856).

The primary hypotheses come from that literature, so v3 is partly a
replication. It uses this project's own gene sets and its own
random-gene-set logic, and asks something the regional analysis couldn't:
where inside the damaged tissue are these genes expressed, and would that
be visible in bulk tissue?

## Question

Within the tissue each disease damages, are its risk genes specifically
expressed in one cell type? How much of that tissue's total signal for
those genes comes from that cell type?

## Data

All from CELLxGENE, raw counts, downloaded to `data/raw/singlecell/`
(gitignored; `src/singlecell.py` lists the URLs).

| Dataset | Tissue | Donors used | Role |
|---|---|---|---|
| Kamath et al. 2022, per-class files (DA neurons, non-DA neurons, astrocytes, microglia, oligodendrocytes, OPCs, endothelial) | SN pars compacta | control donors only (disease = normal) | PD primary |
| Siletti et al. 2023 (Science 382:eadd7046, PMID 37824663), dissection "Midbrain - Substantia Nigra" | SN | all 3 | PD replication; dilution |
| Siletti et al. 2023, medial and lateral entorhinal cortex dissections, pooled | entorhinal cortex | all 3 | AD primary |
| Siletti et al. 2023, caudate body and putamen dissections, pooled | striatum | all 3 | HD exploratory |

The Nurr1-sorted Kamath files aren't used.

## Cell classes

- **Kamath:** one class per file.
- **Siletti**, from `supercluster_term`:
  - Oligodendrocyte → oligodendrocytes
  - Oligodendrocyte precursor and Committed oligodendrocyte precursor → OPCs
  - Astrocyte → astrocytes
  - Microglia supercluster, restricted to `cell_type` = microglial cell → microglia. Other immune cells are dropped.
  - Vascular and Fibroblast → vascular
  - Miscellaneous and Ependymal → dropped
  - every other supercluster → neurons
- **Siletti SN only:** neurons are split into DA neurons and other neurons.
  DA neurons are the neuronal clusters (`cluster_id`) in which at least 20%
  of cells have one or more SLC6A3 (dopamine transporter) UMIs. This rule
  uses a marker gene, never a risk gene.
- **Siletti striatum only:** neurons are split into medium spiny neurons
  (superclusters Medium spiny neuron and Eccentric medium spiny neuron) and
  other neurons.
- Classes with fewer than 50 cells in a dataset are dropped.

## Measures

- **Expression:** each cell's counts are scaled to 10,000 total (CP10K). A
  class's expression of a gene is the mean over its cells, on the linear
  scale.
- **Gene matching:** by Ensembl ID via GENCODE v50, with a gene-symbol
  fallback. Unmatched genes are listed.
- **Universe:** protein-coding genes with nonzero expression in the dataset.
- **Specificity** (EWCE; Skene & Grant 2016, Front Neurosci 10:16,
  PMID 26858593): for gene g and class c,
  s_c(g) = mean_c(g) / sum over classes of mean(g). A set's score for class
  c is the mean of s_c over its genes.

## Tests

**Null.** 10,000 random gene sets of the same size, drawn from the universe.
Each random gene comes from the same expression decile (mean CP10K over all
cells) as the set gene it replaces. One-sided
p = (1 + #{null >= observed}) / (10,000 + 1).

**Primary tests (two), with Benjamini-Hochberg across them:**

1. PD top-10 set, DA neurons, Kamath control donors
2. AD top-10 set, microglia, Siletti entorhinal cortex

**Verdict for each:**

- **supported:** BH p < 0.05 and the target class has the highest set score
  of all classes
- **partial:** BH p < 0.05 but another class scores higher
- **not supported:** otherwise

**Reported with each primary test** (descriptive):

- donor level: the same score per donor (pseudobulk; classes with at least
  20 cells in that donor), against 1,000 matched random sets per donor, and
  the number of donors with a positive z
- leave-one-gene-out p-values

## Secondary and exploratory

Reported unadjusted and labelled.

- every class, for every set, in every dataset
- PD in Siletti SN (replication of test 1)
- PD in oligodendrocytes (Agarwal et al. 2020)
- the all-genome-wide-significant PD and AD sets
- within Kamath DA neurons (control donors), specificity across the 10
  published DA subtypes, and whether SOX6_AGTR1, the subtype lost in PD,
  scores highest (exploratory)
- HD: HTT and onset modifiers in striatal classes (exploratory)

## Positive controls

- **Dopaminergic markers** (the amendment 1 list) in DA neurons. In Siletti
  this is circular, because SLC6A3 defines the class; it's reported with
  that flag.
- **Microglia markers** in microglia (entorhinal): P2RY12, CX3CR1, TMEM119,
  CSF1R, C1QA, C1QB, AIF1, ITGAM, SALL1, HEXB.
- **Striatal projection neuron markers** (amendment 1 list) in medium spiny
  neurons.
- None of the marker genes is in any disease gene set (checked).
- Expected: the target class scores highest, p < 0.05.

## Dilution

Descriptive; Siletti SN, which was dissected without sorting.

- the fraction of SN nuclei that are DA neurons
- for each PD gene and each DA marker:
  - the share of all SN UMIs for that gene that come from DA-neuron nuclei
  - the ratio of its DA-neuron mean to its mean over all SN nuclei

This shows how far a DA-neuron-specific gene could move a bulk SN
measurement.

## Caveats, fixed now

- **Nuclei, not whole cells.** snRNA-seq misses cytoplasmic RNA, while AHBA
  measured tissue RNA.
- **Donor differences.** Siletti has 3 donors. Kamath's control donors are
  mostly aged 75 to 92, and the AHBA donors were 24 to 57.
- **What is measured.** Expression specificity of nearest genes, which is
  not heritability enrichment.
- **Replication.** Positive primary results would largely replicate
  published findings.

Outputs go to `results/v3/`. v3 doesn't change any v2 verdict.

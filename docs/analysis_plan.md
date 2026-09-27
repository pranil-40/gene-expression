# Analysis plan, v2

Frozen 2026-09-24, 19:10 CDT, before any v2 gene-set score, regional
contrast, SNCA rank, or p-value was computed.

## What was already known when this was written

v2 is a documented re-analysis. It is not a preregistered confirmatory
study, because the decisions below were made after seeing:

- the v1 results (commit `220b838`, archived in `archive/v1_exploratory/`):
  PD enrichment in a 6 mm "substantia nigra" mask, FDR p = 0.015; AD
  entorhinal, FDR p = 0.22; HD caudate+putamen, negative contrast
- an external methods review of the v1 draft (2026-09-24)
- Liscovitch & French 2014 (PLoS ONE 9:e115029), which reports SNCA is
  high in AHBA substantia nigra pars compacta (rank 2 of 105 regions) and
  low in pars reticulata (rank 88 of 105)
- the anatomy-only audit in `results/sn_audit/` (no expression values)

Region-level expression matrices for the v2 region definitions had been
built, and checked numerically against abagen's own output, but no gene's
values in them had been looked at.

## Question

Healthy-brain expression can't show why a region degenerates. The v2
question is narrower: do the genes near Parkinson's and Alzheimer's GWAS
loci show higher normalized expression in the region each disease damages
early than in the rest of the brain? Does that hold up against (a) spatial
autocorrelation, (b) comparable random gene sets, (c) dropping any one
donor or any one gene, and (d) other reasonable region definitions?

Huntington's is analysed separately and only as exploratory, because its
gene sources (the causal gene HTT; onset modifiers in mutation carriers)
answer a different question than PD/AD susceptibility loci.

## Data and preprocessing

- AHBA microarray data, 6 donors, processed with abagen 0.1.3 using the
  settings in `src/common.py` (`ABAGEN_KWARGS`): reannotated probes,
  intensity-based filtering at 0.5, differential-stability probe selection,
  SRS sample and gene normalization over matched samples, 2 mm matching
  tolerance, corrected MNI coordinates, no mirroring, **no interpolation of
  empty regions**.
- Left hemisphere plus brainstem only. Only 2 of 6 donors have
  right-hemisphere samples.
- Region value = mean over the samples in that region within a donor, then
  mean over the donors that have samples there.
- Gene universe = all 15,633 genes that survive abagen's filtering. Every
  gene set is scored against this same universe.

## Region definitions

| Name | Use | Substantia nigra definition |
|---|---|---|
| `dk_ontology_sn` | **primary** | AHBA samples whose anatomical label is substantia nigra pars compacta (SNc) or pars reticulata (SNr), left side, as two separate regions |
| `cit168` | sensitivity | SNc/SNr voxels from the CIT168 probabilistic atlas (p(SNc)+p(SNr) >= 0.5) |
| `sn6mm` | exploratory (v1-style) | brainstem voxels within 6 mm of an SN-labelled sample |

The primary SN definition is chosen from the anatomy audit alone.
Coordinate-based definitions don't capture SN tissue in AHBA: only 2 of 52
SN-labelled samples have coordinates inside the CIT168 SN, and the v1-style
mask's samples are mostly pontine nuclei, VTA, red nucleus and other
midbrain tissue (16 of 63 labelled SN). The choice was made before any SN
expression was examined. It was made knowing Liscovitch & French's
label-based SNc result, which is disclosed above.

Targets: PD = `substantia_nigra_pc`. AD = `entorhinal`. HD = `caudate` +
`putamen`. Reference = every other region in the region set. For PD that
includes SNr and brainstem.

## Gene sets

From `src/fetch_gwas_genes.py` (v2 rule): genome-wide significant loci,
clumped by distance and LD, one protein-coding gene per locus.

- Primary: top 10 independent loci by lead p-value (PD: Nalls et al. 2019;
  AD: Bellenguez et al. 2022).
- Sensitivity: all genome-wide significant loci; author-reported genes
  (where the Catalog has them); leave-one-gene-out.
- HD exploratory: HTT alone; GeM-HD onset-modifier genes (hand-entered,
  sourced in `provenance.json`).
- Genes missing from the AHBA universe are dropped and listed, with the
  reason (see `src/gene_probe_audit.py`).

## Statistics

For a gene set G and region r, S(r) = mean over genes in G of normalized
expression. Target contrast:

    T = mean of S over target regions - mean of S over reference regions
    z = T / SD of S over reference regions        (descriptive only)

z describes one map. It is not a Cohen's d, because regions aren't
independent subjects.

Primary tests, one-sided (enrichment):

1. **Spatial null.** 10,000 BrainSMASH surrogates of S, from Euclidean
   distances between region centroids. DK regions use atlas voxel
   centroids. Ontology SNc/SNr use the CIT168 left SNc/SNr centroids.
   Same T on each surrogate; p = (1 + #{T* >= T}) / (B + 1). Variogram fit
   diagnostics are saved.
2. **Gene-set null.** 10,000 random gene sets of the same size drawn from
   the gene universe, same T; p as above.
3. **Correlation-adjusted competitive test** (CAMERA-style, Wu & Smyth
   2012). Per-gene statistic t_g = (x_g(target) - mean x_g(reference)) /
   SD x_g(reference). The set's mean t is compared with the universe's
   mean and SD, inflating the variance by VIF = 1 + (k - 1) * rho, where rho
   is the mean pairwise correlation of the set's genes across regions.
   This covers the problem that random sets are less co-expressed than
   GWAS sets (Fulcher et al. 2021).

Multiple testing: the primary family is {PD, AD} x {spatial null,
gene-set null}, with Benjamini-Hochberg across those 4 p-values. The CAMERA
test, sensitivity analyses and all HD analyses are reported unadjusted and
labelled secondary or exploratory.

Robustness (descriptive, not tested):

- Donor-level T(d) for every donor that has samples in the target and at
  least 10 reference regions; the number of donors with T(d) > 0
- leave-one-donor-out pooled T and spatial p
- leave-one-gene-out pooled T and spatial p
- the same statistics under `cit168` and `sn6mm`, and for the sensitivity
  gene sets

Interpretation, fixed in advance for PD and AD:

- **supported**: both primary nulls significant after BH, T(d) > 0 in at
  least 5 of the donors with coverage, and T > 0 in every leave-one-gene-out
  run
- **partial**: at least one primary null significant, but not all of the
  above
- **not supported**: neither primary null significant

A significant PD result and a non-significant AD result would not be
evidence that the two diseases differ. That would need a direct test,
which isn't planned.

## SNCA checkpoint

A pipeline sanity check, not a hypothesis test. The target, following
Liscovitch & French 2014: raw Allen log2 expression, averaged over all
Allen-mapped SNCA probes; per donor, mean per anatomical structure; left
hemisphere structures sampled in all 6 donors; median across donors; rank.

"Reproduced" = SNc in the top 10% of structures and SNr below the median.
Also reported: the same ranking with abagen's selected probe, and SNc's
rank for SNCA in the primary v2 region matrix.

## External check (secondary)

GTEx median TPM across its 13 brain tissues. log2(TPM + 1), z-scored per
gene across tissues, averaged over the primary gene set. PD: substantia
nigra vs the other 12 tissues. HD: caudate and putamen vs the rest. AD is
not tested, because GTEx has no entorhinal cortex sample. Null: 10,000
random protein-coding gene sets of the same size. Gene sets and
statistics are exactly as above; nothing is re-tuned after seeing GTEx.

## Anything done after this point that isn't listed here

Anything else is labelled exploratory in `docs/revision_log.md` and in the
results tables.

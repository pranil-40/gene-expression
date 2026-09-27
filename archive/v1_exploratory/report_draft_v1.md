# Do Disease Risk Genes Concentrate Where the Disease Actually Strikes?

Pranil

## 1. Research question

Parkinson's kills the substantia nigra. Alzheimer's begins in the entorhinal
cortex. Huntington's destroys the striatum. The genes that raise risk for
each of these diseases are present in every cell in the body, not just the
region that ends up damaged - so why does the damage always land in the same
place?

This project tests whether the brain regions where each disease's GWAS risk
genes are most expressed are the same regions the disease goes on to
destroy, and whether that pattern survives a proper statistical test or just
looks convincing on a bar chart.

## 2. Method

**Gene lists.** Top 10 risk genes per disease, pulled from one large GWAS
study each, picked by PubMed ID before looking at any expression data:
Parkinson's from Nalls et al. 2019 (GCST009325), Alzheimer's from Bellenguez
et al. 2022 (GCST90027158), Huntington's from Moss et al. 2017 (GCST004691).
Genes ranked by association p-value, each lead SNP mapped to a gene through
the GWAS Catalog's own genomic-context annotation. Huntington's disease is
monogenic (caused by CAG repeat expansion in HTT), so there is no
case/control risk GWAS for it the way there is for PD/AD - the Moss et al.
study is a GWAS of disease-progression modifiers in confirmed HTT mutation
carriers, and HTT itself was added to the gene list by hand since a modifier
GWAS wouldn't otherwise surface the causal gene. Full selection rule in
`src/fetch_gwas_genes.py`.

**Expression data.** Regional gene expression from the Allen Human Brain
Atlas (6 postmortem donors), pulled with `abagen` using SRS normalization and
the Desikan-Killiany atlas (83 regions, cortical + subcortical).

**Target region per disease:** Alzheimer's -> entorhinal cortex (direct
atlas match). Huntington's -> striatum, taken as caudate + putamen (direct
match). Parkinson's -> substantia nigra, which the Desikan-Killiany atlas
does not resolve as its own region. Initially used the atlas's "Brain-Stem"
label as a stand-in, but that failed the validation checkpoint outright (see
below) - the actual raw AHBA donor files include real substantia-nigra
tissue samples with their own MNI coordinates, buried inside the generic
brainstem label along with the pons and medulla. Built a corrected region by
taking every brainstem voxel within 6mm of a real substantia-nigra sample
coordinate and relabeling it (`src/build_sn_atlas.py`), pooling 52 real
tissue samples across all 6 donors.

**Validation checkpoint.** Before trusting any of the above, checked that
SNCA (the single strongest Parkinson's GWAS hit, p=4e-154 in Nalls et al.)
shows elevated expression in substantia nigra. It does not - SNCA ranks 38th
of 43 regions, even after building the corrected substantia-nigra region.
This did not resolve into a straightforward pass. See Discussion.

**Statistical test.** For each disease: mean expression of its risk genes in
the target region vs. every other region (`src/regional_comparison.py`),
then a BrainSMASH spatial permutation test (1000 surrogates preserving the
real map's spatial autocorrelation) to get a p-value that accounts for the
fact that nearby brain regions look alike anyway (`src/permutation_test.py`).
Benjamini-Hochberg FDR correction applied across the three diseases.

## 3. Results

| Disease | Target region(s) | Cohen's d | p (spatial null) | p (FDR-corrected) |
|---|---|---|---|---|
| Parkinson's | substantia_nigra | 4.08 | 0.005 | **0.015** |
| Alzheimer's | entorhinal | 1.54 | 0.144 | 0.216 |
| Huntington's | caudate + putamen | -1.95 | 0.986 | 0.986 |

(full table: `results/tables/permutation_results.csv`; figures:
`results/figures/*.png`)

Parkinson's: substantia nigra is the single highest-scoring region out of 43
for the 9-gene PD risk average (`results/figures/parkinsons_regional_heatmap.png`),
and that holds up after accounting for spatial autocorrelation.

Alzheimer's: entorhinal cortex does show elevated risk-gene expression in
the raw comparison (Cohen's d = 1.54, a large effect), but once you build a
spatially-informed null distribution instead of comparing to zero, that
effect is no longer distinguishable from chance (FDR p = 0.22). Nearby
temporal-lobe regions likely share enough of the same expression profile
that "entorhinal is elevated relative to the rest of the brain" isn't as
special as it first looks.

Huntington's: caudate and putamen actually show *lower* average risk-gene
expression than the rest of the brain, the opposite of the hypothesis.

## 4. Discussion

Mixed result, and I think the mixed part is the actual finding here, not a
failure to get a clean answer. PD supports the "risk genes concentrate where
the damage happens" hypothesis under a real statistical test. AD looks
promising in a raw comparison but the effect dissolves once you test it
properly. HD goes the other direction entirely - lower expression in the
target region, not higher.

The Parkinson's validation checkpoint is worth dwelling on, because it
didn't behave the way I expected going in. SNCA - the single strongest PD
GWAS hit by a wide margin - does not show elevated expression in substantia
nigra, even in the properly-resolved region. But the *average* of all 9 PD
risk genes does, and strongly. That's a real distinction: whatever is
driving substantia nigra's vulnerability in PD, it isn't well captured by
looking at SNCA's raw transcript level in healthy tissue. SNCA's relevance
to PD is about protein aggregation dynamics in dopaminergic neurons (Lewy
body pathology), which isn't the same thing as "this region transcribes more
SNCA mRNA than other regions" - and the AHBA donors are neurotypical, not PD
patients, so this is baseline expression in a healthy brain, not diseased
tissue. The composite 9-gene signal picking up the expected pattern even
though the single most disease-relevant gene doesn't is itself worth
sitting with, rather than something to bury by only reporting the multi-gene
result.

For Huntington's, the direction reversal is worth taking seriously rather
than writing off as noise. One read: HD's risk-gene list here comes from a
progression-modifier GWAS (genes that affect how fast the disease progresses
once someone already has the HTT mutation), not a risk GWAS. Those are
answering a different biological question than "why does damage localize to
the striatum" - modifier genes are mostly DNA-repair pathway genes (the
Moss et al. hits are dominated by mismatch-repair-adjacent genes, plausibly
related to somatic CAG-repeat instability), which have no obvious reason to
be regionally enriched in the striatum in the first place. This is arguably
a mismatch between the available data and the question, more than it's
evidence against the underlying hypothesis for HD specifically.

## 5. Limitations

- The substantia nigra region is a coordinate-based reconstruction (6mm
  radius around real tissue-sample coordinates), not a proper anatomical
  segmentation. The radius wasn't systematically justified.
- Single-gene (SNCA) validation didn't pass even after that fix - see
  Discussion for why I don't think that invalidates the multi-gene result,
  but it means the originally-planned "does the known result replicate"
  checkpoint doesn't have a clean answer.
- AHBA is 6 donors, skewed older/male, and is healthy-brain baseline
  expression, not diseased tissue. 5 of 31 risk genes had no probe coverage
  at all (including CR1, a well-established AD gene) and were dropped.
- Huntington's gene list comes from a progression-modifier GWAS in HTT
  carriers, not a case/control risk GWAS - HD is monogenic, so this doesn't
  exist. This may partly explain the HD null/reversed result (see
  Discussion) rather than it being evidence against the hypothesis itself.
- Three of the 31 genes picked up by the automated GWAS-to-gene mapping are
  pseudogene annotations (FNTAP1, PMS2P1, ATG3P1), kept in rather than
  hand-corrected after seeing results.
- Correlational, not causal - genes being expressed where damage happens
  doesn't establish that they cause the damage there.
- Tested one statistic (target-region mean vs. rest-of-brain mean) and one
  direction (enrichment). Didn't check whether the Huntington's result is
  sensitive to gene weighting or a different summary statistic.

## 6. Code

Code and data: (github link once it's up) - see README.md for exact rerun
instructions.

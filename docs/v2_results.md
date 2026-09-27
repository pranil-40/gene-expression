# v2 results

Run 2026-09-24 to 25 under the plan frozen in `docs/analysis_plan.md`.
Every number here comes from `results/`. The tables were generated from the
result CSVs.

## Summary

Under the pre-specified primary analysis, neither the Parkinson's nor the
Alzheimer's GWAS gene set shows enrichment in the region that disease damages
early. Both are "not supported" by the decision rule fixed before running.

Parkinson's enrichment does appear under several alternative choices:
- coordinate-based substantia nigra masks, which the anatomy audit shows are mostly other midbrain and pontine tissue
- the v1 gene list
- turning off abagen's within-sample normalization (a post hoc choice)

Whether a Parkinson's result appears at all depends on those analysis choices. The v1 finding (spatial p = 0.005) doesn't survive the corrections.

The Huntington's gene sets show lower expression in caudate and putamen than elsewhere in all six donors. That is a descriptive observation. The one-sided enrichment test can't establish depletion.

## Primary analysis (pre-specified)

Region set: 41 left-hemisphere Desikan-Killiany regions, brainstem, and
substantia nigra pars compacta (SNc) and pars reticulata (SNr) defined by
AHBA anatomical labels (44 regions). Gene sets: top 10 independent GWS loci,
one protein-coding gene per locus (`data/gene_lists/sets/`).

| Disease | Target | Genes used | Target rank | z | Spatial p | Gene-set p | CAMERA p | BH-adjusted (spatial, gene-set) | Donors with T > 0 | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| Parkinson's | SNc | 9 of 10 (CRHR1 not in AHBA) | 7 of 44 | +1.45 | 0.145 | 0.221 | 0.163 | 0.383, 0.383 | 4 of 6 | not supported |
| Alzheimer's | entorhinal | 9 of 10 (CR1 not in AHBA) | 10 of 44 | +0.64 | 0.383 | 0.336 | 0.312 | 0.383, 0.383 | 4 of 6 | not supported |

z is a descriptive standardized contrast between regions of one map
(target minus the mean of the other regions, divided by their SD). It is not
a Cohen's d.

### Robustness of the primary result

Donor-level contrasts T(d), computed within each donor
(`results/v2/donor_contrasts.csv`):

| Donor | PD samples in SNc | PD T(d) | AD samples in entorhinal | AD T(d) |
|---|---|---|---|---|
| 9861 | 4 | +0.049 | 4 | -0.076 |
| 10021 | 4 | +0.068 | 2 | +0.124 |
| 12876 | 3 | -0.027 | 2 | +0.061 |
| 14380 | 3 | +0.117 | 2 | +0.021 |
| 15496 | 4 | +0.094 | 2 | -0.014 |
| 15697 | 2 | -0.079 | 3 | +0.035 |

Leave-one-out (`results/v2/leave_one_out.csv`, 5,000 surrogates each):

- **PD:** dropping any one gene keeps T above 0. The spatial p ranges from 0.014 (without MCCC1) to 0.474 (without KRTCAP2). Dropping donor 12876 gives p = 0.049; dropping any other donor gives p = 0.065 to 0.255. A single gene or donor can move the PD result across 0.05. That is fragility, and it doesn't count as support.
- **AD:** dropping TREM2 makes T slightly negative. Every leave-one-out spatial p is 0.25 or higher.

## Sensitivity, exploratory and post hoc runs

| Kind | Gene set | SN / region definition | Normalization | Genes used | Target rank | z | Spatial p | Gene-set p | CAMERA p | Donors T>0 |
|---|---|---|---|---|---|---|---|---|---|---|
| primary | alzheimers | dk_ontology_sn | default | 9/10 | 10/44 | +0.64 | 0.383 | 0.336 | 0.312 | 4/6 |
| primary | parkinsons | dk_ontology_sn | default | 9/10 | 7/44 | +1.45 | 0.145 | 0.221 | 0.163 | 4/6 |
| secondary: gene set | alzheimers_all_gws | dk_ontology_sn | default | 70/82 | 6/44 | +0.84 | 0.368 | 0.303 | 0.342 | 5/6 |
| secondary: gene set | parkinsons_all_gws | dk_ontology_sn | default | 79/92 | 37/44 | -0.95 | 0.764 | 0.330 | 0.442 | 3/6 |
| secondary: region definition | alzheimers | cit168 | default | 9/10 | 11/44 | +0.61 | 0.309 | 0.342 | 0.323 | 4/6 |
| secondary: region definition | alzheimers | sn6mm | default | 9/10 | 9/43 | +0.75 | 0.280 | 0.310 | 0.248 | 4/6 |
| secondary: region definition | parkinsons | cit168 | default | 9/10 | 1/44 | +2.77 | 0.023 | 0.229 | 0.145 | 3/3 |
| secondary: region definition | parkinsons | sn6mm | default | 9/10 | 2/43 | +2.15 | 0.058 | 0.159 | 0.043 | 4/5 |
| exploratory: literature-named genes | alzheimers_named | dk_ontology_sn | default | 9/10 | 9/44 | +0.44 | 0.407 | 0.380 | 0.370 | 4/6 |
| exploratory: literature-named genes | parkinsons_named | dk_ontology_sn | default | 10/10 | 26/44 | +0.14 | 0.469 | 0.329 | 0.336 | 3/6 |
| post hoc: normalization | alzheimers | dk_ontology_sn | gene_norm_only | 9/10 | 5/44 | +1.09 | 0.307 | 0.359 | 0.430 | 5/6 |
| post hoc: normalization | huntingtons_htt | dk_ontology_sn | gene_norm_only | 1/1 | 36/44 | -1.14 | 0.899 | 0.792 | 0.772 | 0/6 |
| post hoc: normalization | huntingtons_modifiers | dk_ontology_sn | gene_norm_only | 9/9 | 38/44 | -1.30 | 0.831 | 0.815 | 0.770 | 1/6 |
| post hoc: normalization | parkinsons | dk_ontology_sn | gene_norm_only | 9/10 | 2/44 | +3.23 | 0.023 | 0.076 | 0.153 | 6/6 |
| exploratory: v1 gene list | v1_alzheimers | dk_ontology_sn | default | 8/10 | 7/44 | +0.87 | 0.282 | 0.215 | 0.214 | 5/6 |
| exploratory: v1 gene list | v1_huntingtons | dk_ontology_sn | default | 9/11 | 29/44 | -0.47 | 0.580 | 0.562 | 0.530 | 1/6 |
| exploratory: v1 gene list | v1_parkinsons | dk_ontology_sn | default | 9/10 | 3/44 | +2.32 | 0.040 | 0.130 | 0.126 | 5/6 |
| exploratory: HD | huntingtons_htt | dk_ontology_sn | default | 1/1 | 37/44 | -1.07 | 0.893 | 0.787 | 0.764 | 0/6 |
| exploratory: HD | huntingtons_modifiers | dk_ontology_sn | default | 9/9 | 37/44 | -1.49 | 0.891 | 0.841 | 0.806 | 0/6 |

How to read the PD runs that reach p < 0.05:

- **`cit168`.** The CIT168 SNc region holds 6 samples from 3 donors. Only one of them is labelled SNc; the others include red nucleus, hypothalamus and VTA. Its low p describes that mixed tissue, and the region isn't substantia nigra.
- **`sn6mm`.** This v1-style mask holds 63 samples, 16 of them labelled SN. The rest are mostly pontine nuclei, VTA, central gray and red nucleus.
- **`v1_parkinsons`.** The v1 list counts the 17q21 (MAPT) locus three times (KANSL1, LINC02210-CRHR1, SPPL2C). v2 clumps it to one gene.
- **`gene_norm_only`.** This was run after seeing that abagen's within-sample normalization moves SNCA's SNc rank from 3 to 20. That makes it post hoc, and it can't be used as a primary result. It does show the PD answer depends on a preprocessing default most papers never report changing.
- **None of these runs clears the gene-set null.** The gene-set p is 0.076 or higher in all four, so even where the spatial pattern is strong, a random 9-gene set of comparable size often does as well.

## SNCA checkpoint

`results/snca_checkpoint/summary.txt`. The target was Liscovitch & French
2014 (SNc 2nd of 105 structures, SNr 88th).

| Version | SNc rank | SNr rank | Reproduced? |
|---|---|---|---|
| A: raw log2, all Allen SNCA probes, 103 left structures sampled in all 6 donors | 1 of 103 | 80 of 103 | yes |
| B: raw log2, abagen's selected probe | 1 of 92 | 74 of 92 | yes |
| C: v2 region matrix, abagen default normalization | 20 of 44 | 40 of 44 | n/a (different comparison) |

SNc is above SNr in every donor (by 0.39 to 2.07 log2 units). The post hoc
normalization diagnostic isolates the cause of C: SNc ranks 3rd of 44 with no
normalization or gene-only normalization, and 20th or 21st whenever
within-sample SRS normalization is applied. SNc samples have a wider spread
of expression across genes (median interquartile range 3.30 log2 units)
than cortical (3.16) or other subcortical (3.19) samples. Within-sample
scaling compresses their highly expressed genes more.

The v1 draft explained the failed checkpoint as protein pathology versus
mRNA level. That explanation is withdrawn. The failure came from the
region definition (mostly non-SN tissue) and the normalization.

## External check: GTEx v10

`results/gtex/`. log2(TPM + 1), z-scored per gene across the 13 GTEx brain
tissues.

| Gene set | Target | Genes used | Target rank | T | Gene-set p |
|---|---|---|---|---|---|
| parkinsons | substantia nigra | 10 | 8 of 13 | -0.193 | 0.310 |
| parkinsons_all_gws | substantia nigra | 91 | 9 of 13 | -0.245 | 0.229 |
| huntingtons_htt | caudate + putamen | 1 | 9 of 13 | -0.768 | 0.657 |
| huntingtons_modifiers | caudate + putamen | 9 | 9 of 13 | -0.828 | 0.965 |

GTEx doesn't show PD gene-set enrichment in substantia nigra either.
SNCA itself is below average in GTEx substantia nigra (z = -0.52). GTEx
samples are bulk dissections of the whole nigra from older donors. AHBA's
SNc samples are anatomically labelled pars compacta. The GTEx result is a
check under different conditions and can't count as a direct replication.
AD isn't tested because GTEx has no entorhinal cortex.

## Anatomy audit

`results/sn_audit/anatomy_summary.txt`.

- **52 AHBA samples are labelled substantia nigra**, 37 of them left
  hemisphere. Every donor has 3 or 4 left SNc samples (2 to 4 after
  abagen's anatomy-coordinate consistency check drops two).
- **Their MNI coordinates don't line up with the SN in template space.**
  Only 2 of 52 fall inside the CIT168 SN. Corrected coordinates spread from
  x = -23 to +25 mm and z = -29 to -5 mm.
- **The DK atlas's "brainstem" label doesn't cover most of the SN.** 856 of
  the 915 CIT168 SN voxels are unlabelled in DK.
- **v1's 6 mm mask was mostly other tissue.** It held 63 samples, 16
  labelled SN. The rest were mostly pontine nuclei, VTA, central gray, red
  nucleus and trochlear nucleus.
- **DK regions mix neighbouring tissue too.** The DK "putamen" region holds
  77 samples, only 37 of them labelled putamen. The others include claustrum
  (14), basal nucleus of Meynert and external globus pallidus. The DK
  "caudate" region holds 47 samples, 37 labelled caudate head or body.

## Gene audit

`data/gene_lists/gene_audit.csv`, `data/gene_lists/probe_audit.csv`.

- **Clumping.** PD: 106 GWS variants form 98 independent loci. 8 variants
  merged, including SNCA's secondary signals and the 17q21 (MAPT) region
  variants. AD: 87 variants, 87 loci.
- **Nearest-gene mapping.** At three top-10 loci the nearest protein-coding
  gene isn't the gene usually credited: KRTCAP2 at the GBA1 locus, CRHR1 at
  the MAPT locus, EED at the PICALM locus. Swapping in GBA1, MAPT and PICALM
  (the `*_named` runs) makes both results weaker.
- **CR1 (AD)** has three probes but fails the intensity filter. It is above
  background in 27% of samples, against the 50% threshold. The v1 draft
  called this a probe coverage gap, which was wrong.
- **CRHR1 (PD)** has four Allen probes, but abagen's probe reannotation drops
  all of them.
- **Renamed genes.** GBA1, CYRIB and H2BC13 appear in AHBA under their older
  names. The analysis matches them by Entrez ID.
- **Pseudogenes.** v1's pseudogene picks (FNTAP1, PMS2P1, ATG3P1) can't occur
  under the v2 rule.

## Spatial null diagnostics

`results/v2/variograms/`, fig5. BrainSMASH surrogates are less smooth than the
observed maps at short distances: the observed variogram sits at the lower
edge of the surrogate band below about 15 mm. With 44 regions the fitted null
is probably somewhat too narrow, so the spatial p-values are, if anything,
liberal. That strengthens the "not supported" primary result and weakens the
few sensitivity runs under 0.05.

## Positive controls and power (post hoc)

Added after every result above was known, under
`docs/analysis_plan_amendment_1.md`, written before these were run.
`src/power_check.py`, `results/v2/power/`, fig8. None of it changes a
verdict.

**Positive controls.** Marker genes for cell types that are known to sit in
the target region, run through the primary pipeline unchanged:

| Set | Target | Genes used | Target rank | z | Spatial p | Gene-set p | CAMERA p | Donors with T > 0 |
|---|---|---|---|---|---|---|---|---|
| dopaminergic neuron markers | SNc | 8 of 10 (EN1, LMX1B not in AHBA) | 1 of 44 | +7.54 | 0.0012 | 0.0001 | 0.0021 | 6 of 6 |
| striatal projection neuron markers | caudate + putamen | 9 of 10 (ADORA2A not in AHBA) | 1 of 44 | +7.37 | 0.0002 | 0.0001 | 0.0007 | 6 of 6 |

Both pass every test in every donor. The pipeline finds strong regional
specificity when it is there, including in SNc with only 2 to 4 samples per
donor.

**How large would the primary sets' contrast have had to be?** Shifting only
the target region of the observed PD and AD score maps
(`detection_threshold.csv`):

| Gene set | Observed z | z needed, spatial null | z needed, gene-set null | z needed, both |
|---|---|---|---|---|
| Parkinson's (SNc) | +1.45 | 2.50 | 4.75 | 4.75 |
| Alzheimer's (entorhinal) | +0.64 | 3.75 | 2.25 | 3.75 |

The two nulls bind in different places. SNc is transcriptionally distinct
from the rest of the brain, so random gene sets often show large contrasts
there, and the gene-set null is the hard one for PD: the PD set would need
to rank first by a wide margin (z of about 4.75) before random 9-gene sets
stop matching it. Entorhinal cortex sits among smooth cortical neighbours,
so for AD the spatial null is the hard one. The dopaminergic markers
(z = 7.54) clear both easily.

**Power for random 9-gene sets** (`random_set_power.csv`; nominal one-sided
alpha of 0.05, so upper bounds under the BH correction the primary analysis
used):

| z given to the set | PD gene-set | PD spatial | PD both | AD gene-set | AD spatial | AD both |
|---|---|---|---|---|---|---|
| 1.0 | 0.00 | 0.00 | 0.00 | 0.04 | 0.00 | 0.00 |
| 1.5 | 0.03 | 0.00 | 0.00 | 0.25 | 0.02 | 0.00 |
| 2.0 | 0.13 | 0.12 | 0.05 | 0.55 | 0.28 | 0.20 |
| 2.5 | 0.31 | 0.28 | 0.12 | 0.79 | 0.60 | 0.52 |
| 3.0 | 0.50 | 0.32 | 0.18 | 0.92 | 0.75 | 0.72 |
| 3.5 | 0.64 | 0.32 | 0.22 | 0.98 | 0.87 | 0.85 |
| 4.0 | 0.81 | 0.32 | 0.27 | 1.00 | 0.93 | 0.93 |
| 5.0 | 0.96 | 0.32 | 0.32 | 1.00 | 0.98 | 0.98 |

Gene-set columns use 1,000 draws per level (standard error 0.016 or less).
Spatial and joint columns use 60 draws (standard error up to 0.065).

For AD, both tests together reach 80% power at about z = 3.5. For PD, the
gene-set null reaches 80% at z = 4, but the spatial null levels off at about
a third of draws however large the shift. The likely reason is geometry. The
SNc and SNr centres (CIT168) are 3.2 mm apart and at least 18 mm from every
other region, so the variogram-matched surrogates often place a peak on that
isolated pair, and a surrogate contrast then matches or exceeds the shifted
one. For a small subcortical target with one close neighbour, this spatial
null limits power in a way a cortical target doesn't. The observed PD set
still reached spatial p below 0.05 from z = 2.5 in the threshold analysis,
so the ceiling depends on the shape of the map.

**What this means.** At the observed contrasts (PD z = 1.45, AD z = 0.64),
power for random sets of the same size is close to zero on every test. The
design can detect a gene set that marks its target region about as strongly
as a cell-type marker set does. It can't distinguish a modest enrichment
(z of 1 to 2, a target in the top few regions) from chance. The v2 result is
therefore "no large enrichment", and a modest enrichment of PD or AD risk
genes in their target regions is not ruled out. A study that could detect
one would need more samples in the target regions than the 2 to 4 per donor
AHBA has, or a finer-grained dataset such as single-nucleus RNA-seq.

## What can and can't be claimed

| Claim | Status |
|---|---|
| PD GWAS genes are enriched in SNc in healthy adult brain | Not supported under the pre-specified analysis. Appears only under anatomically incorrect region definitions, the v1 gene list, or a post hoc normalization choice. |
| AD GWAS genes are enriched in entorhinal cortex | Not supported in any run. |
| HD genes are depleted in striatum | Not tested. T(d) is below 0 in all 6 donors, but the planned test was one-sided for enrichment. |
| SNCA is highest in AHBA SNc among left-hemisphere structures | Reproduced (Liscovitch & French 2014), using raw expression. |
| The answer depends on analysis choices | Supported descriptively by fig7 and the table above. |
| The pipeline can detect regional enrichment where it exists | Supported: both cell-type marker sets rank first with every p below 0.003 (post hoc positive controls). |
| PD and AD risk genes have no regional enrichment at all | No. Only large effects were detectable; a modest enrichment isn't ruled out (post hoc power analysis). |
| This explains why the disease damages that region | No. The study measured healthy-brain expression. It didn't measure damage. |
| This could guide clinical screening | No. Nothing individual-level was tested. |

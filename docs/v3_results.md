# v3 results: cell types inside the damaged region

Run 2026-09-27 under `docs/analysis_plan_v3.md`, which was frozen before any
expression value was computed. Every number comes from `results/v3/`.
Figures: `results/v3/figures/fig9_celltypes.png`, `fig10_dilution.png`.

## Summary

- **Alzheimer's.** The AD risk genes are specifically expressed in microglia
  inside entorhinal cortex. The pre-specified test is supported: microglia
  rank first of six cell classes, p = 0.006 (0.011 after correction), and
  the pattern holds in all 3 donors. With all 82 genome-wide significant AD
  loci, the signal is much stronger (z = 9.0).
- **Parkinson's.** The PD risk genes aren't specifically expressed in
  dopamine (DA) neurons, in either SN dataset. DA neurons rank third of
  seven, p = 0.22 in Kamath et al. 2022 and 0.15 in Siletti et al. 2023.
  The shift is small and consistent: all 8 Kamath control donors lean
  positive (z = 0.3 to 1.0), and none reaches significance. SNCA is the one
  strongly DA-leaning gene, and without it the signal disappears.
- **Dilution.** DA neurons are 0.32% of substantia nigra nuclei and carry
  1.5% of its RNA. For DA marker genes, 13% to 88% of the tissue's signal
  comes from those few neurons, which is why bulk AHBA tissue detected them
  (v2 amendment 1, z = 7.5). For the PD risk genes, only 0.3% to 4.7%
  does. Even a real DA-neuron-specific change in these genes would barely
  move a bulk SN measurement.
- **Positive controls.** All four pass, with z from 8 to 21.

## Primary tests

| Disease | Dataset | Target class | Genes used | Target rank | z | p | BH p | Top class | Donors with z > 0 | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| Parkinson's | Kamath SNpc, 8 control donors | DA neurons | 9 of 10 (CRHR1 not in data) | 3 of 7 | +0.70 | 0.222 | 0.222 | OPCs | 8 of 8 | not supported |
| Alzheimer's | Siletti entorhinal (medial + lateral), 3 donors | microglia | 10 of 10 | 1 of 6 | +2.95 | 0.006 | 0.011 | microglia | 3 of 3 | supported |

Score: the mean over the set's genes of each gene's share of expression in
the target class (EWCE specificity). Null: 10,000 random sets, matched gene
by gene on expression decile.

**Which genes carry each result** (`gene_specificity.csv`, fig9):

- **AD:** TREM2 (100% microglial) and MS4A4A (99%) carry most of the signal.
  BIN1 (43%) and EED (37%) lean microglial. The rest point elsewhere: CLU
  to astrocytes, ABCA7 and SPDYE3 to neurons, CR1 and APH1B to vascular
  cells. Leaving out TREM2 or MS4A4A raises p to about 0.05. Microglia stay
  first in every leave-one-gene-out run.
- **PD:** SNCA is 38% DA-neuronal and NSF 27%. LRRK2 (42%) and HIP1R (51%)
  lean toward OPCs, and LRRK2 also toward microglia (32%). Leaving out SNCA
  sends p to 0.50 and DA neurons to sixth place.

## Secondary and exploratory

| Kind | Dataset | Gene set | Target | z | p | Target rank | Top class |
|---|---|---|---|---|---|---|---|
| replication | Siletti SN | PD top 10 | DA neurons | +1.03 | 0.153 | 3 of 7 | OPCs |
| oligodendrocytes (Agarwal et al. 2020) | Kamath SNpc | PD top 10 | oligodendrocytes | +0.07 | 0.420 | 2 of 7 | OPCs |
| oligodendrocytes | Siletti SN | PD top 10 | oligodendrocytes | -0.48 | 0.653 | 4 of 7 | OPCs |
| all GWS loci | Kamath SNpc | PD, 86 genes | DA neurons | +0.89 | 0.184 | 3 of 7 | microglia |
| all GWS loci | Siletti SN | PD, 92 genes | DA neurons | +1.58 | 0.064 | 2 of 7 | microglia |
| all GWS loci | Siletti entorhinal | AD, 82 genes | microglia | +9.01 | 0.0001 | 1 of 6 | microglia |
| DA subtypes | Kamath DA neurons | PD top 10 | SOX6_AGTR1 | +0.47 | 0.287 | 7 of 10 | CALB1_RBP4 |
| DA subtypes | Kamath DA neurons | PD, 86 genes | SOX6_AGTR1 | +1.30 | 0.111 | 7 of 10 | SOX6_DDT |
| HD | Siletti striatum | HTT | MSNs | -0.16 | 0.417 | 1 of 7 | MSNs |
| HD | Siletti striatum | onset modifiers | MSNs | -0.46 | 0.641 | 6 of 7 | microglia |

- The PD genes don't favour the SOX6_AGTR1 subtype that Kamath et al.
  found is lost in PD.
- With all PD loci, microglia score highest. That is a descriptive
  observation, and it wasn't a planned hypothesis.
- The HD genes, the causal gene and DNA-repair onset modifiers, aren't
  specific to medium spiny neurons.
- In Siletti SN, only one donor had the 20 DA neurons needed for a
  donor-level check.

## Positive controls

| Set | Dataset | Target | z | p | Rank |
|---|---|---|---|---|---|
| dopaminergic markers | Kamath SNpc | DA neurons | +20.8 | 0.0001 | 1 of 7 |
| dopaminergic markers (circular: SLC6A3 defines the class) | Siletti SN | DA neurons | +16.7 | 0.0001 | 1 of 7 |
| microglia markers | Siletti entorhinal | microglia | +10.0 | 0.0001 | 1 of 6 |
| striatal projection neuron markers | Siletti striatum | MSNs | +8.2 | 0.0001 | 1 of 7 |

In Siletti SN, the DA rule selected one cluster: 73% of its 192 cells
express SLC6A3, and no other neuronal cluster exceeds 14%
(`da_rule_clusters.csv`).

## Dilution

From `dilution.csv`; Siletti SN, which was dissected without cell sorting.

| | Share of the gene's SN signal from DA-neuron nuclei | DA-neuron mean / SN-wide mean |
|---|---|---|
| DA markers (10 genes) | 13% (NR4A2) to 88% (SLC6A3) | 7x to 178x |
| PD risk genes (10 genes) | 0.3% (KRTCAP2) to 4.7% (SNCA) | 0.2x to 3.1x |

DA neurons are 192 of 59,088 nuclei (0.32%) and carry 1.5% of all UMIs.
Bulk SN tissue reflects DA neurons only for genes those neurons dominate.

## Post hoc check: microglia in bulk tissue

This was added after the results above, to ask whether microglia being
spread through the brain could explain the v2 AD null. The microglia
markers were run through the unchanged v2 regional pipeline
(`bulk_link_microglia_markers.csv`). They are moderately higher in
entorhinal cortex (z = +1.84, rank 2 of 44, spatial p = 0.043, gene-set
p < 0.001) and in SNc (z = +1.77, rank 3 of 44).

So the explanation isn't that microglia are evenly spread. A better
account comes from the gene table above. The AD set's microglial signal
rests mostly on two genes, and the other genes point to astrocytes,
neurons and vascular cells. Averaged over a region, that mixture washes
out.

## What can and can't be claimed

| Claim | Status |
|---|---|
| AD risk genes (top 10 loci) are specifically expressed in microglia in entorhinal cortex | Supported (pre-specified). This replicates published work with this project's gene sets and null. |
| PD risk genes (top 10 loci) are specifically expressed in DA neurons of the SN | Not supported in either dataset. A small, consistent positive lean isn't significant. |
| PD genetic risk isn't linked to DA neurons | No. This test measures expression specificity of nearest genes. Heritability methods over many more genes (Kamath et al. 2022; Agarwal et al. 2020) reach a different kind of conclusion. |
| Bulk regional data can't reflect a DA-neuron-specific signal in these genes | Supported descriptively: at most 4.7% of each PD gene's SN signal comes from DA neurons. |
| Microglia are evenly spread, which explains the v2 AD null | No (post hoc check). Microglia markers are moderately higher in entorhinal cortex. |

## Caveats

- Single-nucleus data miss cytoplasmic RNA. AHBA measured whole-tissue RNA.
- Siletti has 3 donors. Kamath's control donors are mostly aged 75 to 92.
- Nearest-gene mapping (KRTCAP2 at the GBA1 locus, CRHR1 at MAPT, EED at
  PICALM) limits every gene-level analysis in this project.
- The primary hypotheses came from published studies, so the AD result is a
  replication of published work.

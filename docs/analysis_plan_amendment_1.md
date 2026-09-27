# Analysis plan, amendment 1: positive controls and power

Written 2026-09-25, 12:10 CDT. Every v2 result in `docs/v2_results.md` was
already known. Everything below is post hoc. It adds context to the v2
verdicts and doesn't change any of them.

## Why

The frozen plan (`docs/analysis_plan.md`) has no positive control and no
power estimate. A "not supported" verdict only means something if the
design could have detected an effect of a plausible size. Two questions:

1. Does the pipeline find regional enrichment where it is known to exist?
2. How large would a gene set's target contrast have to be before the
   planned tests detect it?

## A. Positive controls

Gene sets with well-known regional specificity, run through exactly the
primary pipeline: `dk_ontology_sn` regions, default normalization, the same
spatial null, gene-set null and CAMERA test, with the same numbers of
surrogates and random sets, plus donor-level contrasts.

| Set | Target | Genes |
|---|---|---|
| `pc_dopaminergic` | substantia nigra pars compacta | TH, SLC6A3, SLC18A2, DDC, KCNJ6, NR4A2, EN1, FOXA2, LMX1B, PITX3 |
| `pc_striatal_msn` | caudate + putamen | PPP1R1B, DRD1, DRD2, PENK, TAC1, ADORA2A, GPR88, RGS9, PDE10A, RASD2 |

The marker choice follows standard usage: Poulin et al. 2020 (Trends
Neurosci 43:155, PMID 32101709) for midbrain dopamine neurons, and Gokce
et al. 2016 (Cell Rep 16:1126, PMID 27425622) for striatal projection
neurons. EN1, LMX1B and ADORA2A aren't in the AHBA gene universe after
abagen's filtering, which was checked before any expression value was
looked at. That leaves 8 and 9 genes.

Expected: each set's target ranks first, and both spatial and gene-set
p < 0.05. If a set fails, the pipeline can't detect even strong regional
specificity, and the v2 null results would say nothing.

There is no Alzheimer's positive control. No established marker set picks
out entorhinal cortex in bulk tissue.

## B. Detection threshold for the actual primary gene sets

For each primary set (PD in SNc, AD in entorhinal), add a constant to the
set's score in the target region so that its standardized contrast equals
z*, where z* runs from 1.0 to 5.0 in steps of 0.25. Everything else in the
map stays as observed. Recompute:

- the spatial p, with 2,000 BrainSMASH surrogates of the shifted map
- the gene-set p, against 10,000 random sets of the same size (unshifted)

Report the smallest z* at which each test, and both together, give
p < 0.05. Only the crossing point matters, so 2,000 surrogates is enough.

## C. Power with random gene sets

k = 9 genes drawn at random from the gene universe. For z* in
{1, 1.5, 2, 2.5, 3, 3.5, 4, 5}, shift each drawn set's target score to z*
as in B. Then:

- gene-set null: 1,000 draws per level, each tested against the 10,000-set null
- spatial null: 60 draws per level, 500 surrogates each (BrainSMASH is the slow step)

Power is the share of draws with p < 0.05. The primary analysis applies
Benjamini-Hochberg across four p-values, so these nominal-level numbers are
upper bounds. The Monte Carlo standard error is reported with each value.

Seeds come from `common.seed_for`. Outputs go to `results/v2/power/`.

## Interpretation, fixed now

- If both positive controls pass, the pipeline can detect strong regional
  enrichment.
- If the z* needed for 80% power is above the observed PD z (+1.45) and AD
  z (+0.64), the v2 results are described as "no large enrichment
  detected". A modest enrichment isn't ruled out.

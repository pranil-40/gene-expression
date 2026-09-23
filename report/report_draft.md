<!--
notes to self: don't fill in results/figures/conclusions until
expression_pipeline.py -> visualize.py have actually run on real data and
validation_checkpoint.py says PASS. see README status section.
-->

# Do Disease Risk Genes Concentrate Where the Disease Actually Strikes?

Pranil
draft skeleton, waiting on real pipeline run

## 1. Research question

(pull from README, that framing already works)

## 2. Method

- gene lists: which studies, how genes were picked (fetch_gwas_genes.py)
- expression data: AHBA/abagen, atlas, normalization settings used
- target region for each disease + the brainstem/SN caveat
- validation checkpoint result (SNCA)
- stats: regional comparison + BrainSMASH permutation, FDR correction

## 3. Results

TBD - heatmaps go here, results table, p-values per disease

## 4. Discussion

TBD - what does the pattern across the three diseases actually mean, and what does it not mean

## 5. Limitations

- no separate SN label in the atlas, brainstem proxy dilutes the PD result
- AHBA is only 6 donors, skewed demographics
- HD gene list is from a progression-modifier GWAS, not a risk GWAS (HD is monogenic)
- 3 of 31 genes are pseudogene annotations that got picked up by the mapping step
- correlation, not causation
- (add whatever else actually turns out to matter once results are in - don't go easy on this section)

## 6. Code

github link here once it's up - see README for how to rerun everything

<!--
STATUS: SKELETON ONLY — do not fill in results, figures, or conclusions
until src/expression_pipeline.py through src/visualize.py have actually
been run on real AHBA data and src/validation_checkpoint.py has PASSED.
See repo README.md "Status" section.
-->

# Neurogenomic Regional Vulnerability: Do Disease Risk Genes Concentrate Where the Disease Strikes?

**Author:** Pranil
**Status:** Draft skeleton — pending real pipeline run

## 1. Research Question

*(See README.md — carry over the framing almost verbatim, it's already tight.)*

## 2. Method

- Gene lists: source studies, selection criteria, exact rule for picking top genes (cite `src/fetch_gwas_genes.py`).
- Expression data: AHBA via abagen, atlas used, normalization/aggregation settings.
- Target-region mapping and why (including the Parkinson's brainstem-proxy caveat).
- Validation checkpoint (SNCA) and its result.
- Statistical test: regional comparison + BrainSMASH spatial permutation, FDR correction.

## 3. Results

*(TBD — insert the three regional heatmaps, the results summary table, and the per-disease p-values once the real run completes.)*

## 4. Discussion

*(TBD — what the pattern across the three diseases actually suggests, and what it doesn't.)*

## 5. Limitations

- Substantia nigra not separately resolvable in the DK+aseg atlas (brainstem proxy).
- AHBA reference: only 6 donors, demographic skew.
- Huntington's gene list comes from a progression-modifier GWAS in HTT carriers, not a risk GWAS (HD is monogenic).
- 3 of 31 mapped genes are pseudogene annotations, kept as-is (see README).
- Correlational, not causal.
- *(Add anything else that turns out to matter once real results are in — don't soften this section.)*

## 6. Data & Code Availability

Code and data: `github.com/<to-fill-in>` — see README.md for exact rerun instructions.

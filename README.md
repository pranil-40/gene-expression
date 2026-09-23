# Neurogenomic Regional Vulnerability

Parkinson's kills the substantia nigra. Alzheimer's begins in the entorhinal
cortex. Huntington's destroys the striatum. The genes that raise risk for
each sit in every cell of the body — so why does the damage always land in
one place?

This project tests whether the brain regions where each disease's GWAS risk
genes are most *expressed* are the same regions the disease goes on to
*destroy*, and whether that overlap survives a spatial null model or is just
a plot that looks convincing.

Original analysis on public data: [Allen Human Brain Atlas](https://human.brain-map.org/)
regional gene expression (via [`abagen`](https://abagen.readthedocs.io/)),
[GWAS Catalog](https://www.ebi.ac.uk/gwas/) risk genes, and
[BrainSMASH](https://brainsmash.readthedocs.io/) spatial permutation testing.

## Method, in order

1. **Gene lists** ([`src/fetch_gwas_genes.py`](src/fetch_gwas_genes.py)) — top
   10 genes per disease from one landmark GWAS study each (picked by sample
   size/recency, not by anything found later in this pipeline), ranked by
   p-value, gene assigned via the GWAS Catalog's own SNP-to-gene genomic
   context. Exact studies and selection rule are documented in the script's
   docstring. Output: [`data/gene_lists/`](data/gene_lists/).
2. **Regional expression** ([`src/expression_pipeline.py`](src/expression_pipeline.py))
   — pulls AHBA microarray expression for those genes across the
   Desikan-Killiany(+subcortical) atlas (~83 regions) via `abagen`.
3. **Validation checkpoint** ([`src/validation_checkpoint.py`](src/validation_checkpoint.py))
   — confirms SNCA (the single strongest PD GWAS hit) is highly expressed
   near the substantia nigra proxy region *before* trusting anything else.
4. **Regional comparison** ([`src/regional_comparison.py`](src/regional_comparison.py))
   — per disease, mean risk-gene expression in the region it destroys vs.
   every other region.
5. **Spatial permutation test** ([`src/permutation_test.py`](src/permutation_test.py))
   — BrainSMASH surrogate maps (N=1000) that preserve spatial autocorrelation,
   used to get an honest p-value instead of naive shuffling; Benjamini-Hochberg
   FDR correction across the 3 diseases.
6. **Figures** ([`src/visualize.py`](src/visualize.py)) — one regional
   expression bar chart per disease plus a results summary table.

## Target-region mapping (decided before running anything, see `src/regional_comparison.py`)

| Disease | Destroyed region | Atlas proxy used |
|---|---|---|
| Parkinson's | Substantia nigra | **Brain-Stem** (DK+aseg has no separate SN label — the brainstem label anatomically contains it, along with pons/medulla; this dilutes any true signal and is a stated limitation) |
| Alzheimer's | Entorhinal cortex | `entorhinal` (direct DK match) |
| Huntington's | Striatum | `caudate` + `putamen` (direct DK+aseg match) |

## Running it

Everything here targets a local venv (`.venv/`) or Google Colab.

```bash
python -m venv .venv
.venv/Scripts/activate   # or `source .venv/bin/activate` on macOS/Linux
pip install -r requirements.txt

python src/fetch_gwas_genes.py        # fast, no large downloads
python src/expression_pipeline.py     # SLOW first run — downloads several GB of AHBA donor data
python src/validation_checkpoint.py   # must PASS before trusting the next two steps
python src/regional_comparison.py
python src/permutation_test.py
python src/visualize.py
```

In Colab: `!pip install -r requirements.txt` in the first cell, then run the
same scripts (or the equivalent notebook cells in `notebooks/`).

## Known limitations (see full report for the honest version)

- The Desikan-Killiany(+aseg) atlas cannot resolve the substantia nigra as
  its own region; the Parkinson's result uses "Brain-Stem" as a coarse proxy.
- The AHBA reference is only 6 donors, all male-skewed and older-adult —
  small and non-representative for a lifespan disease question.
- Huntington's is monogenic; there is no case/control "risk GWAS" for it in
  the PD/AD sense, so its gene list comes from a progression-modifier GWAS in
  confirmed HTT carriers, plus HTT itself added manually (documented in
  `src/fetch_gwas_genes.py`).
- Regional co-localization of risk-gene expression and disease damage is
  correlational, not causal.
- Three of the 31 automatically-mapped genes (FNTAP1, PMS2P1, ATG3P1) are
  pseudogene annotations nearest to their lead SNP per the GWAS Catalog's own
  mapping — kept as-is rather than hand-corrected after the fact (to avoid
  quietly re-tuning the gene list once results are visible), but likely
  contribute little/no real expression signal.

## Status

Scaffolded — pipeline code is complete and runs end-to-end on the fetched
gene lists. The AHBA data pull (`expression_pipeline.py`) has not yet been
run to completion (multi-GB download); results and figures are not yet
generated. See `report/` for the write-up, to be filled in once real results
exist.

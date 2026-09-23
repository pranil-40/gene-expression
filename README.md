# Neurogenomic Regional Vulnerability

Parkinson's kills the substantia nigra. Alzheimer's starts in the entorhinal
cortex. Huntington's destroys the striatum. The genes that raise risk for
each of these diseases are present in every cell in the body, not just the
region that ends up damaged. So why does the damage always land in the same
place?

That's the question this project tries to answer: are the brain regions
where a disease's GWAS risk genes are most expressed the same regions that
disease actually destroys, and does that hold up once you test it properly
instead of just eyeballing a bar chart.

Data sources: gene expression from the [Allen Human Brain Atlas](https://human.brain-map.org/)
(pulled with [abagen](https://abagen.readthedocs.io/)), risk genes from the
[GWAS Catalog](https://www.ebi.ac.uk/gwas/), significance testing with
[BrainSMASH](https://brainsmash.readthedocs.io/) spatial nulls.

## What the pipeline does

1. `src/fetch_gwas_genes.py` - pulls the top 10 risk genes per disease from one
   big GWAS study each (chosen by PMID before looking at anything downstream),
   ranked by p-value, gene assigned to each lead SNP through GWAS Catalog's
   own genomic-context lookup. Writes to `data/gene_lists/`.
2. `src/expression_pipeline.py` - gets AHBA expression for those genes across
   the Desikan-Killiany atlas (~83 regions) via abagen. This is the slow part.
3. `src/validation_checkpoint.py` - checks that SNCA (the #1 PD GWAS hit)
   actually shows up high near the substantia nigra proxy before trusting
   anything past this point.
4. `src/regional_comparison.py` - for each disease, mean risk-gene expression
   in its target region vs. everywhere else.
5. `src/permutation_test.py` - BrainSMASH spatial null (1000 surrogates) so
   the p-values actually mean something, plus BH-FDR correction across the
   three diseases.
6. `src/visualize.py` - the figures.

## Target region for each disease

| Disease | Region destroyed | What I used |
|---|---|---|
| Parkinson's | Substantia nigra | Brain-Stem (DK atlas doesn't break out the SN separately, this is the closest available label and it's not great - see limitations) |
| Alzheimer's | Entorhinal cortex | `entorhinal` (direct match) |
| Huntington's | Striatum | `caudate` + `putamen` |

## Running it

```bash
python -m venv .venv
.venv/Scripts/activate   # source .venv/bin/activate on mac/linux
pip install -r requirements.txt

python src/fetch_gwas_genes.py        # fast
python src/expression_pipeline.py     # slow, downloads several GB first time
python src/validation_checkpoint.py   # should say PASS
python src/regional_comparison.py
python src/permutation_test.py
python src/visualize.py
```

Or just open `notebooks/neurogenomic_regional_vulnerability.ipynb` in Colab -
same steps, `!pip install` instead.

## Limitations

Being upfront about the weak points instead of burying them:

- No separate substantia nigra label in the atlas I used, so the PD result
  is diluted by using the whole brainstem as a stand-in. If PD comes back
  non-significant this is the first thing I'd blame.
- AHBA is only 6 donors and skews older/male. Small reference for a
  lifespan-disease question.
- HD is monogenic, so there isn't really a "risk GWAS" for it like PD/AD
  have - I used a progression-modifier GWAS in confirmed HTT carriers
  instead and added HTT itself by hand (see `fetch_gwas_genes.py`).
- Three of the 31 genes the pipeline picked out (FNTAP1, PMS2P1, ATG3P1) are
  pseudogene annotations, not real protein-coding genes. Left them in rather
  than hand-fixing the list after seeing results - that would be exactly the
  kind of after-the-fact tuning this project is trying to avoid. They
  probably just add noise.
- Correlation, not causation. Genes being expressed where damage happens
  doesn't mean they cause the damage there.

## Status

Code's done and dry-run tested end-to-end on fake data. Real AHBA data pull
is running now - once that's done the rest of the pipeline runs on real
numbers and `report/report_draft.md` gets filled in.

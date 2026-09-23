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

## Result, short version

Ran this on real data. Parkinson's held up - substantia nigra is the #1
region out of 43 for average PD risk-gene expression, and it survives the
spatial permutation test (p=0.005, FDR-corrected p=0.015). Alzheimer's looked
promising on a raw comparison (entorhinal cortex was elevated, Cohen's d=1.5)
but did NOT survive the spatial null (FDR p=0.22) - the effect is
indistinguishable from what you'd expect from spatial autocorrelation alone.
Huntington's actually went the other way: caudate+putamen had *lower*
average risk-gene expression than the rest of the brain, not higher.

So: one disease supports the hypothesis, one doesn't hold up under proper
testing, one contradicts it. See `report/report_draft.md` for the full
writeup and what I think this actually means.

## What the pipeline does

1. `src/fetch_gwas_genes.py` - pulls the top 10 risk genes per disease from one
   big GWAS study each (chosen by PMID before looking at anything downstream),
   ranked by p-value, gene assigned to each lead SNP through GWAS Catalog's
   own genomic-context lookup. Writes to `data/gene_lists/`.
2. `src/build_sn_atlas.py` - the DK atlas's "brainstem" label doesn't separate
   out the substantia nigra, so this carves it out as its own region using
   real AHBA tissue-sample coordinates (see below).
3. `src/expression_pipeline.py` - gets AHBA expression for those genes across
   the modified atlas (~43 regions) via abagen. This is the slow part.
4. `src/validation_checkpoint.py` - checks that SNCA shows up high in
   substantia nigra. It doesn't, actually - see limitations, this turned out
   to be more interesting than a simple pass/fail.
5. `src/regional_comparison.py` - for each disease, mean risk-gene expression
   in its target region vs. everywhere else.
6. `src/permutation_test.py` - BrainSMASH spatial null (1000 surrogates) so
   the p-values actually mean something, plus BH-FDR correction across the
   three diseases.
7. `src/visualize.py` - the figures.

## Target region for each disease

| Disease | Region destroyed | What I used |
|---|---|---|
| Parkinson's | Substantia nigra | `substantia_nigra`, carved out of the DK atlas's brainstem blob (see `build_sn_atlas.py`) using real AHBA tissue-sample coordinates within a 6mm radius |
| Alzheimer's | Entorhinal cortex | `entorhinal` (direct match) |
| Huntington's | Striatum | `caudate` + `putamen` |

## Running it

```bash
python -m venv .venv
.venv/Scripts/activate   # source .venv/bin/activate on mac/linux
pip install -r requirements.txt

python src/fetch_gwas_genes.py        # fast
python src/expression_pipeline.py     # slow first time, downloads a few GB
python src/build_sn_atlas.py          # needs expression_pipeline.py's donor data cached first, then rerun expression_pipeline.py so it picks up the new atlas
python src/validation_checkpoint.py
python src/regional_comparison.py
python src/permutation_test.py
python src/visualize.py
```

(Yeah, `expression_pipeline.py` runs twice - once to get the donor data
cached and build the plain DK atlas result, then `build_sn_atlas.py` needs
that cached data to find real SN sample coordinates, then `expression_pipeline.py`
runs again and picks up the new atlas automatically if it finds the files
`build_sn_atlas.py` wrote. A bit clunky but it means the slow download only
happens once.)

Or just open `notebooks/neurogenomic_regional_vulnerability.ipynb` in Colab -
same steps, `!pip install` instead.

## Limitations

Being upfront about the weak points instead of burying them:

- The substantia nigra region isn't a real anatomical segmentation - it's a
  6mm-radius reconstruction around real tissue-sample coordinates from the
  AHBA donors, carved out of the DK atlas's brainstem label. Better than
  treating all of brainstem as one blob, but still an approximation, and the
  6mm radius was picked without much justification (see `build_sn_atlas.py`).
- SNCA alone does NOT rank high in substantia nigra (rank 38 of 43) even
  after fixing the region-resolution problem. The full 9-gene PD average
  does show up enriched there, and that result holds under the spatial null.
  I think this means single-gene expression level isn't the same thing as
  disease relevance - SNCA's role in PD is about protein aggregation in
  dopaminergic neurons, not necessarily about that region having the highest
  raw SNCA transcript count in a healthy brain. Worth being upfront that the
  validation checkpoint, as originally planned, doesn't cleanly pass.
- AHBA is only 6 donors and skews older/male. Small reference for a
  lifespan-disease question, and 5 of the 31 risk genes weren't in the AHBA
  data at all (probe coverage gaps, not a lookup bug) - notably CR1, an
  actual well-established AD gene.
- HD is monogenic, so there isn't really a "risk GWAS" for it like PD/AD
  have - I used a progression-modifier GWAS in confirmed HTT carriers
  instead and added HTT itself by hand (see `fetch_gwas_genes.py`).
- Three of the 31 genes the pipeline picked out (FNTAP1, PMS2P1, ATG3P1) are
  pseudogene annotations, not real protein-coding genes. Left them in rather
  than hand-fixing the list after seeing results - that would be exactly the
  kind of after-the-fact tuning this project is trying to avoid.
- Correlation, not causation. Genes being expressed where damage happens
  doesn't mean they cause the damage there.
- Only tested one direction (enrichment) with one statistic (target region
  mean minus rest-of-brain mean). Didn't explore whether a different summary
  statistic or gene-weighting scheme changes the Huntington's result.

## Status

Full pipeline run on real data, results and figures in `results/`.
`report/report_draft.md` still needs the results/discussion sections filled
in with real numbers.

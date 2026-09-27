# Revision log: v1 to v2

v1 is commit `220b838` (code) with its outputs archived in
`archive/v1_exploratory/`. v2 was built 2026-09-24 to 25 in response to an
external methods review of the v1 draft. The review's twelve points are
listed below with what changed. Decisions made after the plan was frozen
(`docs/analysis_plan.md`, 2026-09-24 19:10 CDT) are marked **post-plan**.

## Review points

**1. Narrow the question.** The question is now regional enrichment of GWAS
gene sets in healthy tissue (`docs/analysis_plan.md`). Nothing in v2
claims to explain degeneration or to support screening.

**2. SNCA checkpoint.** `src/snca_checkpoint.py` reproduces Liscovitch &
French 2014 exactly as they describe it: SNc ranks 1st of 103 structures and
SNr ranks 80th, and SNc is above SNr in all six donors. The v1 failure had two
causes. The v1 region was mostly non-SN tissue. abagen's within-sample
normalization also moves SNc from 3rd to 20th of 44 (post-plan diagnostic).
The v1 draft's protein-versus-mRNA explanation is withdrawn.

**3. The 6 mm substantia nigra mask.** `src/sn_anatomy_audit.py` looked at
anatomy only, before any expression was examined. Only 2 of 52 SN-labelled
samples have coordinates inside the CIT168 SN. The v1-style mask held 63
samples, 16 labelled SN. The DK "brainstem" label misses most of the SN.

- Primary definition: AHBA anatomical labels (SNc and SNr as separate
  regions).
- Sensitivity: CIT168 atlas (`src/build_atlases.py`).
- Exploratory: the 6 mm mask, rebuilt on the corrected coordinates abagen
  actually matches on (v1 used the uncorrected ones).

**4. Gene selection.** `src/fetch_gwas_genes.py` rebuilt:

- genome-wide significant loci only
- LD-plus-distance clumping (Ensembl, 1000 Genomes EUR)
- one protein-coding gene per locus from GENCODE v50, readthrough genes
  excluded
- a full audit table with every association's decision
  (`data/gene_lists/gene_audit.csv`)
- `src/gene_probe_audit.py`, which separates missing genes into no probe,
  dropped by reannotation, and failed intensity filter, and matches renamed
  genes by Entrez ID

Pseudogene picks can no longer happen, and the 17q21 (MAPT) locus no longer
contributes three genes. Two post-plan items:

- the switch from Catalog genomic contexts to Ensembl positions plus local
  GENCODE, made because the Catalog API rate-limited the script, before any
  v2 gene list existed
- the `*_named` sets, which swap in GBA1, MAPT and PICALM after the audit
  showed nearest-gene mapping misses them

The planned author-reported sensitivity set was dropped because it came out
identical to the primary PD set, and the Catalog has no author genes for the
AD study.

**5. Huntington's separated.** HD is now exploratory only, as two sets:
HTT alone, and the nine GeM-HD onset-modifier genes (hand-entered with
sources). The original v1 HD result is kept in `archive/v1_exploratory/`
and rerun on the v2 pipeline (`v1_huntingtons`).

**6. Gene-set specificity.** Added two tests:

- a random gene-set null (10,000 sets from the same 15,633-gene universe)
- a CAMERA-style competitive test that inflates variance for inter-gene
  correlation (Wu & Smyth 2012), addressing the Fulcher et al. 2021 warning

These are reported separately from the spatial null.

**7. BrainSMASH validation.**

- 10,000 surrogates per test.
- Variogram diagnostics are saved and plotted (fig5). They show the
  surrogates are less smooth than the data at short range, so the spatial
  p-values lean liberal.
- The 83-to-43 region question was a real v1 bug. v1 averaged left and right
  regions, which put every region centroid on the midline (all 43 at
  |x| < 3 mm), so the distance matrix ignored the left-right axis.
- v2 uses the left hemisphere plus brainstem, with true centroids. The
  ontology SNc and SNr regions are placed at the CIT168 left SNc and SNr
  centres.
- v1 also seeded BrainSMASH with Python's `hash()`, which changes between
  runs. v2 uses a fixed CRC32-based seed.

**8. Donor as the unit.** Donor-level contrasts T(d) are computed within each
donor, with sample counts reported. Leave-one-donor-out and
leave-one-gene-out are run for the primary tests. v1's "Cohen's d" is
replaced by a defined descriptive regional contrast z.

**9. Preprocessing.**

- All abagen settings are in `src/common.py`. abagen's own methods text is in
  `data/processed/abagen_methods_report.txt`.
- The rebuilt matrices match abagen's `get_expression_data` to 1e-5
  (`data/processed/abagen_check.txt`).
- Empty regions are no longer filled in (v1 used `missing="interpolate"`).
- Donor demographics are corrected: 5 male and 1 female, ages 24 to 57.
  v1 called them "older."

**10. Sensitivity and external check.** Sensitivity runs cover gene sets
(top 10, all GWS, literature-named), region definitions (ontology, CIT168,
6 mm), normalization (post-plan), and the v1 gene lists. The external check
uses GTEx v10, specified before running.

**11. Claims rewritten.** `docs/v2_results.md` has a table of what can and
can't be claimed. `report/report_draft.md` still describes v1 and carries a
banner saying so.

**12. Reproducibility and related work.** Added:

- `run_all.py`, one command for the whole pipeline
- pinned `requirements.txt`
- `data/gene_lists/provenance.json`, with the Catalog release (2026-09-13),
  GENCODE v50 and the rules used
- `docs/related_work.md`, with 16 PubMed-verified references

## Result of the revision

Neither primary test is supported. PD enrichment appears only under
anatomically incorrect region definitions, the v1 gene list, or a post-plan
normalization choice, and in no case does it clear the gene-set null. See
`docs/v2_results.md`.

## After v2: amendment 1 (2026-09-25)

The review asked whether a null result could be trusted, and the frozen plan
had no way to answer that. `docs/analysis_plan_amendment_1.md` was written
before running `src/power_check.py`, and everything in it is labelled post
hoc:

- positive controls: dopaminergic neuron markers in SNc and striatal
  projection neuron markers in caudate and putamen, through the unchanged
  primary pipeline
- a detection threshold: how high the PD and AD sets' target contrast would
  have had to be before each planned test called it
- power for random 9-gene sets shifted to a given contrast

Both controls pass every test in all six donors. The PD set would have
needed z of about 4.75, against the observed 1.45. The AD set would have
needed about 3.75, against the observed 0.64. The v2 verdicts stand, and
their reading is now "no large enrichment". `analysis.run()` gained an
optional `tgt` argument for the controls; no existing run changes.

## v3: cell types (2026-09-27)

The power analysis showed bulk regional data can only detect large effects,
and published single-cell studies place PD risk in dopamine neurons and AD
risk in microglia. v3 asks the question at that resolution.
`docs/analysis_plan_v3.md` was frozen before any single-nucleus expression
value was computed. `src/singlecell.py` streams public CELLxGENE files
(Kamath et al. 2022; Siletti et al. 2023) and runs the same kind of
random-gene-set test on cell types.

- AD genes in microglia: supported.
- PD genes in dopamine neurons: not supported, in two datasets.
- Dopamine neurons are 0.32% of SN nuclei. At most 4.7% of each PD gene's
  SN signal comes from them.
- A post hoc check showed microglia markers are moderately higher in
  entorhinal cortex in bulk AHBA, so "microglia are evenly spread" doesn't
  explain the v2 AD null. The AD set's mix of cell types does.

Results: `docs/v3_results.md`. New dependency: `h5py`.

## Files

- New code: `src/common.py`, `build_atlases.py`, `expression.py`,
  `sn_anatomy_audit.py`, `snca_checkpoint.py`, `gene_probe_audit.py`,
  `analysis.py`, `gtex_check.py`, `figures.py`, `run_all.py`
- Rewritten: `src/fetch_gwas_genes.py`
- Removed (replaced by the above; the originals are in commit `220b838`):
  `build_sn_atlas.py`, `expression_pipeline.py`, `regional_comparison.py`,
  `permutation_test.py`, `validation_checkpoint.py`, `visualize.py`

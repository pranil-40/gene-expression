"""
Sample-level AHBA expression, and region x gene matrices built from it.

  python src/expression.py extract
      Runs abagen once over every tissue sample (probe reannotation, intensity
      filtering, differential-stability probe selection) with normalization
      switched off, and caches the result per donor in data/cache/. This is
      the slow step. Probe selection doesn't depend on the atlas, so every
      region definition below can reuse the same cache.

  python src/expression.py build
      For each region definition: assign samples to regions with abagen's own
      matching code (2 mm tolerance), keep the matched samples, apply abagen's
      SRS normalization (across genes per sample, then across the kept samples
      per gene within each donor), average samples within region per donor,
      then average donors. Writes donor-level and pooled matrices plus sample
      counts.

  python src/expression.py check
      Runs abagen's get_expression_data end to end on the left_dk atlas and
      confirms the rebuilt matrix matches it. Also saves abagen's own methods
      report.

build also writes dk_ontology_sn__gene_norm_only (abagen's within-sample
normalization switched off), a post hoc sensitivity check. See
results/snca_checkpoint for why.

Region definitions:
  dk_ontology_sn  PRIMARY. left_dk atlas, except samples the AHBA anatomists
                  labelled substantia nigra are assigned by that label to
                  substantia_nigra_pc / substantia_nigra_pr (see the anatomy
                  audit for why coordinates can't be used for the SN).
  cit168          left_cit168 atlas, SN by coordinates in the CIT168 atlas.
  sn6mm           left_sn6mm atlas, v1-style 6 mm mask (exploratory).
  left_dk         plain left_dk, no SN region (used for the abagen check).
"""
from __future__ import annotations

import json
import sys

import abagen
import numpy as np
import pandas as pd
from abagen import correct, images, samples_

from common import ABAGEN_KWARGS, ATLAS_DIR, CACHE, DONORS, PROCESSED, donor_dir

VARIANTS = {
    "dk_ontology_sn": "left_dk",
    "cit168": "left_cit168",
    "sn6mm": "left_sn6mm",
    "left_dk": "left_dk",
}
EXCLUDE_REGIONS = {"substantia_nigra_pc_R", "substantia_nigra_pr_R", "substantia_nigra_6mm_R"}


def atlas_paths(name: str):
    return str(ATLAS_DIR / f"{name}.nii.gz"), pd.read_csv(ATLAS_DIR / f"{name}_info.csv")


def extract():
    CACHE.mkdir(parents=True, exist_ok=True)
    img, info = atlas_paths("left_dk")
    kw = dict(ABAGEN_KWARGS)
    kw.update(sample_norm=None, gene_norm=None, norm_matched=False)
    print("running abagen over all samples (no normalization), this is the slow part...")
    data = abagen.get_expression_data(img, info, region_agg=None, return_donors=True, verbose=0, **kw)
    genes = None
    for uid, df in data.items():
        wells = df.index.get_level_values("well_id").to_numpy()
        if genes is None:
            genes = df.columns.tolist()
        assert df.columns.tolist() == genes
        np.save(CACHE / f"raw_{uid}.npy", df.to_numpy(dtype=np.float32))
        pd.Series(wells, name="well_id").to_csv(CACHE / f"raw_{uid}_wells.csv", index=False)
        print(f"donor {uid}: {df.shape[0]} samples x {df.shape[1]} genes")
    pd.Series(genes, name="gene").to_csv(CACHE / "genes.csv", index=False)
    (CACHE / "extract_settings.json").write_text(json.dumps({**kw, "abagen": abagen.__version__}, indent=2))


def load_raw(uid: str) -> pd.DataFrame:
    genes = pd.read_csv(CACHE / "genes.csv")["gene"].tolist()
    wells = pd.read_csv(CACHE / f"raw_{uid}_wells.csv")["well_id"].to_numpy()
    return pd.DataFrame(np.load(CACHE / f"raw_{uid}.npy"), index=wells, columns=genes)


def donor_annotation(uid: str) -> pd.DataFrame:
    annot = samples_.update_coords(str(donor_dir(uid) / "SampleAnnot.csv"),
                                   corrected_mni=ABAGEN_KWARGS["corrected_mni"])
    return samples_.drop_mismatch_samples(annot, str(donor_dir(uid) / "Ontology.csv"))


def sample_labels(variant: str) -> dict[str, pd.Series]:
    """well_id -> region label for every sample, per donor. '' means unmatched."""
    img, info = atlas_paths(VARIANTS[variant])
    tree = images.check_atlas(img, info)
    id2lab = info.set_index("id")["label"].to_dict()
    out = {}
    for uid in DONORS:
        annot = donor_annotation(uid)
        ids = tree.label_samples(annot, tolerance=ABAGEN_KWARGS["tolerance"])["label"].to_numpy()
        lab = pd.Series([id2lab.get(int(i), "") for i in ids], index=annot["well_id"].to_numpy())
        if variant == "dk_ontology_sn":
            names = annot.set_index("well_id")["structure_name"]
            sn = names[names.str.contains("substantia nigra", case=False, na=False)]
            for well, name in sn.items():
                side = "_R" if "right" in name else ""
                sub = "pc" if "pars compacta" in name else "pr"
                lab.loc[well] = f"substantia_nigra_{sub}{side}"
        out[uid] = lab
    return out


def normalize(raw: pd.DataFrame, sample_norm="srs", gene_norm="srs") -> pd.DataFrame:
    # same calls abagen makes: SRS across genes within each sample, then SRS
    # across the kept samples within each gene
    x = raw
    if sample_norm:
        x = correct.normalize_expression(x.T, norm=sample_norm, ignore_warn=True).T
    if gene_norm:
        x = correct.normalize_expression(x, norm=gene_norm, ignore_warn=True)
    return x


# post-hoc preprocessing sensitivity (not in the frozen plan): abagen's
# within-sample normalization turned out to move SNCA's substantia nigra
# rank from 3rd to 20th, so the gene-set analysis is also rerun without it
NORM_MODES = {"default": ("srs", "srs"), "gene_norm_only": (None, "srs")}


def build(variant: str, save: bool = True, norm_mode: str = "default"):
    sample_norm, gene_norm = NORM_MODES[norm_mode]
    labels = sample_labels(variant)
    donor_means, counts = {}, {}
    for uid in DONORS:
        raw = load_raw(uid)
        lab = labels[uid].reindex(raw.index).fillna("")
        keep = lab.ne("").to_numpy()
        norm = normalize(raw.loc[keep], sample_norm, gene_norm)
        region = lab[keep]
        donor_means[uid] = norm.groupby(region.to_numpy()).mean()
        counts[uid] = region.value_counts()
    counts = pd.DataFrame(counts).fillna(0).astype(int)
    regions = sorted(set().union(*[d.index for d in donor_means.values()]) - EXCLUDE_REGIONS)
    stacked = np.stack([donor_means[u].reindex(regions).to_numpy() for u in DONORS])
    pooled = pd.DataFrame(np.nanmean(stacked, axis=0), index=regions, columns=donor_means[next(iter(DONORS))].columns)
    counts = counts.reindex(regions).fillna(0).astype(int)
    if save:
        name = variant if norm_mode == "default" else f"{variant}__{norm_mode}"
        out = CACHE / name
        out.mkdir(parents=True, exist_ok=True)
        pooled.astype(np.float32).to_pickle(out / "pooled.pkl")
        for uid in DONORS:
            donor_means[uid].reindex(regions).astype(np.float32).to_pickle(out / f"donor_{uid}.pkl")
        counts.to_csv(out / "sample_counts.csv")
        (PROCESSED / variant).mkdir(parents=True, exist_ok=True)
        counts.to_csv(PROCESSED / variant / "sample_counts.csv")
        print(f"{name}: {len(regions)} regions, {pooled.shape[1]} genes, "
              f"{int(counts.values.sum())} samples kept")
    return pooled, donor_means, counts


def load_variant(variant: str):
    out = CACHE / variant
    pooled = pd.read_pickle(out / "pooled.pkl")
    donors = {uid: pd.read_pickle(out / f"donor_{uid}.pkl") for uid in DONORS}
    counts = pd.read_csv(out / "sample_counts.csv", index_col=0)
    return pooled, donors, counts


def check():
    img, info = atlas_paths("left_dk")
    ref, counts, report = abagen.get_expression_data(img, info, region_agg="donors", return_counts=True,
                                                     return_report=True, verbose=0, **ABAGEN_KWARGS)
    ref.index = ref.index.map(info.set_index("id")["label"])
    mine, _, _ = load_variant("left_dk")
    common = ref.index.intersection(mine.index)
    genes = ref.columns.intersection(mine.columns)
    diff = np.nanmax(np.abs(ref.loc[common, genes].to_numpy() - mine.loc[common, genes].to_numpy()))
    msg = (f"abagen get_expression_data vs rebuilt left_dk: {len(common)} regions, {len(genes)} genes, "
           f"max abs difference {diff:.2e}")
    print(msg)
    (PROCESSED / "abagen_check.txt").write_text(msg + "\n")
    (PROCESSED / "abagen_methods_report.txt").write_text(report + "\n", encoding="utf-8")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "build"
    if cmd == "extract":
        extract()
    elif cmd == "build":
        for v in (sys.argv[2:] or VARIANTS):
            build(v)
        build("dk_ontology_sn", norm_mode="gene_norm_only")
    elif cmd == "check":
        check()

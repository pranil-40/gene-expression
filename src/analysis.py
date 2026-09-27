"""
v2 analysis, as laid out in docs/analysis_plan.md.

For each (gene set, region definition, normalization) run:
  - region score map S = mean normalized expression of the set's genes
  - T = mean S over target regions - mean S over the rest, and z = T / SD(rest)
  - spatial null: BrainSMASH surrogate maps of S
  - gene-set null: random gene sets of the same size
  - CAMERA-style competitive test that allows for correlated genes
  - donor-level T for each donor, plus leave-one-gene-out and
    leave-one-donor-out (primary runs only)

python src/analysis.py
-> results/v2/*.csv, results/v2/nulls/*.npz, results/v2/variograms/*.csv, results/v2/scores/*.csv
"""
from __future__ import annotations

import time

import nibabel as nib
import numpy as np
import pandas as pd
from brainsmash.mapgen.base import Base
from brainsmash.mapgen.eval import base_fit
from scipy import stats
from scipy.spatial.distance import pdist, squareform

from build_atlases import CIT168_FILE, CIT168_VOLUMES
from common import (ATLAS_DIR, DONORS, GENE_DIR, RESULTS, ROOT, ahba_symbols, bh_fdr, empirical_p,
                    regional_contrast, seed_for)
from expression import load_variant

OUT = RESULTS / "v2"
N_SURR = 10_000
N_SURR_LOO = 5_000
N_SETS = 10_000
MIN_REF_REGIONS = 10

ATLAS_OF = {"dk_ontology_sn": "left_dk", "cit168": "left_cit168", "sn6mm": "left_sn6mm"}

RUNS = [
    # (kind, gene set, region definition, normalization)
    ("primary", "parkinsons", "dk_ontology_sn", "default"),
    ("primary", "alzheimers", "dk_ontology_sn", "default"),
    ("secondary: gene set", "parkinsons_all_gws", "dk_ontology_sn", "default"),
    ("secondary: gene set", "alzheimers_all_gws", "dk_ontology_sn", "default"),
    # the planned author-reported sensitivity set turned out identical to the
    # primary PD set (the Catalog's author genes are the nearest genes) and
    # the Catalog has no author genes for the AD study, so it isn't run
    ("exploratory: literature-named genes", "parkinsons_named", "dk_ontology_sn", "default"),
    ("exploratory: literature-named genes", "alzheimers_named", "dk_ontology_sn", "default"),
    ("secondary: region definition", "parkinsons", "cit168", "default"),
    ("secondary: region definition", "parkinsons", "sn6mm", "default"),
    ("secondary: region definition", "alzheimers", "cit168", "default"),
    ("secondary: region definition", "alzheimers", "sn6mm", "default"),
    ("exploratory: HD", "huntingtons_htt", "dk_ontology_sn", "default"),
    ("exploratory: HD", "huntingtons_modifiers", "dk_ontology_sn", "default"),
    ("exploratory: v1 gene list", "v1_parkinsons", "dk_ontology_sn", "default"),
    ("exploratory: v1 gene list", "v1_alzheimers", "dk_ontology_sn", "default"),
    ("exploratory: v1 gene list", "v1_huntingtons", "dk_ontology_sn", "default"),
    ("post hoc: normalization", "parkinsons", "dk_ontology_sn", "gene_norm_only"),
    ("post hoc: normalization", "alzheimers", "dk_ontology_sn", "gene_norm_only"),
    ("post hoc: normalization", "huntingtons_htt", "dk_ontology_sn", "gene_norm_only"),
    ("post hoc: normalization", "huntingtons_modifiers", "dk_ontology_sn", "gene_norm_only"),
]


def disease_of(set_name: str) -> str:
    for d in ("parkinsons", "alzheimers", "huntingtons"):
        if d in set_name:
            return d
    raise ValueError(set_name)


def targets(set_name: str, variant: str) -> list[str]:
    d = disease_of(set_name)
    if d == "parkinsons":
        return ["substantia_nigra_6mm"] if variant == "sn6mm" else ["substantia_nigra_pc"]
    if d == "alzheimers":
        return ["entorhinal"]
    return ["caudate", "putamen"]


# Exploratory, added after reading the gene audit (not in the frozen plan).
# At three top-10 loci the nearest protein-coding gene isn't the gene the
# field usually credits: the GBA1 and MAPT loci in PD, the PICALM locus in AD.
# These sets swap only those three.
NAMED_SWAPS = {
    "parkinsons": {"KRTCAP2": "GBA1", "CRHR1": "MAPT"},
    "alzheimers": {"EED": "PICALM"},
}


def load_sets() -> dict[str, list[str]]:
    sets = {f.stem: pd.read_csv(f)["gene"].dropna().tolist() for f in sorted((GENE_DIR / "sets").glob("*.csv"))}
    v1 = pd.read_csv(ROOT / "archive" / "v1_exploratory" / "data" / "gene_lists" / "all_genes.csv")
    for d, g in v1.groupby("disease"):
        sets[f"v1_{d}"] = g["gene"].tolist()
    for d, swaps in NAMED_SWAPS.items():
        if d in sets:
            sets[f"{d}_named"] = [swaps.get(g, g) for g in sets[d]]
    return sets


def centroids(variant: str) -> pd.DataFrame:
    atlas = ATLAS_OF[variant]
    img = nib.load(ATLAS_DIR / f"{atlas}.nii.gz")
    data = np.asarray(img.dataobj)
    info = pd.read_csv(ATLAS_DIR / f"{atlas}_info.csv")
    rows = {}
    for _, r in info.iterrows():
        ijk = np.argwhere(data == r["id"])
        if len(ijk):
            rows[r["label"]] = nib.affines.apply_affine(img.affine, ijk.mean(axis=0))
    if variant == "dk_ontology_sn":
        # ontology SN regions have no voxels; place them at the CIT168 left
        # SNc / SNr centres (probability-weighted)
        cit = nib.load(CIT168_FILE)
        prob = np.asarray(cit.dataobj)
        for label, name in [("substantia_nigra_pc", "SNc"), ("substantia_nigra_pr", "SNr")]:
            vol = prob[..., CIT168_VOLUMES.index(name)]
            ijk = np.argwhere(vol > 0)
            xyz = nib.affines.apply_affine(cit.affine, ijk)
            w = vol[vol > 0] * (xyz[:, 0] < 0)
            rows[label] = (xyz * w[:, None]).sum(axis=0) / w.sum()
    return pd.DataFrame(rows, index=["x", "y", "z"]).T


class Context:
    """Everything that depends only on (region definition, normalization)."""

    def __init__(self, variant: str, norm: str):
        name = variant if norm == "default" else f"{variant}__{norm}"
        pooled, donors, counts = load_variant(name)
        self.regions = pooled.index[pooled.notna().all(axis=1)].tolist()
        self.pooled = pooled.loc[self.regions]
        self.donors = {u: d.reindex(self.regions) for u, d in donors.items()}
        self.counts = counts.reindex(self.regions).fillna(0).astype(int)
        cen = centroids(variant).reindex(self.regions)
        assert cen.notna().all().all(), f"missing centroid in {variant}: {cen[cen.isna().any(axis=1)].index.tolist()}"
        self.centroids = cen
        self.D = squareform(pdist(cen.to_numpy()))
        self.M = self.pooled.to_numpy(dtype=float)
        self.genes = pd.Index(self.pooled.columns)


def contrast_vec(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    # values: (n_maps, n_regions)
    return values[:, mask].mean(axis=1) - values[:, ~mask].mean(axis=1)


def spatial_null(S: np.ndarray, D: np.ndarray, mask: np.ndarray, n: int, seed: int) -> np.ndarray:
    # BrainSMASH defaults (exponential kernel, pv=25, nh=25, no resampling)
    gen = Base(x=S, D=D, seed=seed, n_jobs=8)
    surr = gen(n=n)
    return contrast_vec(np.atleast_2d(surr), mask)


def run(kind: str, set_name: str, genes_req: list[str], ctx: Context, variant: str, norm: str, detail: bool,
        tgt: list[str] | None = None):
    # tgt overrides the disease target (used by the positive controls in power_check.py)
    tag = f"{set_name}|{variant}|{norm}"
    sym = ahba_symbols(list(dict.fromkeys(genes_req)), ctx.genes)
    genes = list(dict.fromkeys(s for s in sym.values() if s))
    missing = [g for g, s in sym.items() if s is None]
    tgt = tgt or targets(set_name, variant)
    mask = np.isin(ctx.regions, tgt)
    assert mask.any(), f"no target region for {tag}"
    S = ctx.pooled[genes].mean(axis=1)
    res = {"kind": kind, "gene_set": set_name, "region_definition": variant, "normalization": norm,
           "n_genes_requested": len(set(genes_req)), "n_genes_used": len(genes),
           "genes_missing_from_ahba": ";".join(missing), "targets": ";".join(tgt)}
    res.update(regional_contrast(S, [r for r in tgt if r in S.index]))

    # spatial null
    null_sp = spatial_null(S.to_numpy(), ctx.D, mask, N_SURR, seed_for(tag + "|spatial"))
    res["p_spatial"] = empirical_p(null_sp, res["T"])

    # gene-set null: T for a random set is the mean over its genes of the
    # per-gene target-minus-reference difference (exact, since there are no
    # missing values after dropping empty regions)
    d_gene = ctx.M[mask].mean(axis=0) - ctx.M[~mask].mean(axis=0)
    rng = np.random.default_rng(seed_for(tag + "|genesets"))
    k = len(genes)
    null_gs = np.array([d_gene[rng.choice(len(d_gene), size=k, replace=False)].mean() for _ in range(N_SETS)])
    res["p_geneset"] = empirical_p(null_gs, res["T"])

    # CAMERA-style competitive test with a variance inflation factor
    sd_ref = ctx.M[~mask].std(axis=0, ddof=1)
    t_gene = np.divide(d_gene, sd_ref, out=np.zeros_like(d_gene), where=sd_ref > 0)
    idx = ctx.genes.get_indexer(genes)
    if k > 1:
        C = np.corrcoef(ctx.M[:, idx].T)
        rho = C[np.triu_indices(k, 1)].mean()
    else:
        rho = 0.0
    vif = 1 + (k - 1) * rho
    z_cam = (t_gene[idx].mean() - t_gene.mean()) / (t_gene.std(ddof=1) * np.sqrt(max(vif, 1e-9) / k))
    res.update({"mean_intergene_corr": rho, "camera_z": z_cam, "p_camera": 1 - stats.norm.cdf(z_cam)})

    # donor-level contrasts
    drows = []
    for uid, Md in ctx.donors.items():
        Sd = Md[genes].mean(axis=1)
        have_t = [r for r in tgt if pd.notna(Sd.get(r))]
        ref = [r for r in ctx.regions if r not in tgt and pd.notna(Sd.get(r))]
        eligible = bool(have_t) and len(ref) >= MIN_REF_REGIONS
        row = {"gene_set": set_name, "region_definition": variant, "normalization": norm, "donor": uid,
               "eligible": eligible, "n_target_samples": int(ctx.counts.loc[tgt, uid].sum()),
               "n_reference_regions": len(ref)}
        if eligible:
            c = regional_contrast(Sd, have_t, ref)
            row.update({"T": c["T"], "z": c["z"]})
        drows.append(row)
    ddf = pd.DataFrame(drows)
    elig = ddf[ddf["eligible"]]
    res["donors_eligible"] = len(elig)
    res["donors_positive"] = int((elig["T"] > 0).sum())

    extra = {"donors": ddf, "null_spatial": null_sp, "null_geneset": null_gs, "S": S}
    if detail:
        loo = []
        for g in genes:
            Sg = ctx.pooled[[x for x in genes if x != g]].mean(axis=1)
            T = regional_contrast(Sg, tgt)["T"]
            p = empirical_p(spatial_null(Sg.to_numpy(), ctx.D, mask, N_SURR_LOO, seed_for(tag + "|loo|" + g)), T)
            loo.append({"gene_set": set_name, "left_out": g, "type": "gene", "T": T, "p_spatial": p})
        for uid in DONORS:
            others = [ctx.donors[u][genes] for u in DONORS if u != uid]
            pooled = pd.concat(others).groupby(level=0).mean().reindex(ctx.regions)
            Sd = pooled.mean(axis=1)
            if Sd.isna().any():
                keep = Sd.notna().to_numpy()
                sub_mask = mask[keep]
                Dd = ctx.D[np.ix_(keep, keep)]
                Sd = Sd[keep]
            else:
                sub_mask, Dd = mask, ctx.D
            T = regional_contrast(Sd, [r for r in tgt if r in Sd.index])["T"]
            p = empirical_p(spatial_null(Sd.to_numpy(), Dd, sub_mask, N_SURR_LOO, seed_for(tag + "|lod|" + uid)), T)
            loo.append({"gene_set": set_name, "left_out": uid, "type": "donor", "T": T, "p_spatial": p})
        extra["loo"] = pd.DataFrame(loo)
        res["loo_gene_all_positive"] = bool((extra["loo"].query("type == 'gene'")["T"] > 0).all())
        emp, u0, surr = base_fit(S.to_numpy(), ctx.D, nsurr=200, return_data=True)
        extra["variogram"] = pd.DataFrame({"distance": u0, "empirical": emp,
                                           "surrogate_mean": surr.mean(axis=0),
                                           "surrogate_p05": np.percentile(surr, 5, axis=0),
                                           "surrogate_p95": np.percentile(surr, 95, axis=0)})
    return res, extra


def verdict(r: pd.Series) -> str:
    sig_sp, sig_gs = r["p_fdr_spatial"] < 0.05, r["p_fdr_geneset"] < 0.05
    if sig_sp and sig_gs and r["donors_positive"] >= 5 and r["loo_gene_all_positive"]:
        return "supported"
    if sig_sp or sig_gs:
        return "partial"
    return "not supported"


def main():
    for sub in ("nulls", "variograms", "scores"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    sets = load_sets()
    contexts, results, donors, loos, scores = {}, [], [], [], {}
    for kind, set_name, variant, norm in RUNS:
        if set_name not in sets:
            print(f"skip {set_name}: no such gene set", flush=True)
            continue
        key = (variant, norm)
        if key not in contexts:
            contexts[key] = Context(variant, norm)
        t0 = time.time()
        res, extra = run(kind, set_name, sets[set_name], contexts[key], variant, norm, detail=(kind == "primary"))
        results.append(res)
        donors.append(extra["donors"])
        scores.setdefault(f"{variant}__{norm}", {})[set_name] = extra["S"]
        stem = f"{set_name}__{variant}__{norm}"
        np.savez_compressed(OUT / "nulls" / f"{stem}.npz", spatial=extra["null_spatial"], geneset=extra["null_geneset"])
        if "loo" in extra:
            loos.append(extra["loo"])
            extra["variogram"].to_csv(OUT / "variograms" / f"{stem}.csv", index=False)
        print(f"{kind:32s} {set_name:24s} {variant:15s} {norm:15s} T={res['T']:+.4f} z={res['z']:+.2f} "
              f"rank={res['best_target_rank']}/{res['n_regions']} p_sp={res['p_spatial']:.4f} "
              f"p_gs={res['p_geneset']:.4f} p_cam={res['p_camera']:.4f} donors+={res['donors_positive']}/"
              f"{res['donors_eligible']} ({time.time() - t0:.0f}s)", flush=True)

    df = pd.DataFrame(results)
    prim = df["kind"] == "primary"
    p = bh_fdr(np.concatenate([df.loc[prim, "p_spatial"], df.loc[prim, "p_geneset"]]))
    n = prim.sum()
    df.loc[prim, "p_fdr_spatial"], df.loc[prim, "p_fdr_geneset"] = p[:n], p[n:]
    df.loc[prim, "verdict"] = df[prim].apply(verdict, axis=1)
    df.to_csv(OUT / "all_results.csv", index=False)
    df[prim].to_csv(OUT / "primary_results.csv", index=False)
    pd.concat(donors).to_csv(OUT / "donor_contrasts.csv", index=False)
    if loos:
        pd.concat(loos).to_csv(OUT / "leave_one_out.csv", index=False)
    for name, s in scores.items():
        pd.DataFrame(s).to_csv(OUT / "scores" / f"{name}.csv")
    for key, ctx in contexts.items():
        ctx.centroids.to_csv(OUT / f"centroids_{key[0]}.csv")
    print("\nprimary:")
    print(df[prim][["gene_set", "T", "z", "best_target_rank", "n_regions", "p_spatial", "p_geneset", "p_camera",
                    "p_fdr_spatial", "p_fdr_geneset", "donors_positive", "donors_eligible",
                    "loo_gene_all_positive", "verdict"]].to_string(index=False))


if __name__ == "__main__":
    main()

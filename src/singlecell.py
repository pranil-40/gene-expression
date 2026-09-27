"""
v3: which cell types inside the damaged region express the risk genes
(docs/analysis_plan_v3.md).

  download    fetch the CELLxGENE h5ad files into data/raw/singlecell/
  aggregate   stream each file and sum CP10K expression and raw UMIs per
              (class, donor); caches to data/cache/singlecell/
  analyze     EWCE-style specificity tests, donor-level checks, positive
              controls, DA subtypes, dilution -> results/v3/

python src/singlecell.py download
python src/singlecell.py aggregate
python src/singlecell.py analyze
"""
from __future__ import annotations

import sys
import urllib.request

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from common import CACHE, GENE_DIR, RAW, RESULTS, bh_fdr, empirical_p, seed_for

SC_RAW = RAW / "singlecell"
SC_CACHE = CACHE / "singlecell"
OUT = RESULTS / "v3"

CXG = "https://datasets.cellxgene.cziscience.com/"
URLS = {
    # Siletti et al. 2023, dissections
    "siletti_SN": "2de1ee2f-d82a-4bea-9288-bae4ea063138",
    "siletti_MEC": "73118fbf-bb19-49c8-bfad-bdf9eb8e103d",
    "siletti_LEC": "b0b0656c-5290-4e40-9886-51d828c5b74c",
    "siletti_Pu": "99a7f50b-40a2-4ac6-8d74-3b1cefe5dcfe",
    "siletti_CaB": "defe49cc-ccee-42ca-9251-a6ec33f11b65",
    # Kamath et al. 2022, one file per cell class
    "kamath_DA": "1bc30a70-d262-4f2e-8f0e-2f4ad89e9534",
    "kamath_nonDA": "d6453205-a7ee-45f5-adcc-231d008551b8",
    "kamath_astro": "ef6083ee-8d24-44d0-8c5b-75c646a430c1",
    "kamath_micro": "5ad1d7bc-ddc5-465f-831a-e86ad73ea60d",
    "kamath_oligo": "d50bd229-0da2-4f1a-89c5-3f5beda53185",
    "kamath_opc": "bcb15673-6a70-4d6f-ad28-b3f42aa70f5e",
    "kamath_endo": "afe7056c-e968-438b-91ab-0fdc58d7b490",
}

KAMATH_CLASS = {"kamath_DA": "DA neurons", "kamath_nonDA": "other neurons", "kamath_astro": "astrocytes",
                "kamath_micro": "microglia", "kamath_oligo": "oligodendrocytes", "kamath_opc": "OPCs",
                "kamath_endo": "endothelial"}

DATASETS = {
    "kamath_sn": list(KAMATH_CLASS),
    "siletti_sn": ["siletti_SN"],
    "siletti_ec": ["siletti_MEC", "siletti_LEC"],
    "siletti_str": ["siletti_CaB", "siletti_Pu"],
}

MIN_CELLS_CLASS = 50
MIN_CELLS_DONOR = 20
DA_RULE_FRACTION = 0.20      # neuronal clusters with >= 20% of cells expressing SLC6A3
N_NULL = 10_000
N_NULL_DONOR = 1_000
N_NULL_LOGO = 2_000
CHUNK = 20_000

MARKERS = {
    "pc_dopaminergic": ["TH", "SLC6A3", "SLC18A2", "DDC", "KCNJ6", "NR4A2", "EN1", "FOXA2", "LMX1B", "PITX3"],
    "pc_microglia": ["P2RY12", "CX3CR1", "TMEM119", "CSF1R", "C1QA", "C1QB", "AIF1", "ITGAM", "SALL1", "HEXB"],
    "pc_striatal_msn": ["PPP1R1B", "DRD1", "DRD2", "PENK", "TAC1", "ADORA2A", "GPR88", "RGS9", "PDE10A", "RASD2"],
}

RUNS = [
    # (kind, dataset, gene set, target class)
    ("primary", "kamath_sn", "parkinsons", "DA neurons"),
    ("primary", "siletti_ec", "alzheimers", "microglia"),
    ("secondary: replication", "siletti_sn", "parkinsons", "DA neurons"),
    ("secondary: oligodendrocytes", "kamath_sn", "parkinsons", "oligodendrocytes"),
    ("secondary: oligodendrocytes", "siletti_sn", "parkinsons", "oligodendrocytes"),
    ("secondary: all GWS loci", "kamath_sn", "parkinsons_all_gws", "DA neurons"),
    ("secondary: all GWS loci", "siletti_sn", "parkinsons_all_gws", "DA neurons"),
    ("secondary: all GWS loci", "siletti_ec", "alzheimers_all_gws", "microglia"),
    ("positive control", "kamath_sn", "pc_dopaminergic", "DA neurons"),
    ("positive control (circular: SLC6A3 defines the class)", "siletti_sn", "pc_dopaminergic", "DA neurons"),
    ("positive control", "siletti_ec", "pc_microglia", "microglia"),
    ("positive control", "siletti_str", "pc_striatal_msn", "MSNs"),
    ("exploratory: HD", "siletti_str", "huntingtons_htt", "MSNs"),
    ("exploratory: HD", "siletti_str", "huntingtons_modifiers", "MSNs"),
    ("exploratory: DA subtypes", "kamath_da_subtypes", "parkinsons", "SOX6_AGTR1"),
    ("exploratory: DA subtypes", "kamath_da_subtypes", "parkinsons_all_gws", "SOX6_AGTR1"),
    ("exploratory: DA subtypes", "kamath_da_subtypes", "pc_dopaminergic", "SOX6_AGTR1"),
]


# ---------------------------------------------------------------- h5ad input

def obs_col(f: h5py.File, name: str) -> pd.Series:
    x = f["obs"][name]
    if isinstance(x, h5py.Group):                 # anndata categorical
        cats = x["categories"][()]
        cats = cats.astype(str) if cats.dtype.kind in "SO" else cats
        return pd.Series(pd.Categorical.from_codes(x["codes"][()], cats)).astype(object)
    v = x[()]
    return pd.Series(v.astype(str) if v.dtype.kind in "SO" else v)


def var_table(f: h5py.File) -> pd.DataFrame:
    def get(name):
        x = f["var"][name]
        if isinstance(x, h5py.Group):
            cats = x["categories"][()].astype(str)
            return cats[x["codes"][()]]
        return x[()].astype(str)
    return pd.DataFrame({"gene_id": get("_index"), "symbol": get("feature_name")})


def stream(path, group_codes: np.ndarray, n_groups: int, gene_col: int | None = None):
    """Sum CP10K expression and raw UMIs per group, one chunk of cells at a time.

    With gene_col set, it instead returns, per group, how many cells have at
    least one UMI of that gene (used for the DA-cluster rule).
    """
    with h5py.File(path, "r") as f:
        X = f["X"]
        n_rows, n_genes = (int(v) for v in X.attrs["shape"])
        indptr = X["indptr"][()]
        if gene_col is None:
            cp10k = np.zeros((n_groups, n_genes))
            umis = np.zeros((n_groups, n_genes))
        else:
            hits = np.zeros(n_groups)
        for start in range(0, n_rows, CHUNK):
            stop = min(start + CHUNK, n_rows)
            g = group_codes[start:stop]
            keep = np.nonzero(g >= 0)[0]
            if not len(keep):
                continue
            lo, hi = indptr[start], indptr[stop]
            m = sp.csr_matrix((X["data"][lo:hi].astype(np.float64), X["indices"][lo:hi], indptr[start:stop + 1] - lo),
                              shape=(stop - start, n_genes))
            G = sp.csr_matrix((np.ones(len(keep)), (g[keep], keep)), shape=(n_groups, stop - start))
            if gene_col is not None:
                hits += G @ (m[:, gene_col].toarray().ravel() > 0).astype(float)
                continue
            total = np.asarray(m.sum(axis=1)).ravel()
            cp10k += (G @ (sp.diags(1e4 / np.maximum(total, 1)) @ m)).toarray()
            umis += (G @ m).toarray()
    return hits if gene_col is not None else (cp10k, umis)


def siletti_class(supercluster: pd.Series, cell_type: pd.Series) -> pd.Series:
    cls = pd.Series("neurons", index=supercluster.index, dtype=object)
    cls[supercluster == "Oligodendrocyte"] = "oligodendrocytes"
    cls[supercluster.isin(["Oligodendrocyte precursor", "Committed oligodendrocyte precursor"])] = "OPCs"
    cls[supercluster == "Astrocyte"] = "astrocytes"
    cls[supercluster == "Microglia"] = None
    cls[(supercluster == "Microglia") & (cell_type == "microglial cell")] = "microglia"
    cls[supercluster.isin(["Vascular", "Fibroblast"])] = "vascular"
    cls[supercluster.isin(["Miscellaneous", "Ependymal"])] = None
    return cls


# ---------------------------------------------------------------- aggregate

def aggregate_files(dataset: str) -> None:
    rows, cp_total, umi_total, genes = [], None, None, None
    notes = []
    for fn in DATASETS[dataset]:
        path = SC_RAW / f"{fn}.h5ad"
        with h5py.File(path, "r") as f:
            var = var_table(f)
            donor = obs_col(f, "donor_id").astype(str)
            if fn.startswith("kamath"):
                keep = obs_col(f, "disease") == "normal"
                cls = pd.Series(np.where(keep, KAMATH_CLASS[fn], None), dtype=object)
                # DA subtypes are kept for the exploratory subtype test; other classes pool their subtypes
                sub = obs_col(f, "author_cell_type").astype(str) if fn == "kamath_DA" else pd.Series("", index=cls.index)
            else:
                cls = siletti_class(obs_col(f, "supercluster_term"), obs_col(f, "cell_type"))
                cluster = obs_col(f, "cluster_id").astype(str)
                sub = pd.Series("", index=cls.index)
                if dataset == "siletti_str":
                    msn = obs_col(f, "supercluster_term").isin(["Medium spiny neuron", "Eccentric medium spiny neuron"])
                    cls[(cls == "neurons") & msn] = "MSNs"
                    cls[cls == "neurons"] = "other neurons"
        if dataset == "siletti_sn":
            # DA rule: neuronal clusters where >= 20% of cells have any SLC6A3 UMI
            neur = cls == "neurons"
            codes, uniq = pd.factorize(cluster.where(neur))
            col = int(np.nonzero(var["symbol"].to_numpy() == "SLC6A3")[0][0])
            hits = stream(path, codes, len(uniq), gene_col=col)
            n = np.bincount(codes[codes >= 0], minlength=len(uniq))
            rule = pd.DataFrame({"cluster_id": uniq, "n_cells": n, "n_slc6a3": hits, "frac_slc6a3": hits / n})
            rule["da_cluster"] = rule["frac_slc6a3"] >= DA_RULE_FRACTION
            OUT.mkdir(parents=True, exist_ok=True)
            rule.sort_values("frac_slc6a3", ascending=False).to_csv(OUT / "da_rule_clusters.csv", index=False)
            da = cluster.isin(rule.loc[rule["da_cluster"], "cluster_id"])
            cls[neur & da] = "DA neurons"
            cls[neur & ~da] = "other neurons"
            notes.append(f"DA clusters: {rule['da_cluster'].sum()} ({int(rule.loc[rule['da_cluster'], 'n_cells'].sum())} cells)")
        key = cls.fillna("").astype(str) + "|" + donor + "|" + sub
        key[cls.isna()] = None
        codes, uniq = pd.factorize(key)
        cp, um = stream(path, codes, len(uniq))
        grp = pd.DataFrame([u.split("|") for u in uniq], columns=["cell_class", "donor", "subtype"])
        grp["n_cells"] = np.bincount(codes[codes >= 0], minlength=len(uniq))
        if genes is None:
            genes = var
        elif not genes["gene_id"].equals(var["gene_id"]):
            # same genes, different order between files: line them up by Ensembl ID
            order = pd.Index(var["gene_id"]).get_indexer(genes["gene_id"])
            if (order < 0).any() or len(var) != len(genes):
                raise ValueError(f"gene sets differ in {fn}")
            cp, um = cp[:, order], um[:, order]
        rows.append(grp)
        cp_total = cp if cp_total is None else np.vstack([cp_total, cp])
        umi_total = um if umi_total is None else np.vstack([umi_total, um])
        print(f"  {fn}: {int(grp['n_cells'].sum())} cells kept in {len(grp)} groups", flush=True)
    groups = pd.concat(rows, ignore_index=True)
    # pool identical (class, donor, subtype) groups coming from different files
    key = groups["cell_class"] + "|" + groups["donor"] + "|" + groups["subtype"]
    codes, uniq = pd.factorize(key)
    G = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
    pooled = groups.groupby(codes).agg(cell_class=("cell_class", "first"), donor=("donor", "first"),
                                       subtype=("subtype", "first"), n_cells=("n_cells", "sum")).reset_index(drop=True)
    SC_CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(SC_CACHE / f"{dataset}.npz", cp10k=G @ cp_total, umis=G @ umi_total)
    pooled.to_csv(SC_CACHE / f"{dataset}_groups.csv", index=False)
    genes.to_csv(SC_CACHE / f"{dataset}_genes.csv", index=False)
    print(f"{dataset}: {int(pooled['n_cells'].sum())} cells, {len(pooled)} groups. " + "; ".join(notes), flush=True)


# ---------------------------------------------------------------- analyze

class Dataset:
    def __init__(self, name: str, by_subtype: bool = False):
        base = "kamath_sn" if by_subtype else name
        z = np.load(SC_CACHE / f"{base}.npz")
        self.groups = pd.read_csv(SC_CACHE / f"{base}_groups.csv", keep_default_na=False)
        self.genes = pd.read_csv(SC_CACHE / f"{base}_genes.csv")
        self.cp10k, self.umis = z["cp10k"], z["umis"]
        if by_subtype:
            keep = (self.groups["cell_class"] == "DA neurons").to_numpy()
            self.groups = self.groups[keep].reset_index(drop=True)
            self.groups["cell_class"] = self.groups["subtype"]
            self.cp10k, self.umis = self.cp10k[keep], self.umis[keep]
        n = self.groups.groupby("cell_class")["n_cells"].sum()
        self.classes = sorted(n[n >= MIN_CELLS_CLASS].index)
        self.class_n = n.loc[self.classes]
        self.means = self._class_means(np.ones(len(self.groups), bool))
        pc = pd.read_csv(RAW / "gencode" / "genes_v50.csv", usecols=["gene_id", "gene", "gene_type"])
        coding = set(pc.loc[pc["gene_type"] == "protein_coding", "gene_id"])
        base_id = self.genes["gene_id"].str.split(".").str[0]
        overall = self.cp10k.sum(axis=0) / self.groups["n_cells"].sum()
        self.universe = np.nonzero(base_id.isin(coding).to_numpy() & (overall > 0))[0]
        self.symbol_to_id = dict(zip(pc["gene"], pc["gene_id"]))
        self.id_index = {g: i for i, g in enumerate(base_id)}
        self.sym_index = {}
        for i, s in enumerate(self.genes["symbol"]):
            self.sym_index.setdefault(s, i)
        dec = pd.qcut(pd.Series(overall[self.universe]).rank(method="first"), 10, labels=False).to_numpy()
        self.decile = pd.Series(dec, index=self.universe)

    def _class_means(self, rows: np.ndarray) -> pd.DataFrame:
        out = {}
        for c in self.classes:
            m = rows & (self.groups["cell_class"] == c).to_numpy()
            n = self.groups.loc[m, "n_cells"].sum()
            if n >= (MIN_CELLS_DONOR if rows.sum() < len(rows) else MIN_CELLS_CLASS):
                out[c] = self.cp10k[m].sum(axis=0) / n
        return pd.DataFrame(out)

    def match(self, symbols: list[str]) -> tuple[list[int], list[str]]:
        idx, missing = [], []
        for s in dict.fromkeys(symbols):
            i = self.id_index.get(self.symbol_to_id.get(s, ""), self.sym_index.get(s))
            if i is None or i not in self.decile.index:
                missing.append(s)
            else:
                idx.append(i)
        return idx, missing


def specificity(means: pd.DataFrame) -> np.ndarray:
    tot = means.to_numpy().sum(axis=1, keepdims=True)
    return np.divide(means.to_numpy(), tot, out=np.zeros_like(means.to_numpy()), where=tot > 0)


def null_scores(spec: np.ndarray, ds: Dataset, genes: list[int], n: int, seed: int) -> np.ndarray:
    """Set scores for random sets whose genes match the set's expression deciles."""
    rng = np.random.default_rng(seed)
    pools = {d: ds.decile.index[ds.decile == d].to_numpy() for d in range(10)}
    need = pd.Series(ds.decile.loc[genes].to_numpy()).value_counts()
    picks = np.concatenate([np.array([rng.choice(pools[d], size=k, replace=False) for _ in range(n)])
                            for d, k in need.items()], axis=1)
    return spec[picks].mean(axis=1)            # (n, n_classes)


def test(ds: Dataset, spec: np.ndarray, classes: list[str], genes: list[int], n: int, seed: int) -> pd.DataFrame:
    obs = spec[genes].mean(axis=0)
    null = null_scores(spec, ds, genes, n, seed)
    rank = pd.Series(obs, index=classes).rank(ascending=False, method="min")
    return pd.DataFrame({"cell_class": classes, "score": obs, "null_mean": null.mean(axis=0),
                         "z": (obs - null.mean(axis=0)) / null.std(axis=0, ddof=1),
                         "p": [empirical_p(null[:, j], obs[j]) for j in range(len(classes))],
                         "rank": rank.to_numpy().astype(int), "n_classes": len(classes)})


def load_gene_sets() -> dict[str, list[str]]:
    sets = {f.stem: pd.read_csv(f)["gene"].dropna().tolist() for f in sorted((GENE_DIR / "sets").glob("*.csv"))}
    sets.update(MARKERS)
    return sets


def analyze() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sets = load_gene_sets()
    cache, long_rows, summary, donor_rows, logo_rows = {}, [], [], [], []
    for kind, dname, set_name, target in RUNS:
        if dname not in cache:
            cache[dname] = Dataset(dname, by_subtype=(dname == "kamath_da_subtypes"))
        ds = cache[dname]
        spec = specificity(ds.means)
        genes, missing = ds.match(sets[set_name])
        tag = f"{dname}|{set_name}"
        res = test(ds, spec, ds.classes, genes, N_NULL, seed_for(tag))
        res.insert(0, "gene_set", set_name)
        res.insert(0, "dataset", dname)
        long_rows.append(res)
        t = res.set_index("cell_class").loc[target]
        row = {"kind": kind, "dataset": dname, "gene_set": set_name, "target": target,
               "n_genes_used": len(genes), "n_genes_requested": len(dict.fromkeys(sets[set_name])),
               "genes_missing": ";".join(missing), "score": t["score"], "null_mean": t["null_mean"], "z": t["z"],
               "p": t["p"], "target_rank": int(t["rank"]), "n_classes": int(t["n_classes"]),
               "top_class": res.sort_values("score", ascending=False)["cell_class"].iloc[0]}
        # donor level
        ti = ds.classes.index(target)
        pos = n_don = 0
        for donor in sorted(ds.groups["donor"].astype(str).unique()):
            m = ds._class_means((ds.groups["donor"].astype(str) == donor).to_numpy())
            if target not in m.columns or m.shape[1] < 3:
                continue
            sd = specificity(m)
            j = list(m.columns).index(target)
            obs = sd[genes, j].mean()
            null = null_scores(sd, ds, genes, N_NULL_DONOR, seed_for(tag + "|donor|" + donor))[:, j]
            z = (obs - null.mean()) / null.std(ddof=1)
            donor_rows.append({"dataset": dname, "gene_set": set_name, "target": target, "donor": donor,
                               "n_target_cells": int(ds.groups.loc[(ds.groups["donor"].astype(str) == donor) &
                                                                   (ds.groups["cell_class"] == target), "n_cells"].sum()),
                               "score": obs, "z": z, "p": empirical_p(null, obs)})
            n_don += 1
            pos += int(z > 0)
        row.update({"donors_eligible": n_don, "donors_positive": pos})
        if kind == "primary":
            for k, g in enumerate(genes):
                rest = [x for x in genes if x != g]
                r = test(ds, spec, ds.classes, rest, N_NULL_LOGO, seed_for(tag + "|logo|" + str(g)))
                logo_rows.append({"dataset": dname, "gene_set": set_name, "left_out": ds.genes["symbol"].iloc[g],
                                  "p": r.set_index("cell_class").loc[target, "p"],
                                  "rank": int(r.set_index("cell_class").loc[target, "rank"])})
        summary.append(row)
        print(f"{kind:45s} {dname:18s} {set_name:22s} {target:16s} genes={len(genes)}/{row['n_genes_requested']} "
              f"z={row['z']:+.2f} p={row['p']:.4f} rank={row['target_rank']}/{row['n_classes']} "
              f"top={row['top_class']} donors+={pos}/{n_don}", flush=True)

    df = pd.DataFrame(summary)
    prim = df["kind"] == "primary"
    df.loc[prim, "p_bh"] = bh_fdr(df.loc[prim, "p"].to_numpy())

    def verdict(r):
        if r["p_bh"] < 0.05 and r["target_rank"] == 1:
            return "supported"
        return "partial" if r["p_bh"] < 0.05 else "not supported"

    df.loc[prim, "verdict"] = df[prim].apply(verdict, axis=1)
    df.to_csv(OUT / "specificity_tests.csv", index=False)
    pd.concat(long_rows).to_csv(OUT / "specificity_all_classes.csv", index=False)
    pd.DataFrame(donor_rows).to_csv(OUT / "donor_level.csv", index=False)
    pd.DataFrame(logo_rows).to_csv(OUT / "leave_one_gene_out.csv", index=False)
    counts = pd.concat([cache[d].class_n.rename("n_cells").to_frame().assign(dataset=d) for d in cache])
    counts.rename_axis("cell_class").reset_index().to_csv(OUT / "class_counts.csv", index=False)
    dilution(cache.get("siletti_sn") or Dataset("siletti_sn"), sets)
    gene_tables(cache, sets)
    bulk_link()
    print("\nprimary:")
    print(df[prim][["dataset", "gene_set", "target", "z", "p", "p_bh", "target_rank", "top_class",
                    "donors_positive", "donors_eligible", "verdict"]].to_string(index=False))


def dilution(ds: Dataset, sets: dict[str, list[str]]) -> None:
    """How much of the SN tissue signal for each gene comes from DA neurons."""
    da = (ds.groups["cell_class"] == "DA neurons").to_numpy()
    n_all, n_da = ds.groups["n_cells"].sum(), ds.groups.loc[da, "n_cells"].sum()
    umi_all, umi_da = ds.umis.sum(axis=0), ds.umis[da].sum(axis=0)
    cp_all, cp_da = ds.cp10k.sum(axis=0) / n_all, ds.cp10k[da].sum(axis=0) / n_da
    rows = []
    for set_name in ["parkinsons", "pc_dopaminergic"]:
        genes, _ = ds.match(sets[set_name])
        for g in genes:
            rows.append({"gene_set": set_name, "gene": ds.genes["symbol"].iloc[g],
                         "share_of_sn_umis_from_da": umi_da[g] / umi_all[g] if umi_all[g] else np.nan,
                         "da_mean_over_sn_mean": cp_da[g] / cp_all[g] if cp_all[g] else np.nan})
    out = pd.DataFrame(rows)
    out.attrs["da_fraction"] = n_da / n_all
    out.to_csv(OUT / "dilution.csv", index=False)
    (OUT / "dilution_summary.txt").write_text(
        f"Siletti SN dissection: {n_da} of {n_all} nuclei are DA neurons ({100 * n_da / n_all:.2f}%).\n"
        f"Share of all SN UMIs from DA-neuron nuclei (all genes): {umi_da.sum() / umi_all.sum():.4f}\n", encoding="utf-8")
    print(f"dilution: DA neurons are {n_da} of {n_all} SN nuclei ({100 * n_da / n_all:.2f}%)")


def gene_tables(cache: dict, sets: dict[str, list[str]]) -> None:
    """Per-gene specificity (rows sum to 1 over classes) for the primary sets."""
    rows = []
    for dname, set_name in [("kamath_sn", "parkinsons"), ("siletti_sn", "parkinsons"), ("siletti_ec", "alzheimers")]:
        ds = cache[dname]
        spec = pd.DataFrame(specificity(ds.means), columns=ds.classes)
        genes, _ = ds.match(sets[set_name])
        t = spec.loc[genes].copy()
        t.insert(0, "gene", ds.genes["symbol"].iloc[genes].to_numpy())
        t.insert(0, "gene_set", set_name)
        t.insert(0, "dataset", dname)
        rows.append(t.melt(id_vars=["dataset", "gene_set", "gene"], var_name="cell_class", value_name="specificity"))
    pd.concat(rows).to_csv(OUT / "gene_specificity.csv", index=False)


def bulk_link() -> None:
    """Post hoc: do microglia markers show a regional contrast in bulk AHBA?

    Added after the v3 results, to check whether microglia being spread
    through the brain could explain the v2 AD null. Uses the v2 primary
    pipeline unchanged.
    """
    from analysis import Context, run
    ctx = Context("dk_ontology_sn", "default")
    rows = []
    for tgt in (["entorhinal"], ["substantia_nigra_pc"]):
        res, _ = run("post hoc: bulk link", "pc_microglia", MARKERS["pc_microglia"], ctx, "dk_ontology_sn", "default",
                     detail=False, tgt=tgt)
        rows.append(res)
    pd.DataFrame(rows).to_csv(OUT / "bulk_link_microglia_markers.csv", index=False)


def download() -> None:
    SC_RAW.mkdir(parents=True, exist_ok=True)
    for name, uid in URLS.items():
        dest = SC_RAW / f"{name}.h5ad"
        if not dest.exists():
            print(f"downloading {name}", flush=True)
            urllib.request.urlretrieve(CXG + uid + ".h5ad", dest)


def main():
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    if step in ("download", "all"):
        download()
    if step in ("aggregate", "all"):
        for d in DATASETS:
            aggregate_files(d)
    if step in ("analyze", "all"):
        analyze()


if __name__ == "__main__":
    main()

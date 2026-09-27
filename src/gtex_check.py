"""
External check in GTEx (secondary, specified in docs/analysis_plan.md before
running): does the same gene set also sit highest in the matching tissue in
an independent RNA-seq dataset?

GTEx v10 median TPM across its 13 brain tissues. Values are log2(TPM + 1),
z-scored per gene across the 13 tissues, then T = mean z over target tissues
minus mean z over the other brain tissues, averaged over the set's genes.
Null: 10,000 random protein-coding gene sets of the same size. PD target =
substantia nigra. HD target = caudate + putamen. AD isn't tested because GTEx
has no entorhinal cortex sample.

python src/gtex_check.py
-> results/gtex/gtex_results.csv, gtex_gene_z.csv
"""
from __future__ import annotations

import gzip
import urllib.request

import numpy as np
import pandas as pd

from analysis import load_sets
from common import RAW, RESULTS, empirical_p, seed_for

URL = ("https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
       "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz")
FILE = RAW / "gtex" / "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
OUT = RESULTS / "gtex"
N_SETS = 10_000

TARGETS = {
    "parkinsons": ["Brain_Substantia_nigra"],
    "huntingtons": ["Brain_Caudate_basal_ganglia", "Brain_Putamen_basal_ganglia"],
}
SETS = ["parkinsons", "parkinsons_all_gws", "huntingtons_htt", "huntingtons_modifiers"]


def load_gtex() -> pd.DataFrame:
    if not FILE.exists():
        FILE.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(URL, FILE)
    with gzip.open(FILE, "rt") as fh:
        fh.readline(); fh.readline()
        df = pd.read_csv(fh, sep="\t")
    brain = [c for c in df.columns if c.startswith("Brain_")]
    df["gene_id"] = df["Name"].str.split(".").str[0]
    return df[["gene_id", "Description"] + brain]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    gtex = load_gtex()
    brain = [c for c in gtex.columns if c.startswith("Brain_")]
    genes = pd.read_csv(RAW / "gencode" / "genes_v50.csv")
    pc = genes[(genes["gene_type"] == "protein_coding") & ~genes["readthrough"]]
    # match on Ensembl ID: GTEx v10 uses older gene symbols (e.g. GBA for GBA1)
    sym_of_id = pc.drop_duplicates("gene_id").set_index("gene_id")["gene"]
    gtex = gtex[gtex["gene_id"].isin(pc["gene_id"])].copy()
    gtex["Description"] = gtex["gene_id"].map(sym_of_id)
    gtex = gtex.drop_duplicates("Description").set_index("Description")
    x = np.log2(gtex[brain] + 1)
    sd = x.std(axis=1, ddof=1)
    x = x[sd > 0]
    z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=1), axis=0)
    sets = load_sets()
    rows, gene_rows = [], []
    for name in SETS:
        disease = next(d for d in TARGETS if d in name)
        tgt = TARGETS[disease]
        ref = [c for c in brain if c not in tgt]
        per_gene = z[tgt].mean(axis=1) - z[ref].mean(axis=1)
        wanted = list(dict.fromkeys(sets[name]))
        used = [g for g in wanted if g in per_gene.index]
        T = per_gene[used].mean()
        rng = np.random.default_rng(seed_for(f"gtex|{name}"))
        vals = per_gene.to_numpy()
        null = np.array([vals[rng.choice(len(vals), size=len(used), replace=False)].mean() for _ in range(N_SETS)])
        rows.append({"gene_set": name, "targets": ";".join(tgt), "n_genes_used": len(used),
                     "genes_missing_from_gtex": ";".join(g for g in wanted if g not in per_gene.index),
                     "T": T, "p_geneset": empirical_p(null, T), "n_tissues": len(brain),
                     "target_rank_of_mean_z": int(pd.Series(z.loc[used].mean(), index=brain)
                                                  .rank(ascending=False)[tgt].min())})
        for g in used:
            gene_rows.append({"gene_set": name, "gene": g, "target_minus_rest_z": per_gene[g],
                              "target_rank": int(z.loc[g].rank(ascending=False)[tgt].min())})
    snca = z.loc["SNCA"].sort_values(ascending=False)
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "gtex_results.csv", index=False)
    pd.DataFrame(gene_rows).to_csv(OUT / "gtex_gene_z.csv", index=False)
    snca.rename("SNCA_z").to_csv(OUT / "gtex_snca_by_tissue.csv")
    print(res.to_string(index=False))
    print("\nSNCA z by GTEx brain tissue:")
    print(snca.round(2).to_string())


if __name__ == "__main__":
    main()

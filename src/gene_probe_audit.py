"""
For every gene in any v1 or v2 gene set, works out whether it's in the AHBA
gene universe and, if it isn't, why:

  no_probe_on_array     no probe for this symbol in Allen's Probes.csv
  dropped_reannotation  Allen has a probe, but abagen's probe reannotation
                        (Arnatkeviciute et al. 2019) doesn't map any probe to it
  failed_ibf            probes exist after reannotation, but none is above
                        background in >= 50% of samples (intensity filter)
  in_universe           survives everything and is in data/cache/genes.csv

python src/gene_probe_audit.py
-> data/gene_lists/probe_audit.csv
"""
from __future__ import annotations

import gzip
from io import StringIO

import numpy as np
import pandas as pd
from abagen import samples_
from importlib.resources import files

from common import ABAGEN_KWARGS, CACHE, DONORS, GENE_DIR, ahba_symbols, donor_dir


def all_genes() -> pd.DataFrame:
    from analysis import load_sets
    rows = [{"gene": g, "set": name} for name, genes in load_sets().items() for g in genes]
    df = pd.DataFrame(rows)
    return df.groupby("gene")["set"].agg(lambda s: ";".join(sorted(set(s)))).reset_index()


def main():
    genes = all_genes()
    first = next(iter(DONORS))
    allen = pd.read_csv(donor_dir(first) / "Probes.csv")
    with gzip.open(files("abagen") / "data" / "reannotated.csv.gz", "rt") as fh:
        reannot = pd.read_csv(StringIO(fh.read()))
    reannot = reannot.merge(allen[["probe_name", "probe_id"]], on="probe_name", how="left").dropna(subset=["entrez_id"])
    universe = set(pd.read_csv(CACHE / "genes.csv")["gene"])
    alias = ahba_symbols(genes["gene"].tolist(), universe)
    genes["ahba_symbol"] = genes["gene"].map(alias)

    wanted = reannot[reannot["gene_symbol"].isin(set(genes["gene"]) | set(genes["ahba_symbol"].dropna()))]
    ids = set(wanted["probe_id"].dropna().astype(int))
    # intensity-based filter, same rule as abagen: fraction of all kept
    # samples (across donors) where the probe is above background
    signal = pd.Series(0.0, index=sorted(ids))
    n = 0
    for uid in DONORS:
        annot = samples_.update_coords(str(donor_dir(uid) / "SampleAnnot.csv"), corrected_mni=True)
        annot = samples_.drop_mismatch_samples(annot, str(donor_dir(uid) / "Ontology.csv"))
        raw_annot = pd.read_csv(donor_dir(uid) / "SampleAnnot.csv")
        keep_cols = raw_annot["well_id"].isin(annot["well_id"]).to_numpy()
        with open(donor_dir(uid) / "PACall.csv") as fh:
            for line in fh:
                pid, rest = line.split(",", 1)
                if int(pid) in ids:
                    vals = np.array(rest.split(","), dtype=float)[keep_cols]
                    signal[int(pid)] += vals.sum()
        n += int(keep_cols.sum())
    frac = signal / n

    rows = []
    for _, g in genes.iterrows():
        sym = g["ahba_symbol"] if pd.notna(g["ahba_symbol"]) else g["gene"]
        a = allen[allen["gene_symbol"] == sym]
        r = reannot[reannot["gene_symbol"] == sym]
        best = frac.reindex(r["probe_id"].dropna().astype(int)).max() if len(r) else np.nan
        if sym in universe:
            status = "in_universe" if sym == g["gene"] else f"in_universe (AHBA name {sym})"
        elif len(r) and best < ABAGEN_KWARGS["ibf_threshold"]:
            status = "failed_ibf"
        elif len(a) and not len(r):
            status = "dropped_reannotation"
        elif not len(a) and not len(r):
            status = "no_probe_on_array"
        else:
            status = "other (check manually)"
        rows.append({"gene": g["gene"], "ahba_symbol": g["ahba_symbol"], "sets": g["set"], "allen_probes": len(a), "reannotated_probes": len(r),
                     "best_probe_frac_above_background": round(best, 3) if pd.notna(best) else None,
                     "status": status})
    out = pd.DataFrame(rows).sort_values(["status", "gene"])
    out.to_csv(GENE_DIR / "probe_audit.csv", index=False)
    print(out[out["status"] != "in_universe"].to_string(index=False))
    print(f"\n{out['status'].str.startswith('in_universe').sum()} of {len(out)} genes are in the AHBA universe")


if __name__ == "__main__":
    main()

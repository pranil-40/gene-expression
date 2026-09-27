"""Shared paths, processing settings and small stats helpers for the v2 pipeline."""
from __future__ import annotations

import zlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
GENE_DIR = DATA / "gene_lists"
ATLAS_DIR = DATA / "atlases"
PROCESSED = DATA / "processed"
CACHE = DATA / "cache"          # big regenerable matrices, gitignored
RAW = DATA / "raw"              # third-party downloads, gitignored
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"

ABAGEN_MICROARRAY = Path.home() / "abagen-data" / "microarray"

# AHBA donor uid -> donor name (abagen/data/donor_info.csv)
DONORS = {
    "9861": "H0351.2001",
    "10021": "H0351.2002",
    "12876": "H0351.1009",
    "14380": "H0351.1012",
    "15496": "H0351.1015",
    "15697": "H0351.1016",
}

# abagen settings for every v2 run. Mostly abagen defaults, written out so the
# methods section can point at one place. missing=None means a region with no
# tissue samples in a donor stays empty for that donor instead of being
# filled in from neighbouring samples (v1 used missing="interpolate").
ABAGEN_KWARGS = dict(
    ibf_threshold=0.5,
    probe_selection="diff_stability",
    donor_probes="aggregate",
    lr_mirror=None,
    missing=None,
    tolerance=2,
    sample_norm="srs",
    gene_norm="srs",
    norm_matched=True,
    norm_structures=False,
    corrected_mni=True,
    reannotated=True,
)

DISEASES = ["parkinsons", "alzheimers"]      # primary comparison
EXPLORATORY = ["huntingtons_htt", "huntingtons_modifiers"]


def donor_dir(uid: str) -> Path:
    return ABAGEN_MICROARRAY / f"normalized_microarray_donor{uid}"


def ahba_symbols(genes: list[str], universe) -> dict[str, str | None]:
    """Current gene symbol -> the symbol AHBA/abagen uses for the same gene.

    abagen's probe annotation dates from 2019, so some genes carry older
    names there (GBA1 is still "GBA"). Genes are matched on Entrez ID via the
    HGNC table. Previous symbols aren't used on their own: HGNC lists
    MIR4697HG as a previous symbol of IGSF9B, but the AHBA probes named
    MIR4697HG carry a different Entrez ID. None means the gene really isn't
    in the AHBA gene universe.
    """
    import gzip
    from importlib.resources import files
    from io import StringIO

    universe = set(universe)
    hgnc = pd.read_csv(RAW / "hgnc" / "hgnc_complete_set.txt", sep="\t", low_memory=False,
                       usecols=["symbol", "entrez_id"])
    with gzip.open(files("abagen") / "data" / "reannotated.csv.gz", "rt") as fh:
        reannot = pd.read_csv(StringIO(fh.read()))
    r = reannot.dropna(subset=["entrez_id"]).copy()
    r["entrez_id"] = r["entrez_id"].astype(int)
    by_entrez = r.groupby("entrez_id")["gene_symbol"].agg(lambda s: sorted(set(s)))
    out = {}
    for g in genes:
        if g in universe:
            out[g] = g
            continue
        rec = hgnc[hgnc["symbol"] == g]
        hit = None
        if len(rec) and pd.notna(rec["entrez_id"].iloc[0]):
            hit = next((s for s in by_entrez.get(int(rec["entrez_id"].iloc[0]), []) if s in universe), None)
        out[g] = hit
    return out


def seed_for(name: str, base: int = 20260924) -> int:
    # hash() on strings changes between Python processes, crc32 doesn't
    return (base + zlib.crc32(name.encode())) % (2**32)


def regional_contrast(scores: pd.Series, target: list[str],
                      reference: list[str] | None = None) -> dict:
    """Target-minus-reference difference in a region-level score map.

    T = mean(score over target regions) - mean(score over reference regions)
    z = T / sd(score over reference regions)

    z is a descriptive standardized contrast between regions of one map. It is
    not a two-group Cohen's d, because regions aren't independent people.
    """
    if reference is None:
        reference = [r for r in scores.index if r not in target]
    t_vals = scores.loc[target].dropna()
    r_vals = scores.loc[reference].dropna()
    T = t_vals.mean() - r_vals.mean()
    sd = r_vals.std(ddof=1)
    ranks = scores.dropna().rank(ascending=False)
    return {
        "T": T,
        "z": T / sd if sd > 0 else np.nan,
        "target_mean": t_vals.mean(),
        "reference_mean": r_vals.mean(),
        "reference_sd": sd,
        "best_target_rank": int(ranks.loc[t_vals.index].min()) if len(t_vals) else np.nan,
        "n_regions": int(scores.notna().sum()),
    }


def empirical_p(null: np.ndarray, observed: float) -> float:
    """One-sided (greater-or-equal) p with the +1 correction, ties count against."""
    null = np.asarray(null)
    null = null[~np.isnan(null)]
    return (1 + np.sum(null >= observed)) / (len(null) + 1)


def bh_fdr(pvals) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = p[order] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out

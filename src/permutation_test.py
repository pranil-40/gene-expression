"""
Step 6: test whether each disease's target-region enrichment (Step 5) holds
up against a spatial null model, not just naive shuffling.

Naive label permutation is the wrong null here: brain gene expression is
spatially autocorrelated (nearby regions look more alike than far-apart
ones), so a naive shuffle overstates significance. BrainSMASH instead
generates surrogate brain maps that preserve the real map's spatial
autocorrelation structure (via its variogram) while randomizing which
anatomical location gets which value — a much harder null to beat.

For each disease:
  1. Take the real per-region gene-score map (mean expression of that
     disease's risk genes, per region — computed in regional_comparison.py).
  2. Generate N_PERMUTATIONS BrainSMASH surrogates of that map using the
     real region x region distance matrix.
  3. For each surrogate, recompute the same statistic as Step 5
     (target-region mean minus all-other-region mean).
  4. Empirical p-value = fraction of surrogate statistics >= the observed
     statistic (one-sided: the hypothesis is enrichment, not depletion).
  5. Correct the 3 resulting p-values for multiple comparisons (Benjamini-
     Hochberg FDR).

Run:  .venv/Scripts/python src/permutation_test.py
Input:  data/processed/disease_gene_scores.csv, data/processed/region_centroids.csv
Output: results/tables/permutation_results.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from brainsmash.mapgen.base import Base
from scipy.spatial.distance import squareform, pdist

from regional_comparison import TARGET_REGION_PATTERNS, match_regions

ROOT = Path(__file__).resolve().parent.parent
SCORES_PATH = ROOT / "data" / "processed" / "disease_gene_scores.csv"
CENTROIDS_PATH = ROOT / "data" / "processed" / "region_centroids.csv"
OUT_PATH = ROOT / "results" / "tables" / "permutation_results.csv"

N_PERMUTATIONS = 1000
SEED = 42


def benjamini_hochberg(pvals: list[float]) -> list[float]:
    pvals = np.asarray(pvals)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity (BH step-up)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n)
    out[order] = adjusted
    return out.tolist()


def target_vs_rest_stat(values: np.ndarray, target_mask: np.ndarray) -> float:
    return values[target_mask].mean() - values[~target_mask].mean()


def main() -> None:
    scores = pd.read_csv(SCORES_PATH, index_col=0)
    centroids = pd.read_csv(CENTROIDS_PATH, index_col=0)

    # Align region ordering across scores and centroids, drop any region
    # missing a centroid (shouldn't happen, but fail loudly if it does).
    common = scores.index.intersection(centroids.index)
    scores = scores.loc[common]
    centroids = centroids.loc[common]

    distance_matrix = squareform(pdist(centroids.to_numpy(), metric="euclidean"))

    results = []
    for disease in scores.columns:
        brain_map = scores[disease].to_numpy()
        target_regions = match_regions(common.tolist(), TARGET_REGION_PATTERNS[disease])
        target_mask = common.isin(target_regions)
        if target_mask.sum() == 0:
            print(f"SKIP {disease}: no target region found")
            continue

        observed = target_vs_rest_stat(brain_map, target_mask)

        rng_seed = SEED + abs(hash(disease)) % 1000
        generator = Base(x=brain_map, D=distance_matrix, seed=rng_seed)
        surrogates = generator(n=N_PERMUTATIONS)  # shape: (N_PERMUTATIONS, n_regions)

        surrogate_stats = np.array(
            [target_vs_rest_stat(row, target_mask) for row in surrogates]
        )
        p_value = (np.sum(surrogate_stats >= observed) + 1) / (N_PERMUTATIONS + 1)

        results.append(
            {
                "disease": disease,
                "n_target_regions": int(target_mask.sum()),
                "observed_statistic": observed,
                "null_mean": surrogate_stats.mean(),
                "null_std": surrogate_stats.std(),
                "p_value_spatial": p_value,
                "n_permutations": N_PERMUTATIONS,
            }
        )

    result_df = pd.DataFrame(results)
    result_df["p_value_fdr_bh"] = benjamini_hochberg(result_df["p_value_spatial"].tolist())

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(OUT_PATH, index=False)
    print(result_df.to_string(index=False))
    print(f"\nWrote -> {OUT_PATH}")


if __name__ == "__main__":
    main()

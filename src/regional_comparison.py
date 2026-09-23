"""
For each disease, compares mean risk-gene expression in the region it
destroys vs. every other region.

Target regions (see README for why): PD -> Brain-Stem (DK atlas doesn't have
a separate substantia nigra label, brainstem is the closest thing available
and it's a real limitation, not ideal). AD -> entorhinal. HD -> caudate +
putamen (= striatum).
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
EXPRESSION_PATH = ROOT / "data" / "processed" / "regional_expression.csv"
GENE_LIST_DIR = ROOT / "data" / "gene_lists"
OUT_PATH = ROOT / "results" / "tables" / "regional_comparison.csv"

TARGET_REGION_PATTERNS = {
    "parkinsons": [r"brain.?stem"],
    "alzheimers": [r"entorhinal"],
    "huntingtons": [r"caudate", r"putamen"],
}


def match_regions(all_regions: list[str], patterns: list[str]) -> list[str]:
    matched = [
        region
        for region in all_regions
        if any(re.search(p, region, re.IGNORECASE) for p in patterns)
    ]
    return matched


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    # a = target region(s), b = everything else. PD and AD only map to one
    # region each so a often has n=1 and the normal pooled-variance formula
    # divides by zero df. When that happens just use b's std instead (basically
    # a z-score of the target against the rest of the brain).
    if len(a) <= 1 or len(b) <= 1:
        denom = b.std(ddof=1) if len(b) > 1 else a.std(ddof=1)
        return (a.mean() - b.mean()) / denom if denom > 0 else np.nan
    pooled_std = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / pooled_std if pooled_std > 0 else np.nan


def main() -> None:
    expression = pd.read_csv(EXPRESSION_PATH, index_col=0)
    all_regions = expression.index.tolist()

    rows = []
    disease_scores = {}
    for disease_path in sorted(GENE_LIST_DIR.glob("*_risk_genes.csv")):
        disease = disease_path.stem.replace("_risk_genes", "")
        genes_df = pd.read_csv(disease_path)
        genes = [g for g in genes_df["gene"].unique() if g in expression.columns]
        if not genes:
            print(f"SKIP {disease}: none of its risk genes found in expression matrix")
            continue

        # mean expression of this disease's risk genes, per region
        gene_score = expression[genes].mean(axis=1)
        disease_scores[disease] = gene_score

        target_regions = match_regions(all_regions, TARGET_REGION_PATTERNS[disease])
        if not target_regions:
            print(f"SKIP {disease}: no atlas region matched target pattern")
            continue

        target_vals = gene_score.loc[target_regions].to_numpy()
        other_vals = gene_score.drop(index=target_regions).to_numpy()

        rows.append(
            {
                "disease": disease,
                "n_genes": len(genes),
                "target_regions": ";".join(target_regions),
                "n_target_regions": len(target_regions),
                "target_mean_expression": target_vals.mean(),
                "other_mean_expression": other_vals.mean(),
                "difference": target_vals.mean() - other_vals.mean(),
                "cohens_d": cohens_d(target_vals, other_vals),
            }
        )

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(OUT_PATH, index=False)
    print(result.to_string(index=False))
    print(f"\nWrote -> {OUT_PATH}")

    # per-disease gene-score maps, needed by permutation_test.py next
    scores_path = ROOT / "data" / "processed" / "disease_gene_scores.csv"
    pd.DataFrame(disease_scores).to_csv(scores_path)
    print(f"Wrote per-region disease gene-score maps -> {scores_path}")


if __name__ == "__main__":
    main()

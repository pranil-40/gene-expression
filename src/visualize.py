"""
Step 7: one brain-region expression heat map per disease, plus the combined
results table (regional comparison + permutation p-values) as a figure.

Run:  .venv/Scripts/python src/visualize.py
Output: results/figures/<disease>_regional_heatmap.png (x3)
        results/figures/results_summary.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from regional_comparison import TARGET_REGION_PATTERNS, match_regions

ROOT = Path(__file__).resolve().parent.parent
SCORES_PATH = ROOT / "data" / "processed" / "disease_gene_scores.csv"
COMPARISON_PATH = ROOT / "results" / "tables" / "regional_comparison.csv"
PERMUTATION_PATH = ROOT / "results" / "tables" / "permutation_results.csv"
FIG_DIR = ROOT / "results" / "figures"


def plot_disease_heatmap(disease: str, gene_score: pd.Series, target_regions: list[str]) -> None:
    ordered = gene_score.sort_values(ascending=False)
    colors = ["#c0392b" if region in target_regions else "#4a6fa5" for region in ordered.index]

    fig, ax = plt.subplots(figsize=(8, max(6, len(ordered) * 0.18)))
    ax.barh(ordered.index, ordered.values, color=colors)
    ax.invert_yaxis()
    ax.set_xlabel("Mean risk-gene expression (AHBA, SRS-normalized)")
    ax.set_title(f"{disease.capitalize()}: regional risk-gene expression\n(red = disease target region)")
    ax.tick_params(axis="y", labelsize=6)
    fig.tight_layout()

    out_path = FIG_DIR / f"{disease}_regional_heatmap.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Wrote -> {out_path}")


def plot_results_summary(comparison: pd.DataFrame, permutation: pd.DataFrame) -> None:
    merged = comparison.merge(permutation, on="disease")
    fig, ax = plt.subplots(figsize=(8, 0.9 + 0.55 * len(merged)))
    ax.axis("off")
    table_cols = ["disease", "target_regions", "cohens_d", "p_value_spatial", "p_value_fdr_bh"]
    display = merged[table_cols].round(4)
    tbl = ax.table(cellText=display.values, colLabels=display.columns, loc="upper center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.8)
    ax.set_title("Regional enrichment vs. spatial null (BrainSMASH, N=1000)", pad=10)
    fig.tight_layout()

    out_path = FIG_DIR / "results_summary.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote -> {out_path}")


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(SCORES_PATH, index_col=0)

    for disease in scores.columns:
        target_regions = match_regions(scores.index.tolist(), TARGET_REGION_PATTERNS[disease])
        plot_disease_heatmap(disease, scores[disease], target_regions)

    comparison = pd.read_csv(COMPARISON_PATH)
    permutation = pd.read_csv(PERMUTATION_PATH)
    plot_results_summary(comparison, permutation)


if __name__ == "__main__":
    main()

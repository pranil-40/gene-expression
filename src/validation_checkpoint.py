"""
Step 4: validation checkpoint. Before trusting anything downstream, confirm
one known, single-gene result reproduces: SNCA (alpha-synuclein, the single
strongest Parkinson's GWAS hit — see data/gene_lists/parkinsons_risk_genes.csv)
should show relatively high expression in/near the substantia nigra region.

This is checked in isolation from the multi-gene risk-score analysis in
regional_comparison.py, specifically so the sanity check doesn't get diluted
or masked by averaging across all 10 PD risk genes.

If this does NOT hold, the pipeline (atlas labeling, gene symbol matching,
donor aggregation) has a bug — fix that before trusting Step 5/6 results for
any disease, per the project's quality bar.

Run:  .venv/Scripts/python src/validation_checkpoint.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from regional_comparison import match_regions

ROOT = Path(__file__).resolve().parent.parent
EXPRESSION_PATH = ROOT / "data" / "processed" / "regional_expression.csv"

CHECK_GENE = "SNCA"
TARGET_PATTERNS = [r"brain.?stem"]  # nearest available proxy — see regional_comparison.py


def main() -> None:
    expression = pd.read_csv(EXPRESSION_PATH, index_col=0)
    if CHECK_GENE not in expression.columns:
        raise SystemExit(f"FAIL: {CHECK_GENE} not found in expression matrix — check gene symbol matching")

    values = expression[CHECK_GENE].sort_values(ascending=False)
    target_regions = match_regions(values.index.tolist(), TARGET_PATTERNS)

    print(f"{CHECK_GENE} expression ranked across {len(values)} regions:")
    print(values.to_string())

    target_rank = [values.index.get_loc(r) + 1 for r in target_regions]
    target_percentile = 100 * (1 - min(target_rank) / len(values)) if target_rank else None

    print(f"\nTarget region(s) matched: {target_regions}")
    print(f"Rank of target region(s) (1 = highest expression): {target_rank}")
    if target_percentile is not None:
        print(f"Best target region is in the top {100 - target_percentile:.1f}% of regions by {CHECK_GENE} expression")

    if not target_regions:
        print("FAIL: no target region matched — check TARGET_PATTERNS / atlas region names")
    elif min(target_rank) <= len(values) // 3:
        print("PASS: target region falls in the top third of regions for SNCA expression — proceed to Step 5.")
    else:
        print("WARNING: target region is not in the top third for SNCA expression — "
              "investigate before trusting downstream results.")


if __name__ == "__main__":
    main()

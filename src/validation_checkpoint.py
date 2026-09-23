"""
Sanity check before trusting anything else: SNCA (alpha-synuclein, the
strongest PD GWAS hit by a mile) should be expressed relatively highly near
the substantia nigra. Checking this one gene alone instead of folding it into
the full 10-gene PD average, otherwise a bug could get averaged out and hide.

If this fails, something's wrong upstream (atlas labels, gene symbol
matching, whatever) -- don't trust regional_comparison.py or
permutation_test.py results until this passes.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from regional_comparison import match_regions

ROOT = Path(__file__).resolve().parent.parent
EXPRESSION_PATH = ROOT / "data" / "processed" / "regional_expression.csv"

CHECK_GENE = "SNCA"
TARGET_PATTERNS = [r"brain.?stem"]  # closest available proxy for substantia nigra


def main() -> None:
    expression = pd.read_csv(EXPRESSION_PATH, index_col=0)
    if CHECK_GENE not in expression.columns:
        raise SystemExit(f"FAIL: {CHECK_GENE} not found in expression matrix, check gene symbol matching")

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
        print("FAIL: no target region matched, check TARGET_PATTERNS / atlas region names")
    elif min(target_rank) <= len(values) // 3:
        print("PASS: target region is in the top third for SNCA expression, good to move on.")
    else:
        print("WARNING: target region is not in the top third for SNCA expression, "
              "look into this before trusting anything downstream.")


if __name__ == "__main__":
    main()

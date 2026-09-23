"""
Pull region-by-region gene expression for the risk-gene lists using abagen
(Allen Human Brain Atlas microarray data).

This step requires an internet connection and downloads several GB of AHBA
donor data on first run (abagen caches it locally after that). Run it in
Colab or a local environment with a good connection and disk space — see
README.md.

Run:  .venv/Scripts/python src/expression_pipeline.py
Output: data/processed/regional_expression.csv
        (rows = brain regions, columns = genes, values = mean microarray
        expression per abagen's donor-aggregation defaults)
"""
from __future__ import annotations

from pathlib import Path

import abagen
import nibabel as nib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
GENE_LIST_PATH = ROOT / "data" / "gene_lists" / "all_genes.csv"
OUT_PATH = ROOT / "data" / "processed" / "regional_expression.csv"
CENTROIDS_OUT_PATH = ROOT / "data" / "processed" / "region_centroids.csv"

# Desikan-Killiany atlas (bundled with abagen's fetch_desikan_killiany) —
# ~83 cortical + subcortical regions, well above the ~30-region target and a
# standard choice so results are comparable to published AHBA work.
ATLAS = "desikan_killiany"


def load_gene_list() -> list[str]:
    df = pd.read_csv(GENE_LIST_PATH)
    genes = sorted(df["gene"].unique().tolist())
    print(f"Loaded {len(genes)} unique genes from {GENE_LIST_PATH.name}: {genes}")
    return genes


def compute_region_centroids(atlas_image: str, atlas_info: pd.DataFrame) -> pd.DataFrame:
    """Mean voxel coordinate (world space, mm) per labeled region, keyed by
    region-name label. The DK atlas gives left/right hemispheres separate
    numeric ids under the *same* label (e.g. two ids both named
    "entorhinal") — averaged together here into one row per label, since
    that's also done for the expression matrix below (see main()), keeping
    both sides aligned on a guaranteed-unique index without relying on
    abagen's internal hemisphere/id bookkeeping.
    """
    img = nib.load(atlas_image)
    data = np.asarray(img.dataobj)
    affine = img.affine

    id_to_label = dict(zip(atlas_info["id"], atlas_info["label"]))
    rows = []
    for region_id, label in id_to_label.items():
        ijk = np.argwhere(data == region_id)
        if ijk.size == 0:
            continue
        centroid_ijk = ijk.mean(axis=0)
        centroid_xyz = nib.affines.apply_affine(affine, centroid_ijk)
        rows.append({"label": label, "x": centroid_xyz[0], "y": centroid_xyz[1], "z": centroid_xyz[2]})

    return pd.DataFrame(rows).groupby("label").mean()


def main() -> None:
    genes = load_gene_list()

    atlas = abagen.fetch_desikan_killiany()
    print("Pulling AHBA expression data via abagen (first run downloads donor "
          "data — this can take a while and several GB of disk)...")
    expression = abagen.get_expression_data(
        atlas["image"],
        atlas["info"],
        gene_norm="srs",
        region_agg="donors",
        missing="interpolate",
    )

    available = [g for g in genes if g in expression.columns]
    missing = sorted(set(genes) - set(available))
    if missing:
        print(f"WARNING: {len(missing)} genes not found in AHBA data: {missing}")

    # Collapse any left/right duplicate-labeled rows (see compute_region_centroids
    # docstring) so the region index is guaranteed unique before saving.
    result = expression[available].groupby(expression.index).mean()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT_PATH)
    print(f"Wrote regional expression matrix -> {OUT_PATH} "
          f"({result.shape[0]} regions x {result.shape[1]} genes)")

    centroids = compute_region_centroids(atlas["image"], atlas["info"])
    centroids = centroids.loc[centroids.index.intersection(result.index)]
    centroids.to_csv(CENTROIDS_OUT_PATH)
    print(f"Wrote region centroids -> {CENTROIDS_OUT_PATH}")


if __name__ == "__main__":
    main()

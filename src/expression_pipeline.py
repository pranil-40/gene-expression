"""
Gets AHBA expression values for the risk genes, region by region, using
abagen. First run downloads a few GB of donor data from the Allen Institute
(abagen caches it after that so it's only slow once). Needs internet.
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

# Desikan-Killiany atlas, comes bundled with abagen. ~83 regions total
# (cortical + subcortical), more than the ~30 the spec mentioned but it's
# the standard atlas people use with abagen so I stuck with it.
ATLAS = "desikan_killiany"


def load_gene_list() -> list[str]:
    df = pd.read_csv(GENE_LIST_PATH)
    genes = sorted(df["gene"].unique().tolist())
    print(f"Loaded {len(genes)} unique genes from {GENE_LIST_PATH.name}: {genes}")
    return genes


def compute_region_centroids(atlas_image: str, atlas_info: pd.DataFrame) -> pd.DataFrame:
    # mean voxel coordinate per region, in mm. Need this later for the
    # BrainSMASH distance matrix. DK atlas has separate L/R ids sharing the
    # same label (e.g. two different ids both called "entorhinal") so this
    # averages them into one row per label -- same thing happens to the
    # expression matrix in main() so the two line up.
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
    print("pulling AHBA expression via abagen, this takes a while on first run...")
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

    # collapse duplicate L/R region labels before saving (see compute_region_centroids)
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

"""
The DK atlas's "brainstem" label is a big undifferentiated blob (midbrain +
pons + medulla), and the validation checkpoint failed hard under it - SNCA
came out almost dead last across all 42 regions, which is backwards from
what the literature says. Checked the raw AHBA donor files
(abagen-data/microarray/*/SampleAnnot.csv) and they actually have real
substantia nigra tissue samples labeled directly in the Allen ontology, with
MNI coordinates - just buried inside "brainstem" along with a lot of other
stuff that dilutes any SN-specific signal.

This carves substantia nigra out as its own region: take every DK atlas
voxel currently labeled "brainstem", and if it's within RADIUS_MM of a real
SN tissue sample coordinate, relabel it to a new id. Leaves every other
region untouched. Saves a modified atlas image + info CSV that
expression_pipeline.py uses instead of the stock DK atlas.

python src/build_sn_atlas.py
"""
from __future__ import annotations

import glob
from pathlib import Path

import abagen
import nibabel as nib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT_IMAGE = ROOT / "data" / "processed" / "atlas_with_sn.nii.gz"
OUT_INFO = ROOT / "data" / "processed" / "atlas_with_sn_info.csv"

BRAINSTEM_LABEL = "brainstem"
SN_LABEL = "substantia_nigra"
RADIUS_MM = 6.0


def pooled_sn_coordinates() -> np.ndarray:
    # abagen caches donor data under ~/abagen-data by default, same place
    # expression_pipeline.py already downloaded it to
    pattern = str(Path.home() / "abagen-data" / "microarray" / "*" / "SampleAnnot.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit("Couldn't find any SampleAnnot.csv under ~/abagen-data - "
                          "run expression_pipeline.py at least once first so the donor data is cached.")

    coords = []
    for f in files:
        df = pd.read_csv(f)
        sn = df[df["structure_name"].str.contains("substantia nigra", case=False, na=False)]
        coords.append(sn[["mni_x", "mni_y", "mni_z"]].to_numpy())
    return np.vstack(coords)


def main() -> None:
    atlas = abagen.fetch_desikan_killiany()
    info = pd.read_csv(atlas["info"])
    img = nib.load(atlas["image"])
    data = np.asarray(img.dataobj).copy()
    affine = img.affine

    brainstem_id = int(info.loc[info["label"] == BRAINSTEM_LABEL, "id"].iloc[0])
    sn_coords = pooled_sn_coordinates()
    print(f"pooled {len(sn_coords)} real substantia nigra tissue samples across all donors")

    brainstem_voxels = np.argwhere(data == brainstem_id)
    inv_affine = np.linalg.inv(affine)
    voxel_mm = nib.affines.apply_affine(affine, brainstem_voxels)

    # distance from every brainstem voxel to the nearest real SN sample
    dists = np.linalg.norm(voxel_mm[:, None, :] - sn_coords[None, :, :], axis=2).min(axis=1)
    is_sn = dists <= RADIUS_MM

    new_id = int(info["id"].max()) + 1
    for voxel in brainstem_voxels[is_sn]:
        data[tuple(voxel)] = new_id

    n_sn_voxels = int(is_sn.sum())
    n_brainstem_left = int((~is_sn).sum())
    print(f"relabeled {n_sn_voxels} of {len(brainstem_voxels)} brainstem voxels as substantia_nigra "
          f"({n_brainstem_left} voxels stay labeled brainstem)")
    if n_sn_voxels == 0:
        raise SystemExit("No voxels got relabeled - RADIUS_MM might be too small, or a coordinate space mismatch.")

    new_info = pd.concat([
        info,
        pd.DataFrame([{"id": new_id, "label": SN_LABEL, "hemisphere": "B",
                        "structure": "subcortex/brainstem"}]),
    ], ignore_index=True)

    OUT_IMAGE.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(data, affine, img.header), OUT_IMAGE)
    new_info.to_csv(OUT_INFO, index=False)
    print(f"wrote {OUT_IMAGE} and {OUT_INFO}")


if __name__ == "__main__":
    main()

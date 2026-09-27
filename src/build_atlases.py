"""
Builds the atlases the v2 pipeline runs abagen on. All of them are the
Desikan-Killiany (DK) atlas restricted to the left hemisphere plus the
brainstem label. Left-only because only 2 of the 6 AHBA donors have
right-hemisphere samples, which is the usual reason AHBA studies drop the
right side.

  left_dk       DK left + brainstem. Base for the ontology-label SN variant
                (substantia nigra samples get pulled out by their AHBA
                anatomical label later, in expression.py).
  left_cit168   left_dk with substantia nigra pars compacta (SNc) and pars
                reticulata (SNr) taken from the CIT168 probabilistic atlas
                (Pauli et al. 2018), an anatomical definition that doesn't
                depend on where the AHBA samples happen to be.
  left_sn6mm    left_dk with the v1-style mask: brainstem voxels within 6 mm
                of an AHBA sample labelled substantia nigra. Kept only as an
                exploratory comparison to v1.

Right-side SN labels are kept as their own regions in left_cit168 and
left_sn6mm so right-side SN samples don't leak into "brainstem"; they get
dropped at the analysis stage.

python src/build_atlases.py
"""
from __future__ import annotations

import urllib.request

import abagen
import nibabel as nib
import numpy as np
import pandas as pd
from abagen import samples_

from common import ATLAS_DIR, RAW, DONORS, donor_dir

CIT168_URL = "https://osf.io/download/w8zq2/"   # CIT168toMNI152-2009c_prob.nii.gz, v1.1.0
CIT168_FILE = RAW / "cit168" / "CIT168toMNI152-2009c_prob.nii.gz"
# volume order in the 4D file, checked by label centroids (see results/atlas_report.txt)
CIT168_VOLUMES = ["Pu", "Ca", "NAC", "EXA", "GPe", "GPi", "SNc", "RN", "SNr",
                  "PBP", "VTA", "VeP", "HN", "HTH", "MN", "STH"]
SN_PROB_THRESHOLD = 0.5     # p(SNc) + p(SNr) >= this counts as substantia nigra
SPHERE_MM = 6.0             # v1 mask radius

BRAINSTEM_ID = 83
LEFT_IDS = range(1, 42)


def load_dk():
    atlas = abagen.fetch_desikan_killiany()
    img = nib.load(atlas["image"])
    info = pd.read_csv(atlas["info"])
    return img, np.asarray(img.dataobj).astype(np.int32), info


def left_only(data: np.ndarray, info: pd.DataFrame):
    keep = set(LEFT_IDS) | {BRAINSTEM_ID}
    out = np.where(np.isin(data, list(keep)), data, 0)
    return out, info[info["id"].isin(keep)].copy()


def voxel_world_coords(img, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ijk = np.argwhere(mask)
    return ijk, nib.affines.apply_affine(img.affine, ijk)


def cit168_on_dk_grid(dk_img) -> dict[str, np.ndarray]:
    if not CIT168_FILE.exists():
        CIT168_FILE.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading CIT168 -> {CIT168_FILE}")
        urllib.request.urlretrieve(CIT168_URL, CIT168_FILE)
    cit = nib.load(CIT168_FILE)
    prob = np.asarray(cit.dataobj)
    # both images are 1 mm MNI grids with identity rotation; only the origin
    # differs, so this is an integer shift, no interpolation
    shift = np.round(np.linalg.inv(cit.affine) @ dk_img.affine[:, 3])[:3].astype(int)
    assert np.allclose(np.abs(np.diag(cit.affine)[:3]), 1) and np.allclose(np.abs(np.diag(dk_img.affine)[:3]), 1)
    out = {}
    for name in ("SNc", "SNr"):
        vol = prob[..., CIT168_VOLUMES.index(name)]
        grid = np.zeros(dk_img.shape, dtype=np.float32)
        i, j, k = np.indices(dk_img.shape)
        ci, cj, ck = i + shift[0], j + shift[1], k + shift[2]
        ok = (ci >= 0) & (ci < vol.shape[0]) & (cj >= 0) & (cj < vol.shape[1]) & (ck >= 0) & (ck < vol.shape[2])
        grid[ok] = vol[ci[ok], cj[ok], ck[ok]]
        out[name] = grid
    return out


def ontology_sn_coords() -> pd.DataFrame:
    rows = []
    for uid in DONORS:
        annot = samples_.update_coords(str(donor_dir(uid) / "SampleAnnot.csv"), corrected_mni=True)
        sn = annot[annot["structure_name"].str.contains("substantia nigra", case=False, na=False)]
        rows.append(sn.assign(donor=uid))
    return pd.concat(rows, ignore_index=True)


def add_regions(info: pd.DataFrame, new: list[tuple[int, str, str]]) -> pd.DataFrame:
    extra = pd.DataFrame([{"id": i, "label": lab, "hemisphere": hemi, "structure": "subcortex/brainstem"}
                          for i, lab, hemi in new])
    return pd.concat([info, extra], ignore_index=True)


def save(name: str, img, data: np.ndarray, info: pd.DataFrame, report: list[str]):
    ATLAS_DIR.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(data.astype(np.int32), img.affine), ATLAS_DIR / f"{name}.nii.gz")
    info.to_csv(ATLAS_DIR / f"{name}_info.csv", index=False)
    counts = {int(i): int((data == i).sum()) for i in info["id"]}
    report.append(f"\n[{name}] {len(info)} regions")
    for _, r in info.iterrows():
        if r["id"] > 40 or r["label"] in ("entorhinal", "caudate", "putamen"):
            report.append(f"  id {r['id']:>3} {r['label']:<28} {r['hemisphere']}  {counts[int(r['id'])]:>6} voxels")


def main():
    report = []
    img, dk, info = load_dk()
    base, base_info = left_only(dk, info)
    save("left_dk", img, base, base_info, report)

    # CIT168 SNc / SNr
    cit = cit168_on_dk_grid(img)
    p_sn = cit["SNc"] + cit["SNr"]
    is_sn = p_sn >= SN_PROB_THRESHOLD
    ijk, xyz = voxel_world_coords(img, is_sn)
    dk_under = dk[is_sn]
    report.append("\nCIT168 SN voxels (p_SNc + p_SNr >= %.2f): %d" % (SN_PROB_THRESHOLD, is_sn.sum()))
    report.append("  DK label underneath them: " +
                  ", ".join(f"{info.set_index('id')['label'].get(k, 'unlabelled')}={v}"
                            for k, v in zip(*np.unique(dk_under, return_counts=True))))
    data = base.copy()
    snc = cit["SNc"][is_sn] >= cit["SNr"][is_sn]
    left = xyz[:, 0] < 0
    # only overwrite brainstem or unlabelled voxels, never another DK region
    writable = np.isin(dk_under, [0, BRAINSTEM_ID])
    for code, sel in [(101, left & snc), (102, left & ~snc), (103, ~left & snc), (104, ~left & ~snc)]:
        sel = sel & writable
        data[tuple(ijk[sel].T)] = code
    cit_info = add_regions(base_info, [(101, "substantia_nigra_pc", "L"), (102, "substantia_nigra_pr", "L"),
                                       (103, "substantia_nigra_pc_R", "R"), (104, "substantia_nigra_pr_R", "R")])
    save("left_cit168", img, data, cit_info, report)

    # v1-style 6 mm spheres around AHBA SN sample coordinates (now using the
    # corrected MNI coordinates abagen matches on; v1 used the original ones)
    sn_coords = ontology_sn_coords()[["mni_x", "mni_y", "mni_z"]].to_numpy()
    bs_ijk, bs_xyz = voxel_world_coords(img, dk == BRAINSTEM_ID)
    d = np.linalg.norm(bs_xyz[:, None, :] - sn_coords[None, :, :], axis=2).min(axis=1)
    near = d <= SPHERE_MM
    data = base.copy()
    data[tuple(bs_ijk[near & (bs_xyz[:, 0] < 0)].T)] = 105
    data[tuple(bs_ijk[near & (bs_xyz[:, 0] >= 0)].T)] = 106
    six_info = add_regions(base_info, [(105, "substantia_nigra_6mm", "L"), (106, "substantia_nigra_6mm_R", "R")])
    save("left_sn6mm", img, data, six_info, report)
    report.append(f"\n6 mm mask built from {len(sn_coords)} AHBA SN sample coordinates "
                  f"({near.sum()} of {len(bs_xyz)} brainstem voxels, {near.sum() / len(bs_xyz):.0%})")

    txt = "\n".join(report)
    (ATLAS_DIR / "atlas_report.txt").write_text(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()

"""
Anatomy-only audit of the substantia nigra (SN) tissue samples. No
expression values are touched here, so the choice of SN definition for v2
can be made before looking at any gene's expression in these samples.

For every AHBA sample labelled substantia nigra it records the original and
corrected MNI coordinates, the region abagen assigns it to under each v2
atlas (same matching code and 2 mm tolerance abagen uses), and the CIT168
SNc/SNr probability at its coordinate. Then for each SN-type region in each
atlas it lists which samples, by their anatomical label, end up inside it.

python src/sn_anatomy_audit.py
-> results/sn_audit/sn_samples.csv, sn_region_composition.csv, anatomy_summary.txt
"""
from __future__ import annotations

import nibabel as nib
import numpy as np
import pandas as pd
from abagen import images, samples_

from build_atlases import CIT168_FILE, CIT168_VOLUMES
from common import ATLAS_DIR, DONORS, RESULTS, donor_dir

OUT = RESULTS / "sn_audit"
ATLASES = ["left_dk", "left_cit168", "left_sn6mm"]
SN_WORDS = "substantia nigra"


def donor_annotation(uid: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(donor_dir(uid) / "SampleAnnot.csv")
    annot = samples_.update_coords(str(donor_dir(uid) / "SampleAnnot.csv"), corrected_mni=True)
    annot = samples_.drop_mismatch_samples(annot, str(donor_dir(uid) / "Ontology.csv"))
    return raw, annot


def subdivision(name: str) -> str:
    if "pars compacta" in name:
        return "SNc"
    if "pars reticulata" in name:
        return "SNr"
    return "SN (unspecified)"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    trees, infos = {}, {}
    for a in ATLASES:
        info = pd.read_csv(ATLAS_DIR / f"{a}_info.csv")
        trees[a] = images.check_atlas(str(ATLAS_DIR / f"{a}.nii.gz"), info)
        infos[a] = info.set_index("id")["label"].to_dict()

    cit = nib.load(CIT168_FILE)
    cprob = np.asarray(cit.dataobj)
    inv = np.linalg.inv(cit.affine)

    def cit_prob(xyz, name):
        ijk = np.round(nib.affines.apply_affine(inv, xyz)).astype(int)
        return cprob[ijk[:, 0], ijk[:, 1], ijk[:, 2], CIT168_VOLUMES.index(name)]

    sn_rows, all_labels = [], []
    for uid in DONORS:
        raw, annot = donor_annotation(uid)
        labels = {a: trees[a].label_samples(annot, tolerance=2)["label"].to_numpy() for a in ATLASES}
        lab_df = annot[["well_id", "structure_name"]].copy()
        lab_df["donor"] = uid
        for a in ATLASES:
            lab_df[a] = [infos[a].get(int(x), "unmatched") if x != 0 else "unmatched" for x in labels[a]]
        all_labels.append(lab_df)

        sn_raw = raw[raw["structure_name"].str.contains(SN_WORDS, case=False, na=False)]
        merged = sn_raw[["well_id", "structure_name", "mni_x", "mni_y", "mni_z"]].merge(
            lab_df.drop(columns="structure_name"), on="well_id", how="left")
        merged = merged.merge(annot[["well_id", "mni_x", "mni_y", "mni_z", "hemisphere"]]
                              .rename(columns={"mni_x": "x_corr", "mni_y": "y_corr", "mni_z": "z_corr"}),
                              on="well_id", how="left")
        merged = merged.rename(columns={"mni_x": "x_orig", "mni_y": "y_orig", "mni_z": "z_orig"})
        merged["donor"] = uid
        merged["kept_by_abagen"] = merged["x_corr"].notna()
        corr = merged[["x_corr", "y_corr", "z_corr"]].to_numpy()
        orig = merged[["x_orig", "y_orig", "z_orig"]].to_numpy()
        xyz = np.where(np.isnan(corr), orig, corr)
        merged["cit168_p_snc"] = cit_prob(xyz, "SNc")
        merged["cit168_p_snr"] = cit_prob(xyz, "SNr")
        merged["subdivision"] = merged["structure_name"].map(subdivision)
        merged["side"] = np.where(merged["structure_name"].str.contains("right"), "R", "L")
        merged["coord_shift_mm"] = np.linalg.norm(
            merged[["x_corr", "y_corr", "z_corr"]].to_numpy() - merged[["x_orig", "y_orig", "z_orig"]].to_numpy(), axis=1)
        sn_rows.append(merged)

    sn = pd.concat(sn_rows, ignore_index=True)
    cols = ["donor", "well_id", "structure_name", "subdivision", "side", "kept_by_abagen",
            "x_orig", "y_orig", "z_orig", "x_corr", "y_corr", "z_corr", "coord_shift_mm",
            "cit168_p_snc", "cit168_p_snr"] + ATLASES
    for a in ATLASES:
        sn[a] = sn[a].fillna("dropped (anatomy/coordinate mismatch)")
    sn[cols].sort_values(["donor", "side", "subdivision"]).to_csv(OUT / "sn_samples.csv", index=False)

    everything = pd.concat(all_labels, ignore_index=True)
    comp_rows = []
    for a in ATLASES:
        sn_regions = [lab for lab in infos[a].values() if lab.startswith("substantia_nigra") or lab == "brainstem"]
        if a == "left_dk":
            sn_regions += ["entorhinal", "caudate", "putamen"]
        for region in sn_regions:
            inside = everything[everything[a] == region]
            for name, n in inside["structure_name"].value_counts().items():
                comp_rows.append({"atlas": a, "region": region, "sample_label": name, "n_samples": n,
                                  "n_donors": inside.loc[inside["structure_name"] == name, "donor"].nunique()})
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(OUT / "sn_region_composition.csv", index=False)

    lines = [f"AHBA samples labelled substantia nigra: {len(sn)} across {sn['donor'].nunique()} donors"]
    lines.append(sn.groupby(["donor", "side", "subdivision"]).size().unstack(fill_value=0).to_string())
    lines.append(f"\nkept after abagen's anatomy/coordinate consistency check: {sn['kept_by_abagen'].sum()}")
    lines.append(f"coordinate correction moved samples by {sn['coord_shift_mm'].median():.1f} mm (median), "
                 f"max {sn['coord_shift_mm'].max():.1f} mm")
    inside = (sn["cit168_p_snc"] + sn["cit168_p_snr"]) >= 0.5
    lines.append(f"SN-labelled samples whose corrected coordinate is inside CIT168 SN (p >= 0.5): {inside.sum()} of {len(sn)}")
    lines.append(f"corrected coordinate ranges for SN-labelled samples: x {sn.x_corr.min():.0f} to {sn.x_corr.max():.0f}, "
                 f"y {sn.y_corr.min():.0f} to {sn.y_corr.max():.0f}, z {sn.z_corr.min():.0f} to {sn.z_corr.max():.0f}")
    for a in ATLASES:
        lines.append(f"\n[{a}] where the SN-labelled samples end up:")
        lines.append(sn.groupby(["side", a]).size().to_string())
    for a in ATLASES:
        sub = comp[(comp["atlas"] == a) & comp["region"].str.startswith("substantia_nigra")]
        for region, g in sub.groupby("region"):
            sn_n = g.loc[g["sample_label"].str.contains(SN_WORDS), "n_samples"].sum()
            lines.append(f"\n[{a}] {region}: {g['n_samples'].sum()} samples, {sn_n} labelled SN; top labels:")
            lines.append(g.sort_values("n_samples", ascending=False).head(8)[["sample_label", "n_samples", "n_donors"]].to_string(index=False))
    for region in ["entorhinal", "caudate", "putamen", "brainstem"]:
        g = comp[(comp["atlas"] == "left_dk") & (comp["region"] == region)]
        lines.append(f"\n[left_dk] {region}: {g['n_samples'].sum()} samples; top anatomical labels:")
        lines.append(g.sort_values("n_samples", ascending=False).head(8)[["sample_label", "n_samples", "n_donors"]].to_string(index=False))
    txt = "\n".join(lines)
    (OUT / "anatomy_summary.txt").write_text(txt + "\n", encoding="utf-8")
    print(txt)


if __name__ == "__main__":
    main()

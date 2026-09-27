"""
Figures for the v2 analysis. Each one answers a question from the review:

  fig1_sn_anatomy       where the AHBA substantia nigra samples actually are,
                        and what each SN definition captures
  fig2_snca_checkpoint  SNCA replication, and how normalization moves it
  fig3_region_scores    risk-gene score in every region, target highlighted
  fig4_gene_heatmap     which individual genes drive each set's score
  fig5_nulls            observed contrast against the spatial and gene-set
                        nulls, plus the BrainSMASH variogram fit
  fig6_donors           the contrast in each of the six donors
  fig7_sensitivity      the same contrast across every alternative choice
  fig8_power            how large an effect the tests could detect (post hoc,
                        docs/analysis_plan_amendment_1.md)
  fig9_celltypes        v3: which cell type expresses each risk gene
  fig10_dilution        v3: how much of the SN tissue signal comes from DA neurons

python src/figures.py
-> results/v2/figures/*.png
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from analysis import load_sets
from build_atlases import CIT168_FILE, CIT168_VOLUMES
from common import RESULTS, ahba_symbols
from expression import load_variant

OUT = RESULTS / "v2" / "figures"
V2 = RESULTS / "v2"

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                                                      "#256abf", "#184f95", "#0d366b"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": BASELINE, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10,
    "axes.titleweight": "semibold", "axes.titlelocation": "left", "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": False, "legend.frameon": False,
})

PRETTY = {"parkinsons": "Parkinson's (top-10 loci)", "alzheimers": "Alzheimer's (top-10 loci)",
          "huntingtons_htt": "Huntington's: HTT", "huntingtons_modifiers": "Huntington's: onset modifiers"}
TARGET_NAME = {"parkinsons": "SN pars compacta", "alzheimers": "entorhinal cortex",
               "huntingtons_htt": "caudate + putamen", "huntingtons_modifiers": "caudate + putamen"}


def hgrid(ax):
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def save(fig, name, out=OUT):
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out / name}.png")


def fig1_sn_anatomy():
    sn = pd.read_csv(RESULTS / "sn_audit" / "sn_samples.csv")
    comp = pd.read_csv(RESULTS / "sn_audit" / "sn_region_composition.csv")
    cit = nib.load(CIT168_FILE)
    prob = np.asarray(cit.dataobj)
    p_sn = prob[..., CIT168_VOLUMES.index("SNc")] + prob[..., CIT168_VOLUMES.index("SNr")]
    proj = (p_sn >= 0.5).any(axis=1)          # collapse the y axis -> x by z footprint
    xs = cit.affine[0, 3] + np.arange(proj.shape[0])
    zs = cit.affine[2, 3] + np.arange(proj.shape[1])

    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw={"width_ratios": [1, 1.25]},
                               layout="constrained")
    a.contourf(xs, zs, proj.T.astype(float), levels=[0.5, 1.5], colors=["#e8e7e2"])
    a.contour(xs, zs, proj.T.astype(float), levels=[0.5], colors=[BASELINE], linewidths=1)
    for sub, color, label in [("SNc", BLUE, "labelled pars compacta"), ("SNr", ORANGE, "labelled pars reticulata")]:
        s = sn[sn["subdivision"] == sub]
        a.scatter(s["x_corr"], s["z_corr"], s=30, color=color, edgecolor=SURFACE, linewidth=1.5, label=label, zorder=3)
    a.set_xlim(-32, 32); a.set_ylim(-36, 2)
    a.set_xlabel("x (mm, MNI; left is negative)"); a.set_ylabel("z (mm, MNI)")
    a.set_title("A  AHBA samples labelled substantia nigra")
    n_in = int(((sn["cit168_p_snc"] + sn["cit168_p_snr"]) >= 0.5).sum())
    a.text(0, -1, f"gray: CIT168 substantia nigra, flattened over y\n"
                  f"in 3D only {n_in} of {len(sn)} labelled samples fall inside it",
           ha="center", va="top", color=MUTED, fontsize=8)
    a.legend(loc="lower left", fontsize=8)

    rows = []
    for atlas, region, name in [("left_cit168", "substantia_nigra_pc", "CIT168 atlas: SNc"),
                                ("left_cit168", "substantia_nigra_pr", "CIT168 atlas: SNr"),
                                ("left_sn6mm", "substantia_nigra_6mm", "v1-style 6 mm spheres")]:
        g = comp[(comp["atlas"] == atlas) & (comp["region"] == region)]
        is_sn = g["sample_label"].str.contains("substantia nigra")
        rows.append((name, int(g.loc[is_sn, "n_samples"].sum()), int(g.loc[~is_sn, "n_samples"].sum())))
    left_sn = sn[(sn["side"] == "L") & sn["kept_by_abagen"]]
    rows.insert(0, ("AHBA anatomical labels: SNc", int((left_sn["subdivision"] == "SNc").sum()), 0))
    names = [r[0] for r in rows][::-1]
    yes = np.array([r[1] for r in rows][::-1]); no = np.array([r[2] for r in rows][::-1])
    y = np.arange(len(rows))
    b.barh(y, yes, height=0.55, color=BLUE, label="labelled substantia nigra")
    b.barh(y, no, left=yes + 0.4, height=0.55, color=BASELINE, label="other tissue")
    for yi, n1, n2 in zip(y, yes, no):
        b.text(n1 + n2 + 1.2, yi, f"{n1} of {n1 + n2}", va="center", color=INK2, fontsize=8)
    b.set_yticks(y, names); b.set_xlabel("tissue samples assigned to the region (left hemisphere)")
    b.set_title("B  What each substantia nigra definition captures")
    b.set_xlim(0, max(yes + no) * 1.25); hgrid(b)
    b.legend(loc="upper right", fontsize=8)
    save(fig, "fig1_sn_anatomy")


def fig2_snca():
    ranks = pd.read_csv(RESULTS / "snca_checkpoint" / "structure_ranks_allen_probes.csv", index_col=0)
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 3.9), gridspec_kw={"width_ratios": [1.6, 1]},
                               layout="constrained")
    v = ranks["median_log2"].sort_values(ascending=False).reset_index(drop=False)
    colors = [BLUE if "pars compacta" in s else ORANGE if "pars reticulata" in s else BASELINE for s in v.iloc[:, 0]]
    sizes = [42 if c != BASELINE else 12 for c in colors]
    a.scatter(np.arange(1, len(v) + 1), v["median_log2"], c=colors, s=sizes, edgecolor=SURFACE, linewidth=1, zorder=3)
    for i, s in enumerate(v.iloc[:, 0]):
        if "substantia nigra" in s:
            a.annotate(f"{'SNc' if 'compacta' in s else 'SNr'}: rank {i + 1} of {len(v)}", (i + 1, v['median_log2'][i]),
                       xytext=(12, 6 if 'compacta' in s else -14), textcoords="offset points", color=INK2, fontsize=8)
    a.set_xlabel("left-hemisphere structures sampled in all 6 donors, ranked")
    a.set_ylabel("SNCA, median log2 expression")
    a.set_title("A  SNCA by structure (raw AHBA values, all SNCA probes)")
    a.grid(axis="y", color=GRID, linewidth=0.8); a.set_axisbelow(True)

    modes = pd.read_csv(RESULTS / "snca_checkpoint" / "normalization_ranks.csv")
    y = np.arange(len(modes))[::-1]
    b.barh(y, modes["snc_rank"], height=0.5, color=[BLUE if r <= 5 else BASELINE for r in modes["snc_rank"]])
    for yi, r, n in zip(y, modes["snc_rank"], modes["n_regions"]):
        b.text(r + 0.6, yi, f"{r} of {n}", va="center", color=INK2, fontsize=8)
    b.set_yticks(y, modes["label"])
    b.set_xlabel("SNc rank for SNCA in the v2 region matrix (1 = highest)")
    b.set_title("B  Rank depends on normalization")
    hgrid(b); b.set_xlim(0, modes["n_regions"].max() + 6)
    save(fig, "fig2_snca_checkpoint")


def fig3_region_scores():
    scores = pd.read_csv(V2 / "scores" / "dk_ontology_sn__default.csv", index_col=0)
    res = pd.read_csv(V2 / "all_results.csv")
    sets = [s for s in ["parkinsons", "alzheimers", "huntingtons_htt", "huntingtons_modifiers"] if s in scores]
    fig, axes = plt.subplots(1, len(sets), figsize=(3.3 * len(sets), 7.5), sharey=False)
    for ax, s in zip(np.atleast_1d(axes), sets):
        r = res[(res["gene_set"] == s) & (res["region_definition"] == "dk_ontology_sn") & (res["normalization"] == "default")].iloc[0]
        tg = r["targets"].split(";")
        v = scores[s].sort_values()
        colors = [BLUE if reg in tg else BASELINE for reg in v.index]
        ax.barh(np.arange(len(v)), v.values, height=0.7, color=colors)
        ax.set_yticks(np.arange(len(v)), [x.replace("_", " ") for x in v.index], fontsize=6.5)
        lo, hi = v.min(), v.max()
        ax.set_xlim(lo - 0.08 * (hi - lo), hi + 0.08 * (hi - lo))
        ax.set_title(f"{PRETTY[s]}\nrank {int(r['best_target_rank'])} of {int(r['n_regions'])}, "
                     f"spatial p = {r['p_spatial']:.3f}", fontsize=9)
        ax.set_xlabel("mean normalized expression")
        hgrid(ax)
    fig.suptitle("Risk-gene score by region (blue = region the disease damages)", x=0.01, ha="left",
                 fontsize=11, fontweight="semibold")
    fig.tight_layout()
    save(fig, "fig3_region_scores")


def fig4_gene_heatmap():
    pooled, _, _ = load_variant("dk_ontology_sn")
    sets = load_sets()
    scores = pd.read_csv(V2 / "scores" / "dk_ontology_sn__default.csv", index_col=0)
    res = pd.read_csv(V2 / "all_results.csv")
    fig, axes = plt.subplots(2, 1, figsize=(11, 9.5), layout="constrained")
    for ax, s in zip(axes, ["parkinsons", "alzheimers"]):
        genes = [g for g in ahba_symbols(sets[s], pooled.columns).values() if g]
        order = scores[s].sort_values(ascending=False).index
        mat = pooled.loc[order, genes].T
        im = ax.imshow(mat.values, aspect="auto", cmap=SEQ, vmin=0, vmax=1)
        ax.set_yticks(range(len(genes)), genes, fontsize=7.5)
        ax.set_xticks(range(len(order)), [o.replace("_", " ") for o in order], rotation=90, fontsize=6.5)
        tg = res[(res["gene_set"] == s) & (res["kind"] == "primary")]["targets"].iloc[0].split(";")
        for t in tg:
            j = list(order).index(t)
            ax.add_patch(plt.Rectangle((j - 0.5, -0.5), 1, len(genes), fill=False, edgecolor=ORANGE, linewidth=2))
        ax.set_title(f"{PRETTY[s]}: each gene across regions, sorted by the set's score "
                     f"(orange box = {TARGET_NAME[s]})", fontsize=9)
    cb = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.01)
    cb.set_label("normalized expression (0 to 1)", color=INK2)
    save(fig, "fig4_gene_heatmap")


def fig5_nulls():
    res = pd.read_csv(V2 / "all_results.csv")
    prim = res[res["kind"] == "primary"]
    fig, axes = plt.subplots(len(prim), 3, figsize=(13, 3.3 * len(prim)), layout="constrained")
    for row, (_, r) in zip(np.atleast_2d(axes), prim.iterrows()):
        stem = f"{r['gene_set']}__{r['region_definition']}__{r['normalization']}"
        nulls = np.load(V2 / "nulls" / f"{stem}.npz")
        short = PRETTY[r["gene_set"]].split(" (")[0]
        for ax, key, title, p in [(row[0], "spatial", "spatial null", r["p_spatial"]),
                                  (row[1], "geneset", "gene-set null", r["p_geneset"])]:
            ax.hist(nulls[key], bins=60, color=BASELINE, edgecolor=SURFACE, linewidth=0.3)
            ax.axvline(r["T"], color=BLUE, linewidth=2)
            ax.text(r["T"], ax.get_ylim()[1] * 0.92, f"  observed\n  p = {p:.4f}", color=INK2, fontsize=8, va="top")
            ax.set_title(f"{short}: {title} (10,000)", fontsize=9)
            ax.set_xlabel("target minus rest (T)")
        vg = pd.read_csv(V2 / "variograms" / f"{stem}.csv")
        ax = row[2]
        ax.fill_between(vg["distance"], vg["surrogate_p05"], vg["surrogate_p95"], color="#cde2fb", label="surrogates, 5th to 95th pct")
        ax.plot(vg["distance"], vg["surrogate_mean"], color=BLUE, linewidth=2, label="surrogate mean")
        ax.scatter(vg["distance"], vg["empirical"], s=16, color=INK, zorder=3, label="observed map")
        ax.set_title(f"{short}: BrainSMASH variogram fit", fontsize=9)
        ax.set_xlabel("distance (mm)"); ax.set_ylabel("semivariance")
        ax.legend(fontsize=7)
    save(fig, "fig5_nulls")


def fig6_donors():
    d = pd.read_csv(V2 / "donor_contrasts.csv")
    res = pd.read_csv(V2 / "all_results.csv")
    sets = ["parkinsons", "alzheimers", "huntingtons_htt", "huntingtons_modifiers"]
    d = d[(d["region_definition"] == "dk_ontology_sn") & (d["normalization"] == "default") & d["gene_set"].isin(sets)]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    for i, s in enumerate(sets):
        g = d[(d["gene_set"] == s) & d["eligible"]]
        jitter = np.linspace(-0.18, 0.18, len(g)) if len(g) > 1 else [0]
        ax.scatter(g["T"], i + np.asarray(jitter), s=40, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3,
                   label="one donor" if i == 0 else None)
        pooled = res[(res["gene_set"] == s) & (res["region_definition"] == "dk_ontology_sn") &
                     (res["normalization"] == "default")]["T"].iloc[0]
        ax.plot([pooled, pooled], [i - 0.3, i + 0.3], color=INK, linewidth=2, label="all donors pooled" if i == 0 else None)
        n_pos, n = int((g["T"] > 0).sum()), len(g)
        ax.annotate(f"{n_pos} of {n} donors above 0", (1.01, i), xycoords=("axes fraction", "data"), va="center",
                    fontsize=8, color=INK2)
    ax.axvline(0, color=BASELINE, linewidth=1)
    ax.set_yticks(range(len(sets)), [PRETTY[s] for s in sets])
    ax.invert_yaxis()
    ax.set_xlabel("target minus rest of brain (T), per donor")
    ax.set_title("Is the pattern there in each person?")
    hgrid(ax); ax.legend(loc="lower right", fontsize=8)
    save(fig, "fig6_donors")


def fig7_sensitivity():
    res = pd.read_csv(V2 / "all_results.csv")
    res["label"] = res.apply(lambda r: f"{r['gene_set']} | {r['region_definition']} | {r['normalization']}", axis=1)
    order = ["primary", "secondary: gene set", "secondary: region definition", "exploratory: literature-named genes",
             "post hoc: normalization", "exploratory: v1 gene list", "exploratory: HD"]
    res["k"] = res["kind"].map({k: i for i, k in enumerate(order)})
    res = res.sort_values(["k", "gene_set"]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.5, 0.28 * len(res) + 1.2))
    y = np.arange(len(res))
    sig = res["p_spatial"] < 0.05
    ax.scatter(res.loc[sig, "z"], y[sig], s=46, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3, label="spatial p < 0.05")
    ax.scatter(res.loc[~sig, "z"], y[~sig], s=46, facecolor=SURFACE, edgecolor=MUTED, linewidth=1.5, zorder=3, label="spatial p >= 0.05")
    ax.axvline(0, color=BASELINE, linewidth=1)
    ax.set_yticks(y, res["label"], fontsize=7)
    ax.invert_yaxis()
    prev = None
    for i, k in enumerate(res["kind"]):
        if k != prev:
            ax.annotate(k, (1.01, i), xycoords=("axes fraction", "data"), fontsize=7.5, color=INK2, va="center")
            prev = k
    ax.set_xlabel("standardized target contrast z (target minus rest, in SDs of the other regions)")
    ax.set_title("Every run, every alternative choice")
    hgrid(ax); ax.legend(loc="lower right", fontsize=8)
    save(fig, "fig7_sensitivity")


def fig8_power():
    power_dir = V2 / "power"
    if not (power_dir / "random_set_power.csv").exists():
        print("skip fig8_power: run src/power_check.py first")
        return
    pw = pd.read_csv(power_dir / "random_set_power.csv")
    thr = pd.read_csv(power_dir / "detection_threshold_summary.csv").set_index("gene_set")
    pc = pd.read_csv(power_dir / "positive_controls.csv").set_index("gene_set")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True, layout="constrained")
    for ax, disease, letter in zip(axes, ["parkinsons", "alzheimers"], "AB"):
        g = pw[pw["disease"] == disease].sort_values("z_shifted")
        for col, color, label in [("geneset", BLUE, "gene-set null"), ("spatial", ORANGE, "spatial null"),
                                  ("both", INK, "both (needed for 'supported')")]:
            p, se = g[f"power_{col}"], g[f"se_{col}"]
            ax.fill_between(g["z_shifted"], (p - se).clip(0), (p + se).clip(upper=1), color=color, alpha=0.12, linewidth=0)
            ax.plot(g["z_shifted"], p, color=color, linewidth=2, marker="o", markersize=4, label=label)
        ax.axhline(0.8, color=BASELINE, linewidth=1, linestyle="--")
        z_obs = thr.loc[disease, "z_observed"]
        ax.axvline(z_obs, color=MUTED, linewidth=1)
        ax.annotate(f"observed\nz = {z_obs:.2f}", (z_obs, 0.97), xytext=(4, 0), textcoords="offset points",
                    va="top", fontsize=8, color=INK2)
        need = thr.loc[disease, "z_needed_both"]
        ax.annotate(f"this gene set would have\nneeded z = {need:.2f} to pass both", (need, 0.05),
                    xytext=(-6, 0), textcoords="offset points", ha="right", fontsize=8, color=INK2)
        ax.plot([need], [0.0], marker="^", markersize=8, color=INK, clip_on=False, zorder=4)
        if disease == "parkinsons" and "pc_dopaminergic" in pc.index:
            ax.annotate(f"dopaminergic marker genes: z = {pc.loc['pc_dopaminergic', 'z']:.1f} →", (1.0, 0.47),
                        xycoords=("axes fraction", "data"), ha="right", fontsize=8, color=INK2)
        ax.set_xlabel("target contrast z the gene set is given")
        ax.set_xticks([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])
        ax.set_title(f"{letter}  {PRETTY[disease].split(' (')[0]}: {TARGET_NAME[disease]}")
        ax.grid(axis="y", color=GRID, linewidth=0.8); ax.set_axisbelow(True)
        ax.set_ylim(0, 1.02)
    axes[0].set_ylabel("share of random 9-gene sets detected (p < 0.05)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=8)
    save(fig, "fig8_power")


def fig9_celltypes():
    v3 = RESULTS / "v3"
    if not (v3 / "gene_specificity.csv").exists():
        print("skip fig9_celltypes: run src/singlecell.py first")
        return
    gs = pd.read_csv(v3 / "gene_specificity.csv")
    tests = pd.read_csv(v3 / "specificity_tests.csv")
    panels = [("kamath_sn", "parkinsons", "DA neurons", "A  Parkinson's genes in substantia nigra (Kamath et al. 2022)"),
              ("siletti_ec", "alzheimers", "microglia", "B  Alzheimer's genes in entorhinal cortex (Siletti et al. 2023)")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), layout="constrained")
    for ax, (dname, sname, target, title) in zip(axes, panels):
        t = gs[(gs["dataset"] == dname) & (gs["gene_set"] == sname)].pivot(index="gene", columns="cell_class",
                                                                           values="specificity")
        cols = [target] + [c for c in t.columns if c != target]
        t = t[cols].sort_values(target, ascending=False)
        im = ax.imshow(t.values, cmap=SEQ, vmin=0, vmax=1, aspect="auto")
        for i in range(t.shape[0]):
            for j in range(t.shape[1]):
                v = t.values[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                        color=SURFACE if v > 0.55 else INK2)
        ax.set_xticks(range(len(cols)), cols, rotation=30, ha="right")
        ax.set_yticks(range(len(t)), t.index)
        ax.add_patch(plt.Rectangle((-0.5, -0.5), 1, len(t), fill=False, edgecolor=ORANGE, linewidth=2))
        r = tests[(tests["dataset"] == dname) & (tests["gene_set"] == sname) & (tests["kind"] == "primary")].iloc[0]
        ax.set_title(f"{title}\n{target}: rank {int(r['target_rank'])} of {int(r['n_classes'])}, "
                     f"p = {r['p']:.3f}, {r['verdict']}", fontsize=9)
        ax.spines[:].set_visible(False)
    cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.01)
    cb.set_label("share of the gene's expression, by cell class (rows sum to 1)", color=INK2)
    save(fig, "fig9_celltypes", out=RESULTS / "v3" / "figures")


def fig10_dilution():
    v3 = RESULTS / "v3"
    if not (v3 / "dilution.csv").exists():
        print("skip fig10_dilution: run src/singlecell.py first")
        return
    d = pd.read_csv(v3 / "dilution.csv")
    summ = (v3 / "dilution_summary.txt").read_text(encoding="utf-8")
    nuc = float(summ.split("(")[1].split("%")[0]) / 100
    umi = float(summ.strip().split(":")[-1])
    d = d.sort_values(["gene_set", "share_of_sn_umis_from_da"], ascending=[False, True]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8, 5.6))
    colors = [BLUE if s == "pc_dopaminergic" else ORANGE for s in d["gene_set"]]
    y = np.arange(len(d))
    ax.hlines(y, 1e-3, d["share_of_sn_umis_from_da"], color=GRID, linewidth=1, zorder=1)
    ax.scatter(d["share_of_sn_umis_from_da"], y, s=46, c=colors, edgecolor=SURFACE, linewidth=1.5, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1.3)
    for x, yy, label in [(nuc, 1.0, f"DA neurons' share of all SN nuclei ({100 * nuc:.2f}%)"),
                         (umi, 3.0, f"DA neurons' share of all SN RNA ({100 * umi:.1f}%)")]:
        ax.axvline(x, color=BASELINE, linewidth=1, linestyle="--")
        ax.text(x * 1.08, yy, label, fontsize=7.5, color=MUTED, va="center")
    ax.set_yticks(y, d["gene"])
    ax.set_xlabel("share of the gene's substantia nigra signal that comes from DA-neuron nuclei (log scale)")
    ax.set_title("Why bulk tissue sees dopamine markers but not dopamine-neuron effects on risk genes")
    ax.scatter([], [], c=BLUE, s=40, label="dopaminergic marker genes")
    ax.scatter([], [], c=ORANGE, s=40, label="Parkinson's risk genes (top 10 loci)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(axis="x", color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    fig.tight_layout()
    save(fig, "fig10_dilution", out=RESULTS / "v3" / "figures")


if __name__ == "__main__":
    fig1_sn_anatomy()
    fig2_snca()
    fig3_region_scores()
    fig4_gene_heatmap()
    fig5_nulls()
    fig6_donors()
    fig7_sensitivity()
    fig8_power()
    fig9_celltypes()
    fig10_dilution()

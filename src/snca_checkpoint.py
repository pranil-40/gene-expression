"""
SNCA checkpoint: tries to reproduce Liscovitch & French 2014 (PLoS ONE
9:e115029), who found SNCA ranked 2nd of 105 left-hemisphere AHBA
structures in substantia nigra pars compacta and 88th in pars reticulata.

Their comparison, as described in the paper: raw Allen log2 expression
averaged over the SNCA probes, mean per anatomical structure within each
donor, left-hemisphere structures sampled in all six donors, median across
donors. That's version A below. B repeats it with the single probe abagen
picks (differential stability). C gives SNc's rank in the v2 region matrix.

"Reproduced" (fixed in docs/analysis_plan.md): SNc in the top 10% of
structures and SNr below the median.

python src/snca_checkpoint.py
-> results/snca_checkpoint/*.csv, summary.txt
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import DONORS, RESULTS, donor_dir
from expression import load_raw, load_variant, normalize, sample_labels

OUT = RESULTS / "snca_checkpoint"
GENE = "SNCA"
SNC = "substantia nigra, pars compacta, left"
SNR = "substantia nigra, pars reticulata, left"


def probe_rows(uid: str, probe_ids: set[int], fname: str = "MicroarrayExpression.csv") -> pd.DataFrame:
    rows = {}
    with open(donor_dir(uid) / fname) as fh:
        for line in fh:
            pid, rest = line.split(",", 1)
            if int(pid) in probe_ids:
                rows[int(pid)] = np.array(rest.split(","), dtype=float)
    return pd.DataFrame(rows).T


def rank_structures(per_sample: dict[str, pd.Series], annot: dict[str, pd.DataFrame]) -> pd.DataFrame:
    per_donor = {}
    for uid, vals in per_sample.items():
        names = annot[uid].set_index("well_id")["structure_name"].reindex(vals.index)
        per_donor[uid] = vals.groupby(names.to_numpy()).mean()
    table = pd.DataFrame(per_donor)
    left = table[table.index.str.contains(", left")]
    complete = left.dropna()
    out = pd.DataFrame({"median_log2": complete.median(axis=1)}).join(complete)
    out["rank"] = out["median_log2"].rank(ascending=False).astype(int)
    return out.sort_values("rank")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    first = next(iter(DONORS))
    probes = pd.read_csv(donor_dir(first) / "Probes.csv")
    snca = probes[probes["gene_symbol"] == GENE]
    ids = set(snca["probe_id"].astype(int))
    annot = {uid: pd.read_csv(donor_dir(uid) / "SampleAnnot.csv") for uid in DONORS}

    # A: all Allen-mapped SNCA probes, raw log2, probe-averaged per sample
    per_probe, per_sample_a = {}, {}
    for uid in DONORS:
        rows = probe_rows(uid, ids)
        rows.columns = annot[uid]["well_id"].to_numpy()
        per_probe[uid] = rows
        per_sample_a[uid] = rows.mean(axis=0)
    rank_a = rank_structures(per_sample_a, annot)
    rank_a.to_csv(OUT / "structure_ranks_allen_probes.csv")

    # per-probe ranks, to see whether the probes agree
    probe_summary = []
    for pid in sorted(ids):
        r = rank_structures({u: per_probe[u].loc[pid] for u in DONORS}, annot)
        name = snca.set_index("probe_id").loc[pid, "probe_name"]
        probe_summary.append({"probe_id": pid, "probe_name": name, "n_structures": len(r),
                              "snc_rank": int(r.loc[SNC, "rank"]) if SNC in r.index else None,
                              "snr_rank": int(r.loc[SNR, "rank"]) if SNR in r.index else None})

    # B: abagen's selected probe (cached, not normalized)
    per_sample_b = {uid: load_raw(uid)[GENE] for uid in DONORS}
    rank_b = rank_structures(per_sample_b, annot)
    rank_b.to_csv(OUT / "structure_ranks_abagen_probe.csv")
    # which Allen probe abagen picked: the one its values correlate with best
    corr = {pid: np.mean([np.corrcoef(per_probe[u].loc[pid, per_sample_b[u].index],
                                       per_sample_b[u].to_numpy())[0, 1] for u in DONORS]) for pid in ids}
    for p in probe_summary:
        p["corr_with_abagen_probe"] = round(corr[p["probe_id"]], 4)
    pd.DataFrame(probe_summary).to_csv(OUT / "probe_summary.csv", index=False)

    # donor-level SNc vs SNr, raw probe-averaged
    donor_rows = []
    for uid in DONORS:
        names = annot[uid].set_index("well_id")["structure_name"]
        v = per_sample_a[uid]
        donor_rows.append({"donor": uid,
                           "snc_n": int((names.reindex(v.index) == SNC).sum()),
                           "snr_n": int((names.reindex(v.index) == SNR).sum()),
                           "snc_mean": v[names.reindex(v.index).eq(SNC).to_numpy()].mean(),
                           "snr_mean": v[names.reindex(v.index).eq(SNR).to_numpy()].mean(),
                           "all_samples_mean": v.mean()})
    donors = pd.DataFrame(donor_rows)
    donors["snc_minus_snr"] = donors["snc_mean"] - donors["snr_mean"]
    donors.to_csv(OUT / "donor_snc_snr.csv", index=False)

    # C: v2 primary region matrix
    pooled, _, _ = load_variant("dk_ontology_sn")
    s = pooled[GENE].sort_values(ascending=False)
    c_rank = int(s.index.get_loc("substantia_nigra_pc") + 1)
    c_rank_r = int(s.index.get_loc("substantia_nigra_pr") + 1)

    # post hoc (not in the plan): C came out much lower than A/B, so check
    # which abagen normalization step moves SNc's rank
    labels = sample_labels("dk_ontology_sn")
    norm_rows = []
    for label, sn_, gn_ in [("no normalization", None, None), ("gene SRS only", None, "srs"),
                            ("sample SRS only", "srs", None), ("both (abagen default)", "srs", "srs")]:
        per_donor = []
        for uid in DONORS:
            raw = load_raw(uid)
            lab = labels[uid].reindex(raw.index).fillna("")
            keep = lab.ne("").to_numpy()
            x = normalize(raw.loc[keep], sn_, gn_)
            per_donor.append(x[GENE].groupby(lab[keep].to_numpy()).mean())
        v = pd.concat(per_donor, axis=1).mean(axis=1)
        v = v.drop([r for r in v.index if r.endswith("_R")]).sort_values(ascending=False)
        norm_rows.append({"label": label, "snc_rank": int(v.index.get_loc("substantia_nigra_pc") + 1),
                          "snr_rank": int(v.index.get_loc("substantia_nigra_pr") + 1), "n_regions": len(v)})
    norm_df = pd.DataFrame(norm_rows)
    norm_df.to_csv(OUT / "normalization_ranks.csv", index=False)

    def verdict(r):
        n = len(r)
        snc, snr = r.loc[SNC, "rank"], r.loc[SNR, "rank"]
        ok = snc <= np.ceil(0.10 * n) and snr > n / 2
        return n, snc, snr, "reproduced" if ok else "NOT reproduced"

    lines = ["SNCA checkpoint (target: Liscovitch & French 2014, SNc 2/105, SNr 88/105)", ""]
    for label, r in [("A  all Allen SNCA probes, raw log2", rank_a), ("B  abagen-selected probe, raw log2", rank_b)]:
        n, snc, snr, v = verdict(r)
        lines.append(f"{label}: {n} left structures in all 6 donors; SNc rank {snc}, SNr rank {snr} -> {v}")
    lines.append(f"C  v2 region matrix (dk_ontology_sn, SRS-normalized): SNc rank {c_rank} of {len(s)}, "
                 f"SNr rank {c_rank_r} of {len(s)}")
    lines.append("\npost hoc: SNc rank for SNCA in the v2 region matrix under each normalization step")
    lines.append(norm_df.to_string(index=False))
    lines.append("\nper-probe (raw log2):")
    lines.append(pd.DataFrame(probe_summary).to_string(index=False))
    lines.append("\ndonor-level SNc vs SNr (raw log2, probe-averaged):")
    lines.append(donors.round(3).to_string(index=False))
    lines.append(f"\ntop 10 structures (A):")
    lines.append(rank_a.head(10)[["median_log2", "rank"]].round(3).to_string())
    txt = "\n".join(lines)
    (OUT / "summary.txt").write_text(txt + "\n", encoding="utf-8")
    print(txt)


if __name__ == "__main__":
    main()

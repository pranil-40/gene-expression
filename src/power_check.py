"""
Positive controls and power, added after the v2 results were known
(docs/analysis_plan_amendment_1.md). None of this changes a v2 verdict.

  A. marker gene sets with known regional specificity, through the primary pipeline
  B. how high the primary PD and AD sets' target contrast would have had to be
     before the planned tests called it
  C. power of the same tests for random 9-gene sets shifted to a given contrast

python src/power_check.py
-> results/v2/power/*.csv
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from analysis import N_SETS, OUT, Context, load_sets, run, spatial_null, targets
from common import ahba_symbols, empirical_p, seed_for

POWER = OUT / "power"
VARIANT, NORM = "dk_ontology_sn", "default"

CONTROLS = {
    # name: (target regions, marker genes)
    "pc_dopaminergic": (["substantia_nigra_pc"],
                        ["TH", "SLC6A3", "SLC18A2", "DDC", "KCNJ6", "NR4A2", "EN1", "FOXA2", "LMX1B", "PITX3"]),
    "pc_striatal_msn": (["caudate", "putamen"],
                        ["PPP1R1B", "DRD1", "DRD2", "PENK", "TAC1", "ADORA2A", "GPR88", "RGS9", "PDE10A", "RASD2"]),
}

Z_GRID_B = np.round(np.arange(1.0, 5.01, 0.25), 2)
N_SURR_B = 2_000

Z_GRID_C = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]
K_C = 9
N_DRAWS_GS = 1_000
N_DRAWS_SP = 60
N_SURR_C = 500

ALPHA = 0.05


def shift_to(S: np.ndarray, mask: np.ndarray, z: float) -> np.ndarray:
    """Copy of S with the target regions moved so the contrast z equals `z`.

    The same constant goes on every target region, so the reference regions,
    and their SD, don't change.
    """
    ref = S[~mask]
    T_now = S[mask].mean() - ref.mean()
    out = S.copy()
    out[mask] += z * ref.std(ddof=1) - T_now
    return out


def geneset_null(ctx: Context, mask: np.ndarray, k: int, seed: int) -> np.ndarray:
    # same construction as analysis.run, so with the same seed it's the same null
    d_gene = ctx.M[mask].mean(axis=0) - ctx.M[~mask].mean(axis=0)
    rng = np.random.default_rng(seed)
    return np.array([d_gene[rng.choice(len(d_gene), size=k, replace=False)].mean() for _ in range(N_SETS)])


def best_rank(S: np.ndarray, mask: np.ndarray) -> int:
    ranks = pd.Series(S).rank(ascending=False).to_numpy()
    return int(ranks[mask].min())


def positive_controls(ctx: Context) -> None:
    rows, donors = [], []
    for name, (tgt, genes) in CONTROLS.items():
        t0 = time.time()
        res, extra = run("positive control", name, genes, ctx, VARIANT, NORM, detail=False, tgt=tgt)
        rows.append(res)
        donors.append(extra["donors"])
        print(f"A  {name:18s} genes={res['n_genes_used']}/{res['n_genes_requested']} z={res['z']:+.2f} "
              f"rank={res['best_target_rank']}/{res['n_regions']} p_sp={res['p_spatial']:.4f} "
              f"p_gs={res['p_geneset']:.4f} p_cam={res['p_camera']:.2g} "
              f"donors+={res['donors_positive']}/{res['donors_eligible']} ({time.time() - t0:.0f}s)", flush=True)
    pd.DataFrame(rows).to_csv(POWER / "positive_controls.csv", index=False)
    pd.concat(donors).to_csv(POWER / "positive_control_donors.csv", index=False)


def detection_threshold(ctx: Context, sets: dict[str, list[str]]) -> None:
    rows = []
    for set_name in ("parkinsons", "alzheimers"):
        mask = np.isin(ctx.regions, targets(set_name, VARIANT))
        sym = ahba_symbols(sets[set_name], ctx.genes)
        genes = list(dict.fromkeys(s for s in sym.values() if s))
        S = ctx.pooled[genes].mean(axis=1).to_numpy()
        sd = S[~mask].std(ddof=1)
        null_gs = geneset_null(ctx, mask, len(genes), seed_for(f"{set_name}|{VARIANT}|{NORM}|genesets"))
        z_obs = (S[mask].mean() - S[~mask].mean()) / sd
        for z in Z_GRID_B:
            Sz = shift_to(S, mask, z)
            T = z * sd
            null_sp = spatial_null(Sz, ctx.D, mask, N_SURR_B, seed_for(f"{set_name}|threshold|{z}"))
            rows.append({"gene_set": set_name, "n_genes": len(genes), "z_observed": z_obs, "z_shifted": z,
                         "T_shifted": T, "target_rank": best_rank(Sz, mask),
                         "p_spatial": empirical_p(null_sp, T), "p_geneset": empirical_p(null_gs, T)})
            r = rows[-1]
            print(f"B  {set_name:12s} z*={z:4.2f} rank={r['target_rank']:2d} p_sp={r['p_spatial']:.4f} "
                  f"p_gs={r['p_geneset']:.4f}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(POWER / "detection_threshold.csv", index=False)

    summary = []
    for set_name, g in df.groupby("gene_set", sort=False):
        g = g.sort_values("z_shifted")
        out = {"gene_set": set_name, "z_observed": g["z_observed"].iloc[0]}
        for col, ok in [("spatial", g["p_spatial"] < ALPHA), ("geneset", g["p_geneset"] < ALPHA),
                        ("both", (g["p_spatial"] < ALPHA) & (g["p_geneset"] < ALPHA))]:
            # smallest z* from which p stays under ALPHA for every larger z* on the grid
            stays = ok.astype(int)[::-1].cummin()[::-1].astype(bool)
            out[f"z_needed_{col}"] = g.loc[stays, "z_shifted"].min() if stays.any() else np.nan
        summary.append(out)
    pd.DataFrame(summary).to_csv(POWER / "detection_threshold_summary.csv", index=False)
    print(pd.DataFrame(summary).to_string(index=False), flush=True)


def random_set_power(ctx: Context) -> None:
    rows = []
    n_genes = ctx.M.shape[1]
    for disease in ("parkinsons", "alzheimers"):
        tgt = targets(disease, VARIANT)
        mask = np.isin(ctx.regions, tgt)
        null_gs = geneset_null(ctx, mask, K_C, seed_for(f"power|{disease}|null"))
        rng = np.random.default_rng(seed_for(f"power|{disease}|draws"))
        maps = [ctx.M[:, rng.choice(n_genes, size=K_C, replace=False)].mean(axis=1) for _ in range(N_DRAWS_GS)]
        sds = np.array([S[~mask].std(ddof=1) for S in maps])
        for z in Z_GRID_C:
            t0 = time.time()
            p_gs = np.array([empirical_p(null_gs, z * sd) for sd in sds])
            p_sp = np.array([empirical_p(spatial_null(shift_to(S, mask, z), ctx.D, mask, N_SURR_C,
                                                      seed_for(f"power|{disease}|{z}|{i}")), z * sds[i])
                             for i, S in enumerate(maps[:N_DRAWS_SP])])
            hit_gs = p_gs < ALPHA
            hit_sp = p_sp < ALPHA
            hit_both = hit_gs[:N_DRAWS_SP] & hit_sp
            row = {"disease": disease, "target": ";".join(tgt), "k": K_C, "z_shifted": z,
                   "n_draws_geneset": len(p_gs), "n_draws_spatial": len(p_sp)}
            for col, hits in [("geneset", hit_gs), ("spatial", hit_sp), ("both", hit_both)]:
                p = hits.mean()
                row[f"power_{col}"] = p
                row[f"se_{col}"] = np.sqrt(p * (1 - p) / len(hits))
            rows.append(row)
            print(f"C  {disease:12s} z*={z:3.1f} power: geneset={row['power_geneset']:.2f} "
                  f"spatial={row['power_spatial']:.2f} both={row['power_both']:.2f} ({time.time() - t0:.0f}s)",
                  flush=True)
    pd.DataFrame(rows).to_csv(POWER / "random_set_power.csv", index=False)


def main():
    POWER.mkdir(parents=True, exist_ok=True)
    ctx = Context(VARIANT, NORM)
    positive_controls(ctx)
    detection_threshold(ctx, load_sets())
    random_set_power(ctx)


if __name__ == "__main__":
    main()

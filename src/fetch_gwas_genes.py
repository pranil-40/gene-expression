"""
Builds the v2 risk-gene lists and a traceable audit table covering every
association considered.

Changes from v1 (v1 code: commit 220b838; v1 lists: archive/v1_exploratory):
  - Only genome-wide significant associations (p < 5e-8) count for PD/AD.
  - Associations are grouped into independent loci first, so one locus can't
    contribute several genes. Two lead variants are the same locus if they
    sit on the same chromosome within 1 Mb and are in LD (r^2 >= 0.1,
    1000 Genomes EUR, Ensembl REST). If the LD lookup fails for a pair, the
    fallback is distance only (same locus if within 250 kb), and the audit
    table says which rule decided.
  - Each locus gets one gene from GENCODE v50: a protein-coding gene whose
    body contains the lead variant, otherwise the nearest protein-coding gene
    within 1 Mb. Readthrough genes (GENCODE tag) are skipped. Ties go to an
    author-reported gene, then alphabetical. Pseudogenes and lncRNAs can't be
    picked any more (v1 filtered them by name pattern and missed some).
  - Primary set = top 10 loci by lead p-value. Sensitivity sets = all
    genome-wide significant loci, and author-reported genes where the
    Catalog has them.
  - Huntington's is split off (see HD_SETS): HTT on its own, and the onset
    modifier genes named by the GeM-HD consortium. Those aren't curated in the
    Catalog, so they're typed in by hand with the source recorded.

Gene annotation used to come from the Catalog's per-variant records, but that
endpoint rate-limited this script to a crawl (2026-09-24), so positions now
come from one batched Ensembl request and genes from a pinned local GENCODE
file. That switch was made before any v2 gene list existed.

python src/fetch_gwas_genes.py
-> data/gene_lists/gene_audit.csv, data/gene_lists/sets/*.csv, data/gene_lists/provenance.json
"""
from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import json
import re
import time
import urllib.request

import numpy as np
import pandas as pd
import requests

from common import CACHE, GENE_DIR, RAW

API = "https://www.ebi.ac.uk/gwas/rest/api"
ENSEMBL = "https://rest.ensembl.org"
GENCODE_URL = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/"
               "gencode.v50.basic.annotation.gtf.gz")
GENCODE_GTF = RAW / "gencode" / "gencode.v50.basic.annotation.gtf.gz"
API_CACHE = CACHE / "api"

GWS = 5e-8
WINDOW_BP = 1_000_000
FALLBACK_BP = 250_000
LD_R2 = 0.1
LD_POP = "1000GENOMES:phase_3:EUR"
TOP_N = 10

STUDIES = {
    "parkinsons": {"accession": "GCST009325", "pubmed_id": "31701892",
                   "citation": "Nalls et al. 2019, Lancet Neurology 18:1091-1102"},
    "alzheimers": {"accession": "GCST90027158", "pubmed_id": "35379992",
                   "citation": "Bellenguez et al. 2022, Nature Genetics 54:412-436"},
}

HD_SETS = {
    "huntingtons_htt": {
        "genes": ["HTT"],
        "source": "HTT CAG-repeat expansion causes HD (Huntington's Disease Collaborative Research Group 1993, Cell 72:971-983)",
    },
    "huntingtons_modifiers": {
        # one gene per genome-wide significant onset-modifier locus; where the
        # paper names two candidates (MSH3/DHFR, RRM2B/UBR5) the one listed in
        # the 2022 follow-up summary is used
        "genes": ["FAN1", "MSH3", "MLH1", "PMS1", "PMS2", "LIG1", "RRM2B", "TCERG1", "CCDC82"],
        "source": ("GeM-HD Consortium 2019, Cell 178:887-900 (Table 1, age-at-onset modifier loci); "
                   "gene choices as summarized in Lee et al. 2022, Am J Hum Genet 109:885-899"),
    },
}

session = requests.Session()
_last = {"t": 0.0}


class ApiError(RuntimeError):
    pass


def _cache_path(key: str):
    return API_CACHE / f"{hashlib.sha1(key.encode()).hexdigest()}.json"


def request_json(method: str, url: str, params=None, body=None, tries=8, gap=1.0):
    """HTTP with an on-disk cache and patient retries. Raises instead of
    returning nothing, so throttling can't silently drop records."""
    path = _cache_path(method + url + json.dumps(params or {}, sort_keys=True) + json.dumps(body or {}, sort_keys=True))
    if path.exists():
        return json.loads(path.read_text())
    delay = 10
    for _ in range(tries):
        wait = gap - (time.time() - _last["t"])
        if wait > 0:
            time.sleep(wait)
        _last["t"] = time.time()
        try:
            r = session.request(method, url, params=params, json=body, timeout=90,
                                headers={"Accept": "application/json", "Content-Type": "application/json"})
        except requests.RequestException:
            time.sleep(delay); delay = min(delay * 2, 300); continue
        if r.status_code == 200:
            data = r.json()
            API_CACHE.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(data))
            return data
        if r.status_code in (429, 500, 502, 503, 504):
            ra = r.headers.get("Retry-After")
            time.sleep(max(delay, int(ra)) if ra and ra.isdigit() else delay)
            delay = min(delay * 2, 300)
            continue
        raise ApiError(f"{r.status_code} from {url}")
    raise ApiError(f"gave up after {tries} tries: {url}")


def load_gencode() -> pd.DataFrame:
    table = GENCODE_GTF.with_name("genes_v50.csv")
    if table.exists():
        return pd.read_csv(table)
    if not GENCODE_GTF.exists():
        GENCODE_GTF.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(GENCODE_URL, GENCODE_GTF)
    rows = []
    with gzip.open(GENCODE_GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] != "gene":
                continue
            a = f[8]
            rows.append({"chr": f[0].replace("chr", ""), "start": int(f[3]), "end": int(f[4]),
                         "gene_id": re.search(r'gene_id "([^"]+)"', a).group(1).split(".")[0],
                         "gene": re.search(r'gene_name "([^"]+)"', a).group(1),
                         "gene_type": re.search(r'gene_type "([^"]+)"', a).group(1),
                         "readthrough": 'tag "readthrough_gene"' in a})
    g = pd.DataFrame(rows)
    g.to_csv(table, index=False)
    return g


def associations(accession: str) -> list[dict]:
    data = request_json("GET", f"{API}/studies/{accession}/associations", gap=2.0)
    rows = []
    for a in data["_embedded"]["associations"]:
        snps = [al["riskAlleleName"].split("-")[0] for l in a["loci"] for al in l["strongestRiskAlleles"]]
        rsid = next((s for s in snps if s.startswith("rs")), None)
        rows.append({
            "association_id": a["_links"]["self"]["href"].rsplit("/", 1)[-1],
            "lead_snp": rsid,
            "pvalue": a.get("pvalue"),
            "author_reported_genes": ";".join(sorted({g["geneName"] for l in a["loci"]
                                                     for g in l.get("authorReportedGenes", [])})),
        })
    return rows


def ncbi_position(rsid: str) -> tuple[str, int] | None:
    """GRCh38 position from NCBI dbSNP, used when Ensembl is unavailable."""
    try:
        d = request_json("GET", f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid[2:]}", tries=5, gap=0.4)
    except ApiError:
        return None
    for pl in d.get("primary_snapshot_data", {}).get("placements_with_allele", []):
        traits = pl.get("placement_annot", {}).get("seq_id_traits_by_assembly", [])
        if not pl.get("is_ptlp") or not any(t.get("assembly_name", "").startswith("GRCh38") for t in traits):
            continue
        m = re.match(r"NC_0+(\d+)\.", pl["seq_id"])
        if not m:
            continue
        chrom = {"23": "X", "24": "Y"}.get(m.group(1), m.group(1))
        spdi = pl["alleles"][0]["allele"]["spdi"]
        return chrom, int(spdi["position"]) + 1
    return None


def positions(rsids: list[str]) -> dict[str, tuple[str, int]]:
    out = {}
    ids = sorted(set(rsids))
    for i in range(0, len(ids), 150):
        try:
            data = request_json("POST", f"{ENSEMBL}/variation/homo_sapiens", body={"ids": ids[i:i + 150]},
                                tries=3, gap=1.0)
        except ApiError:
            continue
        for rsid, rec in data.items():
            maps = [m for m in rec.get("mappings", []) if m.get("assembly_name") == "GRCh38"
                    and re.fullmatch(r"[0-9XY]+", str(m.get("seq_region_name")))]
            if maps:
                out[rsid] = (str(maps[0]["seq_region_name"]), int(maps[0]["start"]))
    missing = [r for r in ids if r not in out]
    if missing:
        print(f"  {len(missing)} positions not from Ensembl, trying NCBI dbSNP", flush=True)
    for rsid in missing:
        p = ncbi_position(rsid)
        if p:
            out[rsid] = p
    return out


def ld_r2(a: str, b: str) -> float | None:
    try:
        res = request_json("GET", f"{ENSEMBL}/ld/human/pairwise/{a}/{b}",
                           params={"population_name": LD_POP}, tries=5, gap=0.5)
    except ApiError:
        return None
    if not res:
        return 0.0           # Ensembl returns nothing when the pair isn't in LD
    return float(res[0]["r2"])


def clump(df: pd.DataFrame) -> pd.DataFrame:
    df = df.sort_values("pvalue").reset_index(drop=True)
    df["locus_id"], df["merged_by"] = None, ""
    leads = []
    for i, row in df.iterrows():
        home = None
        for j in leads:
            lead = df.loc[j]
            dist = abs(row["pos_grch38"] - lead["pos_grch38"])
            if row["chr"] != lead["chr"] or dist > WINDOW_BP:
                continue
            r2 = ld_r2(lead["lead_snp"], row["lead_snp"])
            if r2 is None:
                if dist <= FALLBACK_BP:
                    home = j
                    df.loc[i, "merged_by"] = f"distance fallback ({dist / 1000:.0f} kb from {lead['lead_snp']}, LD lookup failed)"
                    break
                continue
            if r2 >= LD_R2:
                home = j
                df.loc[i, "merged_by"] = f"r2={r2:.2f} with {lead['lead_snp']}"
                break
        if home is None:
            leads.append(i)
            df.loc[i, "locus_id"] = f"chr{row['chr']}:{int(row['pos_grch38'])}"
        else:
            df.loc[i, "locus_id"] = df.loc[home, "locus_id"]
    df["is_locus_lead"] = df.index.isin(leads)
    return df


def pick_gene(chrom: str, pos: int, author: set[str], genes: pd.DataFrame):
    pc = genes[(genes["chr"] == chrom) & (genes["gene_type"] == "protein_coding") & ~genes["readthrough"]]
    near = pc[(pc["end"] >= pos - WINDOW_BP) & (pc["start"] <= pos + WINDOW_BP)].copy()
    if near.empty:
        return None, "no protein-coding gene within 1 Mb", None, "", ""
    near["dist"] = np.where((near["start"] <= pos) & (near["end"] >= pos), 0,
                            np.minimum(np.abs(near["start"] - pos), np.abs(near["end"] - pos)))
    best = near["dist"].min()
    tied = sorted(near.loc[near["dist"] == best, "gene"])
    pref = [g for g in tied if g in author] or tied
    gene = pref[0]
    basis = "lead variant inside gene" if best == 0 else "nearest protein-coding gene"
    others = ";".join(sorted(set(near.loc[near["dist"] == 0, "gene"]) - {gene}))
    gid = near.loc[near["gene"] == gene, "gene_id"].iloc[0]
    return gene, basis, int(best), others, gid


def build_disease(disease: str, meta: dict, genes: pd.DataFrame) -> pd.DataFrame:
    print(f"{disease}: associations for {meta['accession']}", flush=True)
    df = pd.DataFrame([r for r in associations(meta["accession"]) if r["lead_snp"]])
    pos = positions(df["lead_snp"].tolist())
    df["chr"] = df["lead_snp"].map(lambda s: pos.get(s, (None, None))[0])
    df["pos_grch38"] = df["lead_snp"].map(lambda s: pos.get(s, (None, None))[1])
    df["genome_wide_significant"] = df["pvalue"] < GWS
    ok = df["genome_wide_significant"] & df["chr"].notna()
    print(f"  {ok.sum()} genome-wide significant variants with a GRCh38 position; clumping", flush=True)
    gws = clump(df[ok].copy())
    other = df[~ok].copy()
    other["locus_id"], other["is_locus_lead"], other["merged_by"] = None, False, ""
    df = pd.concat([gws, other], ignore_index=True)

    picks = []
    for _, r in df.iterrows():
        if pd.isna(r["chr"]):
            picks.append((None, "no GRCh38 position", None, "", ""))
        else:
            picks.append(pick_gene(r["chr"], int(r["pos_grch38"]), set(r["author_reported_genes"].split(";")), genes))
    df["mapped_gene"], df["mapping_basis"], df["distance_bp"], df["other_genes_containing_variant"], \
        df["mapped_gene_ensembl"] = zip(*picks)
    gtype = genes.drop_duplicates("gene").set_index("gene")["gene_type"]
    df["author_gene_types"] = df["author_reported_genes"].map(
        lambda s: ";".join(f"{g}({gtype.get(g, 'not in GENCODE')})" for g in s.split(";") if g))

    leads = df[df["is_locus_lead"]].sort_values("pvalue")
    df["locus_rank"] = df["locus_id"].map({lid: k + 1 for k, lid in enumerate(leads["locus_id"])})
    df["disease"], df["study_accession"], df["pubmed_id"], df["citation"] = \
        disease, meta["accession"], meta["pubmed_id"], meta["citation"]
    return df


def assign_sets(df: pd.DataFrame) -> pd.DataFrame:
    df["in_primary_top10"] = False
    df["in_all_gws"] = False
    df["decision"] = ""
    for disease, g in df.groupby("disease"):
        seen, n = set(), 0
        for idx, row in g[g["is_locus_lead"]].sort_values("pvalue").iterrows():
            gene = row["mapped_gene"]
            if gene is None:
                df.loc[idx, "decision"] = "excluded: " + row["mapping_basis"]
                continue
            if gene in seen:
                df.loc[idx, "decision"] = f"excluded: {gene} already taken by a stronger locus"
                continue
            seen.add(gene)
            n += 1
            df.loc[idx, "in_all_gws"] = True
            df.loc[idx, "in_primary_top10"] = n <= TOP_N
            df.loc[idx, "decision"] = "included (top 10)" if n <= TOP_N else "included (all-GWS sensitivity set only)"
        for idx in g.index[~g["is_locus_lead"]]:
            if not df.loc[idx, "genome_wide_significant"]:
                df.loc[idx, "decision"] = "excluded: p >= 5e-8"
            elif pd.isna(df.loc[idx, "chr"]):
                df.loc[idx, "decision"] = "excluded: no GRCh38 position"
            else:
                df.loc[idx, "decision"] = f"excluded: same locus as a stronger variant ({df.loc[idx, 'merged_by']})"
    return df


def write_sets(df: pd.DataFrame, genes: pd.DataFrame):
    out = GENE_DIR / "sets"
    out.mkdir(parents=True, exist_ok=True)
    cols = ["gene", "locus_rank", "locus_id", "lead_snp", "pvalue", "mapping_basis", "source"]
    pcg = set(genes.loc[(genes["gene_type"] == "protein_coding") & ~genes["readthrough"], "gene"])
    for disease, g in df.groupby("disease"):
        cite = g["citation"].iloc[0]
        for suffix, mask in [("", g["in_primary_top10"]), ("_all_gws", g["in_all_gws"])]:
            s = g[mask].sort_values("locus_rank").rename(columns={"mapped_gene": "gene"})
            s["source"] = cite
            s[cols].to_csv(out / f"{disease}{suffix}.csv", index=False)
        auth = []
        for _, r in g[g["in_primary_top10"]].sort_values("locus_rank").iterrows():
            cands = [x for x in r["author_reported_genes"].split(";") if x in pcg]
            if cands:
                pick = r["mapped_gene"] if r["mapped_gene"] in cands else cands[0]
                auth.append({"gene": pick, "locus_rank": r["locus_rank"], "locus_id": r["locus_id"],
                             "lead_snp": r["lead_snp"], "pvalue": r["pvalue"],
                             "mapping_basis": "author-reported gene", "source": cite})
        if len(auth) >= 5:
            pd.DataFrame(auth).drop_duplicates("gene")[cols].to_csv(out / f"{disease}_author.csv", index=False)
    for name, meta in HD_SETS.items():
        pd.DataFrame({"gene": meta["genes"], "locus_rank": range(1, len(meta["genes"]) + 1), "locus_id": "",
                      "lead_snp": "", "pvalue": None, "mapping_basis": "named by source paper",
                      "source": meta["source"]})[cols].to_csv(out / f"{name}.csv", index=False)


def main():
    GENE_DIR.mkdir(parents=True, exist_ok=True)
    genes = load_gencode()
    try:
        stats = request_json("GET", "https://www.ebi.ac.uk/gwas/api/search/stats", tries=3)
    except ApiError:
        stats = {}
    df = pd.concat([build_disease(d, m, genes) for d, m in STUDIES.items()], ignore_index=True)
    df = assign_sets(df)
    order = ["disease", "study_accession", "pubmed_id", "association_id", "lead_snp", "pvalue",
             "genome_wide_significant", "chr", "pos_grch38", "locus_id", "is_locus_lead", "locus_rank",
             "merged_by", "author_reported_genes", "author_gene_types", "mapped_gene", "mapped_gene_ensembl",
             "mapping_basis", "distance_bp", "other_genes_containing_variant", "in_primary_top10",
             "in_all_gws", "decision", "citation"]
    df = df[order].sort_values(["disease", "pvalue"])
    df.to_csv(GENE_DIR / "gene_audit.csv", index=False)
    write_sets(df, genes)
    prov = {
        "accessed": dt.date.today().isoformat(),
        "gwas_catalog_release": stats,
        "gene_annotation": "GENCODE v50 basic (release date 2026-06-09), protein_coding, readthrough genes excluded",
        "variant_positions": "Ensembl REST /variation, GRCh38",
        "ld": f"Ensembl REST /ld pairwise, {LD_POP}",
        "rules": {"gws_p": GWS, "window_bp": WINDOW_BP, "ld_r2": LD_R2, "ld_fallback_bp": FALLBACK_BP, "top_n": TOP_N},
        "studies": STUDIES,
        "hand_entered_sets": HD_SETS,
    }
    (GENE_DIR / "provenance.json").write_text(json.dumps(prov, indent=2))
    for disease, g in df.groupby("disease"):
        top = g[g["in_primary_top10"]].sort_values("locus_rank")
        print(f"\n{disease}: {g['is_locus_lead'].sum()} independent GWS loci; top 10 genes:")
        print(top[["locus_rank", "mapped_gene", "lead_snp", "pvalue", "mapping_basis", "distance_bp"]].to_string(index=False))
        merges = g[g["merged_by"].astype(str).str.len() > 0]
        print(f"  {len(merges)} variants merged into stronger loci; "
              f"{merges['merged_by'].str.contains('fallback').sum()} by distance fallback")


if __name__ == "__main__":
    main()

"""
Pull top disease-associated genes from the GWAS Catalog REST API.

Selection criteria (fixed here, before any downstream expression analysis —
documented so the choice can't be quietly re-tuned after seeing results):

  1. One landmark GWAS study per disease, chosen for sample size / recency and
     identified by its PubMed ID, independent of anything found later in this
     project:
       - Parkinson's:  Nalls et al. 2019, Lancet Neurology, PMID 31701892
                        (GWAS Catalog accession GCST009325)
       - Alzheimer's:  Bellenguez et al. 2022, Nat Genet, PMID 35379992
                        (GCST90027158)
       - Huntington's: Moss et al. 2017, Lancet Neurology, PMID 28642124 —
                        a GWAS of disease-progression modifiers *within*
                        confirmed HTT mutation carriers (GCST004691). HD is
                        monogenic, so there is no case/control risk GWAS in
                        the PD/AD sense; this is the closest equivalent
                        (genetic modifiers of clinical severity/onset), and
                        HTT itself is added manually as the causal gene since
                        no risk-GWAS pipeline would otherwise surface it.
  2. Associations ranked by ascending reported p-value.
  3. Each lead SNP mapped to a gene via the GWAS Catalog's own genomic-context
     annotation: prefer a gene whose body the SNP falls inside
     (isIntergenic == False); otherwise take the nearest flanking gene.
  4. Top TOP_N unique genes kept per disease, in p-value rank order.

Run:  .venv/Scripts/python src/fetch_gwas_genes.py
Output: data/gene_lists/<disease>_risk_genes.csv, data/gene_lists/all_genes.csv
"""
from __future__ import annotations

import csv
import re
import time
from pathlib import Path

import requests

# Non-protein-coding / placeholder annotations to skip when mapping a SNP to
# a gene (uncharacterized loci, pseudogenes, antisense/lncRNA, micro-RNAs) —
# we want the nearest *named* gene a wet-lab reader would recognize.
_SKIP_GENE_PATTERN = re.compile(
    r"^(LOC\d+|LINC\d+|MIR\d+.*|RNU\d+.*|RNA5S.*|SNOR.*|.*-AS\d*)$"
)

API = "https://www.ebi.ac.uk/gwas/rest/api"
TOP_N = 10
OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "gene_lists"

STUDIES = {
    "parkinsons": {
        "accession": "GCST009325",
        "pubmed_id": "31701892",
        "citation": "Nalls et al. 2019, Lancet Neurology",
        "manual_add": [],
    },
    "alzheimers": {
        "accession": "GCST90027158",
        "pubmed_id": "35379992",
        "citation": "Bellenguez et al. 2022, Nature Genetics",
        "manual_add": [],
    },
    "huntingtons": {
        "accession": "GCST004691",
        "pubmed_id": "28642124",
        "citation": "Moss et al. 2017, Lancet Neurology (progression-modifier GWAS in HTT carriers)",
        # HTT is the monogenic cause of HD and would not appear in a modifier
        # GWAS; included manually and flagged as such in the output.
        "manual_add": ["HTT"],
    },
}


def _get_with_retry(url: str, timeout: int, max_retries: int = 5) -> requests.Response:
    delay = 3
    for attempt in range(max_retries):
        resp = requests.get(url, timeout=timeout)
        if resp.status_code != 429:
            return resp
        time.sleep(delay)
        delay *= 2
    resp.raise_for_status()
    return resp


def fetch_associations(accession: str) -> list[dict]:
    resp = _get_with_retry(f"{API}/studies/{accession}/associations", timeout=60)
    resp.raise_for_status()
    return resp.json()["_embedded"]["associations"]


def lead_snp_rsid(assoc: dict) -> str | None:
    for locus in assoc.get("loci", []):
        for allele in locus.get("strongestRiskAlleles", []):
            name = allele.get("riskAlleleName", "")
            rsid = name.split("-")[0]
            if rsid.startswith("rs"):
                return rsid
    return None


def mapped_gene_for_snp(rsid: str, cache: dict) -> str | None:
    if rsid in cache:
        return cache[rsid]
    resp = _get_with_retry(f"{API}/singleNucleotidePolymorphisms/{rsid}", timeout=30)
    if resp.status_code != 200:
        cache[rsid] = None
        return None
    contexts = resp.json().get("genomicContexts", [])
    gene = None
    best_distance = None
    for ctx in contexts:
        gene_name = ctx.get("gene", {}).get("geneName")
        if not gene_name or _SKIP_GENE_PATTERN.match(gene_name):
            continue  # skip uncharacterized/pseudogene/antisense annotations
        if not ctx.get("isIntergenic", True):
            gene = gene_name
            break
        distance = ctx.get("distance")
        if distance is not None and (best_distance is None or distance < best_distance):
            best_distance, gene = distance, gene_name
    cache[rsid] = gene
    return gene


def top_genes_for_study(accession: str, snp_cache: dict) -> list[tuple[str, float, str]]:
    associations = fetch_associations(accession)
    associations.sort(key=lambda a: a.get("pvalue") if a.get("pvalue") is not None else 1.0)

    ranked: list[tuple[str, float, str]] = []
    seen_genes: set[str] = set()
    for assoc in associations:
        if len(seen_genes) >= TOP_N:
            break
        rsid = lead_snp_rsid(assoc)
        if not rsid:
            continue
        gene = mapped_gene_for_snp(rsid, snp_cache)
        time.sleep(0.3)  # be polite to the free public API
        if not gene or gene in seen_genes:
            continue
        seen_genes.add(gene)
        ranked.append((gene, assoc.get("pvalue"), rsid))
    return ranked


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    snp_cache: dict = {}
    all_rows = []

    for disease, meta in STUDIES.items():
        print(f"Fetching {disease} ({meta['accession']}, PMID {meta['pubmed_id']})...")
        ranked = top_genes_for_study(meta["accession"], snp_cache)
        rows = [
            {
                "disease": disease,
                "gene": gene,
                "pvalue": pvalue,
                "lead_snp": rsid,
                "source": "gwas_catalog",
                "study_accession": meta["accession"],
                "pubmed_id": meta["pubmed_id"],
                "citation": meta["citation"],
            }
            for gene, pvalue, rsid in ranked
        ]
        for gene in meta["manual_add"]:
            rows.append(
                {
                    "disease": disease,
                    "gene": gene,
                    "pvalue": None,
                    "lead_snp": None,
                    "source": "manual_causal_gene",
                    "study_accession": meta["accession"],
                    "pubmed_id": meta["pubmed_id"],
                    "citation": meta["citation"],
                }
            )

        out_path = OUT_DIR / f"{disease}_risk_genes.csv"
        with out_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        print(f"  -> {out_path} ({len(rows)} genes: {[r['gene'] for r in rows]})")
        all_rows.extend(rows)

    all_path = OUT_DIR / "all_genes.csv"
    with all_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_rows[0].keys())
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nWrote combined list -> {all_path} ({len(all_rows)} rows total)")


if __name__ == "__main__":
    main()

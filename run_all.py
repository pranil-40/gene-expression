"""
Reruns the whole v2 pipeline in order.

    python run_all.py            # everything
    python run_all.py --from analysis

Needs internet for the first run (AHBA ~4 GB via abagen, CIT168 atlas,
GENCODE v50, GTEx v10 medians, GWAS Catalog/Ensembl/NCBI lookups). After
that, everything reads from data/raw and data/cache.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PY = sys.executable

STEPS = [
    ("atlases", ["src/build_atlases.py"]),
    ("extract", ["src/expression.py", "extract"]),
    ("build", ["src/expression.py", "build"]),
    ("check", ["src/expression.py", "check"]),
    ("anatomy", ["src/sn_anatomy_audit.py"]),
    ("genes", ["src/fetch_gwas_genes.py"]),
    ("probes", ["src/gene_probe_audit.py"]),
    ("snca", ["src/snca_checkpoint.py"]),
    ("analysis", ["src/analysis.py"]),
    ("power", ["src/power_check.py"]),       # post hoc, docs/analysis_plan_amendment_1.md
    ("gtex", ["src/gtex_check.py"]),
    ("figures", ["src/figures.py"]),
]


def main():
    names = [s[0] for s in STEPS]
    start = 0
    if "--from" in sys.argv:
        start = names.index(sys.argv[sys.argv.index("--from") + 1])
    for name, cmd in STEPS[start:]:
        if name == "extract" and (ROOT / "data" / "cache" / "genes.csv").exists():
            print("== extract: cache already there, skipping (delete data/cache to redo)")
            continue
        print(f"== {name}: {' '.join(cmd)}", flush=True)
        subprocess.run([PY, *cmd], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()

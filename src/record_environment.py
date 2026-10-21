"""Record the software environment, data hashes and RNG seeds used for the analysis.

Writes ENVIRONMENT.txt and requirements.lock.txt in the project root.
"""
import hashlib
import importlib.metadata as metadata
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ["numpy", "pandas", "scipy", "scikit-learn", "matplotlib", "seaborn", "biopython", "ete3",
            "networkx", "pyyaml", "requests", "tqdm", "pyfaidx"]
HASHED_FILES = [
    "raw_data/zoonomia/241-mammalian-2020v2.phast-242.nh",
    "raw_data/zoonomia/zoonomia_mammal_assemblies.tsv",
    "raw_data/toga_inventory/overview.table.tsv",
    "raw_data/toga_gene_status/download_manifest.tsv",
    "analysis/gene_status_matrix.tsv",
    "results/agentic_triage/ranked_co_loss_candidates.tsv",
    "results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv",
    "results/phylo_null/phylo_null_pairs.tsv",
    "results/phylogenetic_validation_unlinked/pagel_dependent_vs_independent.tsv",
    "results/matched_pair_null/matched_pair_null.tsv",
    "results/triage/candidate_triage.tsv",
]


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()


versions = {}
for name in PACKAGES:
    try:
        versions[name] = metadata.version(name)
    except metadata.PackageNotFoundError:
        pass

(ROOT / "requirements.lock.txt").write_text("".join(f"{k}=={v}\n" for k, v in versions.items()))

lines = [
    f"recorded_utc: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
    f"python: {platform.python_version()}",
    f"platform: {platform.platform()}",
    f"git_commit: {git('rev-parse', 'HEAD')}",
    f"git_uncommitted_files: {len(git('status', '--porcelain').splitlines())}",
    "",
    "packages (see requirements.lock.txt):",
    *[f"  {k}=={v}" for k, v in versions.items()],
    "",
    "sha256 of key inputs and outputs:",
]
for rel in HASHED_FILES:
    path = ROOT / rel
    digest = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else "MISSING"
    lines.append(f"  {digest}  {rel}")
lines += [
    "",
    "random seeds:",
    "  src/phylo_null_coloss.py: numpy default_rng(20261014)",
    "  src/phylogenetic_candidate_test.py: random starts default_rng(20261021 + parameter_count)",
    "",
    "run settings:",
    "  src/phylo_null_coloss.py 100000 (BURDEN_AWARE=1 by default)",
    "  src/phylogenetic_candidate_test.py --starts 4 --maxiter 200 (unlinked pairs, one job per pair)",
    "",
    "data source: TOGA human-reference loss_summ_data for 162 Zoonomia assemblies (exact accession match);",
    "  genome, alignment and ortholog files were never downloaded or needed.",
]
(ROOT / "ENVIRONMENT.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))

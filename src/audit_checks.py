"""Integrity and bias checks on the gene-status matrix and the raw screen (writes results/audit/).

1. Rebuild the matrix from the raw TOGA files and compare it with analysis/gene_status_matrix.tsv.
2. Verify SHA-256 of every raw status file, the tree, the roster and the TOGA catalog.
3. Correlate per-species status counts with assembly contiguity (contig/scaffold N50).
4. Rank the left-right module pairs inside the full raw screen (same scoring as agentic_candidate_triage.py).
"""
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "audit"
OUT.mkdir(parents=True, exist_ok=True)
LR_GENES = ["DAND5", "CIROP", "CIROZ", "PKD1L1", "MMP21"]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


matrix = pd.read_csv(ROOT / "analysis/gene_status_matrix.tsv", sep="\t", index_col=0)
manifest = pd.read_csv(ROOT / "raw_data/toga_gene_status/download_manifest.tsv", sep="\t")

# 1-2. Hashes and independent rebuild
mismatched, columns = [], {}
for record in manifest.itertuples():
    path = ROOT / record.local_file
    if sha256(path) != record.sha256:
        mismatched.append(record.zoonomia_tree_tip)
    genes = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            entity, identifier, status = line.rstrip("\n").split("\t")
            if entity == "GENE":
                genes[identifier] = status
    columns[record.zoonomia_tree_tip] = genes
rebuilt = pd.DataFrame(columns).sort_index()[matrix.columns]

assembly_manifest = json.loads((ROOT / "orthologs/zoonomia_toga_assembly_manifest.summary.json").read_text())
integrity = {
    "raw_status_files": len(manifest),
    "raw_status_sha256_mismatches": mismatched,
    "rebuilt_matrix_identical": bool(rebuilt.equals(matrix)),
    "matrix_shape": list(matrix.shape),
    "matrix_sha256": sha256(ROOT / "analysis/gene_status_matrix.tsv"),
    "tree_sha256_matches_manifest": sha256(ROOT / "raw_data/zoonomia/241-mammalian-2020v2.phast-242.nh") == assembly_manifest["tree_sha256"],
    "roster_sha256_matches_manifest": sha256(ROOT / "raw_data/zoonomia/zoonomia_mammal_assemblies.tsv") == assembly_manifest["zoonomia_roster_sha256"],
    "toga_catalog_sha256_matches_manifest": sha256(ROOT / "raw_data/toga_inventory/overview.table.tsv") == assembly_manifest["toga_catalog_sha256"],
    "status_counts": {k: int(v) for k, v in matrix.stack().value_counts().items()},
}
(OUT / "integrity.json").write_text(json.dumps(integrity, indent=2) + "\n")
print(json.dumps(integrity, indent=2))

# 3. Assembly quality vs status counts
eligible = pd.read_csv(ROOT / "orthologs/zoonomia_toga_eligible_assemblies.tsv", sep="\t").set_index("zoonomia_tree_tip")
burden = pd.DataFrame({status: (matrix == status).sum() for status in ["I", "L", "UL", "M", "PM"]})
burden["contig_n50"] = pd.to_numeric(eligible.loc[burden.index, "selected_contig_n50_bp"], errors="coerce")
burden["scaffold_n50"] = pd.to_numeric(eligible.loc[burden.index, "selected_scaffold_n50_bp"], errors="coerce")
burden.to_csv(OUT / "species_status_counts_with_n50.tsv", sep="\t")
rows = []
for status, metric in itertools.product(["I", "L", "UL", "M", "PM"], ["contig_n50", "scaffold_n50"]):
    result = stats.spearmanr(burden[status], burden[metric], nan_policy="omit")
    rows.append(dict(status=status, assembly_metric=metric, spearman_rho=result.statistic, p_value=result.pvalue))
pd.DataFrame(rows).to_csv(OUT / "status_vs_assembly_quality.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).round(3).to_string(index=False))

# 4. Rank of LR module pairs in the raw screen (blocked matrix products, 7,204 eligible genes)
lost = (matrix == "L").values.astype(np.float32)
comparable = matrix.isin(["L", "I"]).values.astype(np.float32)
eligible_genes = lost.sum(axis=1) >= 2
ids = matrix.index[eligible_genes]
lost, comparable = lost[eligible_genes], comparable[eligible_genes]
symbols = pd.read_csv(ROOT / "raw_data/ensembl/ensembl_id_symbol_map.tsv", sep="\t").drop_duplicates("Gene name")
symbol_to_id = dict(zip(symbols["Gene name"], symbols["Gene stable ID"]))
position = {g: i for i, g in enumerate(ids)}


def pair_excess(a, b):
    ia, ib = position[symbol_to_id[a]], position[symbol_to_id[b]]
    both = comparable[ia] * comparable[ib]
    shared = (lost[ia] * lost[ib] * both).sum()
    return shared - (lost[ia] * both).sum() * (lost[ib] * both).sum() / both.sum(), shared


targets = {f"{a}-{b}": pair_excess(a, b) for a, b in itertools.combinations(LR_GENES, 2)}
better = {k: 0 for k in targets}
total_pairs = 0
block = 1000
for start in range(0, len(ids), block):
    stop = min(start + block, len(ids))
    units = comparable[start:stop] @ comparable.T
    shared = lost[start:stop] @ lost.T
    with np.errstate(divide="ignore", invalid="ignore"):
        excess = shared - (lost[start:stop] @ comparable.T) * (comparable[start:stop] @ lost.T) / units
    valid = (units >= 20) & (shared >= 2) & (np.arange(len(ids))[None, :] > np.arange(start, stop)[:, None])
    excess = np.where(valid, excess, -np.inf)
    total_pairs += int(valid.sum())
    for key, (value, _) in targets.items():
        better[key] += int((excess > value + 1e-6).sum())
rank = pd.DataFrame(
    [dict(pair=k, shared_L=int(s), excess=round(float(e), 2), rank_in_screen=better[k] + 1, pairs_in_screen=total_pairs)
     for k, (e, s) in targets.items()]
)
rank.to_csv(OUT / "lr_module_rank_in_screen.tsv", sep="\t", index=False)
print(rank.to_string(index=False))

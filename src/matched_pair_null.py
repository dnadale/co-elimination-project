"""Empirical null for co-loss: compare each pair with all gene pairs that have matching loss counts.

For a pair (A, B), every gene whose total L count is within 10% (min 2) of A's, paired with every
gene matched to B, forms the comparison set. Shared L is counted within the pair's own comparable
species (both I/L calls). This keeps clade structure, species loss burden and assembly effects of the
real data without a model. It is a filter for gene-family blocks, not a genome-wide significance test:
the tested pairs were selected for high excess, so they sit at the top of their own strata.

n_equal_or_more counts the pair itself, so 1 means no other matched pair shares as many losses.
"""
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "matched_pair_null"
OUT.mkdir(parents=True, exist_ok=True)
LR_GENES = ["DAND5", "CIROP", "CIROZ", "PKD1L1", "MMP21"]

df = pd.read_csv(ROOT / "analysis/gene_status_matrix.tsv", sep="\t", index_col=0)
L = (df == "L").values.astype(np.float32)
C = df.isin(["L", "I"]).values.astype(np.float32)
n_lost = L.sum(axis=1)
pos = {g: i for i, g in enumerate(df.index)}

symbols = pd.read_csv(ROOT / "raw_data/ensembl/ensembl_id_symbol_map.tsv", sep="\t").drop_duplicates("Gene name")
symbol_to_id = dict(zip(symbols["Gene name"], symbols["Gene stable ID"]))

unlinked = pd.read_csv(
    ROOT / "results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv", sep="\t", keep_default_na=False
)
tests = [
    ("candidate", r.gene_a_symbol or r.gene_a, r.gene_b_symbol or r.gene_b, r.gene_a, r.gene_b)
    for r in unlinked.itertuples()
]
tests += [("LR_control", a, b, symbol_to_id[a], symbol_to_id[b]) for a, b in itertools.combinations(LR_GENES, 2)]


def matched(index):
    target = n_lost[index]
    return np.where(np.abs(n_lost - target) <= max(2, 0.1 * target))[0]


rows = []
for kind, label_a, label_b, gene_a, gene_b in tests:
    ia, ib = pos[gene_a], pos[gene_b]
    comparable = C[ia] * C[ib]
    observed = float((L[ia] * L[ib] * comparable).sum())
    set_a, set_b = matched(ia), matched(ib)
    shared = (L[set_a] * comparable) @ (L[set_b] * comparable).T
    distinct = set_a[:, None] != set_b[None, :]
    n_pairs = int(distinct.sum())
    n_equal_or_more = int(((shared >= observed) & distinct).sum())
    rows.append(dict(
        kind=kind, gene_a=label_a, gene_b=label_b, observed_shared_L=int(observed),
        n_matched_a=len(set_a), n_matched_b=len(set_b), n_matched_pairs=n_pairs,
        n_equal_or_more=n_equal_or_more, fraction=n_equal_or_more / n_pairs,
    ))

result = pd.DataFrame(rows)
result.to_csv(OUT / "matched_pair_null.tsv", sep="\t", index=False)
print(result.to_string(index=False))

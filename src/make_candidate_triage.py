"""Combine the vetting, tree-null, Pagel and matched-pair results for the 19 unlinked pairs.

Tier rules (exploratory triage, not a significance claim; all pairs were selected from 7.4M scored pairs):
  block  : >= 10 other gene pairs with matched loss counts share at least as many losses, i.e. the pair
           belongs to a group of genes lost in the same species (gene-family or lineage effect).
  lead   : not a block, q < 0.05 in both the tree null and the Pagel test.
  weak   : not a block, fails at least one of those two tests.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "triage"
OUT.mkdir(parents=True, exist_ok=True)
BLOCK_MIN_OTHER_PAIRS = 10
Q_CUTOFF = 0.05
MANUAL_NOTES = {"FAM200A-FAM200B": "paralogous pair; ortholog assignment artifact not excluded"}

read = lambda path: pd.read_csv(ROOT / path, sep="\t", keep_default_na=False)
pairs = read("results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv")
vetting = read("results/vetting/unlinked_pairs_vetting.tsv")
tree_null = read("results/phylo_null/phylo_null_pairs.tsv").query("kind == 'candidate'").reset_index(drop=True)
pagel = read("results/phylogenetic_validation_unlinked/pagel_dependent_vs_independent.tsv")
matched = read("results/matched_pair_null/matched_pair_null.tsv").query("kind == 'candidate'").reset_index(drop=True)

assert len(pairs) == len(vetting) == len(tree_null) == len(pagel) == len(matched) == 19
assert (pairs.gene_a.values == pagel.gene_a.values).all() and (pairs.gene_b.values == pagel.gene_b.values).all()

label = lambda symbol, gene_id: symbol or gene_id
table = pd.DataFrame({
    "pair": [f"{label(a, ia)}-{label(b, ib)}" for a, b, ia, ib in zip(pairs.gene_a_symbol, pairs.gene_b_symbol, pairs.gene_a, pairs.gene_b)],
    "shared_L": pairs.shared_loss_count.astype(int),
    "n_clades_with_both_L": vetting.n_clades_with_both_L.astype(int),
    "tree_null_mean": tree_null.tree_null_mean.astype(float),
    "q_tree_null": tree_null.q_tree_BH.astype(float).round(4),
    "pagel_lrt": pagel.lrt_statistic.astype(float).round(1),
    "q_pagel": pagel.q_bh_within_selected_shortlist.astype(float).round(4),
    "other_matched_pairs_as_extreme": matched.n_equal_or_more.astype(int) - 1,
})
is_block = table.other_matched_pairs_as_extreme >= BLOCK_MIN_OTHER_PAIRS
passes = (table.q_tree_null < Q_CUTOFF) & (table.q_pagel < Q_CUTOFF)
table["tier"] = np.where(is_block, "block", np.where(passes, "lead", "weak"))
table["note"] = table.pair.map(MANUAL_NOTES).fillna("")
table = table.sort_values(["tier", "q_pagel"], key=lambda s: s.map({"lead": 0, "weak": 1, "block": 2}) if s.name == "tier" else s)
table.to_csv(OUT / "candidate_triage.tsv", sep="\t", index=False)
print(table.to_string(index=False))

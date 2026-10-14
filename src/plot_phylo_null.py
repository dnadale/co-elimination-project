"""Merge vetting table with phylogeny-aware null results; write a TSV and a chart."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phylo_null"

pairs = pd.read_csv(ROOT / "results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv", sep="\t")
vet = pd.read_csv(ROOT / "results/vetting/unlinked_pairs_vetting.tsv", sep="\t")
null = pd.read_csv(OUT / "phylo_null_pairs.tsv", sep="\t")
cand = null[null.kind == "candidate"].reset_index(drop=True)
assert len(cand) == len(vet) == len(pairs)

label = [f"{a if isinstance(a, str) else i[:15]}-{b if isinstance(b, str) else j[:15]}"
         for a, b, i, j in zip(pairs.gene_a_symbol, pairs.gene_b_symbol, pairs.gene_a, pairs.gene_b)]
vet.insert(0, "pair", label)
tab = vet.drop(columns=["clades_majority_lost"]).copy()
tab["naive_expected"] = cand["naive_expected"]
tab["tree_null_mean"] = cand["tree_null_mean"]
tab["p_tree"] = cand["p_tree"]
tab["q_tree_BH"] = cand["q_tree_BH"].round(3)
tab["survives_tree_correction"] = cand["q_tree_BH"] < 0.05
tab.to_csv(OUT / "unlinked_pairs_vetting_with_tree_null.tsv", sep="\t", index=False)
print(tab[["pair", "both_L", "tree_null_mean", "q_tree_BH", "survives_tree_correction"]].to_string(index=False))

d = tab.sort_values("q_tree_BH").reset_index(drop=True)
fig, axes = plt.subplots(1, 2, figsize=(14, 7.5), gridspec_kw={"width_ratios": [3, 1.3]})
ax = axes[0]
y = np.arange(len(d))[::-1]
ax.barh(y + 0.2, d.both_L, 0.4, color=np.where(d.survives_tree_correction, "#2a9d8f", "#9aa5b1"), label="Observed shared losses")
ax.barh(y - 0.2, d.tree_null_mean, 0.4, color="#e76f51", label="Expected by chance on the tree (simulated mean)")
ax.set_yticks(y)
ax.set_yticklabels(d.pair)
ax.set_xlabel("Species where both genes are lost")
ax.set_title("Observed vs tree-simulated shared losses\n(teal = survives correction, grey = does not)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2, fontsize=8, frameon=False)
ax2 = axes[1]
ax2.barh(y, -np.log10(d.q_tree_BH.clip(lower=1e-4)), color=np.where(d.survives_tree_correction, "#2a9d8f", "#9aa5b1"))
ax2.axvline(-np.log10(0.05), color="k", ls="--", lw=1)
ax2.text(-np.log10(0.05), len(d) - 0.3, " q = 0.05", fontsize=8, va="bottom")
ax2.set_yticks(y)
ax2.set_yticklabels([])
ax2.set_xlabel("-log10 q (tree-corrected)")
ax2.set_title("Significance")
fig.suptitle("19 unlinked co-loss pairs after phylogeny-aware correction", fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / "unlinked_pairs_tree_corrected.png", dpi=160)

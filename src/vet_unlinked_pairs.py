"""Vet the 19 unlinked co-loss pairs: artifact/clade check + STRING interaction check."""
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from plot_lr_module import GENUS2CLADE  # noqa: E402

OUT = ROOT / "results" / "vetting"
OUT.mkdir(parents=True, exist_ok=True)

pairs = pd.read_csv(ROOT / "results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv", sep="\t")
mat = pd.read_csv(ROOT / "analysis/gene_status_matrix.tsv", sep="\t", index_col="gene_id")
clade = pd.Series({s: GENUS2CLADE.get(s.split("_")[0], "Other") for s in mat.columns})


def string_score(a, b):
    if not (isinstance(a, str) and isinstance(b, str)):
        return None
    url = "https://string-db.org/api/json/network?" + urllib.parse.urlencode(
        {"identifiers": f"{a}\r{b}", "species": 9606, "caller_identity": "co_elim"})
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            hits = json.load(r)
    except Exception:
        return None
    return max((h["score"] for h in hits), default=0.0)


rows = []
for p in pairs.itertuples():
    a, b = mat.loc[p.gene_a], mat.loc[p.gene_b]
    both_L = (a == "L") & (b == "L")
    either_UL = ((a == "UL") | (b == "UL")) & ~both_L
    only_a = (a == "L") & (b != "L")
    only_b = (b == "L") & (a != "L")
    cl = clade[both_L[both_L].index].value_counts()
    # clades where the pair is lost in >=50% of that clade's species
    sizes = clade.value_counts()
    full = [c for c, n in cl.items() if n / sizes[c] >= 0.5]
    rows.append(dict(
        gene_a=p.gene_a_symbol, gene_b=p.gene_b_symbol,
        both_L=int(both_L.sum()), only_a_L=int(only_a.sum()), only_b_L=int(only_b.sum()),
        uncertain_UL_in_either=int(either_UL.sum()),
        n_clades_with_both_L=len(cl), n_clades_majority_lost=len(full),
        top_clade=cl.index[0] if len(cl) else "", top_clade_share=round(cl.iloc[0] / both_L.sum(), 2) if len(cl) else 0,
        clades_majority_lost=";".join(full),
        string_score=string_score(p.gene_a_symbol, p.gene_b_symbol),
    ))

df = pd.DataFrame(rows)


def verdict(r):
    if r.string_score and r.string_score >= 0.4:
        return "functional link (STRING)"
    if r.n_clades_majority_lost >= 3:
        return "multi-clade, no known link"
    if r.n_clades_majority_lost >= 2:
        return "2 clades, no known link"
    return "single clade (likely lineage-driven)"


df["verdict"] = df.apply(verdict, axis=1)
df.to_csv(OUT / "unlinked_pairs_vetting.tsv", sep="\t", index=False)
print(df.drop(columns=["clades_majority_lost"]).to_string())

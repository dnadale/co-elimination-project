"""Phylogeny-aware null for co-loss: simulate each gene independently on the real tree.

For each gene, fit a 2-state (intact/lost) continuous-time Markov model on the pruned Zoonomia tree
by maximum likelihood, then simulate N independent histories per gene. The ML rate alone does not
reproduce a gene's observed loss frequency, so simulations are kept only if the gene's simulated
loss count among the pair's comparable species is within 10% (min 2) of its observed count.
A pair's p-value is the share of retained simulation pairs in which two independent genes share at
least as many losses as observed. This keeps the tree structure (related species share state) that
the plain Fisher test ignores. Only I and L calls are used; all other statuses are ambiguous.
Branch lengths are also scaled by each lineage's mean relative loss burden (set BURDEN_AWARE=0 to disable).
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.optimize import minimize
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phylo_null"
OUT.mkdir(parents=True, exist_ok=True)
N_SIM = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
MIN_ACCEPTED = 200
rng = np.random.default_rng(20261014)

mat = pd.read_csv(ROOT / "analysis/gene_status_matrix.tsv", sep="\t", index_col="gene_id")
species = list(mat.columns)

tree = Phylo.read(str(ROOT / "raw_data/zoonomia/241-mammalian-2020v2.phast-242.nh"), "newick")
keep = set(species)
missing = keep - {t.name for t in tree.get_terminals()}
assert not missing, f"species not in tree: {sorted(missing)[:5]}"
for tip in [t for t in tree.get_terminals() if t.name not in keep]:
    tree.prune(tip)

nodes = list(tree.find_clades(order="preorder"))
idx = {id(n): i for i, n in enumerate(nodes)}
parent = np.full(len(nodes), -1)
blen = np.zeros(len(nodes))
for n in nodes:
    for c in n.clades:
        parent[idx[id(c)]] = idx[id(n)]
        blen[idx[id(c)]] = max(c.branch_length or 0.0, 1e-6)
tip_node = {n.name: idx[id(n)] for n in nodes if n.is_terminal()}
tip_idx = np.array([tip_node[s] for s in species])
children = [[idx[id(c)] for c in n.clades] for n in nodes]
postorder = list(range(len(nodes)))[::-1]  # preorder reversed
scale = np.median(blen[blen > 0])
blen = blen / scale

# Lineage rate heterogeneity: genes share a per-branch rate factor equal to the mean relative loss
# burden of the descendant species, so high-burden clades (e.g. cetaceans) lose many genes together
# even when the genes are independent.
burden = (mat == "L").sum(axis=0).values.astype(float)
tip_factor = burden / np.median(burden)
tip_sum = np.zeros(len(nodes))
tip_cnt = np.zeros(len(nodes))
for s_i, n in enumerate(tip_idx):
    tip_sum[n], tip_cnt[n] = tip_factor[s_i], 1
for n in postorder:
    if children[n]:
        tip_sum[n] = sum(tip_sum[c] for c in children[n])
        tip_cnt[n] = sum(tip_cnt[c] for c in children[n])
if os.environ.get("BURDEN_AWARE", "1") == "1":
    blen = blen * (tip_sum / tip_cnt)
print(f"tree: {len(species)} tips, {len(nodes)} nodes", flush=True)


def trans(a, b, t):
    # a: intact->lost rate, b: lost->intact rate
    return expm(np.array([[-a, a], [b, -b]]) * t)


def neg_loglik(theta, obs):
    a, b = np.exp(theta)
    L = np.ones((len(nodes), 2))
    for s_i, n in enumerate(tip_idx):
        v = obs[s_i]
        if v == 0:
            L[n] = [1, 0]
        elif v == 1:
            L[n] = [0, 1]
    cache = {}
    for n in postorder:
        if not children[n]:
            continue
        L[n] = 1.0
        for c in children[n]:
            key = round(blen[c], 8)
            if key not in cache:
                cache[key] = trans(a, b, blen[c])
            L[n] *= cache[key] @ L[c]
    # root fixed to intact: the reference genes are present in the mammalian ancestor
    return -np.log(max(L[0][0], 1e-300))


def fit_and_simulate(obs):
    best = min((minimize(neg_loglik, x0, args=(obs,), method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-3})
                for x0 in ([0.0, 0.0], [-1.0, 1.0])), key=lambda r: r.fun)
    a, b = np.exp(best.x)
    P = {}
    state = np.zeros((len(nodes), N_SIM), dtype=bool)
    for n in range(1, len(nodes)):
        key = round(blen[n], 8)
        if key not in P:
            P[key] = trans(a, b, blen[n])[:, 1]
        p_lost = np.where(state[parent[n]], P[key][1], P[key][0])
        state[n] = rng.random(N_SIM) < p_lost
    return state[tip_idx], a, b


def encode(row):
    out = np.full(len(species), -1)
    out[(row == "I").values] = 0
    out[(row == "L").values] = 1
    return out


def sym2id():
    s = pd.read_csv(ROOT / "raw_data/ensembl/ensembl_id_symbol_map.tsv", sep="\t").drop_duplicates("Gene stable ID")
    return {n: g for g, n in zip(s["Gene stable ID"], s["Gene name"])}


base2full = {g.split(".")[0]: g for g in mat.index}
name2base = sym2id()


def gid(symbol_or_id):
    b = name2base.get(symbol_or_id, symbol_or_id).split(".")[0]
    return base2full[b]


pairs = pd.read_csv(ROOT / "results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv", sep="\t", keep_default_na=False)
tests = [("candidate", r.gene_a, r.gene_b, r.gene_a_symbol, r.gene_b_symbol) for r in pairs.itertuples()]
lr = ["DAND5", "CIROP", "CIROZ", "PKD1L1", "MMP21"]
for i in range(len(lr)):
    for j in range(i + 1, len(lr)):
        tests.append(("LR_control", gid(lr[i]), gid(lr[j]), lr[i], lr[j]))

cache = {}


def get(g):
    if g not in cache:
        obs = encode(mat.loc[g])
        cache[g] = (obs,) + fit_and_simulate(obs)
        print(f"fit {g}: lost tips={int((obs == 1).sum())} a={cache[g][2]:.3g} b={cache[g][3]:.3g}", flush=True)
    return cache[g]


rows = []
for kind, ga, gb, sa, sb in tests:
    oa, sima, a_a, b_a = get(ga)
    ob, simb, a_b, b_b = get(gb)
    comp = (oa >= 0) & (ob >= 0)
    obs_shared = int(((oa == 1) & (ob == 1) & comp).sum())

    def accepted(obs, sim):
        target = int(((obs == 1) & comp).sum())
        counts = (sim & comp[:, None]).sum(axis=0)
        return np.where(np.abs(counts - target) <= max(2, 0.1 * target))[0]

    ia, ib = accepted(oa, sima), accepted(ob, simb)
    k = min(len(ia), len(ib))
    # naive (tree-blind) expectation from observed frequencies, for comparison
    naive = (oa[comp] == 1).mean() * (ob[comp] == 1).mean() * comp.sum()
    if k < MIN_ACCEPTED:
        null_mean, null_p95, null_max, p_tree = float("nan"), -1, -1, float("nan")
    else:
        rng.shuffle(ia)
        rng.shuffle(ib)
        null = (sima[:, ia[:k]] & simb[:, ib[:k]] & comp[:, None]).sum(axis=0)
        null_mean, null_p95, null_max = float(null.mean()), int(np.percentile(null, 95)), int(null.max())
        p_tree = (1 + int((null >= obs_shared).sum())) / (k + 1)
    rows.append(dict(
        kind=kind, gene_a=sa, gene_b=sb, comparable_species=int(comp.sum()), observed_shared_L=obs_shared,
        naive_expected=round(float(naive), 1), tree_null_mean=round(null_mean, 1),
        tree_null_p95=null_p95, tree_null_max=null_max, n_accepted_sims=k, p_tree=p_tree,
    ))
    print(rows[-1], flush=True)

df = pd.DataFrame(rows)
tested = df["p_tree"].notna()
m = int(tested.sum())
srt = df[tested].sort_values("p_tree", ascending=False)
bh = srt["p_tree"] * m / (m - np.arange(len(srt)))
df["q_tree_BH"] = bh.cummin().clip(upper=1).reindex(df.index)
df["survives_q05"] = df["q_tree_BH"] < 0.05
df.to_csv(OUT / "phylo_null_pairs.tsv", sep="\t", index=False)
print(df[["kind", "gene_a", "gene_b", "observed_shared_L", "naive_expected", "tree_null_mean", "p_tree", "q_tree_BH"]].to_string())

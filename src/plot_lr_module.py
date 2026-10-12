"""Figures and tables for the left-right (LR) module analysis.

Inputs : analysis/gene_status_matrix.tsv, raw_data/ensembl/ensembl_id_symbol_map.tsv,
         raw_data/zoonomia/241-mammalian-2020v2.phast-242.nh, results/agentic_triage/*.tsv
Outputs: results/lr_module/ (figures + tables)
"""
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "lr_module"
OUT.mkdir(parents=True, exist_ok=True)

MODULE = ["DAND5", "CIROP", "CIROZ", "PKD1L1", "MMP21"]

# Broad clade per genus (assigned manually from standard mammal taxonomy).
CLADES = {
    "Chiroptera": "Anoura Artibeus Carollia Craseonycteris Eptesicus Hipposideros Lasiurus Macroglossus "
                  "Megaderma Micronycteris Miniopterus Mormoops Murina Myotis Noctilio Pipistrellus "
                  "Pteronotus Pteropus Rhinolophus Tadarida Tonatia",
    "Cetacea": "Balaenoptera Eschrichtius Eubalaena Inia Kogia Lipotes Mesoplodon Neophocaena Phocoena "
               "Platanista Tursiops Ziphius",
    "Ruminantia/other Artiodactyla": "Ammotragus Beatragus Bison Elaphurus Giraffa Hemitragus Odocoileus "
                                     "Ovis Pantholops Rangifer Saiga Vicugna",
    "Carnivora": "Canis Enhydra Felis Helogale Leptonychotes Mellivora Mungos Odobenus Panthera "
                 "Paradoxurus Pteronura Spilogale Suricata Ursus Vulpes",
    "Perissodactyla": "Ceratotherium Dicerorhinus Equus Tapirus",
    "Primates": "Alouatta Ateles Callicebus Cebus Cercocebus Cercopithecus Chlorocebus Colobus "
                "Daubentonia Erythrocebus Eulemur Indri Lemur Macaca Mandrillus Mirza Nycticebus "
                "Otolemur Pithecia Propithecus Pygathrix Rhinopithecus Saguinus Saimiri Semnopithecus",
    "Rodentia": "Acomys Allactaga Aplodontia Cavia Chinchilla Cricetomys Ctenodactylus Ctenomys Dasyprocta "
                "Dinomys Dipodomys Dolichotis Ellobius Glis Graphiurus Heterocephalus Hydrochoerus Hystrix "
                "Ictidomys Jaculus Marmota Mesocricetus Microtus Mus Muscardinus Myocastor Nannospalax "
                "Octodon Ondatra Perognathus Petromus Rattus Sigmodon Spermophilus Thryonomys Xerus Zapus",
    "Lagomorpha": "Lepus Ochotona",
    "Eulipotyphla": "Condylura Erinaceus Scalopus Solenodon Sorex Uropsilus",
    "Afrotheria": "Echinops Heterohyrax Microgale Trichechus",
    "Xenarthra": "Choloepus Dasypus Myrmecophaga Tamandua Tolypeutes",
    "Pholidota": "Manis",
    "Scandentia": "Tupaia",
}
GENUS2CLADE = {g: c for c, gs in CLADES.items() for g in gs.split()}
CLADE_COLORS = dict(zip(CLADES, plt.cm.tab20(np.linspace(0, 1, len(CLADES)))))

STATUS_ORDER = ["I", "PI", "PG", "PM", "M", "UL", "L"]
STATUS_COLORS = ["#4c9f70", "#a8d5ba", "#f2e394", "#f4b183", "#9e9e9e", "#8e7cc3", "#c0392b"]
STATUS_LABELS = {"I": "Intact", "PI": "Partially intact", "PG": "Probably gene", "PM": "Probably missing",
                 "M": "Missing", "UL": "Uncertain loss", "L": "Lost"}


def load():
    mat = pd.read_csv(ROOT / "analysis/gene_status_matrix.tsv", sep="\t", index_col=0)
    m = pd.read_csv(ROOT / "raw_data/ensembl/ensembl_id_symbol_map.tsv", sep="\t", dtype=str).fillna("")
    m.columns = [c.strip().lower() for c in m.columns]
    m["gene stable id"] = m["gene stable id"].str.split(".").str[0]
    name = dict(zip(m["gene stable id"], m["gene name"]))
    sym2gid = {v: k for k, v in name.items() if v}
    rows = {}
    for s in MODULE:
        g = [x for x in mat.index if x.split(".")[0] == sym2gid[s]]
        rows[s] = mat.loc[g[0]]
    return mat, name, pd.DataFrame(rows)


def clade_of(sp):
    return GENUS2CLADE[sp.split("_")[0]]


def fig_heatmap(status, species_order, path, title):
    codes = status.loc[species_order].apply(lambda c: c.map({s: i for i, s in enumerate(STATUS_ORDER)}))
    fig, ax = plt.subplots(figsize=(7, 30))
    ax.imshow(codes.values, aspect="auto", cmap=ListedColormap(STATUS_COLORS), vmin=0, vmax=len(STATUS_ORDER) - 1)
    ax.set_xticks(range(len(MODULE)), MODULE)
    ax.set_yticks(range(len(species_order)), species_order, fontsize=6)
    for lab, sp in zip(ax.get_yticklabels(), species_order):
        lab.set_color(CLADE_COLORS[clade_of(sp)])
    ax.xaxis.tick_top()
    ax.legend(handles=[Patch(color=c, label=STATUS_LABELS[s]) for s, c in zip(STATUS_ORDER, STATUS_COLORS)],
              loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8)
    ax.set_title(title, pad=30)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def tree_layout(tree):
    ypos = {}
    for i, tip in enumerate(tree.get_terminals()):
        ypos[tip] = i
    for clade in tree.get_nonterminals(order="postorder"):
        ypos[clade] = np.mean([ypos[c] for c in clade])
    depths = tree.depths()
    return ypos, depths


def fig_tree(status, path):
    tree = Phylo.read(str(ROOT / "raw_data/zoonomia/241-mammalian-2020v2.phast-242.nh"), "newick")
    keep = set(status.index)
    for tip in [t for t in tree.get_terminals() if t.name not in keep]:
        tree.prune(tip)
    ypos, depths = tree_layout(tree)
    maxd = max(depths.values())
    fig, ax = plt.subplots(figsize=(10, 28))

    def draw(clade):
        for ch in clade.clades:
            ax.plot([depths[clade], depths[ch]], [ypos[ch]] * 2, color="#555", lw=0.6)
        if clade.clades:
            ys = [ypos[c] for c in clade.clades]
            ax.plot([depths[clade]] * 2, [min(ys), max(ys)], color="#555", lw=0.6)
        for ch in clade.clades:
            draw(ch)

    draw(tree.root)
    for tip in tree.get_terminals():
        y, sp = ypos[tip], tip.name
        lost = status.loc[sp, "CIROZ"] == "L"
        ax.text(maxd * 1.14 + 0.002, y, sp, va="center", fontsize=6, color=CLADE_COLORS[clade_of(sp)])
        for j, g in enumerate(MODULE):
            s = status.loc[sp, g]
            ax.add_patch(plt.Rectangle((maxd * 1.01 + j * maxd * 0.025, y - 0.45), maxd * 0.024, 0.9,
                                       color=STATUS_COLORS[STATUS_ORDER.index(s)]))
        if lost:
            ax.plot(maxd * 1.005, y, "o", color="#c0392b", ms=3)
    for j, g in enumerate(MODULE):
        ax.text(maxd * 1.01 + j * maxd * 0.025 + maxd * 0.012, -1.5, g, rotation=90, fontsize=7, ha="center")
    ax.set_ylim(len(tree.get_terminals()), -4)
    ax.axis("off")
    ax.legend(handles=[Patch(color=c, label=STATUS_LABELS[s]) for s, c in zip(STATUS_ORDER, STATUS_COLORS)],
              loc="lower left", fontsize=8)
    ax.set_title("LR module status on the Zoonomia tree (red dot = CIROZ lost)")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def pairwise(status):
    lost = status == "L"
    res = []
    for a, b in combinations(MODULE, 2):
        x, y = lost[a], lost[b]
        both, oa, ob, neither = int((x & y).sum()), int((x & ~y).sum()), int((~x & y).sum()), int((~x & ~y).sum())
        or_, p = fisher_exact([[both, oa], [ob, neither]], alternative="greater")
        res.append({"gene_a": a, "gene_b": b, "both_lost": both, "only_a": oa, "only_b": ob, "neither": neither,
                    "odds_ratio": or_, "fisher_p": p})
    df = pd.DataFrame(res).sort_values("fisher_p").reset_index(drop=True)
    n = len(df)
    ranked_p = df["fisher_p"].values
    bh = np.minimum.accumulate((ranked_p * n / np.arange(1, n + 1))[::-1])[::-1]
    df["q_bh"] = np.clip(bh, 0, 1)
    return df


def fig_pair_heatmap(df, path):
    m = pd.DataFrame(np.nan, index=MODULE, columns=MODULE)
    for _, r in df.iterrows():
        v = -np.log10(r["fisher_p"])
        m.loc[r["gene_a"], r["gene_b"]] = v
        m.loc[r["gene_b"], r["gene_a"]] = v
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(m.values, cmap="Reds")
    ax.set_xticks(range(5), MODULE)
    ax.set_yticks(range(5), MODULE)
    for i in range(5):
        for j in range(5):
            if i != j:
                ax.text(j, i, f"{m.values[i, j]:.1f}", ha="center", va="center", fontsize=9)
    fig.colorbar(im, label="-log10(Fisher p), one-sided")
    ax.set_title("Pairwise co-loss enrichment in the LR module")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_pair_bar(df, path):
    d = df.sort_values("both_lost")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.barh(d["gene_a"] + "–" + d["gene_b"], d["both_lost"], color="#c0392b")
    ax.set_xlabel("Species where both genes are lost")
    ax.set_title("Co-loss counts per pair (162 species)")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_forest(df, path):
    d = df.sort_values("odds_ratio").copy()
    finite_max = d.loc[np.isfinite(d["odds_ratio"]), "odds_ratio"].max()
    d["plot_or"] = d["odds_ratio"].replace(np.inf, finite_max * 3)
    fig, ax = plt.subplots(figsize=(7, 4))
    labels = d["gene_a"] + "–" + d["gene_b"]
    ax.scatter(d["plot_or"], labels, color=np.where(np.isinf(d["odds_ratio"]), "#8e7cc3", "#c0392b"))
    ax.set_xscale("log")
    ax.axvline(1, color="k", lw=0.8, ls="--")
    ax.set_xlabel("Odds ratio (log scale)")
    ax.set_title("Odds ratios of co-loss (purple = infinite OR, zero cell; plotted at 3x max finite)")
    ax.title.set_fontsize(8)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_gene_counts(status, path):
    counts = status.apply(lambda c: c.value_counts()).fillna(0).reindex(STATUS_ORDER).T
    counts.to_csv(OUT / "lr_module_status_counts_per_gene.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(8, 4))
    left = np.zeros(len(counts))
    for s, c in zip(STATUS_ORDER, STATUS_COLORS):
        ax.barh(counts.index, counts[s], left=left, color=c, label=STATUS_LABELS[s])
        left += counts[s].values
    ax.set_xlabel("Species (of 162)")
    ax.legend(fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.set_title("TOGA status per module gene")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_clade_bar(status, path):
    df = pd.DataFrame({"clade": [clade_of(s) for s in status.index], "lost": status["CIROZ"] == "L"})
    t = df.groupby("clade")["lost"].agg(lost="sum", total="count")
    t["fraction_lost"] = t["lost"] / t["total"]
    t = t.sort_values("fraction_lost")
    t.to_csv(OUT / "CIROZ_loss_by_clade.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(t.index, t["fraction_lost"], color=[CLADE_COLORS[c] for c in t.index])
    for i, (c, r) in enumerate(t.iterrows()):
        ax.text(r["fraction_lost"] + 0.01, i, f"{int(r['lost'])}/{int(r['total'])}", va="center", fontsize=8)
    ax.set_xlim(0, 1.15)
    ax.set_xlabel("Fraction of species with CIROZ lost")
    ax.set_title("CIROZ loss by clade")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_network(name, path):
    d = pd.read_csv(ROOT / "results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv", sep="\t")
    sym = lambda g: name.get(g.split(".")[0], "") or g.split(".")[0]
    g = nx.Graph()
    for _, r in d.iterrows():
        g.add_edge(sym(r["gene_a"]), sym(r["gene_b"]), w=r["shared_loss_count"])
    fig, ax = plt.subplots(figsize=(9, 7))
    pos = nx.spring_layout(g, seed=3, k=0.9)
    nx.draw_networkx(g, pos, ax=ax, node_size=[300 + 250 * g.degree(n) for n in g],
                     node_color=["#c0392b" if n == "CIROZ" else "#9ecae1" for n in g], font_size=8,
                     width=[g[u][v]["w"] / 15 for u, v in g.edges])
    ax.set_title("19 unlinked co-loss candidate pairs (edge width = shared losses; CIROZ red)")
    ax.axis("off")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_linked_vs_unlinked(path):
    base = ROOT / "results/agentic_triage"
    top = pd.read_csv(base / "ranked_co_loss_candidates.tsv", sep="\t")
    core = pd.read_csv(base / "ranked_co_loss_candidates_no_histone.tsv", sep="\t")
    unl = pd.read_csv(base / "ranked_co_loss_candidates_unlinked.tsv", sep="\t")
    parts = {"Histone–histone": len(top) - len(core), "Linked (<1 Mb)": len(core) - len(unl),
             "Unlinked": len(unl)}
    pd.Series(parts, name="pairs").to_csv(OUT / "candidate_pair_breakdown.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(parts.keys(), parts.values(), color=["#bbb", "#9ecae1", "#c0392b"])
    for i, v in enumerate(parts.values()):
        ax.text(i, v + 0.5, str(v), ha="center")
    ax.set_ylabel("Top-100 candidate pairs")
    ax.set_title("What the top 100 co-loss pairs are made of")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    mat, name, status = load()
    clades = pd.Series({s: clade_of(s) for s in status.index}, name="clade")
    order = sorted(status.index, key=lambda s: (clades[s], s))

    table = status.join(clades)
    table["genes_lost"] = (status == "L").sum(axis=1)
    table.index.name = "species"
    table.to_csv(OUT / "lr_module_species_status.tsv", sep="\t")
    table[table["genes_lost"] >= 4].sort_values("genes_lost", ascending=False).to_csv(
        OUT / "species_lost_4plus_of_5.tsv", sep="\t")
    status[status["PKD1L1"] == "UL"].join(clades).to_csv(OUT / "PKD1L1_uncertain_loss_species.tsv", sep="\t")

    pw = pairwise(status)
    pw.to_csv(OUT / "lr_module_pairwise_fisher_bh.tsv", sep="\t", index=False)

    fig_heatmap(status, order, OUT / "01_module_status_heatmap.png",
                "LR module status per species (grouped by clade)")
    fig_tree(status, OUT / "02_tree_with_module_status.png")
    fig_pair_heatmap(pw, OUT / "03_pairwise_coloss_heatmap.png")
    fig_pair_bar(pw, OUT / "04_pair_coloss_counts.png")
    fig_forest(pw, OUT / "05_odds_ratio_forest.png")
    fig_gene_counts(status, OUT / "06_status_counts_per_gene.png")
    fig_clade_bar(status, OUT / "07_CIROZ_loss_by_clade.png")
    fig_network(name, OUT / "08_unlinked_pairs_network.png")
    fig_linked_vs_unlinked(OUT / "09_candidate_pair_breakdown.png")
    print("done:", sorted(p.name for p in OUT.iterdir()))


if __name__ == "__main__":
    main()

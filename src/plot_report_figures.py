"""Presentation figures for each main result file. Output: results/figures/

1 matrix + QC overview        (analysis/gene_status_matrix.tsv + .qc.json)
2 top-100 candidate pairs     (results/agentic_triage/ranked_co_loss_candidates*.tsv)
3 19 unlinked pairs           (ranked_co_loss_candidates_unlinked.tsv)
4 LR module Fisher results    (results/lr_module_pairwise_fisher.tsv)
5 CIROZ loss species          (results/CIROZ_lost_species.txt)
Also bundles all figures into results/figures/all_figures.pdf
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from plot_lr_module import GENUS2CLADE, STATUS_COLORS, STATUS_LABELS, STATUS_ORDER, clade_of  # noqa: E402

OUT = ROOT / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 12, "axes.titleweight": "bold"})
RED, BLUE, GREY = "#c0392b", "#2c7fb8", "#9a9a9a"
FIGS = []


def save(fig, name):
    fig.savefig(OUT / name, dpi=180, bbox_inches="tight")
    FIGS.append(fig)


def sym(df, col):
    return df[f"{col}_symbol"].where(df[f"{col}_symbol"].fillna("") != "", df[col].str.split(".").str[0])


def fig_matrix_qc():
    q = json.load(open(ROOT / "analysis/gene_status_matrix.qc.json"))
    tot = pd.Series(q["combined_status_counts"]).reindex(STATUS_ORDER)
    sp = pd.DataFrame([{"species": r["species"], **r["status_counts"]} for r in q["species_qc"]]).set_index("species")
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), gridspec_kw={"width_ratios": [1.1, 1.2, 1.1]})
    ax = axs[0]
    ax.barh([STATUS_LABELS[s] for s in tot.index][::-1], tot.values[::-1], color=STATUS_COLORS[::-1])
    for i, v in enumerate(tot.values[::-1]):
        ax.text(v, i, f" {v:,} ({v / tot.sum():.1%})", va="center", fontsize=8)
    ax.set_xlim(0, tot.max() * 1.45)
    ax.set_xlabel("Gene-species calls")
    ax.set_title("Status of all 3,153,168 calls")
    ax = axs[1]
    ax.hist(sp["L"], bins=30, color=RED, edgecolor="white")
    ax.axvline(sp["L"].median(), color="k", ls="--", lw=1)
    ax.text(sp["L"].median(), ax.get_ylim()[1] * 0.95, f"  median {sp['L'].median():.0f}", va="top")
    ax.set_xlabel("Genes called Lost per species")
    ax.set_ylabel("Species")
    ax.set_title("Loss burden per species")
    ax = axs[2]
    ax.axis("off")
    rows = [("Species", f"{q['species_count']}"), ("Genes", f"{q['gene_count']:,}"),
            ("Matrix cells", f"{q['matrix_cell_count']:,}"), ("Missing cells", f"{q['missing_cells']}"),
            ("Per-species source files", f"{len(q['species_qc'])} (SHA256 recorded)"),
            ("Ancestry correction", "not yet applied")]
    for i, (k, v) in enumerate(rows):
        ax.text(0, 0.92 - i * 0.14, k, fontweight="bold")
        ax.text(0.55, 0.92 - i * 0.14, v)
    ax.set_title("QC checks")
    fig.suptitle("Gene status matrix: overview and QC", fontsize=14, fontweight="bold", y=1.03)
    save(fig, "1_matrix_qc_overview.png")
    sp.assign(loss_total=sp["L"]).sort_values("L", ascending=False)[["L", "UL", "M"]].to_csv(
        OUT / "species_loss_burden.tsv", sep="\t")


def fig_top100():
    base = ROOT / "results/agentic_triage"
    d = pd.read_csv(base / "ranked_co_loss_candidates_symbolic.tsv", sep="\t")
    core = pd.read_csv(base / "ranked_co_loss_candidates_no_histone.tsv", sep="\t")
    unl = pd.read_csv(base / "ranked_co_loss_candidates_unlinked.tsv", sep="\t")
    key = lambda x: x["gene_a"] + "|" + x["gene_b"]
    cat = np.where(key(d).isin(key(unl)), "Unlinked", np.where(key(d).isin(key(core)), "Linked (<1 Mb)", "Histone–histone"))
    d["category"] = cat
    colors = {"Histone–histone": GREY, "Linked (<1 Mb)": BLUE, "Unlinked": RED}
    fig, axs = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1.1, 1]})
    ax = axs[0]
    for c, col in colors.items():
        s = d[d["category"] == c]
        ax.scatter(s["shared_loss_count"], s["enrichment_naive"], s=45, color=col, alpha=0.85,
                   edgecolor="white", label=f"{c} (n={len(s)})")
    ax.set_xlabel("Species with both genes lost")
    ax.set_ylabel("Naive enrichment (observed / expected)")
    ax.legend(frameon=False)
    ax.set_title("Top 100 candidate co-loss pairs")
    ax = axs[1]
    top = d.head(20).iloc[::-1]
    lab = sym(top, "gene_a") + " – " + sym(top, "gene_b")
    ax.barh(lab, top["shared_loss_count"], color=[colors[c] for c in top["category"]])
    ax.set_xlabel("Species with both genes lost")
    ax.set_title("Top 20 by rank")
    ax.tick_params(axis="y", labelsize=8)
    fig.suptitle("Phenotype-blind co-loss screen: 25.9 million pairs scored, top 100 reported",
                 fontsize=14, fontweight="bold", y=1.0)
    save(fig, "2_top100_candidates.png")
    d.to_csv(OUT / "top100_with_category.tsv", sep="\t", index=False)


def fig_unlinked():
    u = pd.read_csv(ROOT / "results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv", sep="\t")
    u["pair"] = sym(u, "gene_a") + " – " + sym(u, "gene_b")
    u = u.sort_values("shared_loss_count")
    fig, axs = plt.subplots(1, 2, figsize=(14, 6.5), sharey=True)
    hit = u["pair"].str.contains("CIROZ")
    axs[0].barh(u["pair"], u["shared_loss_count"], color=np.where(hit, RED, BLUE))
    axs[0].set_xlabel("Species with both genes lost")
    axs[0].set_title("Shared losses")
    axs[1].barh(u["pair"], u["enrichment_naive"], color=np.where(hit, RED, BLUE))
    axs[1].axvline(1, color="k", ls="--", lw=1)
    axs[1].set_xlabel("Naive enrichment (1 = chance)")
    axs[1].set_title("Enrichment over expectation")
    fig.suptitle("19 physically unlinked co-loss pairs (red = involves CIROZ)", fontsize=14,
                 fontweight="bold", y=1.0)
    save(fig, "3_unlinked_pairs.png")


def fig_fisher():
    f = pd.read_csv(ROOT / "results/lr_module/lr_module_pairwise_fisher_bh.tsv", sep="\t")
    f["pair"] = f["gene_a"] + " – " + f["gene_b"]
    f = f.sort_values("fisher_p", ascending=False)
    fig, axs = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    axs[0].barh(f["pair"], -np.log10(f["fisher_p"]), color=RED)
    axs[0].axvline(-np.log10(0.05), color="k", ls="--", lw=1)
    axs[0].set_xlabel("-log10(Fisher p), one-sided")
    axs[0].set_title("Significance")
    w = 0.38
    y = np.arange(len(f))
    axs[1].barh(y - w / 2, f["both_lost"], w, color=RED, label="both lost")
    axs[1].barh(y + w / 2, f["only_a"] + f["only_b"], w, color=GREY, label="only one lost")
    axs[1].legend(frameon=False)
    axs[1].set_xlabel("Species (of 162)")
    axs[1].set_title("Co-loss counts")
    for i, (_, r) in enumerate(f.iterrows()):
        o = "inf" if np.isinf(r["odds_ratio"]) else f"{r['odds_ratio']:.0f}"
        axs[1].text(max(r["both_lost"], r["only_a"] + r["only_b"]) + 1, i, f"OR {o}", va="center", fontsize=8)
    axs[1].set_xlim(0, axs[1].get_xlim()[1] * 1.12)
    fig.suptitle("Left–right module: DAND5, CIROP, CIROZ, PKD1L1, MMP21 (all 10 pairs, BH q < 0.001)",
                 fontsize=14, fontweight="bold", y=1.0)
    save(fig, "4_LR_module_fisher.png")


def fig_ciroz():
    sp = pd.read_csv(ROOT / "results/lr_module/lr_module_species_status.tsv", sep="\t", index_col=0)
    lost = set((ROOT / "results/CIROZ_lost_species.txt").read_text().split())
    sp["ciroz_lost"] = sp.index.isin(lost)
    t = sp.groupby("clade")["ciroz_lost"].agg(lost="sum", total="count")
    t["kept"] = t["total"] - t["lost"]
    t = t.sort_values("total")
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5), gridspec_kw={"width_ratios": [1.3, 1]})
    axs[0].barh(t.index, t["lost"], color=RED, label="CIROZ lost")
    axs[0].barh(t.index, t["kept"], left=t["lost"], color="#cfe3d6", label="not lost")
    for i, (_, r) in enumerate(t.iterrows()):
        axs[0].text(r["total"] + 0.5, i, f"{int(r['lost'])}/{int(r['total'])}", va="center", fontsize=9)
    axs[0].set_xlabel("Species")
    axs[0].legend(frameon=False, loc="lower right")
    axs[0].set_title("CIROZ status by clade")
    ax = axs[1]
    n_other = sp.loc[sp["ciroz_lost"], ["DAND5", "CIROP", "PKD1L1", "MMP21"]].eq("L").sum(axis=1)
    ax.hist(n_other, bins=np.arange(-0.5, 5.5, 1), color=RED, edgecolor="white")
    ax.set_xlabel("Other module genes also lost (of 4)")
    ax.set_ylabel("CIROZ-lost species")
    ax.set_title(f"Co-loss within the {len(lost)} CIROZ-lost species")
    fig.suptitle("CIROZ is lost in 91 of 162 species, concentrated in a few clades", fontsize=14,
                 fontweight="bold", y=1.0)
    save(fig, "5_CIROZ_losses.png")


def main():
    fig_matrix_qc()
    fig_top100()
    fig_unlinked()
    fig_fisher()
    fig_ciroz()
    with PdfPages(OUT / "all_figures.pdf") as pdf:
        for f in FIGS:
            pdf.savefig(f, bbox_inches="tight")
    print(sorted(p.name for p in OUT.iterdir()))


if __name__ == "__main__":
    main()

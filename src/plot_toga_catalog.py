"""Plot assembly and taxon counts from the TOGA reference catalog.

This is a resource-coverage figure only. It is not a summary of Zoonomia-matched
species or gene-loss results.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "raw_data" / "toga_inventory" / "overview.table.tsv"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "toga_catalog_coverage.png"


def summarize_catalog(catalog: pd.DataFrame) -> pd.DataFrame:
    required_columns = {
        "Species Taxonomy ID",
        "Taxonomic Lineage",
        "Assembly name",
    }
    missing_columns = required_columns - set(catalog.columns)
    if missing_columns:
        raise ValueError(f"Catalog is missing columns: {sorted(missing_columns)}")
    if catalog.empty:
        raise ValueError("Catalog contains no assembly records.")

    lineage = catalog["Taxonomic Lineage"].fillna("").astype(str)
    is_placental = lineage.str.contains(r"(?:^|;\s*)Eutheria(?:;|$)", regex=True)
    summary = pd.DataFrame(
        {
            "catalog_group": is_placental.map(
                {
                    True: "Placental (Eutheria)",
                    False: "Non-Eutherian entries",
                }
            ),
            "taxon_id": catalog["Species Taxonomy ID"],
            "assembly": catalog["Assembly name"],
        }
    )
    return (
        summary.groupby("catalog_group", sort=False)
        .agg(
            assembly_records=("assembly", "count"),
            distinct_taxa=("taxon_id", "nunique"),
        )
        .reindex(["Placental (Eutheria)", "Non-Eutherian entries"])
        .fillna(0)
        .astype(int)
        .reset_index()
    )


def make_chart(summary: pd.DataFrame, output_path: Path) -> None:
    figure, axis = plt.subplots(figsize=(9.2, 5.4), facecolor="#f7f7f2")
    axis.set_facecolor("#f7f7f2")

    positions = range(len(summary))
    bar_height = 0.34
    assembly_bars = axis.barh(
        [position - bar_height / 2 for position in positions],
        summary["assembly_records"],
        height=bar_height,
        color="#176b68",
        label="Assembly records",
    )
    taxon_bars = axis.barh(
        [position + bar_height / 2 for position in positions],
        summary["distinct_taxa"],
        height=bar_height,
        color="#d77a35",
        label="Distinct taxonomy IDs",
    )

    axis.set_yticks(list(positions), summary["catalog_group"])
    axis.invert_yaxis()
    axis.set_xlabel("Count in TOGA human-reference catalog")
    axis.set_title(
        "TOGA mammal catalog coverage",
        loc="left",
        fontsize=17,
        fontweight="bold",
        color="#183b3b",
        pad=16,
    )
    axis.text(
        0,
        1.025,
        "Assembly records exceed taxa because some taxa have multiple assemblies",
        transform=axis.transAxes,
        fontsize=9.5,
        color="#526365",
    )
    axis.grid(axis="x", color="#c8d0ca", linewidth=0.8, alpha=0.65)
    axis.set_axisbelow(True)
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.spines["bottom"].set_color("#aab4af")
    axis.tick_params(axis="y", length=0, pad=10)
    axis.legend(frameon=False, loc="lower right")

    for bars in (assembly_bars, taxon_bars):
        axis.bar_label(bars, padding=4, fontsize=10, color="#263b3b")

    figure.text(
        0.09,
        0.025,
        "Archive inventory only. Includes non-Zoonomia taxa and alternate assemblies; not gene-loss results.",
        fontsize=8.5,
        color="#526365",
    )
    figure.tight_layout(rect=(0.04, 0.08, 0.98, 0.94))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    catalog = pd.read_csv(args.input, sep="\t", dtype={"Species Taxonomy ID": str})
    summary = summarize_catalog(catalog)
    make_chart(summary, args.output)
    summary_path = args.output.with_suffix(".tsv")
    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"Assembly records: {len(catalog):,}")
    print(f"Distinct taxa: {catalog['Species Taxonomy ID'].nunique():,}")
    print(f"Placental assembly records: {int(summary.loc[summary['catalog_group'] == 'Placental (Eutheria)', 'assembly_records'].iloc[0]):,}")
    print(f"Chart: {args.output}")
    print(f"Summary: {summary_path}")
    print("Interpretation: archive coverage only; not a Zoonomia-matched set or biological result.")


if __name__ == "__main__":
    main()
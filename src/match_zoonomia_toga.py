"""Match official Zoonomia mammal-tree tips to available TOGA assemblies.

Matching uses normalized scientific names only. Multiple TOGA assemblies for
one species are all retained in the manifest; none is selected automatically.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TREE = PROJECT_ROOT / "raw_data" / "zoonomia" / "241-mammalian-2020v2.phast-242.nh"
DEFAULT_CATALOG = PROJECT_ROOT / "raw_data" / "toga_inventory" / "overview.table.tsv"
DEFAULT_ROSTER = PROJECT_ROOT / "raw_data" / "zoonomia" / "zoonomia_mammal_assemblies.tsv"
DEFAULT_OUTPUT = PROJECT_ROOT / "orthologs" / "zoonomia_toga_assembly_manifest.tsv"
DEFAULT_ELIGIBLE_OUTPUT = PROJECT_ROOT / "orthologs" / "zoonomia_toga_eligible_assemblies.tsv"
ASSEMBLY_ACCESSION = re.compile(r"^GC[AF]_\d+\.\d+$")


def parse_newick_tips(text: str) -> list[str]:
    """Parse leaf labels from a Newick tree, ignoring internal node labels."""
    position = 0
    tips: list[str] = []

    def skip_space() -> None:
        nonlocal position
        while position < len(text) and text[position].isspace():
            position += 1

    def parse_label() -> str:
        nonlocal position
        skip_space()
        start = position
        while position < len(text) and text[position] not in ":,();":
            position += 1
        return text[start:position].strip().strip("'").strip('"')

    def skip_branch_length() -> None:
        nonlocal position
        skip_space()
        if position < len(text) and text[position] == ":":
            position += 1
            while position < len(text) and text[position] not in ",();":
                position += 1

    def parse_subtree() -> None:
        nonlocal position
        skip_space()
        if position >= len(text):
            raise ValueError("Unexpected end of Newick tree.")

        if text[position] == "(":
            position += 1
            while True:
                parse_subtree()
                skip_space()
                if position >= len(text):
                    raise ValueError("Unclosed parenthesis in Newick tree.")
                if text[position] == ",":
                    position += 1
                    continue
                if text[position] == ")":
                    position += 1
                    break
                raise ValueError(f"Unexpected Newick character at offset {position}.")
            parse_label()
            skip_branch_length()
            return

        label = parse_label()
        if not label:
            raise ValueError(f"Empty Newick tip at offset {position}.")
        tips.append(label)
        skip_branch_length()

    parse_subtree()
    skip_space()
    if position < len(text) and text[position] == ";":
        position += 1
    skip_space()
    if position != len(text):
        raise ValueError(f"Unexpected trailing Newick content at offset {position}.")
    if len(tips) != len(set(tips)):
        raise ValueError("Newick tree contains duplicate terminal labels.")
    return tips


def normalize_name(name: str) -> str:
    return " ".join(name.replace("_", " ").strip().split()).casefold()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_toga_catalog(path: Path) -> tuple[dict[str, list[dict[str, str]]], dict[str, list[dict[str, str]]]]:
    by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_accession: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Species",
            "Common name",
            "Species Taxonomy ID",
            "Assembly name",
            "NCBI accession / source",
            "contig N50 (bp)",
            "scaffold N50 (bp)",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"TOGA catalog is missing columns: {sorted(missing)}")
        for row in reader:
            by_species[normalize_name(row["Species"])].append(row)
            accession = row["NCBI accession / source"].strip()
            if ASSEMBLY_ACCESSION.fullmatch(accession):
                by_accession[accession].append(row)
    return by_species, by_accession


def read_zoonomia_roster(path: Path) -> tuple[dict[str, dict[str, str]], list[dict[str, str]]]:
    roster: list[dict[str, str]] = []
    by_name: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"scientific_name", "assembly_accession", "assembly_url"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Zoonomia roster is missing columns: {sorted(missing)}")
        for row in reader:
            normalized = normalize_name(row["scientific_name"])
            if normalized in by_name:
                raise ValueError(f"Duplicate Zoonomia roster species: {row['scientific_name']}")
            if not ASSEMBLY_ACCESSION.fullmatch(row["assembly_accession"]):
                raise ValueError(f"Invalid NCBI assembly accession for {row['scientific_name']}")
            by_name[normalized] = row
            roster.append(row)
    return by_name, roster


def find_roster_match(
    scientific_name: str,
    roster_by_name: dict[str, dict[str, str]],
    roster: list[dict[str, str]],
) -> tuple[dict[str, str] | None, str]:
    exact = roster_by_name.get(normalize_name(scientific_name))
    if exact:
        return exact, "exact_scientific_name"

    normalized = normalize_name(scientific_name)
    subspecies_matches = [
        row
        for row in roster
        if normalize_name(row["scientific_name"]).startswith(f"{normalized} ")
    ]
    if len(subspecies_matches) == 1:
        return subspecies_matches[0], "unique_binomial_to_subspecies_match"
    return None, "no_unique_roster_match"


def create_manifest(
    tree_path: Path,
    catalog_path: Path,
    roster_path: Path,
    output_path: Path,
    eligible_output_path: Path,
) -> dict[str, int | str]:
    tips = parse_newick_tips(tree_path.read_text(encoding="utf-8"))
    by_species, by_accession = read_toga_catalog(catalog_path)
    roster_by_name, roster = read_zoonomia_roster(roster_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    counts: dict[str, int] = defaultdict(int)
    eligible_rows: list[dict[str, str]] = []
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "zoonomia_tree_tip",
            "roster_scientific_name",
            "roster_match_method",
            "zoonomia_assembly_accession",
            "qc_status",
            "toga_candidate_count",
            "candidate_species_names",
            "candidate_assembly_names",
            "candidate_accessions_or_sources",
            "candidate_contig_n50_bp",
            "candidate_scaffold_n50_bp",
            "selected_assembly",
            "selected_accession",
            "selected_contig_n50_bp",
            "selected_scaffold_n50_bp",
            "selection_note",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for tip in sorted(tips):
            scientific_name = tip.replace("_", " ")
            roster_record, roster_method = find_roster_match(scientific_name, roster_by_name, roster)
            roster_name = roster_record["scientific_name"] if roster_record else ""
            zoonomia_accession = roster_record["assembly_accession"] if roster_record else ""
            exact_accession_matches = by_accession.get(zoonomia_accession, []) if roster_record else []
            same_species = by_species.get(normalize_name(roster_name), []) if roster_record else []
            candidates_by_identity = {
                (item["Assembly name"], item["NCBI accession / source"]): item
                for item in [*same_species, *exact_accession_matches]
            }
            candidates = list(candidates_by_identity.values())

            if normalize_name(scientific_name) == "homo sapiens":
                qc_status = "reference_species_not_query"
                selected = None
                note = "Human is the TOGA reference; there is no query-species loss summary for it."
            elif roster_record is None:
                qc_status = "excluded_no_unique_zoonomia_roster_match"
                selected = None
                note = "No exact or unique binomial/subspecies match in the official roster."
            elif len(exact_accession_matches) == 1:
                selected = exact_accession_matches[0]
                qc_status = (
                    "eligible_exact_accession_binomial_name_match"
                    if roster_method == "unique_binomial_to_subspecies_match"
                    else "eligible_exact_accession"
                )
                note = (
                    "Selected because the TOGA NCBI assembly accession exactly matches the Zoonomia roster accession."
                    if roster_method == "exact_scientific_name"
                    else "Assembly accession matches exactly; tree tip maps uniquely to this roster subspecies by binomial."
                )
            elif len(exact_accession_matches) > 1:
                qc_status = "excluded_duplicate_exact_accession_records"
                selected = None
                note = "More than one TOGA record has the Zoonomia assembly accession; needs duplicate-record QC."
            elif same_species:
                valid_accessions = [
                    row["NCBI accession / source"].strip()
                    for row in same_species
                    if ASSEMBLY_ACCESSION.fullmatch(row["NCBI accession / source"].strip())
                ]
                qc_status = (
                    "excluded_alternative_assembly_only"
                    if valid_accessions
                    else "excluded_unverified_assembly_source"
                )
                selected = None
                note = (
                    "TOGA has data for this taxon, but not for the exact Zoonomia assembly accession; "
                    "do not mix assembly versions in the primary matrix."
                    if valid_accessions
                    else "TOGA source does not provide a verifiable NCBI assembly accession."
                )
            else:
                qc_status = "excluded_exact_assembly_not_in_toga"
                selected = None
                note = "No TOGA record for the exact Zoonomia assembly accession or the matched taxon."

            counts[qc_status] += 1
            row = {
                "zoonomia_tree_tip": tip,
                "roster_scientific_name": roster_name,
                "roster_match_method": roster_method,
                "zoonomia_assembly_accession": zoonomia_accession,
                "qc_status": qc_status,
                "toga_candidate_count": str(len(candidates)),
                "candidate_species_names": ";".join(sorted({item["Species"] for item in candidates})),
                "candidate_assembly_names": ";".join(sorted({item["Assembly name"] for item in candidates})),
                "candidate_accessions_or_sources": ";".join(
                    sorted({item["NCBI accession / source"] for item in candidates})
                ),
                "candidate_contig_n50_bp": ";".join(
                    sorted({item["contig N50 (bp)"] for item in candidates})
                ),
                "candidate_scaffold_n50_bp": ";".join(
                    sorted({item["scaffold N50 (bp)"] for item in candidates})
                ),
                "selected_assembly": selected["Assembly name"] if selected else "",
                "selected_accession": selected["NCBI accession / source"] if selected else "",
                "selected_contig_n50_bp": selected["contig N50 (bp)"] if selected else "",
                "selected_scaffold_n50_bp": selected["scaffold N50 (bp)"] if selected else "",
                "selection_note": note,
            }
            writer.writerow(row)
            if selected:
                eligible_rows.append(row)

    eligible_output_path.parent.mkdir(parents=True, exist_ok=True)
    with eligible_output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(eligible_rows)

    summary: dict[str, int | str] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "tree_source": "https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2.phast-242.nh",
        "tree_file": str(tree_path.resolve()),
        "tree_sha256": sha256_file(tree_path),
        "zoonomia_roster_file": str(roster_path.resolve()),
        "zoonomia_roster_sha256": sha256_file(roster_path),
        "toga_catalog_source": "https://genome.senckenberg.de/download/TOGA/human_hg38_reference/overview.table.tsv",
        "toga_catalog_file": str(catalog_path.resolve()),
        "toga_catalog_sha256": sha256_file(catalog_path),
        "tree_tip_count": len(tips),
        "zoonomia_roster_entries": len(roster),
        "assembly_match_rule": "Select only TOGA rows whose NCBI assembly accession exactly matches the official Zoonomia roster accession.",
        "name_match_rule": "Exact scientific name, or a unique binomial-to-subspecies mapping flagged explicitly.",
        "eligible_assembly_count": len(eligible_rows),
        "qc_status_counts": dict(sorted(counts.items())),
    }
    summary_path = output_path.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, default=DEFAULT_TREE)
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--roster", type=Path, default=DEFAULT_ROSTER)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--eligible-output", type=Path, default=DEFAULT_ELIGIBLE_OUTPUT)
    args = parser.parse_args()

    for path in (args.tree, args.catalog, args.roster):
        if not path.is_file():
            raise FileNotFoundError(f"Required input does not exist: {path}")

    summary = create_manifest(
        args.tree,
        args.catalog,
        args.roster,
        args.output,
        args.eligible_output,
    )
    print(f"Tree tips: {summary['tree_tip_count']}")
    print(f"Eligible exact-assembly records: {summary['eligible_assembly_count']}")
    for status, count in summary["qc_status_counts"].items():
        print(f"{status}: {count}")
    print(f"Manifest: {args.output}")
    print(f"Eligible assemblies: {args.eligible_output}")
    print(f"Summary: {args.output.with_suffix('.summary.json')}")


if __name__ == "__main__":
    main()
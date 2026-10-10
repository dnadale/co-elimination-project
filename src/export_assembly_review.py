"""Export assembly records excluded by accession-based quality control."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "orthologs" / "zoonomia_toga_assembly_manifest.tsv"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "zoonomia_toga_assembly_review.tsv"
EXCLUDED_PREFIX = "excluded_"
OUTPUT_FIELDS = [
    "zoonomia_tree_tip",
    "roster_scientific_name",
    "roster_match_method",
    "zoonomia_assembly_accession",
    "qc_status",
    "toga_candidate_count",
    "candidate_assembly_names",
    "candidate_accessions_or_sources",
    "candidate_contig_n50_bp",
    "candidate_scaffold_n50_bp",
    "selection_note",
]


def export_review(input_path: Path, output_path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        required = set(OUTPUT_FIELDS) | {"qc_status"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Match manifest is missing columns: {sorted(missing)}")

        with output_path.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.DictWriter(
                destination,
                fieldnames=OUTPUT_FIELDS,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            for row in reader:
                status = row["qc_status"]
                if status.startswith(EXCLUDED_PREFIX):
                    writer.writerow(row)
                    counts[status] = counts.get(status, 0) + 1

    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(f"Assembly match manifest does not exist: {args.input}")

    counts = export_review(args.input, args.output)
    print(f"Excluded species rows: {sum(counts.values())}")
    for status, count in sorted(counts.items()):
        print(f"{status}: {count}")
    print(f"Review table: {args.output}")


if __name__ == "__main__":
    main()
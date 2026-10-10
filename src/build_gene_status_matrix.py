"""Build a rectangular gene-by-species matrix from validated TOGA summaries."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DOWNLOAD_MANIFEST = (
    PROJECT_ROOT / "raw_data" / "toga_gene_status" / "download_manifest.tsv"
)
DEFAULT_MATRIX = PROJECT_ROOT / "analysis" / "gene_status_matrix.tsv"
DEFAULT_QC = PROJECT_ROOT / "analysis" / "gene_status_matrix.qc.json"
EXPECTED_ENTITIES = {"GENE", "TRANSCRIPT", "PROJECTION"}
EXPECTED_STATUSES = {"I", "PI", "UL", "M", "PM", "L", "PG"}
STATUS_MEANINGS = {
    "I": "intact",
    "PI": "partially_intact",
    "UL": "uncertain_loss",
    "M": "missing",
    "PM": "partially_missing",
    "L": "lost",
    "PG": "paralogous_projection",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_download_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "zoonomia_tree_tip",
            "local_file",
            "bytes",
            "sha256",
            "gene_rows",
            "download_status",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Download manifest is missing columns: {sorted(missing)}")
        records = list(reader)

    if not records:
        raise ValueError("Download manifest contains no files.")
    tips = [record["zoonomia_tree_tip"] for record in records]
    if len(tips) != len(set(tips)):
        raise ValueError("Download manifest contains duplicate species/tree-tip rows.")
    failed = [record["zoonomia_tree_tip"] for record in records if record["download_status"] == "FAILED"]
    if failed:
        raise ValueError(f"Download manifest contains failures: {failed[:10]}")
    return sorted(records, key=lambda record: record["zoonomia_tree_tip"])


def read_gene_status_file(
    path: Path,
    tip: str,
    expected_sha256: str,
    expected_bytes: int,
    expected_gene_rows: int,
) -> tuple[dict[str, str], Counter[str], dict[str, int]]:
    if not path.is_file():
        raise FileNotFoundError(f"Status file for {tip} is missing: {path}")
    actual_bytes = path.stat().st_size
    if actual_bytes != expected_bytes:
        raise ValueError(f"Size mismatch for {tip}: expected {expected_bytes}, found {actual_bytes}.")
    actual_sha256 = sha256_file(path)
    if actual_sha256 != expected_sha256:
        raise ValueError(f"SHA-256 mismatch for {tip}.")

    gene_status: dict[str, str] = {}
    status_counts: Counter[str] = Counter()
    entity_counts: Counter[str] = Counter()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 3:
                raise ValueError(f"{tip}: expected 3 columns at line {line_number}.")
            entity, identifier, status = fields
            if entity not in EXPECTED_ENTITIES:
                raise ValueError(f"{tip}: unknown entity {entity!r} at line {line_number}.")
            if status not in EXPECTED_STATUSES:
                raise ValueError(f"{tip}: unknown status {status!r} at line {line_number}.")
            entity_counts[entity] += 1
            if entity == "GENE":
                if not identifier:
                    raise ValueError(f"{tip}: empty gene identifier at line {line_number}.")
                if identifier in gene_status:
                    raise ValueError(f"{tip}: duplicate GENE identifier {identifier!r}.")
                gene_status[identifier] = status
                status_counts[status] += 1

    if len(gene_status) != expected_gene_rows:
        raise ValueError(
            f"{tip}: expected {expected_gene_rows} GENE rows, found {len(gene_status)}."
        )
    return gene_status, status_counts, dict(entity_counts)


def build_matrix(download_manifest: Path, matrix_path: Path, qc_path: Path) -> dict[str, object]:
    records = read_download_manifest(download_manifest)
    species = [record["zoonomia_tree_tip"] for record in records]
    all_gene_status: dict[str, list[str]] = {}
    species_qc = []
    combined_status_counts: Counter[str] = Counter()
    reference_gene_ids: set[str] | None = None

    for species_index, record in enumerate(records):
        file_path = (PROJECT_ROOT / record["local_file"]).resolve()
        gene_status, status_counts, entity_counts = read_gene_status_file(
            file_path,
            record["zoonomia_tree_tip"],
            record["sha256"],
            int(record["bytes"]),
            int(record["gene_rows"]),
        )
        current_gene_ids = set(gene_status)
        if reference_gene_ids is None:
            reference_gene_ids = current_gene_ids
            for gene_id in sorted(reference_gene_ids):
                all_gene_status[gene_id] = [""] * len(records)
        elif current_gene_ids != reference_gene_ids:
            missing = sorted(reference_gene_ids - current_gene_ids)
            unexpected = sorted(current_gene_ids - reference_gene_ids)
            raise ValueError(
                f"Gene-ID set differs for {record['zoonomia_tree_tip']}: "
                f"missing={missing[:5]}, unexpected={unexpected[:5]}"
            )

        for gene_id, status in gene_status.items():
            all_gene_status[gene_id][species_index] = status
        combined_status_counts.update(status_counts)
        species_qc.append(
            {
                "species": record["zoonomia_tree_tip"],
                "gene_rows": len(gene_status),
                "entity_rows": entity_counts,
                "status_counts": dict(sorted(status_counts.items())),
                "source_file": record["local_file"],
                "source_sha256": record["sha256"],
            }
        )

    if not all_gene_status:
        raise ValueError("No gene records were loaded.")
    if any(not status for statuses in all_gene_status.values() for status in statuses):
        raise ValueError("Matrix contains empty statuses after species reconciliation.")

    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    with matrix_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *species])
        for gene_id in sorted(all_gene_status):
            writer.writerow([gene_id, *all_gene_status[gene_id]])

    matrix_qc = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_download_manifest": str(download_manifest.resolve()),
        "input_manifest_sha256": sha256_file(download_manifest),
        "matrix_file": str(matrix_path.resolve()),
        "matrix_sha256": sha256_file(matrix_path),
        "species_count": len(species),
        "gene_count": len(all_gene_status),
        "matrix_cell_count": len(all_gene_status) * len(species),
        "missing_cells": 0,
        "status_meanings": STATUS_MEANINGS,
        "combined_status_counts": dict(sorted(combined_status_counts.items())),
        "species_qc": species_qc,
        "scientific_status": "validated pairwise human-reference TOGA gene classifications; ancestry correction not yet applied",
    }
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.write_text(json.dumps(matrix_qc, indent=2) + "\n", encoding="utf-8")
    return matrix_qc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-manifest", type=Path, default=DEFAULT_DOWNLOAD_MANIFEST)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--qc-report", type=Path, default=DEFAULT_QC)
    args = parser.parse_args()
    if not args.download_manifest.is_file():
        raise FileNotFoundError(f"Download manifest does not exist: {args.download_manifest}")

    qc = build_matrix(args.download_manifest, args.matrix, args.qc_report)
    print(f"Species: {qc['species_count']:,}")
    print(f"Genes: {qc['gene_count']:,}")
    print(f"Status cells: {qc['matrix_cell_count']:,}")
    print(f"Matrix: {args.matrix}")
    print(f"QC report: {args.qc_report}")


if __name__ == "__main__":
    main()
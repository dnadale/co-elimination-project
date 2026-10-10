"""Download and validate exact-assembly TOGA loss-summary files.

Files are saved per Zoonomia tree tip. A rerun skips previously downloaded files
only when compressed size and the complete gzip/TSV contents pass QC.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "orthologs" / "zoonomia_toga_file_availability.tsv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "raw_data" / "toga_gene_status"
ALLOWED_ENTITIES = {"GENE", "TRANSCRIPT", "PROJECTION"}
ALLOWED_STATUSES = {"I", "PI", "UL", "M", "PM", "L", "PG"}
SAFE_TIP = re.compile(r"^[A-Za-z0-9_.-]+$")
HEADERS = {"User-Agent": "co-elimination-project reproducible data acquisition/1.0"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_gzip_tsv(path: Path) -> dict[str, int]:
    counts = {"GENE": 0, "TRANSCRIPT": 0, "PROJECTION": 0}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 3:
                raise ValueError(f"Expected 3 tab-separated fields at row {line_number}.")
            entity, identifier, status = fields
            if entity not in ALLOWED_ENTITIES:
                raise ValueError(f"Unknown entity {entity!r} at row {line_number}.")
            if not identifier:
                raise ValueError(f"Empty identifier at row {line_number}.")
            if status not in ALLOWED_STATUSES:
                raise ValueError(f"Unknown status {status!r} at row {line_number}.")
            counts[entity] += 1

    if counts["GENE"] == 0:
        raise ValueError("File has no GENE-level classifications.")
    return counts


def download_one(record: dict[str, str], output_dir: Path, retries: int) -> dict[str, str | int]:
    tip = record["zoonomia_tree_tip"]
    if not SAFE_TIP.fullmatch(tip):
        raise ValueError(f"Unsafe tree-tip label for filename: {tip!r}")
    if record["availability_qc_status"] != "status_file_verified":
        raise ValueError(f"Refusing unverified source for {tip}.")
    url = record["status_file_url"]
    if not url.startswith("https://genome.senckenberg.de/download/TOGA/"):
        raise ValueError(f"Unexpected TOGA URL for {tip}: {url}")

    expected_size = int(record["status_file_bytes"])
    output_path = output_dir / f"{tip}.loss_summ_data.tsv.gz"
    if output_path.is_file() and output_path.stat().st_size == expected_size:
        try:
            entity_counts = validate_gzip_tsv(output_path)
            return {
                "zoonomia_tree_tip": tip,
                "url": url,
                "local_file": str(output_path.relative_to(PROJECT_ROOT)),
                "http_status": 200,
                "bytes": output_path.stat().st_size,
                "sha256": sha256_file(output_path),
                "gene_rows": entity_counts["GENE"],
                "transcript_rows": entity_counts["TRANSCRIPT"],
                "projection_rows": entity_counts["PROJECTION"],
                "download_status": "SKIPPED_ALREADY_VALID",
                "error": "",
            }
        except (OSError, EOFError, UnicodeDecodeError, ValueError):
            pass

    last_error = "download failed"
    for attempt in range(1, retries + 1):
        temporary_path: Path | None = None
        try:
            request = Request(url, headers=HEADERS)
            with urlopen(request, timeout=120) as response:
                if response.status != 200:
                    raise RuntimeError(f"Unexpected HTTP status: {response.status}")
                fd, temp_name = tempfile.mkstemp(
                    prefix=f".{tip}.", suffix=".part", dir=output_dir
                )
                temporary_path = Path(temp_name)
                with os.fdopen(fd, "wb") as target:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        target.write(chunk)
            actual_size = temporary_path.stat().st_size
            if actual_size != expected_size:
                raise ValueError(f"Expected {expected_size} bytes; received {actual_size}.")
            entity_counts = validate_gzip_tsv(temporary_path)
            os.replace(temporary_path, output_path)
            temporary_path = None
            return {
                "zoonomia_tree_tip": tip,
                "url": url,
                "local_file": str(output_path.relative_to(PROJECT_ROOT)),
                "http_status": 200,
                "bytes": actual_size,
                "sha256": sha256_file(output_path),
                "gene_rows": entity_counts["GENE"],
                "transcript_rows": entity_counts["TRANSCRIPT"],
                "projection_rows": entity_counts["PROJECTION"],
                "download_status": "DOWNLOADED_AND_VALIDATED",
                "error": "",
            }
        except (HTTPError, URLError, TimeoutError, OSError, EOFError, UnicodeDecodeError, ValueError, RuntimeError) as error:
            last_error = f"attempt {attempt}/{retries}: {error}"
        finally:
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)

    return {
        "zoonomia_tree_tip": tip,
        "url": url,
        "local_file": "",
        "http_status": "FAILED",
        "bytes": 0,
        "sha256": "",
        "gene_rows": 0,
        "transcript_rows": 0,
        "projection_rows": 0,
        "download_status": "FAILED",
        "error": last_error,
    }


def read_verified_sources(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"zoonomia_tree_tip", "availability_qc_status", "status_file_url", "status_file_bytes"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Availability manifest is missing columns: {sorted(missing)}")
        return [row for row in reader if row["availability_qc_status"] == "status_file_verified"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--limit", type=int, help="Optional small count for a smoke test")
    args = parser.parse_args()
    if args.workers < 1 or args.retries < 1:
        raise ValueError("--workers and --retries must be positive.")
    if not args.input.is_file():
        raise FileNotFoundError(f"Availability manifest is missing: {args.input}")

    records = read_verified_sources(args.input)
    records.sort(key=lambda row: row["zoonomia_tree_tip"])
    if args.limit is not None:
        if args.limit < 1:
            raise ValueError("--limit must be positive.")
        records = records[: args.limit]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(download_one, row, args.output_dir, args.retries) for row in records]
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"{result['download_status']}: {result['zoonomia_tree_tip']}")
    results.sort(key=lambda row: str(row["zoonomia_tree_tip"]))

    manifest_path = args.output_dir / "download_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    failures = sum(result["download_status"] == "FAILED" for result in results)
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_availability_manifest": str(args.input.resolve()),
        "requested_count": len(records),
        "downloaded_and_validated": sum(
            result["download_status"] == "DOWNLOADED_AND_VALIDATED" for result in results
        ),
        "skipped_already_valid": sum(
            result["download_status"] == "SKIPPED_ALREADY_VALID" for result in results
        ),
        "failures": failures,
        "total_bytes": sum(int(result["bytes"]) for result in results),
        "manifest": str(manifest_path.resolve()),
    }
    summary_path = args.output_dir / "download_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Validated downloads: {summary['downloaded_and_validated']}")
    print(f"Already valid: {summary['skipped_already_valid']}")
    print(f"Failures: {failures}")
    print(f"Download manifest: {manifest_path}")
    print(f"Summary: {summary_path}")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
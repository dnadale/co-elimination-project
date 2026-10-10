"""Verify selected Zoonomia-matched TOGA gene-status summaries exist."""

from __future__ import annotations

import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "orthologs" / "zoonomia_toga_eligible_assemblies.tsv"
DEFAULT_OUTPUT = PROJECT_ROOT / "orthologs" / "zoonomia_toga_file_availability.tsv"
REFERENCE_URL = "https://genome.senckenberg.de/download/TOGA/human_hg38_reference/"
ORDER_DIRECTORIES = [
    "Afrotheria",
    "Carnivora",
    "Chiroptera",
    "Dermoptera",
    "Eulipotyphla",
    "Lagomorpha",
    "Metatheria",
    "Perissodactyla",
    "Pholidota",
    "Primates",
    "Prototheria",
    "Rodentia",
    "Ruminantia",
    "Scandentia",
    "Suina",
    "Tylopoda",
    "Whippomorpha",
    "Xenarthra",
]
HEADERS = {"User-Agent": "co-elimination-project reproducibility QC/1.0"}


class DirectoryLinks(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.directories: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        href = dict(attrs).get("href")
        if href and href.endswith("/") and href not in {"../", "/"}:
            self.directories.append(href)


def fetch_bytes(url: str, method: str = "GET") -> tuple[int, dict[str, str], bytes]:
    request = Request(url, headers=HEADERS, method=method)
    with urlopen(request, timeout=60) as response:
        return response.status, dict(response.headers.items()), response.read() if method == "GET" else b""


def list_order_directories(order: str) -> tuple[str, list[tuple[str, str]]]:
    order_url = urljoin(REFERENCE_URL, f"{order}/")
    status, _, content = fetch_bytes(order_url)
    parser = DirectoryLinks()
    parser.feed(content.decode("utf-8", "replace"))
    candidates = []
    for link in parser.directories:
        directory_name = link.rstrip("/").rsplit("/", 1)[-1]
        if "__" in directory_name:
            assembly = directory_name.rsplit("__", 1)[-1]
            directory_url = urljoin(order_url, link)
            candidates.append((assembly, directory_url))
    return order, candidates


def build_assembly_directory_map(max_workers: int) -> tuple[dict[str, list[tuple[str, str]]], dict[str, str]]:
    by_assembly: dict[str, list[tuple[str, str]]] = {}
    failures: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(list_order_directories, order): order for order in ORDER_DIRECTORIES}
        for future in as_completed(futures):
            order = futures[future]
            try:
                _, records = future.result()
                for assembly, directory_url in records:
                    by_assembly.setdefault(assembly, []).append((order, directory_url))
            except Exception as error:
                failures[order] = str(error)
    return by_assembly, failures


def check_status_file(record: dict[str, str], candidates: list[tuple[str, str]]) -> dict[str, str]:
    assembly = record["selected_assembly"]
    file_results = []
    for order, directory_url in candidates:
        url = urljoin(directory_url, "loss_summ_data.tsv.gz")
        try:
            status, headers, _ = fetch_bytes(url, method="HEAD")
            file_results.append((status, headers.get("Content-Length", ""), order, url))
        except HTTPError as error:
            file_results.append((error.code, "", order, url))
        except (URLError, TimeoutError) as error:
            file_results.append(("ERROR", str(error), order, url))

    successes = [result for result in file_results if result[0] == 200]
    if len(successes) == 1:
        status, byte_count, order, url = successes[0]
        qc_status = "status_file_verified"
        note = "Exactly one matching assembly directory contains the gene-status summary."
    elif len(successes) > 1:
        status, byte_count, order, url = successes[0]
        qc_status = "duplicate_status_files"
        note = "Multiple archive directories contain summaries for this assembly label."
    elif file_results:
        status, byte_count, order, url = file_results[0]
        qc_status = "status_file_unavailable"
        note = "Assembly directory was located, but its gene-status summary was not available."
    else:
        status, byte_count, order, url = "NOT_FOUND", "", "", ""
        qc_status = "assembly_directory_not_found"
        note = "No matching assembly directory was found in the searched TOGA order indexes."

    return {
        **record,
        "archive_order": order,
        "status_file_http_status": str(status),
        "status_file_bytes": byte_count,
        "status_file_url": url,
        "availability_qc_status": qc_status,
        "availability_qc_note": note,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be positive.")
    if not args.input.is_file():
        raise FileNotFoundError(f"Eligible assembly manifest is missing: {args.input}")

    with args.input.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        records = list(reader)
    if not records:
        raise ValueError("Eligible assembly manifest has no records.")

    directory_map, index_failures = build_assembly_directory_map(args.workers)
    output_rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                check_status_file,
                record,
                directory_map.get(record["selected_assembly"], []),
            )
            for record in records
        ]
        for future in as_completed(futures):
            output_rows.append(future.result())
    output_rows.sort(key=lambda row: row["zoonomia_tree_tip"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output_rows[0])
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    counts: dict[str, int] = {}
    for row in output_rows:
        status = row["availability_qc_status"]
        counts[status] = counts.get(status, 0) + 1
    summary = {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "eligible_input_count": len(records),
        "order_indexes_checked": len(ORDER_DIRECTORIES),
        "order_index_failures": index_failures,
        "availability_qc_counts": counts,
        "verified_status_file_count": counts.get("status_file_verified", 0),
        "manifest": str(args.output.resolve()),
    }
    summary_path = args.output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    for status, count in sorted(counts.items()):
        print(f"{status}: {count}")
    print(f"Availability manifest: {args.output}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
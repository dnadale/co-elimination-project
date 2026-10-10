"""Save the official Zoonomia mammal roster and its assembly accessions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.request import Request, urlopen


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PAGE_API_URL = "https://zoonomiaproject.org/wp-json/wp/v2/pages/418"
DEFAULT_RAW_JSON = PROJECT_ROOT / "raw_data" / "zoonomia" / "zoonomia_mammal_list_view.json"
DEFAULT_ROSTER_TSV = PROJECT_ROOT / "raw_data" / "zoonomia" / "zoonomia_mammal_assemblies.tsv"
ASSEMBLY_ACCESSION = re.compile(r"^GC[AF]_\d+\.\d+$")


class ZoonomiaSpeciesParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.depth = 0
        self.active = False
        self.in_species_heading = False
        self.current: dict[str, object] | None = None
        self.records: list[dict[str, str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        classes = set((attributes.get("class") or "").split())
        if tag == "div" and "tree" in classes and not self.active:
            self.active = True
            self.depth = 1
            self.current = {"scientific_name": "", "assembly_accession": ""}
            return

        if not self.active or self.current is None:
            return
        if tag == "div":
            self.depth += 1
        if tag == "h3" and "species" in classes:
            self.in_species_heading = True
        if tag == "a":
            href = attributes.get("href") or ""
            if "ncbi.nlm.nih.gov/assembly/" in href:
                accession = href.rstrip("/").rsplit("/", 1)[-1]
                if ASSEMBLY_ACCESSION.fullmatch(accession):
                    self.current["assembly_accession"] = accession

    def handle_data(self, data: str) -> None:
        if self.active and self.in_species_heading and self.current is not None:
            current_name = str(self.current["scientific_name"])
            self.current["scientific_name"] = f"{current_name} {data.strip()}".strip()

    def handle_endtag(self, tag: str) -> None:
        if not self.active or self.current is None:
            return
        if tag == "h3":
            self.in_species_heading = False
        if tag == "div":
            self.depth -= 1
            if self.depth == 0:
                name = " ".join(str(self.current["scientific_name"]).split())
                accession = str(self.current["assembly_accession"])
                if name and accession:
                    self.records.append(
                        {
                            "scientific_name": name,
                            "assembly_accession": accession,
                            "assembly_url": f"https://www.ncbi.nlm.nih.gov/assembly/{accession}",
                        }
                    )
                self.current = None
                self.active = False


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-json", type=Path, default=DEFAULT_RAW_JSON)
    parser.add_argument("--output", type=Path, default=DEFAULT_ROSTER_TSV)
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Fetch the latest official page instead of reusing the saved API response",
    )
    args = parser.parse_args()

    if args.refresh or not args.raw_json.exists():
        request = Request(
            PAGE_API_URL,
            headers={"User-Agent": "co-elimination-project metadata inventory/1.0"},
        )
        with urlopen(request, timeout=60) as response:
            raw_json = response.read()
        args.raw_json.parent.mkdir(parents=True, exist_ok=True)
        args.raw_json.write_bytes(raw_json)
    else:
        raw_json = args.raw_json.read_bytes()

    page = json.loads(raw_json)
    content = page.get("content", {}).get("rendered", "")
    if not content:
        raise ValueError("Official Zoonomia page response has no rendered page content.")

    species_parser = ZoonomiaSpeciesParser()
    species_parser.feed(content)
    records = species_parser.records
    names = [record["scientific_name"].casefold() for record in records]
    if not records or len(names) != len(set(names)):
        raise ValueError("Roster extraction returned no records or duplicate scientific names.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["scientific_name", "assembly_accession", "assembly_url"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(sorted(records, key=lambda row: row["scientific_name"].casefold()))

    summary = {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": PAGE_API_URL,
        "source_page_url": "https://zoonomiaproject.org/the-mammal-tree-list-view/",
        "source_page_id": page.get("id"),
        "source_sha256": sha256_bytes(raw_json),
        "roster_count": len(records),
        "output": str(args.output.resolve()),
    }
    summary_path = args.output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"Official species/assembly records: {len(records)}")
    print(f"Saved roster: {args.output}")
    print(f"Saved provenance: {summary_path}")


if __name__ == "__main__":
    main()
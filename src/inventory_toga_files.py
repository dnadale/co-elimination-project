"""Inventory TOGA documentation and selected public directory listings."""

from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin

import pandas as pd
import requests

PROJECT_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_DIR / "raw_data" / "toga_inventory"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
BASE_URL = "https://genome.senckenberg.de/download/TOGA/"
HEADERS = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120 Safari/537.36"}


class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.links.append(href)


def main():
    session = requests.Session()
    session.headers.update(HEADERS)

    documentation = {
        "README.txt": "README.txt",
        "README.download.txt": "README.download.txt",
        "README.multipleCodonAlignments.txt": "README.multipleCodonAlignments.txt",
    }
    for filename, relative_url in documentation.items():
        response = session.get(urljoin(BASE_URL, relative_url), timeout=(20, 90))
        response.raise_for_status()
        (OUTPUT_DIR / filename).write_bytes(response.content)
        print(f"Saved {filename}: {len(response.content):,} bytes")

    directories = {
        "main": BASE_URL,
        "human_reference": urljoin(BASE_URL, "human_hg38_reference/"),
        "mouse_reference": urljoin(BASE_URL, "mouse_mm10_reference/"),
    }
    rows = []

    for directory_name, directory_url in directories.items():
        response = session.get(directory_url, timeout=(20, 90))
        response.raise_for_status()
        parser = LinkParser()
        parser.feed(response.text)

        for href in sorted(set(parser.links)):
            if href in {"../", "/", "#"}:
                continue
            name = href.rstrip("/").split("/")[-1]
            lower_name = name.lower()
            rows.append(
                {
                    "directory": directory_name,
                    "name": name,
                    "url": urljoin(directory_url, href),
                    "is_directory": href.endswith("/"),
                    "likely_relevant": any(
                        keyword in lower_name
                        for keyword in [
                            "loss", "ortholog", "annotation", "gtf", "bed",
                            "codon", "gene", "isoform", "mutation",
                        ]
                    ),
                }
            )

    inventory = pd.DataFrame(rows)
    inventory["checked_at_utc"] = datetime.now(timezone.utc).isoformat()
    output_file = OUTPUT_DIR / "toga_available_files.tsv"
    inventory.to_csv(output_file, sep="\t", index=False)
    print(f"TOGA inventory saved to: {output_file}")
    print(f"Listed entries: {len(inventory):,}")


if __name__ == "__main__":
    main()

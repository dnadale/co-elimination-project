"""Download small official Zoonomia and TOGA source pages.

This script deliberately does not download genome assemblies or whole-genome
alignments. It records source URLs, HTTP status, sizes, timestamps, and hashes.
"""

from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import importlib
import time

import pandas as pd
import requests

PROJECT_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_DIR / "raw_data" / "source_inventory"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SOURCES = {
    "zoonomia_data_page": [
        "https://zoonomiaproject.org/the-data/",
        "https://zoonomiaproject.org/the-data",
        "https://zoonomiaproject.org/",
    ],
    "zoonomia_alignment_page": [
        "https://cglgenomics.ucsc.edu/data/cactus/",
        "http://cglgenomics.ucsc.edu/data/cactus/",
    ],
    "toga_download_index": [
        "https://genome.senckenberg.de/download/TOGA/",
        "http://genome.senckenberg.de/download/TOGA/",
    ],
    "toga_repository": ["https://github.com/hillerlab/TOGA2"],
}

HEADERS = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120 Safari/537.36",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.9",
}


def check_packages():
    packages = {
        "numpy": "numpy",
        "pandas": "pandas",
        "scipy": "scipy",
        "scikit-learn": "sklearn",
        "biopython": "Bio",
        "pyyaml": "yaml",
        "requests": "requests",
    }
    print("Checking Python packages:")
    for display_name, import_name in packages.items():
        try:
            importlib.import_module(import_name)
            print(f"  OK      {display_name}")
        except ImportError:
            print(f"  MISSING {display_name}")


def main():
    check_packages()
    session = requests.Session()
    session.headers.update(HEADERS)
    rows = []

    for source_name, candidate_urls in SOURCES.items():
        output_file = OUTPUT_DIR / f"{source_name}.html"
        success = False

        for url in candidate_urls:
            try:
                response = session.get(url, timeout=(20, 90))
                response.raise_for_status()
                content = response.content
                if len(content) > 50_000_000:
                    raise RuntimeError("response exceeded 50 MB metadata limit")

                output_file.write_bytes(content)
                rows.append(
                    {
                        "source_name": source_name,
                        "url": response.url,
                        "status": response.status_code,
                        "bytes": len(content),
                        "local_file": str(output_file.relative_to(PROJECT_DIR)),
                        "sha256": sha256(content).hexdigest(),
                        "error": "",
                    }
                )
                print(f"SUCCESS: {source_name} ({len(content):,} bytes)")
                success = True
                break
            except Exception as error:
                print(f"FAILED: {source_name} via {url}: {error}")
                time.sleep(2)

        if not success:
            rows.append(
                {
                    "source_name": source_name,
                    "url": "",
                    "status": "FAILED",
                    "bytes": 0,
                    "local_file": "",
                    "sha256": "",
                    "error": "All candidate URLs failed",
                }
            )

    manifest = pd.DataFrame(rows)
    manifest["checked_at_utc"] = datetime.now(timezone.utc).isoformat()
    manifest_path = PROJECT_DIR / "raw_data" / "source_inventory_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    print(f"Manifest saved to: {manifest_path}")


if __name__ == "__main__":
    main()

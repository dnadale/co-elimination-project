#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"

TREE="raw_data/zoonomia/241-mammalian-2020v2.phast-242.nh"
CATALOG="raw_data/toga_inventory/overview.table.tsv"
ROSTER="raw_data/zoonomia/zoonomia_mammal_assemblies.tsv"

for input_file in "$TREE" "$CATALOG"; do
    if [[ ! -s "$input_file" ]]; then
        printf 'Required input is missing or empty: %s\n' "$input_file" >&2
        exit 1
    fi
done

PYTHON="${PYTHON:-python3}"
if ! command -v "$PYTHON" >/dev/null 2>&1; then
    printf 'Python executable not found: %s\n' "$PYTHON" >&2
    exit 1
fi

"$PYTHON" -m py_compile src/match_zoonomia_toga.py
"$PYTHON" -m py_compile src/inventory_zoonomia_roster.py
"$PYTHON" -m py_compile src/export_assembly_review.py
"$PYTHON" -m py_compile src/qc_toga_availability.py

if [[ ! -s "$ROSTER" ]]; then
    "$PYTHON" src/inventory_zoonomia_roster.py
fi

"$PYTHON" src/match_zoonomia_toga.py \
    --tree "$TREE" \
    --catalog "$CATALOG" \
    --roster "$ROSTER" \
    --output orthologs/zoonomia_toga_assembly_manifest.tsv \
    --eligible-output orthologs/zoonomia_toga_eligible_assemblies.tsv

"$PYTHON" src/export_assembly_review.py
"$PYTHON" src/qc_toga_availability.py

printf '\nAudit and eligible manifests are saved under orthologs/. Exclusions are in results/.\n'
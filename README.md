# co_elimination_project

Comparative genomics project for coordinated gene co-elimination / phenotype-agnostic pathway discovery.

## Current status

The project environment is set up and the official Zoonomia and TOGA metadata have been inventoried. Public Zoonomia resources include the 240-species comparative alignment, phylogeny, conservation scores, assemblies, and TOGA annotations. The next step is to select a manageable set of mammalian assemblies and download the relevant precomputed TOGA status and orthology files.

No large genome or whole-genome alignment files have been downloaded yet. See [PROGRESS_2026-09-20.md](PROGRESS_2026-09-20.md) for the current milestone record.

## Directory layout
- data/            : project data
- raw_data/        : raw input files
- genomes/         : genome FASTA files
- alignments/      : multiple sequence alignments
- orthologs/       : ortholog mapping outputs
- analysis/        : intermediate analysis outputs
- results/         : final tables, figures, summaries
- logs/            : log files
- notebooks/       : notebooks for exploration
- src/             : analysis scripts

## Data inventory

- `raw_data/source_inventory/`: official Zoonomia and TOGA source pages
- `raw_data/source_inventory_manifest.tsv`: source URLs, status codes, sizes, timestamps, and checksums
- `raw_data/toga_inventory/`: TOGA documentation and available-file inventory

## Reproducibility

- `notebooks/Project1.ipynb`: permanent notebook entry point
- `src/inventory_sources.py`: checks core packages and inventories official source pages
- `src/inventory_toga_files.py`: records TOGA documentation and public directory listings
- `requirements.txt`: baseline Python dependencies
- `.gitignore`: excludes large genomes, alignments, compressed result files, caches, and secrets

Run the inventory scripts from the project root with:

```bash
python src/inventory_sources.py
python src/inventory_toga_files.py
```

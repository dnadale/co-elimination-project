# Project Goal

## Main goal

Identify genes or biological pathways that are repeatedly lost together across genomes, using comparative genomics to discover coordinated gene co-elimination patterns without relying on a predefined phenotype.

## Research question

Which genes or pathways show statistically meaningful patterns of joint presence and absence across a carefully selected set of species, and what might those patterns reveal about functional dependencies or biological adaptation?

## Objectives

1. Assemble and document a reproducible set of genomes and relevant raw data.
2. Identify comparable genes or orthologous groups across the selected species.
3. Build a presence/absence matrix for genes or ortholog groups.
4. Detect pairs or groups of genes that are lost together more often than expected by chance.
5. Evaluate candidate relationships with statistical tests and appropriate controls for shared ancestry.
6. Summarize the strongest candidate pathways and provide tables, figures, and reproducible analysis notebooks.

## Planned workflow

1. Define the species set and inclusion criteria.
2. Collect genome and annotation data in `raw_data/` and `genomes/`.
3. Standardize and quality-check the inputs.
4. Infer ortholog groups and save mappings in `orthologs/`.
5. Generate alignments or other sequence comparisons in `alignments/` when needed.
6. Construct presence/absence data and intermediate files in `analysis/`.
7. Test co-elimination patterns and produce final outputs in `results/`.
8. Record commands, parameters, and errors in `logs/` and keep the workflow documented in `notebooks/` and `src/`.

## Expected outputs

- A documented genome and species inventory.
- Ortholog-group assignments and quality-control summaries.
- A gene or ortholog-group presence/absence matrix.
- Ranked candidate co-elimination pairs or pathways.
- Statistical results and visualizations.
- A reproducible notebook and supporting scripts.

## Success criteria

The project will be considered successful when the analysis produces a reproducible, documented set of statistically supported co-elimination candidates, along with enough information to inspect the underlying genomes, ortholog assignments, assumptions, and limitations.

## Current status

The project environment and directory structure are set up. Data collection, ortholog inference, co-elimination testing, and biological interpretation are next.

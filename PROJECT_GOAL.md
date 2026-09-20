# Project Goal

## Main goal

Identify genes or biological pathways that are repeatedly lost together across placental mammal genomes, using Zoonomia-scale comparative genomics to discover coordinated gene co-elimination patterns without relying on a predefined phenotype.

## Research question

Which genes or pathways show statistically meaningful patterns of joint retention, disruption, and loss across approximately 240 placental mammal species spanning about 100 million years of evolution, and what might those patterns reveal about functional dependencies or biological adaptation?

## Objectives

1. Assemble and document a reproducible pilot set and, eventually, a Zoonomia-scale set of placental mammal genomes.
2. Identify comparable genes or orthologous groups across the selected species.
3. Use TOGA and/or CESAR, where feasible, to identify coding-sequence disruptions and gene-inactivation evidence.
4. Build a phylogenetically informed matrix of gene retention, disruption, and loss.
5. Detect pairs or groups of genes that are lost together more often than expected by chance.
6. Use phyloP and other conservation measures as supporting evidence when interpreting candidate modules.
7. Evaluate candidate relationships with statistical tests and appropriate controls for shared ancestry.
8. Summarize the strongest candidate pathways and provide tables, figures, and reproducible analysis notebooks.

## Planned workflow

1. Define the pilot species set and inclusion criteria, then plan expansion to the full Zoonomia resource.
2. Collect genome, annotation, alignment, and phylogenetic data in `raw_data/` and `genomes/`.
3. Standardize and quality-check the inputs.
4. Run or obtain TOGA/CESAR gene-status calls where feasible.
5. Infer ortholog groups and save mappings in `orthologs/`.
6. Generate alignments or other sequence comparisons in `alignments/` when needed.
7. Construct phylogenetically informed retention/loss data and intermediate files in `analysis/`.
8. Test co-elimination patterns and produce final outputs in `results/`.
9. Record commands, parameters, software versions, and errors in `logs/` and keep the workflow documented in `notebooks/` and `src/`.

## Expected outputs

- A documented Zoonomia genome and species inventory.
- Ortholog-group assignments and quality-control summaries.
- TOGA/CESAR-informed gene-status calls where available.
- A phylogenetically informed gene or ortholog-group retention/loss matrix.
- Ranked candidate co-elimination pairs or pathways.
- Statistical results, conservation summaries, and visualizations.
- A reproducible notebook and supporting scripts.

## Success criteria

The project will be considered successful when the analysis produces a reproducible, documented set of statistically supported co-elimination candidates, along with enough information to inspect the underlying genomes, ortholog assignments, assumptions, and limitations.

## Scope boundary for the master's project

The minimum thesis deliverable is a reproducible pilot consisting of Zoonomia data, orthology and gene-status evidence, a retention/loss matrix, phylogenetically informed co-elimination scoring, and constraint or functional interpretation. Expression, alternative splicing, TACIT enhancer prediction, and experimental validation are valuable extensions but should not be required for completion unless compatible data and pipelines are already available.

The known CIROP/DAND5 left-right patterning module may be used as a positive-control case. Recovering a known evolutionary profile before searching for new modules will provide an important pipeline validation step.

## Current status

The project environment and directory structure are set up. Official Zoonomia and TOGA metadata have been successfully accessed and recorded locally. The TOGA archive confirms that reusable orthology, gene-status, mutation, annotation, and alignment outputs exist. No large genome or alignment files have been downloaded yet. The next step is to select a manageable set of mammalian assemblies and retrieve the relevant precomputed TOGA files before constructing the retention/loss matrix.

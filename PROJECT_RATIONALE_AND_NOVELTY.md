# Project rationale and novelty summary

## Core idea
This project is a phenotype-agnostic screen for coordinated gene co-elimination across placental mammals. The goal is to identify gene pairs or modules that show unusually frequent shared loss or retention patterns across a large comparative dataset, using phylogenetically informed tests to distinguish genuine evolutionary coupling from shared ancestry.

## Why this is interesting
Complex biological systems are often encoded by interacting gene networks. If a pathway or developmental process becomes dispensable or is rewired in a lineage, genes within that module may be lost or retained together across species. By scanning across placental mammals without first assuming a phenotype, the project aims to discover biologically meaningful modules that might not have been obvious from a trait-first study.

## Relation to the known left-right module
The left-right patterning module involving CIROP/DAND5 and related genes is a motivating positive-control case. A successful pipeline should recover that type of known evolutionary footprint before exploring less-characterized candidate modules. This is not the entire project goal; it is a benchmark for the analytic framework.

## Why wet-lab validation matters
The computational screen produces ranked candidate modules, not final biological proof. Candidate genes must be followed up with CRISPR perturbation, functional assays, or orthogonal genetic evidence to determine whether the evolutionary signal reflects a true functional dependency.

## Novelty
This project is novel in its explicit phenotype-agnostic, Zoonomia-scale, phylogenetically informed co-elimination framework. It is not simply a follow-up on one disease gene or one known pathway. The combination of broad comparative screening, explicit assembly QC and metadata control, and ranked hypothesis generation is the defining novelty.

## Expected impact
The expected impact is not to solve most diseases directly, but to generate a set of testable biological hypotheses: conserved gene dependencies, pathway modules, and candidate genes that deserve functional follow-up in model systems. It is most valuable as a discovery and prioritization framework rather than a final clinical solution.

# Bruno project brief: coordinated gene co-elimination across placental mammals

## Scope
This project is a phenotype-agnostic comparative genomics screen for genes or pathways that are repeatedly lost or retained together across placental mammal genomes.

The goal is not to predefine a phenotype. Instead, the workflow is to build a broad, phylogenetically informed gene-status matrix and then identify candidate co-elimination modules that may reflect shared biological function, dependency, or evolutionary constraint.

## Central biological motivation
A motivating positive-control case is a conserved left-right patterning module involving CIROP/DAND5 and related genes (for example MMP21, PKD1L1, and CIROZ). The project may test whether this known module is recoverable from the data before exploring previously uncharacterized gene modules.

## Project logic
1. Assemble and quality-check a feasible Zoonomia-scale mammal dataset.
2. Build a gene retention/loss matrix across species.
3. Detect pairs or modules with unusually frequent shared loss.
4. Use phylogeny-aware statistical tests to guard against false positives from shared ancestry.
5. Rank candidate modules and interpret them as biological hypotheses, not as final proof.

## Important caveat
The computationally detected gene-loss patterns are evidence for candidate modules, not definitive proof of a biological mechanism. They should be treated as ranked hypotheses requiring follow-up validation.

## Alignment with current work
The current project is consistent with this brief:
- exact-accession QC on the Zoonomia/TOGA data
- construction of a validated mammal genome set
- gene-status matrix generation
- raw shared-loss candidate ranking
- phylogenetic validation of ranked candidates

The work is therefore aligned with Bruno’s project brief, although the present stage remains exploratory and not yet a finished final biological interpretation.

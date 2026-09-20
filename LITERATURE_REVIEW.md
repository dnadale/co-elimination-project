# Literature Review Plan

This reading list supports the master’s project on phenotype-agnostic coordinated gene co-elimination across placental mammals.

## Priority 1: Core project concept

### Szenker-Ravi et al. (2022)
**Discovery of a genetic module essential for assigning left-right asymmetry in humans and ancestral vertebrates**

Journal: Nature Genetics

Why it matters: This is the closest example of phenotype-agnostic gene-module discovery. It demonstrates how shared evolutionary loss profiles can group known and uncharacterized genes before the biological function is fully understood.

Questions to answer:
- What species and genomic data were used?
- How were shared gene-loss profiles defined?
- How were candidate genes grouped?
- Which parts of the computational workflow are reproducible?

### Szenker-Ravi et al. (2025)
**CIROZ is dispensable in ancestral vertebrates but essential for left-right patterning in humans**

Journal: American Journal of Human Genetics

Why it matters: This paper shows how a computationally discovered candidate can be followed up with clinical and experimental validation.

Questions to answer:
- What evidence connected CIROZ to the phenotype?
- Which validation steps are feasible for this master’s project?
- How were limitations of the original computational prediction addressed?

## Priority 2: Zoonomia and comparative genomics

### Zoonomia Consortium (2020)
**A comparative genomics multitool for scientific discovery and conservation**

Journal: Nature

Why it matters: This is the foundation for understanding the Zoonomia dataset, mammalian genome alignments, species coverage, and available comparative resources.

Questions to answer:
- Which genome assemblies and alignments are available?
- What species tree and annotations are provided?
- Which processed files can be used instead of regenerated?
- What storage and compute requirements should be expected on IBEX?

### Armstrong et al. (2020)
**The role of long-term conservation in the evolution of mammalian genes**

Journal: Nature

Why it matters: Provides useful examples of how broad mammalian comparative data can reveal conserved and rapidly evolving genomic regions.

## Priority 3: Gene-loss and orthology tools

### Kirilenko et al. (2023)
**Integrating gene annotation with orthology inference at scale: TOGA**

Journal: Genome Biology

Why it matters: TOGA may provide the gene-status information needed to distinguish intact genes, disrupted genes, and likely gene losses across genomes.

Questions to answer:
- What inputs does TOGA require?
- What outputs indicate a likely gene loss?
- Are relevant precomputed results available?
- How does TOGA handle incomplete assemblies and uncertain annotations?

### CESAR 2.0
Search terms: **CESAR 2.0 coding exon annotation genome alignments**

Why it matters: CESAR can help identify disrupted coding exons, frameshifts, premature stop codons, splice-site changes, and exon deletions in genome alignments.

Questions to answer:
- How does CESAR complement TOGA?
- Can it be run on the chosen pilot dataset within the project timeline?
- Are precomputed CESAR results available?

### Emms and Kelly (2019)
**OrthoFinder: phylogenetic orthology inference for comparative genomics**

Journal: Genome Biology

Why it matters: Provides an alternative if suitable precomputed ortholog assignments are unavailable.

Questions to answer:
- Can OrthoFinder handle the selected pilot species?
- What input sequence and annotation quality is required?
- Is it necessary if TOGA already supplies orthology relationships?

## Priority 4: Evolutionary conservation and statistics

### Pollard et al. (2010)
**Detection of nonneutral substitution rates on mammalian phylogenies**

Journal: Genome Research

Why it matters: Foundational paper for phyloP and the detection of conserved or rapidly evolving sequence regions across a phylogeny.

Project use: phyloP should be used as supporting evidence for evolutionary constraint or acceleration, not as the primary definition of gene loss.

### Hubisz et al. (2011)
**Mapping gene regulatory elements by phylogenetic footprinting**

Journal: Molecular Biology and Evolution

Why it matters: Provides background on the PHAST framework and phylogenetic footprinting methods related to phyloP.

### Pagel (1994)
**Detecting correlated evolution on phylogenies**

Journal: Proceedings of the Royal Society B

Why it matters: Foundational background for testing whether two traits or gene states evolve dependently while accounting for phylogenetic relationships.

## Paper-extraction template

For every paper, record:

- Full citation and DOI
- Dataset and number of species
- Required input files
- Software and versions
- Definition of gene loss
- Orthology method
- Statistical method
- Phylogenetic controls
- Available code
- Available processed data
- What can be reproduced on IBEX
- What could become a contribution of this project

## Recommended reading order

1. Szenker-Ravi et al. (2022)
2. Zoonomia Consortium (2020)
3. Kirilenko et al. (2023), TOGA
4. Pollard et al. (2010), phyloP
5. Szenker-Ravi et al. (2025)
6. CESAR 2.0 paper
7. Emms and Kelly (2019), OrthoFinder
8. Pagel (1994)
9. Hubisz et al. (2011)

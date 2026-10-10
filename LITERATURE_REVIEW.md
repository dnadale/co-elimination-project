# Literature Review Plan

This reading list supports the master’s project on phenotype-agnostic coordinated gene co-elimination across placental mammals.

The broader summer reading provides the evolutionary rationale for the project: gene duplication creates redundancy, pleiotropy and developmental importance impose constraint, and relaxed selection or regulatory rewiring can produce lineage-specific innovation. These topics should inform interpretation of the full feasible Zoonomia-scale analysis.

The official Zoonomia and TOGA metadata have now been inventoried. This confirms that the literature review can be connected to reusable public annotations and gene-status outputs rather than assuming that all computational preprocessing must be recreated.

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

Journal: Science, 380(6643), eabn3107
DOI: [10.1126/science.abn3107](https://doi.org/10.1126/science.abn3107)

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
- Can it be run for the full selected species set within the project timeline?
- Are precomputed CESAR results available?

### Emms and Kelly (2019)
**OrthoFinder: phylogenetic orthology inference for comparative genomics**

Journal: Genome Biology

Why it matters: Provides an alternative if suitable precomputed ortholog assignments are unavailable.

Questions to answer:
- Can OrthoFinder handle the selected full-scale species set?
- What input sequence and annotation quality is required?
- Is it necessary if TOGA already supplies orthology relationships?

## Supporting summer reading

The following topics are relevant for the literature review and discussion, but are not all required pipeline components:

- Cañestro et al. (2013): gene duplication, paralogs, orthologs, synteny, and the evolutionary fates of duplicate genes.
- Wang et al. (2020): transcriptome/translatome differences and evolutionary buffering.
- Mazin et al. (2021): alternative splicing, exon-level constraint, and organ development.
- Cooper (2024): limitations of single-gene explanations for complex evolutionary traits.
- Cardoso-Moreira and colleagues: developmental expression trajectories, pleiotropy, tissue-specific constraint, and evolutionary timing.
- Fukushima and Pollock (2020): organ-specific expression evolution and phylogenetic expression models.
- Reviews and case studies on convergent evolution, lineage-specific innovation, ageing, and antagonistic pleiotropy: useful conceptual context, but outside the core co-elimination analysis.

These readings support a network-level interpretation of co-elimination and caution against assuming that similar phenotypes always use identical genetic mechanisms.

## TACIT and regulatory follow-up

TACIT is relevant as a possible downstream method for predicting enhancer activity and associating regulatory evolution with an explicitly labeled phenotype. It should not replace the phenotype-agnostic retention/loss screen because TACIT requires phenotype-associated training or comparison data.

A sensible order is:

1. Discover candidate modules from gene retention, disruption, and loss profiles.
2. Overlay phyloP/phastCons and other constraint evidence.
3. Use TACIT for enhancer-focused follow-up if a suitable phenotype and training dataset are available.

## Positive-control strategy

The CIROP/DAND5/MMP21/PKD1L1/C1orf127 or CIROZ left-right patterning material can motivate a positive-control analysis. The pipeline should first test whether a known module or a subset of its evolutionary profile is recoverable. Any exact gene functions, clade-specific loss claims, disease percentages, and recent-paper claims must be verified against the primary sources before being reported as results.

## Verification requirement

The summer study guide is a useful synthesis but is not itself a primary source. Before citing exact numbers, software versions, species counts, gene functions, or recent publications, check the original papers, supplementary files, and official code or data repositories.

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

## Priority 5: AI-assisted genome mining and primary data resources

### Anthropic (2026)
**Claude discovers a novel enzyme system with CRISPR-like repeats**

Source: [Anthropic Science announcement](https://www.anthropic.com/news/claude-discovers-novel-enzyme-system) and its linked [technical report](https://www-cdn.anthropic.com/22573675ada52a8ca8a97a1a4b4326b2f208a071.pdf), published September 23, 2026.

Why it matters: This is a case study in agent-assisted genome mining. The report describes Claude agents surveying reverse-transcriptase families, prioritizing candidate systems, writing evidence summaries, and supporting expert review. ART means array-associated reverse transcriptases and is a biological system reported mainly in bacteriophages; it is not a released algorithm or plug-in for mammalian comparative genomics, and its biological function remains under investigation.

Project use: Adapt the general workflow only as an auditable candidate-triage layer after deterministic gene-status and phylogenetic analyses. AI-generated hypotheses must remain separate from measured data and require reproducible checks and expert review.

### Official Zoonomia and TOGA data resources

- [Zoonomia data portal](https://zoonomiaproject.org/the-data/): official access point for assemblies, the 240-species Cactus alignment, conservation scores, phylogeny, and TOGA annotations.
- [TOGA public archive](https://genome.senckenberg.de/download/TOGA/): precomputed reference/query annotations and gene-status summaries. Consult archive README files for schemas, classifications, and pairwise-comparison limitations.
- [TOGA source repository](https://github.com/hillerlab/TOGA2): software and implementation documentation; confirm the exact version associated with any output before citing or reproducing it.

These are primary data/software resources rather than substitutes for peer-reviewed papers. Record access dates, assembly identifiers, file checksums, and exact software versions when using them.

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
5. Anthropic (2026), ART genome-mining case study and technical report
6. Szenker-Ravi et al. (2025)
7. CESAR 2.0 paper
8. Emms and Kelly (2019), OrthoFinder
9. Pagel (1994)
10. Hubisz et al. (2011)

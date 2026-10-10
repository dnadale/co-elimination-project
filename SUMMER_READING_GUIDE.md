# Summer Reading Guide

This document records the summer reading assigned for the project and identifies how each topic informs the master’s research direction.

## Central synthesis

The reading supports a layered view of evolutionary change:

1. Gene duplication creates redundancy and opportunities for neofunctionalization or subfunctionalization.
2. Pleiotropy and developmental importance impose strong constraint on some genes and tissues.
3. Relaxed selection permits faster sequence, expression, and splicing evolution in other contexts.
4. Regulatory rewiring and lineage-specific innovations can alter how conserved biological systems operate.
5. Comparative genomics can detect these patterns across species and generate hypotheses for functional validation.

This supports the project’s central premise: coordinated gene loss can reveal biological modules, but gene loss must be interpreted alongside orthology, phylogeny, constraint, regulatory context, and data quality.

## Reading themes and project relevance

### Gene duplication and vertebrate diversification

Cañestro et al. (2013), *Impact of gene gains, losses and duplication modes on the origin and diversification of vertebrates*, provides background on two rounds of whole-genome duplication, paralogs, orthologs, synteny, and the possible fates of duplicated genes.

Project relevance:

- Helps distinguish orthology from paralogy before constructing a retention/loss matrix.
- Provides an evolutionary explanation for redundancy and later gene loss.
- Supports interpreting co-elimination as module-level evolution rather than isolated gene events.

### Transcriptome and translatome evolution

Wang et al. (2020), *Transcriptome and translatome co-evolution in mammals*, emphasizes that RNA abundance and translated protein output can evolve at different rates and that translational regulation can buffer sequence or transcript-level change.

Project relevance:

- Supports treating gene presence or coding integrity as only one layer of biological function.
- Provides context for why a retained gene may not have identical activity across species.
- Makes expression and translatome data useful downstream validation, not prerequisites for the initial screen.

### Alternative splicing and organ development

Mazin et al. (2021), *Alternative splicing during mammalian organ development*, highlights constitutive exons, cassette exons, developmental splicing changes, and tissue-specific evolutionary constraints.

Project relevance:

- Helps interpret exon-level disruption detected by CESAR or related tools.
- Suggests that loss of an exon or splice site should be evaluated in transcript and developmental context.
- Provides a possible extension after the core coding gene-loss screen.

### Complexity and network-level explanations

Cooper (2024), *The case against simplistic genetic explanations of evolution*, cautions against treating complex traits as the result of a single gene and emphasizes interconnected, small-effect, and developmental networks.

Project relevance:

- Strongly supports a module-based co-elimination design.
- Supports cautious interpretation: computational co-loss is evidence of a candidate relationship, not proof of a single-gene mechanism.
- Reinforces the need for phylogenetic controls and independent validation.

### Organ evolution, developmental constraint, and expression divergence

Cardoso-Moreira and colleagues, including the 2019 mammalian organ-development expression study and the later organ-evolution synthesis, provide context for pleiotropy, developmental timing, tissue-specific constraint, and the contrast between highly constrained organs and more permissive evolutionary settings.

Project relevance:

- Useful for downstream prioritization of candidate modules by tissue and developmental role.
- Supports interpreting expression divergence as complementary evidence rather than direct evidence of gene elimination.
- The 2019 multi-organ, multi-species dataset may be useful for validation if accessible and compatible with the full analysis.

### Preadaptive expression pathways

Fukushima and Pollock (2020), *Amalgamated cross-species transcriptomes reveal organ-specific propensity in gene expression evolution*, examines directional expression shifts, organ-to-organ expression patterns, duplication, and phylogenetic models such as Ornstein-Uhlenbeck processes.

Project relevance:

- Provides possible methods for a later expression-based extension.
- Supports the idea that evolutionary changes are constrained by ancestral architecture.
- Should not be part of the minimum three-to-four-month deliverable unless suitable processed data are already available.

### Zoonomia and evolutionary constraint

The Zoonomia literature provides the comparative mammal resource, multi-species alignments, phylogeny, and conservation framework for the project. The supplied reading emphasizes constrained coding and noncoding sequence, ultraconserved regions, accelerated regions, and unannotated constrained regions.

Project relevance:

- Zoonomia remains the primary comparative resource.
- phyloP and phastCons can provide supporting constraint features.
- Constraint scores should help interpret candidate genes or regulatory regions, not define gene loss by themselves.
- Exact species counts, assemblies, score definitions, and percentages must be checked against the original publication and data release.

Primary resources to consult while writing:

- [Zoonomia data portal](https://zoonomiaproject.org/the-data/) for assembly, alignment, conservation-score, phylogeny, and annotation access.
- [TOGA archive](https://genome.senckenberg.de/download/TOGA/) for per-assembly outputs and archive documentation. Its gene-status classifications are reference-versus-query comparisons, so missingness and reference-specific genes must not be interpreted as confirmed lineage losses without phylogenetic checks.
- Kirilenko et al. (2023), *Integrating gene annotation with orthology inference at scale*, Science 380(6643), eabn3107, [doi:10.1126/science.abn3107](https://doi.org/10.1126/science.abn3107), for TOGA methods and evaluation.

### AI-assisted genome mining: ART case study

Anthropic's September 2026 announcement, ["Claude discovers a novel enzyme system with CRISPR-like repeats"](https://www.anthropic.com/news/claude-discovers-novel-enzyme-system), and its linked [technical report](https://www-cdn.anthropic.com/22573675ada52a8ca8a97a1a4b4326b2f208a071.pdf) describe an AI-assisted search that identified candidate array-associated reverse transcriptases (ARTs). ART is a newly reported biological system, not an algorithm or released software package; its function is still being investigated.

Project relevance:

- Offers a methodological example of broad sequence screening, candidate prioritization, evidence summaries, expert review, and experimental follow-up.
- Supports using AI to help triage ranked co-elimination candidates and organize evidence, after the matrix and statistical tests are generated reproducibly.
- Does not establish that AI review improves accuracy by itself. Candidate status must come from the underlying data and validated analyses; AI suggestions need independent checking and an audit trail.

### TACIT and phenotype-linked enhancer inference

TACIT (Transcriptional Activity Conservation and Inference Toolkit) uses sequence-based prediction of open-chromatin or enhancer activity across species and phylogenetically informed association with phenotypes.

Project relevance:

- TACIT is valuable for a later regulatory follow-up, especially for left-right patterning or other explicitly labeled phenotypes.
- It is not a direct replacement for the phenotype-agnostic co-elimination screen because it requires phenotype labels or phenotype-associated comparisons.
- A practical order is: discover candidate modules by gene retention/loss, overlay constraint, then use TACIT to investigate candidate enhancers.
- Phylogenetic regression and permutation testing are important safeguards against false associations.

### Left-right asymmetry and the CIROP/DAND5 example

The CIROP, DAND5, MMP21, PKD1L1, and C1orf127/CIROZ material provides a motivating positive-control case for a flow-dependent left-right patterning module. The proposed biological chain is broadly represented as:

`motile cilia -> leftward flow -> flow sensing -> CIROP/DAND5 regulation -> Nodal -> Lefty/Pitx2 -> organ asymmetry`

Project relevance:

- The module can serve as a positive-control case if the relevant ortholog and gene-status data are available.
- A positive-control test can assess whether known shared evolutionary profiles are recoverable before searching for novel modules.
- Claims about exact gene functions, clade-specific losses, disease percentages, and recent papers must be verified against primary sources before being used as conclusions.

### Convergent evolution and lineage-specific innovation

The pregnancy, GSDMC, KEAP1 recoding, and related readings illustrate that similar phenotypes can arise through different molecular routes and that lineage-specific innovations can involve coding, regulatory, or translational changes.

Project relevance:

- Prevents assuming that phenotypic similarity implies identical genetic mechanisms.
- Supports searching for repeated module-level patterns while allowing lineage-specific rewiring.
- Provides conceptual background for classifying retained, disrupted, rewired, and lineage-specific genes.
- These topics are contextual rather than required computational components of the core co-elimination analysis.

### Ageing, pleiotropy, and evolutionary trade-offs

The ageing material emphasizes mutation accumulation, antagonistic pleiotropy, selection shadows, conserved pathways, and lineage-specific longevity mechanisms.

Project relevance:

- Useful background for understanding why gene effects depend on life stage, tissue, and ecological context.
- Supports caution when interpreting relaxed constraint or gene loss as simple biological improvement.
- The ageing and archaic-aware CRISPR topics are outside the current co-elimination thesis scope unless the project is explicitly expanded.

## Methods to prioritize

### Core methods

- Zoonomia assemblies, alignments, annotations, and phylogeny
- Orthology and paralogy quality control
- TOGA and/or CESAR gene-status assessment
- Presence, disruption, and loss matrix construction
- Phylogenetically informed co-elimination scoring
- Null models or phylogenetic permutation tests
- Basic candidate-module clustering and visualization

### Supporting methods

- phyloP and phastCons constraint overlays
- Ancestral-state reconstruction
- Expression or translatome comparison where processed data are available
- Pathway and protein-interaction annotation
- Literature and clinical-genetics validation

### Future extensions

- TACIT enhancer prediction and phenotype association
- Cross-species ATAC-seq or H3K27ac integration
- Alternative-splicing and exon-level evolutionary analysis
- Deep-learning regulatory models such as Enformer or Basenji
- Experimental validation through CRISPR, reporter assays, or model organisms

## Reading-to-thesis decision

The summer reading is worth incorporating, but not every method belongs in the first deliverable. The thesis should use the reading to justify the biological model and prioritize a reproducible computational core. The minimum defensible project is:

`Zoonomia -> orthology/gene status -> retention/loss matrix -> phylogenetic co-elimination -> constraint and functional interpretation`

TACIT, expression, splicing, ageing, and experimental validation should be presented as extensions unless the lab already has compatible data and pipelines.

## Verification note

This guide is a synthesis of assigned reading notes, not a substitute for the original publications. Before citing numerical results, software versions, gene functions, dates, or claims from recent papers, verify them against the primary paper, supplementary material, and official data or code repository.

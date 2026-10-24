# Study guide: comparative genomics and coordinated gene loss

## Purpose

This guide is a learning path for the project on phenotype-agnostic discovery of
coordinated gene co-elimination across placental mammals. It connects the
biological ideas, comparative-genomics data, statistical challenges, and
current repository outputs. It is intended to help you read the literature
critically and explain what the analysis can—and cannot—conclude.

The project’s broad question is:

> Do some genes show unusually coordinated histories of retention and loss
> across mammals, beyond what shared ancestry and technical artifacts would
> lead us to expect?

The goal is to generate and prioritize biological hypotheses. A pattern of
co-loss is not, by itself, proof that two genes interact, belong to the same
pathway, or were eliminated for the same adaptive reason.

## How to use this guide

Work through the sections in order. For each paper, use the **paper notes
template** near the end and record what the authors actually tested, not just
the headline result. Treat the questions and exercises as checkpoints: if you
cannot answer them from a source, mark the point as unresolved rather than
guessing.

## Learning objectives

By the end, you should be able to:

1. Explain why gene loss can be informative about changing biological systems.
2. Distinguish orthology, gene absence, gene disruption, and uncertainty.
3. Describe how a gene-by-species status matrix is made and quality-checked.
4. Explain why species are not statistically independent observations.
5. Compare a raw co-loss screen with phylogenetically informed tests.
6. Interpret candidate rankings without overstating exploratory or
   post-selection evidence.
7. Propose independent evidence and validation for a candidate module.
8. Reproduce the project’s data lineage from source metadata to reported
   candidate.

## 1. Biological foundations

### 1.1 Gene gain, duplication, and loss

Genes are gained, duplicated, modified, and lost over evolutionary time.
Duplication can create redundancy: after one copy changes or disappears, the
other may retain enough function for the organism. Conversely, a gene with
important dosage-sensitive or broadly pleiotropic functions may be difficult
to lose. These are tendencies, not rules. A gene’s fate depends on its
biological context, the organism, and the particular evolutionary changes.

Gene loss can therefore raise useful questions. Did a function become
unnecessary, was it replaced by another mechanism, or did another part of the
system change first? Comparative evidence can identify lineages and genes for
which those questions are worth pursuing. It usually cannot answer the
mechanistic question alone.

### 1.2 Why might genes be lost together?

Several explanations can produce similar patterns:

- **Functional dependence:** two genes contribute to the same process, and
  loss of one makes the other less useful.
- **Pathway or module change:** an entire biological process changes, so
  multiple components become dispensable.
- **Shared ancestry:** one loss in a common ancestor is inherited by many
  descendant species. This can look like many independent co-losses if
  ancestry is ignored.
- **Physical linkage:** nearby genes may be removed or misassembled together.
- **Gene-family complexity:** paralogs, repetitive sequence, or annotation
  ambiguity can create correlated calls.
- **Technical artifact:** assembly gaps, low sequence quality, or reference
  bias can make several genes appear absent or disrupted.
- **Chance and selection effects:** scanning millions of gene pairs will
  produce striking-looking patterns even when no meaningful relationship
  exists.

The analysis must compare these explanations. “Co-elimination” is a useful
screening term, not a causal conclusion.

### 1.3 Phenotype-first and phenotype-agnostic reasoning

A phenotype-first analysis starts with a trait and asks which genes may explain
it. A phenotype-agnostic screen instead looks for unusual similarities in
evolutionary histories before assigning a trait. The latter can suggest
unexpected connections, but it also increases the need for careful statistical
controls and independent follow-up: candidates were selected because they
looked unusual in the same data used to evaluate them.

## 2. Comparative-genomics concepts

### 2.1 Orthologs and paralogs

**Orthologs** are genes in different species descended from a gene in their
last common ancestor. **Paralogs** arise through gene duplication. A one-to-one
ortholog is comparatively straightforward to track, but gene families often
have duplications, losses, and ambiguous mappings. A missing ortholog call can
mean true loss, an unrecognized homolog, a mapping limitation, or missing
sequence.

Before comparing gene status across species, ask:

- What is the reference gene set?
- How were query genes matched to reference genes?
- How are duplications and one-to-many mappings represented?
- Does the data distinguish loss from an unresolved or unannotated case?

### 2.2 Gene inactivation and confidence

Loss of function can be supported by coding-sequence disruptions such as
frameshifts, premature stop codons, splice-site disruption, or deleted exons.
An intact-looking sequence is not necessarily proof of function, and an
apparent absence is not necessarily a true evolutionary loss. The strength of
the conclusion depends on the evidence, assembly quality, annotation, and
classification scheme.

Keep at least three conceptual states separate:

1. **Retained/intact evidence**
2. **Disrupted or likely lost evidence**
3. **Unknown, uncertain, or not comparable**

Do not silently convert unknown or uncertain states to loss or retention.
Check the relevant tool’s documentation and output schema for the exact status
definitions before analysis.

### 2.3 Reference-versus-query data and ancestral inference

TOGA results in this project are derived from reference/query comparisons.
A query marked as lost relative to a reference is not automatically proof that
the gene was lost on that query’s terminal branch. The reference may itself
differ from the ancestral state, and calls can be affected by assembly or
annotation quality. Inferring when a loss occurred requires comparisons among
species and an explicit evolutionary model on a species tree.

This distinction matters throughout the project: observed pairwise status
calls are input data; inferred ancestral or independent loss events are
analysis results and require additional assumptions.

### 2.4 Zoonomia and TOGA in this project

The repository documents Zoonomia as the source of mammalian species,
assemblies, alignments, and phylogenetic resources, and TOGA as a source of
precomputed orthology and gene-status outputs. Prefer validated, documented
processed outputs when suitable, while recording assembly identifiers,
source files, checksums, and software versions.

The current documented matrix contains **162 species and 19,464 genes**, with
**3,153,168 cells and no missing cells in the matrix**. This describes the
current project artifact, not all possible Zoonomia species or a claim that
every gene call is biologically certain. See
[`QUALITY_CONTROL_SUMMARY.md`](QUALITY_CONTROL_SUMMARY.md) and
[`analysis/gene_status_matrix.qc.json`](analysis/gene_status_matrix.qc.json).

## 3. From source data to a co-loss candidate

The project’s conceptual pipeline is:

1. **Select and document species and assemblies.** Prevent mixing assembly
   versions and record exclusions.
2. **Collect comparable gene-status outputs.** Verify file identity, format,
   and status definitions.
3. **Build a gene-by-species matrix.** Preserve uncertainty and comparability
   information rather than treating every cell as equally reliable.
4. **Screen gene pairs or modules.** Identify patterns that warrant further
   investigation.
5. **Account for phylogeny and selection.** Ask whether the pattern remains
   unusual after shared ancestry and the scale of the search are considered.
6. **Inspect technical and biological context.** Review assembly quality,
   gene family, genomic proximity, annotations, and independent evidence.
7. **Report a ranked hypothesis with limitations.** Keep measured data,
   model-based inference, and interpretation distinct.

### 3.1 What a status matrix represents

Rows are genes (or defined ortholog groups); columns are species; cells hold
the evidence category for that gene in that species. A matrix is a convenient
summary, but it hides important details unless it is accompanied by
provenance, confidence, and comparability metadata. Always retain links back
to the source calls.

The project’s matrix QC checks include dimensions, missingness, and records of
the source files and their SHA256 checksums. These checks support traceability
and data integrity; they do not establish that all biological calls are
correct.

### 3.2 Raw pairwise co-loss screens

A simple screen compares the number of species in which both genes are called
lost with a baseline expectation. The project’s candidate tables report
quantities such as shared-loss count, individual loss counts, an expected
shared-loss count under a naive calculation, excess shared losses, and
enrichment.

These are ranking or screening statistics. A naive expectation that treats
species as independent can be misleading because mammals share evolutionary
history. It can also be affected by:

- unequal status coverage between genes;
- uncertainty and non-comparable calls;
- clade-specific patterns;
- a high background loss rate for one or both genes;
- linked genes, paralogs, or repetitive families; and
- searching millions of pairs and reporting only the strongest.

A large enrichment or shared-loss count is therefore a reason to investigate,
not a standalone significance claim.

### 3.3 Multiple testing and post-selection inference

When millions of pairs are screened, some will rank highly by chance or due
to technical structure. Multiple-testing correction can address a defined
family of tests under its assumptions, but it does not repair confounding,
poor calls, or a misspecified null model.

If the same dataset is used first to choose top candidates and then to test
those candidates, the second-stage p-values are **post-selection**. They can
look more convincing than evidence from a candidate set specified in advance.
Be explicit about this limitation. Stronger confirmation can come from a
pre-specified test, independent taxa/data, robust sensitivity analyses, and
orthogonal biological evidence.

### 3.4 Phylogenetic dependence

Closely related species tend to share traits because of common ancestry. For
example, one ancestral loss inherited by many descendants should not be
counted as dozens of independent evolutionary events. Phylogenetic methods
model character changes on a tree or otherwise account for shared ancestry.

Pagel’s correlated-evolution framework is a useful foundation for discrete
traits: it compares models in which two characters evolve independently and
dependently along a phylogeny. Interpretation still depends on the tree, state
coding, model fit, and data quality. A phylogenetic test is not automatically
immune to post-selection, systematic call error, or an incorrect species tree.

### 3.5 Positive controls and negative controls

A **positive control** asks whether the analysis can recover a known or
expected signal. In this project, the left-right patterning module
(MMP21, DAND5, PKD1L1, CIROP, and CIROZ) is a motivating example. Because the
genes are chosen from prior knowledge, recovery validates aspects of the
pipeline but is not a new discovery.

Negative or null controls can help reveal artifacts. Examples include
comparing candidates with genes matched for loss frequency and call
completeness, testing whether physical linkage explains a pattern, and
permuting data in a way that preserves relevant structure. The null must
reflect the question being asked; arbitrary shuffling can destroy the very
phylogenetic or technical structure that should be controlled.

## 4. Functional interpretation and validation

After statistical screening, interpret candidates using evidence that is as
independent as possible from the ranking:

- pathway and ontology annotations;
- protein interactions and known complex membership;
- expression in relevant tissues or developmental stages;
- human genetics and clinical phenotypes;
- model-organism experiments;
- gene-family and genomic-neighborhood context;
- independent sequence or assembly checks; and
- conservation or acceleration evidence.

Conservation scores such as phyloP can support interpretation of sequence
constraint or acceleration. They do not, on their own, establish gene loss or
functional elimination. Similarly, expression correlation or interaction
database membership does not prove that two genes were co-eliminated for a
shared reason.

The strongest next step for a computational candidate is a test that can
distinguish competing explanations. Depending on the candidate, this might be
manual inspection of underlying sequence evidence, an independent annotation
resource, a pre-specified comparison in additional taxa, or experimental
perturbation. Experimental validation is a future follow-up, not a prerequisite
for honestly reporting a computational screen; the claims must match the
evidence.

## 5. Current project status and interpretation

According to the repository’s October 2026 status and QC documents:

- The matrix and candidate screening artifacts have been generated and
  checked for documented dimensions and provenance.
- A broad exploratory screen scored approximately 25.9 million gene pairs.
- The left-right module has been examined as a positive control.
- A shortlist of unlinked pairs is available for biological and technical
  review.
- Phylogenetic validation has been run on selected candidates, but these
  results are exploratory and post-selection.

The status report notes that many top raw candidates are histone-histone or
physically linked pairs, while unlinked candidates require further review.
It also records uncertainty around some gene calls and the possibility that
some apparent losses reflect technical artifacts. Read
[`writeups/STATUS_REPORT_2026-10-12.md`](writeups/STATUS_REPORT_2026-10-12.md)
and [`QUALITY_CONTROL_SUMMARY.md`](QUALITY_CONTROL_SUMMARY.md) before
describing project findings.

The current evidence supports a **reproducible discovery-stage screen**. It
does not yet establish novel functional modules or causal mechanisms.
Candidate-specific claims should be revisited as the analyses and QC advance.

## 6. Guided reading sequence

Use the repository’s [`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) as the
starting bibliography. The reading list below follows a learning order. Verify
full bibliographic details, claims, species counts, and methods against the
primary paper and its supplements before citing them.

### A. Learn the motivating discovery paradigm

**Szenker-Ravi et al. (2021), “Discovery of a genetic module essential for
assigning left-right asymmetry in humans and ancestral vertebrates,” listed in
the literature review as published in *Nature Genetics*.**

Focus on how evolutionary profiles were used to identify a module, what data
and candidate-selection steps were involved, and what was known before
follow-up. Separate the discovery logic from later evidence about biological
function.

**Questions:** What counted as a shared evolutionary profile? How was the
candidate set constructed? Which observations were independent validation?
Which parts of this project’s design are similar, and which are not?

### B. Understand follow-up from computation to biology

**Szenker-Ravi et al. (2025), “CIROZ is dispensable in ancestral vertebrates
but essential for left-right patterning in humans,” listed in the literature
review as published in AJHG.**

Trace the evidence chain from computational candidate to biological claim.
Identify which experiments or clinical evidence address function and which
claims still depend on evolutionary inference.

**Questions:** What was the original prediction? What evidence tested it?
Which experiments would be feasible here, and what would they establish?

### C. Understand the comparative resource

**Zoonomia Consortium (2020), “A comparative genomics multitool for scientific
discovery and conservation,” listed in the literature review as published in
*Nature*.**

Read for the resource design: species, assemblies, alignments, phylogeny, and
the distinctions between raw and processed data. Compare the paper with the
actual metadata and files used in this repository.

**Questions:** Which species and assemblies are represented in this project?
What did the project actually download or reuse? Which resources remain
unavailable locally?

### D. Learn gene-status inference and orthology

**Kirilenko et al. (2023), “Integrating gene annotation with orthology
inference at scale: TOGA,” *Science* 380(6643), eabn3107,
DOI: 10.1126/science.abn3107.**

Study how orthology and gene annotation are combined, what output categories
mean, and how alignment or assembly limitations affect calls. Read the
appropriate TOGA documentation alongside the paper; do not infer status
semantics from a filename alone.

**CESAR 2.0** is relevant for coding-exon annotation in genome alignments.
Compare its role with TOGA and determine whether the current project uses
CESAR outputs or only treats them as a possible extension.

**Emms and Kelly (2019), OrthoFinder**, is an alternative orthology resource
to understand if precomputed mappings are unavailable or unsuitable. Do not
assume a second orthology pipeline is automatically better; first define the
mapping question and compare outputs.

**Questions:** How are orthology and gene loss distinguished? Which files
contain status evidence, mutations, or annotations? How does the pipeline
represent uncertainty?

### E. Learn evolutionary constraint and correlated evolution

**Pollard et al. (2010), “Detection of nonneutral substitution rates on
mammalian phylogenies,” listed in the literature review as published in
*Genome Research*.**

Use this as background for phyloP and conservation/acceleration. Keep
nucleotide-level evolutionary-rate evidence separate from gene-inactivation
calls.

**Pagel (1994), “Detecting correlated evolution on phylogenies: a general
method for the comparative analysis of discrete characters,” listed in the
literature review as published in *Proceedings of the Royal Society B*.**

Focus on the independent versus dependent character models and why a
phylogeny-aware comparison differs from counting matching species states.

**Questions:** What is the unit of evolutionary evidence? What assumptions
does the model make? How would candidate selection affect interpretation of
a test result?

### F. Broaden context only after the core methods

The literature review lists Cañestro et al. (2013), Wang et al. (2020), and
Mazin et al. (2021) as supporting reading on gene gain/loss and duplication,
mammalian transcriptome/translatome evolution, and alternative splicing.
Use these to enrich biological interpretation only after confirming the exact
paper and relevance. These topics are context, not substitutes for testing
gene-status calls and phylogenetic co-loss.

TACIT is a possible phenotype-linked regulatory follow-up when a suitable
phenotype and training data exist. It is not a replacement for this
phenotype-agnostic gene-status screen.

## 7. Practical exercises using this repository

Do these exercises without changing the source data. Keep notes on the exact
file, column, and assumption used.

### Exercise 1: Reconstruct data provenance

1. Read the QC summary and matrix QC JSON.
2. Record the species count, gene count, missing-cell count, and source-file
   checksum information.
3. Trace one species from the matrix QC report back to its TOGA summary and
   assembly identifier.
4. Explain what the checksums verify and what they cannot verify.

### Exercise 2: Inspect one candidate carefully

1. Choose a candidate from
   `results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv`.
2. Record its shared-loss count, individual loss counts, comparable species,
   and ranking statistic.
3. Check whether either gene has uncertain or incomplete calls and inspect
   gene symbols/biotypes.
4. Ask whether the genes are physically linked, family-related, or affected by
   the same assembly-quality issue.
5. Write two competing explanations for the pattern before proposing a
   biological interpretation.

### Exercise 3: Contrast raw and phylogenetic evidence

1. Pick a candidate with both a screening result and a phylogenetic-validation
   result.
2. Identify the species tree and coding of the gene states used.
3. Explain how inherited clade-level states affect the count of independent
   evidence.
4. State whether the phylogenetic result is confirmatory or post-selection.
5. List what additional data or analysis would make the inference more
   independent.

### Exercise 4: Evaluate the positive control

Read the left-right module status report and results tables. Separate:

- prior biological knowledge used to choose the control;
- observed status patterns in the project data;
- statistical tests applied;
- phylogenetic interpretation; and
- experimental or clinical evidence from the literature.

Explain why successful recovery is useful pipeline evidence but not a novel
discovery by itself.

### Exercise 5: Design a more defensible follow-up

Write a short proposal for one candidate that includes:

- a specific, falsifiable biological hypothesis;
- the main technical alternative explanation;
- a test that distinguishes these explanations;
- a plan to account for phylogeny and multiple testing;
- an independent data source or validation step; and
- a cautious statement of what a positive result would mean.

## 8. Paper notes template

Copy this checklist into your own notes for each source:

- Full citation and DOI or stable source link:
- Primary paper, review, or data/software resource:
- Research question:
- Species, assemblies, and sample size:
- Input data and provenance:
- Orthology method:
- Definition of intact, disrupted, lost, and uncertain:
- Statistical model and null hypothesis:
- How phylogeny was handled:
- Multiple-testing or candidate-selection procedure:
- Main result, in the authors’ own scope:
- Independent validation:
- Limitations acknowledged by the authors:
- Code, data, and supplement locations:
- Relevance to this project:
- Claim that must be checked before citing:

## 9. Glossary

- **Assembly:** A reconstructed genome sequence for an individual species or
  sample; assembly quality affects the ability to detect genes.
- **Co-loss / co-elimination:** A pattern in which two or more genes have
  evidence of loss in overlapping evolutionary lineages. The pattern alone
  does not establish causation or functional dependence.
- **Comparability:** Whether a gene’s status can be meaningfully evaluated in
  a species, given coverage, orthology, and data quality.
- **Gene-status matrix:** A table of gene states across species, ideally
  accompanied by confidence and provenance.
- **Homology:** Evolutionary relatedness between sequences; it does not by
  itself specify orthology or equivalent function.
- **Ortholog:** A homologous gene relationship arising through speciation.
- **Paralog:** A homologous gene relationship arising through gene
  duplication.
- **Phylogenetic non-independence:** Dependence among species observations
  caused by shared ancestry.
- **Pseudogenization:** Accumulation of changes that disrupt a gene’s
  protein-coding function; evidence and confidence vary.
- **Post-selection inference:** Statistical evaluation performed after
  selecting candidates based on the same or related data.
- **Positive control:** A known or expected case used to check whether a
  method can recover a signal; it is not an independent discovery.
- **TOGA:** Tool for alignment-based orthology and gene annotation; consult
  the paper and output documentation for exact status interpretation.
- **Zoonomia:** A comparative genomics resource for mammals, including
  genomes and comparative data described in the project literature.

## 10. Reporting principles

When presenting results:

1. State the actual species and gene counts for the analyzed dataset, not the
   resource’s headline scope.
2. Distinguish source calls from inferred ancestral events.
3. Report uncertain and non-comparable calls explicitly.
4. Give the null model, candidate-selection process, and multiple-testing
   scope.
5. Identify physical linkage, gene-family complexity, and assembly quality as
   potential alternatives.
6. Label exploratory and post-selection tests as such.
7. Separate statistical association, functional interpretation, and causal
   claims.
8. Include source versions, assembly identifiers, checksums, and reproducible
   analysis details where available.
9. Verify exact literature claims against primary sources before citing.

## Repository reading map

- Project question and intended scope: [`PROJECT_GOAL.md`](PROJECT_GOAL.md)
- Conceptual framework and methods: [`RESEARCH_FOCUS.md`](RESEARCH_FOCUS.md)
- Reading bibliography and data links: [`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md)
- Verified matrix and candidate counts/caveats:
  [`QUALITY_CONTROL_SUMMARY.md`](QUALITY_CONTROL_SUMMARY.md)
- Current interpretation and next steps:
  [`writeups/STATUS_REPORT_2026-10-12.md`](writeups/STATUS_REPORT_2026-10-12.md)
- Matrix and QC artifact: `analysis/gene_status_matrix.tsv` and
  `analysis/gene_status_matrix.qc.json`
- Ranked screening outputs: `results/agentic_triage/`
- Phylogenetic validation outputs: `results/phylogenetic_validation/`

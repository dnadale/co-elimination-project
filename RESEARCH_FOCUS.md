# Coordinated Gene Co-Elimination as a Phenotype-Agnostic Paradigm for Pathway Discovery

## Project concept

Traditional comparative genomics often follows a phenotype-first approach: researchers begin with a known trait, such as loss of vision or adaptation to an aquatic environment, and search for genes lost in the species showing that trait. This project takes a genotype-first, phenotype-agnostic approach.

Because complex biological processes are governed by interacting protein networks, the evolutionary obsolescence of a trait may be accompanied by coordinated decay of an entire gene module. By clustering genes according to their evolutionary retention and loss profiles, we aim to discover functionally linked genetic networks without first specifying the phenotype they control. These candidate modules can then be reverse-engineered using clinical, expression, conservation, and model-organism data.

## Current progress

Some of the project's initial deliverables are complete, including a quality-controlled gene-status matrix and an initial shared-loss screen. Candidate genes and pairs remain preliminary, hypothesis-generating candidates, not validated discoveries.

Gene-status interpretation follows the classification framework from TOGA:

- **Intact (I):** classified as intact.
- **Partially intact (PI):** classified as partially intact.
- **Uncertain loss (UL):** possible loss with uncertain status.
- **Missing (M):** classified as missing.
- **Partially missing (PM):** missing classification where the overlapping chain covers less than 50% of the gene's CDS.
- **Lost (L):** classified as lost.
- **Paralogous projection (PG):** classified as a paralogous projection.

The validated matrix contains **162 species** and **19,464 genes** (**3,153,168 gene–species cells**), with **0 missing cells** (every position is populated). Across those entries, the TOGA status distribution is **74.73% intact**, **11.05% partially intact**, **4.54% uncertain loss**, **4.53% lost**, **2.94% missing**, **1.07% partially missing**, and **1.14% paralogous projection**. The QC report records SHA256 hashes for the source TOGA files for reproducibility. Zero missing cells means every matrix position has a recorded classification; M and PM remain explicit TOGA classifications.

The initial screen scored **25,945,206 gene pairs** for shared loss and retained the **top 100 candidate pairs**. This shortlist includes **53 histone–histone pairs**, with the remaining pairs comprising **28 physically linked** and **19 unlinked pairs** for follow-up. These counts describe a screening shortlist, not experimentally or independently validated findings.

## Primary research focus

The initial focus is the Zoonomia comparative genomics resource, which provides whole-genome comparisons across approximately 240 placental mammal species spanning roughly 100 million years of evolution. The exact species count and available annotations will be confirmed during data acquisition.

## Methodology

1. **Large-scale comparative genomics**
   - Use Zoonomia whole-genome alignments and associated annotations where available.
   - Establish a quality-controlled set of placental mammal genomes and a species phylogeny.

2. **Systematic gene-inactivation screening**
   - Use tools such as TOGA and CESAR to identify gene losses and disrupted coding sequences.
   - Track premature stop codons, frameshifting insertions or deletions, disrupted splice sites, missing exons, and other evidence of pseudogenization.
   - Distinguish confirmed losses from missing data or low-confidence predictions.

3. **Presence/absence and evolutionary-profile construction**
   - Build a genome-wide matrix describing gene retention, functional disruption, and loss across species or ancestral lineages.
   - Incorporate phylogenetic information so that shared ancestry is considered when evaluating correlated losses.

4. **Unsupervised clustering and statistical testing**
   - Identify gene pairs or modules with unusually correlated evolutionary histories.
   - Compare observed co-elimination scores with appropriate null models and phylogenetically informed controls.
   - Rank candidate modules by replication, confidence, and evidence quality.

5. **Functional contextualization**
   - Examine pathway annotations, protein interactions, tissue-specific expression, and conservation or constraint patterns.
   - Use phyloP scores as supporting evidence for evolutionary constraint or accelerated evolution, not as the primary definition of gene loss.

6. **Optional regulatory extension**
   - If a candidate phenotype and suitable training data are available, use TACIT to predict conserved or lineage-specific enhancer activity.
   - Treat TACIT as a phenotype-linked follow-up method rather than a replacement for the phenotype-agnostic gene-loss screen.
   - Use phylogenetic regression and permutation testing to reduce false associations.

7. **Reverse functional validation**
   - Compare retained human orthologs with clinical genetics databases.
   - Where appropriate, evaluate candidate genes using model organisms or existing experimental evidence.

## Example motivating case

The left-right asymmetry module involving MMP21, DAND5, PKD1L1, CIROP, and CIROZ illustrates the type of discovery this framework is intended to support. A shared, unusual evolutionary footprint can group known and previously uncharacterized genes before the relevant phenotype is fully understood. Subsequent clinical genetics and experimental work can then connect the computationally discovered module to left-right patterning and heterotaxy.

This example is motivating context rather than a result of the current project.

It can also serve as a positive-control case for the pilot: first test whether a known left-right patterning module can be recovered from the available orthology and gene-status data, then search for previously uncharacterized modules.

## Role of phyloP

phyloP scores are relevant as an auxiliary feature. They can help quantify conservation or lineage-specific acceleration at nucleotide sites and can support interpretation of candidate genes. They should not replace explicit gene-inactivation evidence from tools such as TOGA or CESAR, because conservation scores alone do not reliably establish that a gene has been eliminated.

## Main deliverables

- A documented Zoonomia species and genome inventory.
- Quality-controlled TOGA/CESAR gene-status data where computational access permits.
- A phylogenetically aware gene retention and loss matrix.
- Ranked co-eliminated gene modules with statistical support.
- Functional annotations and conservation summaries for candidate modules.
- Reproducible notebooks, scripts, parameters, logs, and final figures.

## Confirmed public data resources

The official Zoonomia and TOGA resources have now been checked. Zoonomia provides the 240-species Cactus alignment, a mammalian phylogeny, conservation scores, genome assemblies, and associated comparative data. The TOGA archive provides reusable annotations and gene-status outputs for reference/query assembly comparisons, including `loss_summ_data.tsv.gz`, `orthologsClassification.tsv.gz`, `geneInactivatingMutations.tsv.gz`, BED12/GTF annotations, and codon or protein alignments.

The project will use these precomputed resources where possible rather than rebuilding the entire annotation pipeline. The exact set of species and assembly directories still needs to be selected, and the pairwise nature of TOGA calls must be accounted for in ancestral-state and co-elimination analyses.

## Important feasibility note

The full Zoonomia-scale analysis will require substantial storage, compute time, software setup, and careful validation. The project should begin by downloading and validating a manageable subset of the existing TOGA outputs before scaling to the complete resource.

For a three-to-four-month master's project, the core workflow should remain:

`Zoonomia -> orthology/gene status -> retention/loss matrix -> phylogenetic co-elimination -> constraint and functional interpretation`

Expression, translatome, alternative-splicing, TACIT, ageing, and experimental-validation analyses should be treated as supporting context or future work unless the lab already has compatible processed data and pipelines.

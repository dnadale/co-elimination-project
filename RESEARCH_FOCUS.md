# Coordinated Gene Co-Elimination as a Phenotype-Agnostic Paradigm for Pathway Discovery

## Project concept

Traditional comparative genomics often follows a phenotype-first approach: researchers begin with a known trait, such as loss of vision or adaptation to an aquatic environment, and search for genes lost in the species showing that trait. This project takes a genotype-first, phenotype-agnostic approach.

Because complex biological processes are governed by interacting protein networks, the evolutionary obsolescence of a trait may be accompanied by coordinated decay of an entire gene module. By clustering genes according to their evolutionary retention and loss profiles, we aim to discover functionally linked genetic networks without first specifying the phenotype they control. These candidate modules can then be reverse-engineered using clinical, expression, conservation, and model-organism data.

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

6. **Reverse functional validation**
   - Compare retained human orthologs with clinical genetics databases.
   - Where appropriate, evaluate candidate genes using model organisms or existing experimental evidence.

## Example motivating case

The left-right asymmetry module involving MMP21, DAND5, PKD1L1, CIROP, and CIROZ illustrates the type of discovery this framework is intended to support. A shared, unusual evolutionary footprint can group known and previously uncharacterized genes before the relevant phenotype is fully understood. Subsequent clinical genetics and experimental work can then connect the computationally discovered module to left-right patterning and heterotaxy.

This example is motivating context rather than a result of the current project.

## Role of phyloP

phyloP scores are relevant as an auxiliary feature. They can help quantify conservation or lineage-specific acceleration at nucleotide sites and can support interpretation of candidate genes. They should not replace explicit gene-inactivation evidence from tools such as TOGA or CESAR, because conservation scores alone do not reliably establish that a gene has been eliminated.

## Main deliverables

- A documented Zoonomia species and genome inventory.
- Quality-controlled TOGA/CESAR gene-status data where computational access permits.
- A phylogenetically aware gene retention and loss matrix.
- Ranked co-eliminated gene modules with statistical support.
- Functional annotations and conservation summaries for candidate modules.
- Reproducible notebooks, scripts, parameters, logs, and final figures.

## Important feasibility note

The full Zoonomia-scale analysis will require substantial storage, compute time, software setup, and careful validation. The project should begin with a small pilot subset of mammals and a limited set of genes or annotations before scaling to the complete resource.

If you use these data, please cite our manuscript:
"Integrating gene annotation with orthology inference at scale"
Bogdan M. Kirilenko, Chetan Munegowda, Ekaterina Osipova, David Jebb, Virag Sharma, Moritz Blumer, Ariadna Morales,
Alexis-Walid Ahmed, Dimitrios-Georgios Kontopoulos, Leon Hilgers, Zoonomia Consortium, and Michael Hiller
Science, 380(6643), eabn3107, 2023



The subdirectories contain data that we generated with
- human (hg38) assembly
- mouse (mm10 assembly
- chicken (galGal6) assembly
and other species as reference.

[[
The directory MammalianDNAZooAssemblies/ contains the 2bit files of all mammalian assemblies we obtained from the DNA zoo consortium.
The BatAnnotationEVM/ directory contains the input evidences for the bat genome annotation test (Figure 3C).
]]

For human, mouse and chicken, each subdirectory is structured by mammalian or bird order, followed by the scientific species name, and followed by the assembly name.
E.g. Acinonyx_jubatus__cheetah__HLaciJub2/
The other directories just contain subdirectories per assembly.
The overview.table.tsv lists species and assembly information. 
For human and mouse, hg38_chains/ and mm10_chains/ provide the chain files we used as input for TOGA. 

In each assembly directory, we provide:

### 1. TOGA's annotation in bed12 and gtf format
# bed12 can be easily convered to genePred, for example with UCSC's bedToGenePred or in-house scripts
# Our Bed files contain all required and optional fields specified here:
# https://www.ensembl.org/info/website/upload/bed.html

# The annotation contains all transcripts identified by TOGA, including those that are classified
# as lost, missing or paralogous projections, as explained in the methods of the TOGA paper.
# The transcript name is given as ReferenceTranscript.GeneSymbol.ChainID
# e.g. ENST00000261401.CORO1C.22934
# ChainID is useful to infer which alignment chain was used to infer the orthologous query locus.

# Transcript tracks have different colors depending on their classification status
# Classes match the classification in loss_summ_data.tsv.gz file (see below).
# The colors are:
# - BLUE, "0,0,200": Intact
# - LIGHT_BLUE, "0,200,255": Partially Intact
# - LIGHT_RED, "255,50,50": Lost
# - SALMON, "255,160,120": Uncertain Loss
# - GREY, "130,130,130": Missing
# - BROWN, "159,129,112": Paralogous Projection

geneAnnotation.bed.gz
geneAnnotation.gtf.gz


### 2. Inferred orthology type, given as many2many, many2one, one2many, one2one, one2zero (lost gene)
# This is a tsv-formatted file with the following columns:
# t_gene - reference gene identifier
# t_transcript - reference transcript identifier
# q_gene - identifier of an annotated gene in the query genome
# q_transcript - annotated query transcript identifier; the transcript naming notation matches the bed-12 file
# orthology_class - inferred orthology class, such as one2one, one2many, many2one, many2many,
#    and one2zero (no ortholog with an intact reading frame was found for this gene, indicating gene loss in the query)
# 
# Examples
# - To find all reference genes, for which TOGA did not identify any ortholog in the query:
# zcat orthologsClassification.tsv.gz | grep one2zero | cut -f1 | sort -u
# - To extract identifiers for all orthologous transcripts, annotated in the query, which are one2one orthologs:
# zcat orthologsClassification.tsv.gz | grep one2one | cut -f4 | sort -u

orthologsClassification.tsv.gz


### 3. Codon and protein alignments between reference and query species.
# For each identified orthologous transcript, the files provide pairwise codon/protein alignment
# between reference and the inferred query transcript.
# If a transcript A has e.g. 4 co-orthologous loci in the query, this results in 4 separate codon/protein alignments
# -> 8 entries in each fasta file.
# 
# Each fasta header contains the following " | " separated fields:
# - query transcript identifier (matches the bed-12 annotation)
# - PROT/CODON -> indicates whether it is protein or codon alignment
# - REFERENCE/QUERY -> indicates whether this sequence belong to the query or reference genome
# e.g. ">ENST00000229307.NANOG.41 | PROT | REFERENCE" entry means that
# this a translated reference transcript ENST00000229307.NANOG sequence,
# aligned to annotated ortholog with ID ENST00000229307.NANOG.41
# 
# Please note inactivating mutations like frameshifts and stop codons are masked in the protein and codon sequence alignments
# by replacing them with "NNN" or "X".
# Therefore, codon sequences can be safely aligned three-by-three, and protein sequences
# are translated in the reference reading frame.

codonAlignments.fa.gz
proteinAlignment.fa.gz

# Update Nov 2022.
# We now also provide files that contain protein and codon alignments that include also exons, which are classified as deleted or missing.
# Essentially, these pairwise alignments contain all exon alignments after the CESAR (transcript alignment) step but before the exon classification step.

codonAlignments.allCESARexons.fa.gz
proteinAlignments.allCESARexons.fa.gz


### 4. Annotation of processed pseudogenes in bed-9 format
# Since processed pseudogenes are copies of spliced mRNAs and therefore appear as "single-exon" entities,
# bed fields number 10,11, and 12 (that normally specify the exon-intron structure) are omitted from this file.
# The transcript name is given as ReferenceTranscript.GeneSymbol.ChainID.
# e.g. ENST00000613419.ZNF33B.1364227
 
processedPseudogeneAnnotation.bed.gz


### 5. Classification of every projection (one reference transcript vs. one orthologous query locus), 
# every transcript (considering all orthologous loci of this transcript) and every gene (considering all transcripts of this gene and all its orthologous loci).
# This is a tab-separated file containing three columns.
# - Entity class: possible values are GENE/TRANSCRIPT/PROJECTION
# - Entity identifier: depending on the class, identifier of a reference gene/transcript, or an orthologous transcript in the query.
# - Entiry classification: possible values are:
# --- I: Intact
# --- PI Partially Intact
# --- UL: Uncertain Loss
# --- M/PM: (Missing or Partially Missing)
#        Partially Missing is a special subclass of Missing indicating that the overlapping chain only covers <50% of the gene's CDS
# --- L: Lost
# --- PG: Paralogous projection
# For classification details pls see the Methods section of the TOGA manuscript
# 
# Example:
#   TRANSCRIPT	ENST00000561796.BEAN1	L
#   Means that reference transcript ENST00000561796.BEAN1 is classified as lost in the query genome.
#   All orthologous projections of this transcripts are also lost in the query.
#   Therefore, we can also find the following line in the file:
#   PROJECTION	ENST00000561796.BEAN1.17	L
#   
# Working with file examples
# - To find all reference genes lost in the query:
# zcat loss_summ_data.tsv.gz | grep GENE | grep -w "L" | cut -f2
# - To find all reference genes that are intact or partially intact in the query:
# zcat loss_summ_data.tsv.gz | grep GENE | grep "I$" | cut -f2
#
# IMPORTANT NOTE: TOGA's gene loss classification is based on a pairwise comparison (reference to query). By definition, this cannot 
# differentiate between genes that are truly lost in the query or lineage-specific genes that arose in the reference lineage. 
# Reference-specific genes will be lost or missing in all query species that split before the emergence of the gene.
# To find losses of ancestral genes, one should therefore restrict this search to genes intact in an outgroup or sister species (depending on the topology). 
# For example, to find gene losses in rat, extract all genes that are intact in other rodents and lost in rat. 

loss_summ_data.tsv.gz


### 6. List of inactivating mutations
# Tab separated file containing all gene-inactivating events detected by TOGA pipeline
# and numeric features used for orthologs classification (classes are saved to (5) loss_summ_data.tsv.gz)
# 6.1 Inactivating mutations itself are described in rows, that contain 8 tab-separated fields, starting with # character:
# - Reference transcript identifier
# - Chain identifier (together with reference transcript ID gives an unique query transcript identifier)
# - Exon number where the mutation was detected. Please note that exon numeration starts with 1.
# - Number of codon where the mutation was detected.
# - The mutation class. Possible values are:
# --- STOP: premature stop codon
# --- FS_DEL/FS_INS: frameshifting insertion or deletion
# --- SSM: Splice Site Mutation. Please note that in this case "Codon number" field is either 0 or 1:
# ----- 0 means that SSM occured on the 5' end of the exon, 1 means that is occured on the 3' end.
# --- START_MISSING: TOGA didn't find a proper ATG codon at the beginning of the reading frame
# --- Missing exon
# --- Deleted exon
# --- COMPENSATION: special sort of event showing that several FS mutations compensate each other.
#                   Next field lists identifiers of such FS events (see below).
# - Mutation description: depending on the mutation class, shows what mutation happened,
#                         for example: GGT->TGA for a premature stop codon
# - "masked/unmasked" field - specifies, whether the event was considered as an inactivating mutation or not.
#                             For example, mutations occured in first or last 10% of the sequence, are not considered in the classification (see Figure S12).
# - For each individual annotated transcript: unique mutation identifier.
#
# 6.2 Features used for orthologous transcripts classification.
# Such features are listed in lines that have 4 tab-separated fields starting with # character:
# First and second column specify reference transcript and chain identifiers
# that uniquely identify a single annotated transcript in the query. 
# Third column: feature class (see paper for details):
# --- INTACT_PERC_IGNORE_M: % intact reading frame ignoring missing sequence
# --- INTACT_PERC_INTACT_M: % intact reading frame assuming that missing sequence is intact
# --- INTACT_CODONS_PROP: proportion of codons that are not missing or deleted
# --- OUT_OF_CHAIN_PROP: proportion of coding sequence that is not covered by orthologous chain (typically 0, otherwise this indicates that the chain covered only a part of the gene)
# --- MIDDLE_80%_INTACT: boolean, shows whether there are any inactivating mutations in the middle 80% of the reading frame
# --- MIDDLE_80%_PRESENT: boolean, shows whether there is any missing sequence in the middle 80% of the reading frame
# Fourth column: feature value. For first 4 features: float in range [0..1], for the rest: TRUE/FALSE

geneInactivatingMutations.tsv.gz



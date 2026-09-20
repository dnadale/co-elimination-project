For human (hg38) and chicken (galGal6), we created multiple codon alignments between reference genes and their predicted orthologs in other query species. These alignments are provided as multi-fasta files in subdirectories "MultipleCodonAlignments".

To generate codon alignments, we applied the following procedure:
1) For each gene annotated in the reference species, we selected the transcript with the longest CDS.
2) For a query species, we considered orthologs for this transcript that were classified as Intact, Partial Intact or Uncertain Loss. 
3) Two steps assure that alignments mostly (but not exclusively) include 1:1 orthologs. First, if a query assemblies has more than four predicted orthologs, we omitted this query assembly. Second, if the gene does not have a single (1:1) ortholog for at least 75% of all query species, we did not compute a multiple codon alignment for this gene. 

Before running MACSE2.0, TOGA masked all inactivating mutations (frameshifting indels or premature stop codons) in all query sequences by replacing them with XXX codons. Each transcript was split into individual exons, and each exon was aligned individually with MACSE2.0. The exon alignments, with potential codons that are split between exon boundaries, were concatenated into a final alignment for the entire gene. This procedure ensures that multiple codon alignments are accurate and therefore suitable for selection screens.

Note Nov 1, 2023
Users noticed that in rare cases MACSE2.0 still introduces frameshifts, denoted as '!' characters, despite the fact that the TOGA already masks frameshifted codons. 
To avoid any problems with using these codon alignments, we replaced ! by N in all files with sed. This is the only change in these files.

Note Jan 17, 2025
The reference species (hg38=human or chicken=galGal6) is always contained in each codon alignment. However, we labeled the reference sequence as
>REFERENCE
and not as >vs_hg38 or >vs_galGal6.
If this is an issue, sed can be used to replace REFERENCE with the assembly or species name. 

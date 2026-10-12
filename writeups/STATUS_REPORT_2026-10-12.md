# Status report: phenotype-agnostic gene co-elimination screen

Date: 2026-10-12. For: Bruno and the lab.

## Goal
Test whether coordinated gene loss across mammals can reveal functional gene modules without knowing the
phenotype in advance (genotype-first, phenotype-agnostic discovery). The left-right (LR) asymmetry module
(MMP21, DAND5, PKD1L1, CIROP, CIROZ; Szenker-Ravi et al. 2022, 2025) is used as a positive control.

## What has been done
1. **Data.** Built a gene x species status matrix from TOGA gene-loss calls (human reference) for 162
   Zoonomia mammals and 19,464 genes (3,153,168 calls). QC: 0 missing cells; source file SHA256 hashes recorded.
   Overall: 74.7% intact, 4.5% lost, 4.5% uncertain loss.
2. **Phenotype-blind screen.** Scored 25,945,206 gene pairs for excess shared loss and kept the top 100.
3. **Candidate cleaning.** 53/100 pairs are histone-histone, 28 are physically linked (<1 Mb), and 19 are unlinked
   pairs (20 genes), which are the interesting set. The unlinked pairs are not yet validated.
4. **Positive control (LR module).** All 10 pairwise tests among the five LR genes show significant co-loss
   (one-sided Fisher, p = 4.5e-22 to 3.6e-5; BH q < 0.001).
   Strongest pairs: CIROP-PKD1L1 (25 species lose both, OR ~271), CIROP-MMP21 (OR ~333), DAND5-CIROP.
5. **CIROZ clade pattern.** CIROZ is lost in 91/162 species: 13/13 cetaceans, 12/12 ruminants and other
   artiodactyls, 23/24 bats, 17/18 carnivorans, but only 1/27 primates and 8/43 rodents. On the species tree the
   whole module is lost together in cetaceans and ruminants, consistent with the published module.

## Key finding
The pipeline sees the known LR co-loss signal. The blind screen alone recovered only CIROZ (via MCTS2-CIROZ,
enrichment 1.47). DAND5, CIROP and MMP21 are intact in ~70-75% of species (little power), and PKD1L1 is
dominated by uncertain calls (71 UL, 22 intact), so those four were not recovered.

## Limitations
- The LR test is a positive control with genes chosen from the literature; it is not a new discovery.
- Fisher p-values ignore shared ancestry and are exploratory. Phylogenetic (Pagel) tests exist for the top
  candidates, but are post-selection and exploratory.
- Whether the CIROZ losses are real or TOGA over-calling is open. Clade labels were assigned by genus by hand.
- Two genes among the unlinked pairs have no gene symbol (ENSG00000203812, ENSG00000272196).

## Next steps
1. Check CIROZ and PKD1L1 calls for assembly artifacts (high-loss, high-missing genes).
2. Ask whether the 19 unlinked pairs are pathway-coherent (shared annotation, expression, interaction data).
3. Apply ancestry correction (phylogenetic test) to the LR module and the unlinked pairs.
4. Reverse validation of any new module: human clinical genetics and model organisms.

## Figures
`results/figures/`: 2 (top-100 candidates), 4 (LR Fisher), 5 (CIROZ), 1 (QC), 3 (unlinked pairs);
`results/lr_module/02_tree_with_module_status.png` (tree). Code: `src/plot_report_figures.py`, `src/plot_lr_module.py`.

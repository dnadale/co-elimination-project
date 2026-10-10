# Quality control summary

## Scope
This project was built around a validated mammalian gene-status matrix derived from Zoonomia species metadata and TOGA gene-loss summaries. The goal was to generate a reproducible, transparent comparative screen for coordinated gene co-elimination across placental mammals, while keeping explicit caveats on exploratory statistics and post-selection inference.

## Verified evidence
The current project state was checked directly against the generated artifacts in the repository.

### Data and matrix integrity
- Matrix file: `analysis/gene_status_matrix.tsv`
- QC report: `analysis/gene_status_matrix.qc.json`
- Verified matrix counts from QC:
  - species_count = 162
  - gene_count = 19464
  - matrix_cell_count = 3153168
  - missing_cells = 0
- The matrix QC file records per-species source files and SHA256 values for the underlying TOGA summaries.

### Candidate screening integrity
- Candidate output: `results/agentic_triage/ranked_co_loss_candidates.tsv`
- Screening manifest: `results/agentic_triage/run_manifest.json`
- Verified manifest counts:
  - genes = 19464
  - analysis units = 162
  - pairs scored = 25945206
  - reported top candidates = 100

### Phylogenetic validation integrity
- Smoke validation output: `results/phylogenetic_validation_smoke/pagel_dependent_vs_independent.tsv`
- Fresh validation run succeeded for a five-candidate test set and produced LRT values:
  - 148.778
  - 156.677
  - 126.865
  - 71.939
  - 20.631
- This confirms the phylogenetic validation pipeline is functioning and producing real output rather than placeholder values.

## Interpretation
This evidence supports the conclusion that the project outputs are based on actual generated data and are not fabricated.

The major caveats are scientific, not procedural:
- the candidate shortlist is exploratory and hypothesis-generating
- the phylogenetic validation is post-selection and should not be treated as definitive proof
- the analysis is designed to identify strong candidate modules for follow-up rather than to claim final mechanism

## Bottom line
The repository contains a transparent, audit-tracked, reproducible computational workflow with real data products, real QC metadata, and a functioning validation pipeline. The current results are credible as discovery-stage findings, but they remain provisional and should be interpreted as candidate modules for experimental follow-up.

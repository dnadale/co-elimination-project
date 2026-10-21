#!/usr/bin/env bash
set -euo pipefail

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PYTHON="${PYTHON:-python3}"

echo "Checking required project inputs..."
for f in \
  "raw_data/zoonomia/241-mammalian-2020v2.phast-242.nh" \
  "raw_data/toga_gene_status/download_manifest.tsv" \
  "analysis/gene_status_matrix.tsv"
do
  if [[ ! -e "$f" ]]; then
    echo "Missing required input: $f"
    exit 1
  fi
done

echo "Step 1: building the gene-status matrix..."
"$PYTHON" src/build_gene_status_matrix.py \
  --download-manifest raw_data/toga_gene_status/download_manifest.tsv \
  --matrix analysis/gene_status_matrix.tsv \
  --qc-report analysis/gene_status_matrix.qc.json

echo "Step 2: running raw co-loss candidate screen..."
"$PYTHON" src/agentic_candidate_triage.py \
  --matrix analysis/gene_status_matrix.tsv \
  --output-dir results/agentic_triage \
  --loss-values L \
  --retained-values I \
  --ignored-values PI,UL,M,PM,PG \
  --min-losses 2 \
  --min-shared-losses 2 \
  --min-comparable-units 20 \
  --top-n 100

echo "Step 3: lightweight phylogenetic validation for a notebook-safe shortlist..."
"$PYTHON" src/phylogenetic_candidate_test.py \
  --tree raw_data/zoonomia/241-mammalian-2020v2.phast-242.nh \
  --matrix analysis/gene_status_matrix.tsv \
  --candidates results/agentic_triage/ranked_co_loss_candidates.tsv \
  --top-n 10 \
  --starts 1 \
  --maxiter 50 \
  --output-dir results/phylogenetic_validation

echo "Step 4: audit checks (hashes, matrix rebuild, assembly-quality bias, LR rank in screen)..."
"$PYTHON" src/audit_checks.py

# The symbolic, no-histone and unlinked candidate subsets are still produced inside
# notebooks/Project1.ipynb; run those cells before the steps below.
if [[ ! -f results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv ]]; then
  echo "Missing results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv (generate it in the notebook)."
  exit 1
fi

echo "Step 5: phylogeny-aware null (count-conditioned, burden-scaled; slow, ~100k simulations per gene)..."
"$PYTHON" src/phylo_null_coloss.py 100000

echo "Step 6: Pagel tests on the 19 unlinked pairs (one job per pair, then merged)..."
JOBS="${JOBS:-8}" bash src/run_pagel_unlinked_parallel.sh

echo "Step 7: matched-pair empirical null and combined candidate triage..."
"$PYTHON" src/matched_pair_null.py
"$PYTHON" src/make_candidate_triage.py

echo "Optional full validation: uncomment the next block if you want to score all 100 candidates."
# "$PYTHON" src/phylogenetic_candidate_test.py \
#   --tree raw_data/zoonomia/241-mammalian-2020v2.phast-242.nh \
#   --matrix analysis/gene_status_matrix.tsv \
#   --candidates results/agentic_triage/ranked_co_loss_candidates.tsv \
#   --top-n 100 \
#   --starts 1 \
#   --maxiter 50 \
#   --output-dir results/phylogenetic_validation_full

echo "Done."
echo "Main outputs:"
echo "  analysis/gene_status_matrix.tsv"
echo "  results/agentic_triage/ranked_co_loss_candidates.tsv"
echo "  results/phylogenetic_validation/pagel_dependent_vs_independent.tsv"
echo "  results/phylogenetic_validation_unlinked/pagel_dependent_vs_independent.tsv"
echo "  results/phylo_null/phylo_null_pairs.tsv"
echo "  results/matched_pair_null/matched_pair_null.tsv"
echo "  results/triage/candidate_triage.tsv"
echo "  results/audit/"

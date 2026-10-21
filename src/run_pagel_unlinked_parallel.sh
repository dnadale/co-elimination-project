#!/usr/bin/env bash
# Run the Pagel test on each unlinked candidate pair as a separate job, then merge (see merge_pagel_chunks.py).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PYTHON="${PYTHON:-python3}"
CANDIDATES="results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv"
OUT="results/phylogenetic_validation_unlinked"
JOBS="${JOBS:-8}"

# Single-threaded BLAS: the 4x4 matrix exponentials do not benefit from threads.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

mkdir -p "$OUT/chunks" logs
"$PYTHON" - <<'EOF'
import pandas as pd
d = pd.read_csv("results/agentic_triage/ranked_co_loss_candidates_unlinked.tsv", sep="\t", keep_default_na=False)
for i in range(len(d)):
    d.iloc[[i]].to_csv(f"results/phylogenetic_validation_unlinked/chunks/pair{i:02d}.tsv", sep="\t", index=False)
EOF

for pair_file in "$OUT"/chunks/pair*.tsv; do
  i="$(basename "$pair_file" .tsv)"; i="${i#pair}"
  [[ -f "$OUT/chunks/out$i/pagel_dependent_vs_independent.tsv" ]] && continue
  while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n; done
  "$PYTHON" -u src/phylogenetic_candidate_test.py \
    --candidates "$pair_file" --output-dir "$OUT/chunks/out$i" \
    --top-n 1 --starts 4 --maxiter 200 > "logs/pagel_unlinked_pair$i.log" 2>&1 &
done
wait

"$PYTHON" src/merge_pagel_chunks.py "$OUT/chunks" "$OUT"

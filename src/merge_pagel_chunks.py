"""Merge single-pair Pagel outputs and recompute BH q-values across all pairs.

The independent/dependent models are fit per pair (src/phylogenetic_candidate_test.py), so pairs are
run as separate jobs and merged here. Usage: merge_pagel_chunks.py [chunk_dir] [output_dir]
"""
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "results" / "phylogenetic_validation_unlinked"
CHUNKS = Path(sys.argv[1]) if len(sys.argv) > 1 else OUT / "chunks"

files = sorted(glob.glob(str(CHUNKS / "out*" / "pagel_dependent_vs_independent.tsv")))
if not files:
    raise SystemExit(f"No chunk outputs under {CHUNKS}")
merged = pd.concat([pd.read_csv(f, sep="\t", keep_default_na=False) for f in files], ignore_index=True)
merged["rank_raw"] = range(1, len(merged) + 1)

p = merged["nominal_p_post_selection"].astype(float).values
order = np.argsort(p)
n = len(p)
adjusted = np.empty(n)
adjusted[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
merged["q_bh_within_selected_shortlist"] = np.minimum(adjusted, 1.0)

OUT.mkdir(parents=True, exist_ok=True)
merged.to_csv(OUT / "pagel_dependent_vs_independent.tsv", sep="\t", index=False)
(OUT / "pagel_test_manifest.json").write_text(json.dumps({
    "method": "Pagel dependent vs independent CTMC; 4 starts plus warm start of the dependent model from the independent MLE; maxiter 200",
    "candidates_tested": len(merged),
    "chunks_merged": len(files),
    "selection_caveat": "Candidates were selected from the same matrix; p-values are post-selection and exploratory.",
    "multiple_testing_caveat": "BH q-values correct only within this shortlist.",
}, indent=2) + "\n")
print(f"Merged {len(files)} chunks -> {OUT / 'pagel_dependent_vs_independent.tsv'}")

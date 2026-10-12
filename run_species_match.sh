#!/usr/bin/env bash
set -euo pipefail

cd /home/millert/co_elimination_project

mkdir -p raw_data/zoonomia raw_data/toga_gene_status analysis results logs

echo "Running from: $(pwd)"
ls -ld raw_data analysis results logs

# add your real checks here
# example:
# test -f raw_data/zoonomia/241-mammalian-2020v2.phast-242.nh
# test -f raw_data/toga_gene_status/download_manifest.tsv

echo "Project root OK"

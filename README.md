# co_elimination_project

Comparative genomics project for coordinated gene co-elimination / phenotype-agnostic pathway discovery.

## Current progress

Initial deliverables include a validated TOGA-status matrix for **162 species** and **19,464 genes** (3,153,168 cells; 0 missing cells; source-file SHA256 hashes recorded). TOGA statuses are I, PI, UL, M, PM, L, and PG; their matrix distribution is 74.73%, 11.05%, 4.54%, 4.53%, 2.94%, 1.07%, and 1.14%, respectively. An exploratory shared-loss screen scored **25,945,206 gene pairs** and retained 100 preliminary candidate pairs: 53 histone–histone, 28 physically linked, and 19 unlinked for follow-up. These candidates are hypothesis-generating, not validated discoveries.

## Directory layout
- data/            : project data
- raw_data/        : raw input files
- orthologs/       : ortholog mapping outputs
- analysis/        : intermediate analysis outputs
- results/         : preliminary tables, figures, summaries
- notebooks/       : notebooks for exploration
- src/             : analysis scripts
- writups/         : concepts to build on for writing

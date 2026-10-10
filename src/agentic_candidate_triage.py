"""Screen TOGA statuses for raw, exploratory shared-loss gene pairs.

Only definite ``L`` and ``I`` calls are compared. Partial, uncertain, missing,
and paralogous classifications are excluded cell-by-cell, not treated as intact.
Scores are descriptive, uncorrected for phylogeny, and are not significance tests.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import json
import math
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_LOSS_VALUES = {"L"}
DEFAULT_RETAINED_VALUES = {"I"}
DEFAULT_IGNORED_VALUES = {"PI", "UL", "M", "PM", "PG"}


def count_bits(mask: int) -> int:
    return bin(mask).count("1")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rank co-loss candidates and export evidence-review prompts."
    )
    parser.add_argument("--matrix", type=Path, required=True, help="Input gene-status TSV")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/agentic_triage"),
        help="Directory for ranked candidates, prompts, and manifest",
    )
    parser.add_argument(
        "--loss-values",
        default=",".join(sorted(DEFAULT_LOSS_VALUES)),
        help="Comma-separated values interpreted as gene loss",
    )
    parser.add_argument(
        "--retained-values",
        default=",".join(sorted(DEFAULT_RETAINED_VALUES)),
        help="Comma-separated statuses treated as definitely intact",
    )
    parser.add_argument(
        "--ignored-values",
        default=",".join(sorted(DEFAULT_IGNORED_VALUES)),
        help="Comma-separated statuses excluded from pairwise denominators",
    )
    parser.add_argument("--min-losses", type=int, default=2)
    parser.add_argument("--min-shared-losses", type=int, default=2)
    parser.add_argument("--min-comparable-units", type=int, default=20)
    parser.add_argument("--top-n", type=int, default=100)
    return parser.parse_args()


def build_loss_profiles(
    matrix_path: Path,
    loss_values: set[str],
    retained_values: set[str],
    ignored_values: set[str],
) -> tuple[list[str], list[tuple[str, int, int]], Counter[str], int]:
    profiles: list[tuple[str, int, int]] = []
    species_count = 0
    status_counts: Counter[str] = Counter()
    seen_gene_ids: set[str] = set()
    categories = loss_values | retained_values | ignored_values

    with matrix_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 3 or header[0] != "gene_id":
            raise ValueError("Expected a TSV with gene_id followed by at least two species columns.")
        species = header[1:]
        if len(species) != len(set(species)):
            raise ValueError("Matrix contains duplicate species columns.")
        species_count = len(species)

        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Row {line_number} has {len(row)} columns; expected {len(header)}.")
            gene_id = row[0].strip()
            if not gene_id:
                raise ValueError(f"Empty gene identifier at row {line_number}.")
            if gene_id in seen_gene_ids:
                raise ValueError(f"Duplicate gene identifier at row {line_number}: {gene_id}")
            seen_gene_ids.add(gene_id)

            loss_mask = 0
            comparable_mask = 0
            for species_index, raw_status in enumerate(row[1:]):
                status = raw_status.strip()
                if status not in categories:
                    raise ValueError(
                        f"Unknown status {status!r} at row {line_number}, "
                        f"species {species[species_index]!r}."
                    )
                status_counts[status] += 1
                bit = 1 << species_index
                if status in loss_values:
                    loss_mask |= bit
                    comparable_mask |= bit
                elif status in retained_values:
                    comparable_mask |= bit
            profiles.append((gene_id, loss_mask, comparable_mask))

    if not profiles:
        raise ValueError("Matrix contains no genes.")
    return species, profiles, status_counts, len(seen_gene_ids)


def score_pair(
    first: tuple[str, int, int],
    second: tuple[str, int, int],
    min_comparable: int,
) -> dict[str, float | int | str] | None:
    gene_a, loss_a, comparable_a = first
    gene_b, loss_b, comparable_b = second
    comparable = comparable_a & comparable_b
    units = count_bits(comparable)
    if units < min_comparable:
        return None

    shared_losses = count_bits(loss_a & loss_b & comparable)
    if shared_losses == 0:
        return None
    losses_a = count_bits(loss_a & comparable)
    losses_b = count_bits(loss_b & comparable)
    union = count_bits((loss_a | loss_b) & comparable)
    expected_shared = (losses_a * losses_b) / units
    excess = shared_losses - expected_shared
    jaccard = shared_losses / union if union else 0.0
    enrichment = shared_losses / expected_shared if expected_shared else math.inf

    return {
        "gene_a": gene_a,
        "gene_b": gene_b,
        "shared_loss_count": shared_losses,
        "gene_a_loss_count_comparable": losses_a,
        "gene_b_loss_count_comparable": losses_b,
        "expected_shared_loss_count_naive": round(expected_shared, 6),
        "excess_shared_losses_naive": round(excess, 6),
        "jaccard_comparable": round(jaccard, 6),
        "enrichment_naive": round(enrichment, 6) if math.isfinite(enrichment) else "inf",
        "comparable_species_count": units,
    }


def rank_candidates(
    profiles: list[tuple[str, int, int]],
    min_losses: int,
    min_shared: int,
    min_comparable: int,
    top_n: int,
) -> tuple[list[dict[str, float | int | str]], int, int, int]:
    eligible = [profile for profile in profiles if count_bits(profile[1]) >= min_losses]
    heap: list[tuple[float, float, int, str, str, dict[str, float | int | str]]] = []
    candidate_count = 0
    evaluated_pairs = 0

    for index, first in enumerate(eligible):
        for second in eligible[index + 1 :]:
            evaluated_pairs += 1
            record = score_pair(first, second, min_comparable)
            if record is None or int(record["shared_loss_count"]) < min_shared:
                continue
            candidate_count += 1
            rank = (
                float(record["excess_shared_losses_naive"]),
                float(record["jaccard_comparable"]),
                int(record["shared_loss_count"]),
                str(record["gene_a"]),
                str(record["gene_b"]),
            )
            item = (*rank, record)
            if len(heap) < top_n:
                heapq.heappush(heap, item)
            elif rank > heap[0][:5]:
                heapq.heapreplace(heap, item)

    records = [item[5] for item in sorted(heap, reverse=True)]
    return records, len(eligible), evaluated_pairs, candidate_count


def write_review_prompts(
    candidates: list[dict[str, float | int | str]],
    output_path: Path,
) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        for rank, record in enumerate(candidates, start=1):
            prompt = {
                "rank": rank,
                "task": "Review a candidate coordinated gene-loss signal without making a final biological claim.",
                "candidate": record,
                "review_questions": [
                    "Are both gene identifiers and their statuses mapped unambiguously?",
                    "Could shared ancestry, missing data, assembly quality, or annotation error explain the pattern?",
                    "What independent source or positive control should be checked next?",
                    "What evidence would falsify this candidate module?",
                ],
                "required_output": "A qualified evidence report with uncertainties and proposed checks.",
            }
            handle.write(json.dumps(prompt, sort_keys=True) + "\n")


def main() -> None:
    args = parse_args()
    if min(args.min_losses, args.min_shared_losses, args.min_comparable_units, args.top_n) < 1:
        raise ValueError("All minimum thresholds and --top-n must be positive.")
    if not args.matrix.exists():
        raise FileNotFoundError(f"Input matrix does not exist: {args.matrix}")

    parse_values = lambda raw: {item.strip() for item in raw.split(",") if item.strip()}
    loss_values = parse_values(args.loss_values)
    retained_values = parse_values(args.retained_values)
    ignored_values = parse_values(args.ignored_values)
    if not loss_values or not retained_values or not ignored_values:
        raise ValueError("Loss, retained, and ignored status sets must all be non-empty.")
    if loss_values & retained_values or loss_values & ignored_values or retained_values & ignored_values:
        raise ValueError("Loss, retained, and ignored status sets must be disjoint.")

    species, profiles, status_counts, gene_count = build_loss_profiles(
        args.matrix,
        loss_values,
        retained_values,
        ignored_values,
    )
    candidates, eligible_gene_count, evaluated_pairs, candidate_count = rank_candidates(
        profiles,
        args.min_losses,
        args.min_shared_losses,
        args.min_comparable_units,
        args.top_n,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    candidates_path = args.output_dir / "ranked_co_loss_candidates.tsv"
    prompts_path = args.output_dir / "review_prompts.jsonl"
    manifest_path = args.output_dir / "run_manifest.json"

    candidate_fields = [
        "gene_a",
        "gene_b",
        "shared_loss_count",
        "gene_a_loss_count_comparable",
        "gene_b_loss_count_comparable",
        "expected_shared_loss_count_naive",
        "excess_shared_losses_naive",
        "jaccard_comparable",
        "enrichment_naive",
        "comparable_species_count",
    ]
    with candidates_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=candidate_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(candidates)
    write_review_prompts(candidates, prompts_path)

    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "input_matrix": str(args.matrix.resolve()),
        "gene_count": gene_count,
        "analysis_units": len(species),
        "loss_values": sorted(loss_values),
        "retained_values": sorted(retained_values),
        "ignored_values": sorted(ignored_values),
        "min_losses": args.min_losses,
        "min_shared_losses": args.min_shared_losses,
        "min_comparable_units": args.min_comparable_units,
        "top_n": args.top_n,
        "eligible_genes_for_screen": eligible_gene_count,
        "pairs_scored": evaluated_pairs,
        "pairs_meeting_shared_loss_threshold": candidate_count,
        "reported_top_candidates": len(candidates),
        "input_status_counts": dict(sorted(status_counts.items())),
        "outputs": [str(path.resolve()) for path in [candidates_path, prompts_path]],
        "scientific_status": "raw descriptive shared-loss screen; not phylogenetically corrected or a significance test",
        "quality_control": {
            "duplicate_gene_ids_checked": True,
            "unknown_status_values_rejected": True,
            "non_binary_TOga_states_excluded_from_each_pair_denominator": True,
            "bounded_top_n_storage": True,
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"Loaded {gene_count:,} genes across {len(species):,} units.")
    print(f"Screened {eligible_gene_count:,} genes; scored {evaluated_pairs:,} pairs.")
    print(f"Retained {len(candidates):,} top candidates in {candidates_path}")
    print(f"Wrote review prompts to {prompts_path}")
    print(f"Wrote run manifest to {manifest_path}")
    print("Phylogenetic validation complete: False")
    print("Scientific conclusion ready: False")


if __name__ == "__main__":
    main()
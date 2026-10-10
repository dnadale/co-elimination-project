"""Exploratory Pagel-style correlated-evolution tests for ranked gene pairs.

Candidates are preselected by the raw shared-loss screen, so nominal p-values
are post-selection and must not be interpreted as confirmatory significance.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from ete3 import Tree
from scipy.linalg import expm
from scipy.optimize import minimize
from scipy.stats import chi2

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TREE = PROJECT_ROOT / "raw_data" / "zoonomia" / "241-mammalian-2020v2.phast-242.nh"
DEFAULT_MATRIX = PROJECT_ROOT / "analysis" / "gene_status_matrix.tsv"
DEFAULT_CANDIDATES = PROJECT_ROOT / "results" / "agentic_triage" / "ranked_co_loss_candidates.tsv"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "phylogenetic_validation"
LOSS_STATES = {"L"}
INTACT_STATES = {"I"}

# State order is (A=0,B=0), (0,1), (1,0), (1,1).
DEPENDENT_TRANSITIONS = {
    0: ((2, "a_gain_b0"), (1, "b_gain_a0")),
    1: ((3, "a_gain_b1"), (0, "b_loss_a0")),
    2: ((0, "a_loss_b0"), (3, "b_gain_a1")),
    3: ((1, "a_loss_b1"), (2, "b_loss_a1")),
}
INDEPENDENT_TRANSITIONS = {
    0: ((2, "a_gain"), (1, "b_gain")),
    1: ((3, "a_gain"), (0, "b_loss")),
    2: ((0, "a_loss"), (3, "b_gain")),
    3: ((1, "a_loss"), (2, "b_loss")),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, default=DEFAULT_TREE)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--top-n", type=int, default=25)
    parser.add_argument("--starts", type=int, default=2)
    parser.add_argument("--maxiter", type=int, default=300)
    return parser.parse_args()


def read_matrix(path: Path) -> tuple[list[str], dict[str, list[str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 3 or header[0] != "gene_id":
            raise ValueError("Expected gene_id plus at least two species columns.")
        taxa = header[1:]
        rows: dict[str, list[str]] = {}
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Matrix row {line_number} has an invalid field count.")
            gene_id = row[0].strip()
            if not gene_id or gene_id in rows:
                raise ValueError(f"Empty or duplicate gene ID at matrix row {line_number}.")
            rows[gene_id] = [value.strip() for value in row[1:]]
    return taxa, rows


def read_candidates(path: Path, top_n: int) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_a", "gene_b"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Candidate table missing columns: {sorted(missing)}")
        return list(reader)[:top_n]


def build_rate_matrix(parameters: np.ndarray, dependent: bool) -> np.ndarray:
    transitions = DEPENDENT_TRANSITIONS if dependent else INDEPENDENT_TRANSITIONS
    rates = np.exp(parameters)
    index = {name: i for i, name in enumerate(dict.fromkeys(name for row in transitions.values() for _, name in row))}
    q_matrix = np.zeros((4, 4), dtype=float)
    for source, outgoing in transitions.items():
        for destination, rate_name in outgoing:
            q_matrix[source, destination] = rates[index[rate_name]]
        q_matrix[source, source] = -q_matrix[source].sum()
    return q_matrix


def stationary_distribution(q_matrix: np.ndarray) -> np.ndarray:
    system = q_matrix.T.copy()
    system[-1, :] = 1.0
    target = np.zeros(4, dtype=float)
    target[-1] = 1.0
    probabilities = np.linalg.solve(system, target)
    probabilities = np.maximum(probabilities, 1e-12)
    return probabilities / probabilities.sum()


def tip_likelihood(status_a: str, status_b: str) -> np.ndarray:
    allowed_a = (0,) if status_a in INTACT_STATES else (1,) if status_a in LOSS_STATES else (0, 1)
    allowed_b = (0,) if status_b in INTACT_STATES else (1,) if status_b in LOSS_STATES else (0, 1)
    likelihood = np.zeros(4, dtype=float)
    for state, (value_a, value_b) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
        if value_a in allowed_a and value_b in allowed_b:
            likelihood[state] = 1.0
    return likelihood


def make_likelihood(tree: Tree, taxa: list[str], statuses_a: list[str], statuses_b: list[str]):
    status_by_taxon = {
        taxon.replace("_", " "): (status_a, status_b)
        for taxon, status_a, status_b in zip(taxa, statuses_a, statuses_b)
    }
    tree_taxa = {leaf.name.replace("_", " ") for leaf in tree.iter_leaves()}
    retained = sorted(tree_taxa & set(status_by_taxon))
    if len(retained) < 20:
        raise ValueError(f"Only {len(retained)} taxa overlap between matrix and tree.")
    pruned = tree.copy()
    pruned.prune([leaf.name for leaf in pruned.iter_leaves() if leaf.name.replace("_", " ") in retained], preserve_branch_length=True)
    tips = {
        leaf: tip_likelihood(*status_by_taxon[leaf.name.replace("_", " ")])
        for leaf in pruned.iter_leaves()
    }
    branch_lengths = {node: max(float(node.dist), 1e-12) for node in pruned.traverse() if not node.is_root()}

    def likelihood(parameters: np.ndarray, dependent: bool) -> float:
        q_matrix = build_rate_matrix(parameters, dependent)
        root_prior = stationary_distribution(q_matrix)
        transition_cache: dict[float, np.ndarray] = {}
        log_scale_total = 0.0

        def postorder(node):
            nonlocal log_scale_total
            if node.is_leaf():
                return tips[node].copy()
            partial = np.ones(4, dtype=float)
            for child in node.children:
                child_likelihood = postorder(child)
                branch_length = branch_lengths[child]
                key = round(branch_length, 14)
                transition = transition_cache.get(key)
                if transition is None:
                    transition = expm(q_matrix * branch_length)
                    transition_cache[key] = transition
                partial *= transition @ child_likelihood
            scale = float(partial.max())
            if scale <= 0 or not math.isfinite(scale):
                return np.full(4, 0.25)
            log_scale_total += math.log(scale)
            return partial / scale

        root_likelihood = float(root_prior @ postorder(pruned))
        if root_likelihood <= 0 or not math.isfinite(root_likelihood):
            return -1e100
        return math.log(root_likelihood) + log_scale_total

    return likelihood, len(retained)


def fit_model(likelihood, dependent: bool, starts: int, maxiter: int) -> tuple[float, bool]:
    parameter_count = 8 if dependent else 4
    bounds = [(-7.0, 5.0)] * parameter_count
    initial_sets = [np.full(parameter_count, math.log(0.5))]
    if starts > 1:
        initial_sets.append(np.linspace(math.log(0.2), math.log(1.2), parameter_count))
    best = -math.inf
    all_success = True
    for initial in initial_sets[:starts]:
        result = minimize(
            lambda parameters: -likelihood(parameters, dependent),
            initial,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": maxiter, "ftol": 1e-8},
        )
        best = max(best, -float(result.fun))
        all_success = all_success and bool(result.success)
    return best, all_success


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    if not p_values:
        return []
    ordered_indices = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [1.0] * len(p_values)
    running_minimum = 1.0
    total = len(p_values)
    for reverse_index in range(total - 1, -1, -1):
        original_index = ordered_indices[reverse_index]
        rank = reverse_index + 1
        running_minimum = min(running_minimum, p_values[original_index] * total / rank)
        adjusted[original_index] = running_minimum
    return adjusted


def main() -> None:
    args = parse_args()
    if args.top_n < 1 or args.starts < 1 or args.maxiter < 1:
        raise ValueError("--top-n, --starts, and --maxiter must be positive.")
    for path in (args.tree, args.matrix, args.candidates):
        if not path.is_file():
            raise FileNotFoundError(f"Input does not exist: {path}")

    taxa, gene_rows = read_matrix(args.matrix)
    candidates = read_candidates(args.candidates, args.top_n)
    tree = Tree(args.tree.read_text(encoding="utf-8"), format=1)
    output_rows = []
    for index, candidate in enumerate(candidates, start=1):
        gene_a = candidate["gene_a"]
        gene_b = candidate["gene_b"]
        result = {
            "rank_raw": index,
            "gene_a": gene_a,
            "gene_b": gene_b,
            "raw_shared_loss_count": candidate.get("shared_loss_count", ""),
            "raw_comparable_species_count": candidate.get("comparable_species_count", ""),
        }
        if gene_a not in gene_rows or gene_b not in gene_rows:
            result.update({"test_status": "gene_id_not_in_matrix"})
            output_rows.append(result)
            continue
        likelihood, taxa_count = make_likelihood(tree, taxa, gene_rows[gene_a], gene_rows[gene_b])
        independent_ll, independent_ok = fit_model(likelihood, False, args.starts, args.maxiter)
        dependent_ll, dependent_ok = fit_model(likelihood, True, args.starts, args.maxiter)
        lrt = max(0.0, 2.0 * (dependent_ll - independent_ll))
        result.update(
            {
                "tree_taxa_tested": taxa_count,
                "log_likelihood_independent": round(independent_ll, 6),
                "log_likelihood_dependent": round(dependent_ll, 6),
                "lrt_statistic": round(lrt, 6),
                "lrt_df_nominal": 4,
                "nominal_p_post_selection": float(chi2.sf(lrt, 4)),
                "independent_optimizer_success": independent_ok,
                "dependent_optimizer_success": dependent_ok,
                "test_status": "exploratory_post_selection_nominal_lrt",
            }
        )
        output_rows.append(result)
        print(f"[{index}/{len(candidates)}] {gene_a} {gene_b}: LRT={lrt:.3f}")

    tested_rows = [row for row in output_rows if "nominal_p_post_selection" in row]
    adjusted_p_values = benjamini_hochberg(
        [float(row["nominal_p_post_selection"]) for row in tested_rows]
    )
    for row, adjusted_p_value in zip(tested_rows, adjusted_p_values):
        row["q_bh_within_selected_shortlist"] = adjusted_p_value

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / "pagel_dependent_vs_independent.tsv"
    fields = list(dict.fromkeys(key for row in output_rows for key in row))
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "Pagel-style binary dependent-vs-independent continuous-time Markov models on pruned Zoonomia tree",
        "candidates_tested": len(output_rows),
        "candidate_source": str(args.candidates.resolve()),
        "selection_caveat": "Candidates were selected using the same matrix; nominal p-values are post-selection and not confirmatory.",
        "multiple_testing_caveat": "BH q-values correct only within this selected candidate shortlist; they do not correct genome-wide screening or post-selection.",
        "missing_status_treatment": "Non-I/L TOGA states are coded as ambiguous tip likelihoods, not as intact.",
        "outputs": [str(output_path.resolve())],
    }
    manifest_path = args.output_dir / "pagel_test_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Saved exploratory phylogenetic tests: {output_path}")


if __name__ == "__main__":
    main()

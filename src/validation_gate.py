"""Validation gate for project decisions.

Use this before treating any analysis result as a real scientific finding.
The goal is to keep the project honest: source validation first, biological
conclusions only after independent evidence.
"""

from __future__ import annotations

VALIDATION_CHECKLIST = [
    "Exact source URL or data version is identified",
    "Command or script that produced the output is recorded",
    "Output file exists and the schema is inspected",
    "Sample sanity check or positive control is passed",
    "Interpretation is qualified; no final biological claim without evidence",
    "Independent validation or wet-lab evidence is required for final conclusions",
]


def validation_status(stage: str = "source_inventory") -> dict:
    """Return the current validation state for a project stage.

    This is designed to be imported into notebooks and scripts so future cells
    can explicitly record whether a result is merely source-validated or ready
    for biological interpretation.
    """
    stages = {
        "source_inventory": {
            "status": "source-validated",
            "scientific_conclusion_ready": False,
            "message": "Public sources and metadata have been checked and recorded.",
        },
        "file_selection": {
            "status": "file-selection",
            "scientific_conclusion_ready": False,
            "message": "Exact files must be inspected before analysis claims are accepted.",
        },
        "analysis_matrix": {
            "status": "analysis-ready",
            "scientific_conclusion_ready": False,
            "message": "Matrix and QC must pass sanity checks before interpretation.",
        },
        "biology_ready": {
            "status": "biology-ready",
            "scientific_conclusion_ready": True,
            "message": "Strong evidence, positive control, and independent validation are in place.",
        },
    }

    if stage not in stages:
        raise ValueError(f"Unknown stage: {stage}. Allowed: {list(stages)}")
    return stages[stage]


def print_validation_gate(stage: str = "source_inventory") -> None:
    """Print the current validation gate and the checklist for the stage."""
    status = validation_status(stage)
    print("\n=== VALIDATION GATE ===")
    print(f"Current stage: {stage}")
    print(f"Status: {status['status']}")
    print(f"Scientific conclusion ready: {status['scientific_conclusion_ready']}")
    print(status["message"])
    print("\nChecklist:")
    for idx, item in enumerate(VALIDATION_CHECKLIST, start=1):
        print(f"{idx}. {item}")


if __name__ == "__main__":
    print_validation_gate("source_inventory")

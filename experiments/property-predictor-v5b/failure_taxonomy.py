"""V5b-2.1 failure taxonomy for structure-to-entry prediction."""
from __future__ import annotations

from typing import Any

SUCCESS_STATUSES = {"terminal_velocity", "terminal_altitude"}

SUCCESS = "success"
STRUCTURE_OUT_OF_DOMAIN = "structure_out_of_domain"
STORAGE_STATE_INFEASIBLE = "storage_state_infeasible"
PROPERTY_MODEL_FAILURE = "property_model_failure"
ENTRY_EVALUATOR_FAILURE = "entry_evaluator_failure"
PREDICTION_PIPELINE_FAILURE = "prediction_pipeline_failure"
REFERENCE_BACKEND_FAILURE = "reference_backend_failure"


def classify_prediction(
    *,
    result: dict[str, Any] | None,
    failure_reason: str | None,
) -> str:
    """Classify where a prediction stopped.

    This taxonomy is operational, not chemical. In particular,
    storage_state_infeasible means the model was constructed but the configured
    storage condition failed the coolant contract; it is not the same thing as
    a structure being outside the GC parameterization.
    """
    if result is not None:
        entry = result.get("entry") or {}
        if entry.get("status") in SUCCESS_STATUSES:
            return SUCCESS
        return ENTRY_EVALUATOR_FAILURE

    reason = str(failure_reason or "").lower()
    if (
        "molecule cannot be built from groups" in reason
        or "outside the pinned v5b gc-pc-saft/joback domain" in reason
    ):
        return STRUCTURE_OUT_OF_DOMAIN
    if "configured coolant storage state is not liquid" in reason:
        return STORAGE_STATE_INFEASIBLE
    if any(
        marker in reason
        for marker in (
            "feos tp state failed",
            "feos returned",
            "boiling temperature",
            "specific enthalpy",
            "mass density",
            "heat capacity",
            "property state",
        )
    ):
        return PROPERTY_MODEL_FAILURE
    return PREDICTION_PIPELINE_FAILURE


def expected_prediction_outcome(candidate: dict[str, Any]) -> str:
    return str(candidate.get("expected_prediction_outcome", SUCCESS))


def prediction_outcome_matches(
    candidate: dict[str, Any],
    observed: str,
) -> bool:
    return expected_prediction_outcome(candidate) == observed


def classify_reference_failure(failure_reason: str | None) -> str:
    return REFERENCE_BACKEND_FAILURE
